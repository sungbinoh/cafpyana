#!/usr/bin/env python3
"""Cut flows (cuts.CHAINS: our preselection, the MicroBooNE reference
selection, or both) by truth category.

Input: configs/sp_mc*.py outputs (wc + mcnu). Candidates join to the truth on
(file, __ntuple, entry, tmatch_idx); unmatched ones are cosmic. Event
cuts (cuts.EVENT_CUTS, e.g. wc_valid from wcevt) come first and are rows
of true neutrinos, before "all candidates".
  efficiency  true signal neutrinos with >= 1 passing candidate matched to
              them / all true signal neutrinos
  purity      signal candidates / passing candidates
eff_err and purity_err are binomial, from the MC counts.

syst_unc() gives the fractional uncertainty on a selected rate from the
systematic weights (wgt / slimwgt tables, make_weights.py), per bin:
  multisim (univ_<i>)   rms of the universes about the nominal rate
  +/- sigma unisim      max(|ps1 - cv|, |ms1 - cv|)
  morph unisim          |morph - nominal|
summed in quadrature over the knobs of each group (flux, GENIE, Geant4).

  python analysis_village/SinglePhotonWC/cutflow.py -i <dir or .df files> \
      [-c presel uboone presel+uboone] [--uboone-bdt <xml dir>]
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from analysis_village.SinglePhotonWC import cuts, dfio  # noqa: E402

NU_KEY = ["file", "__ntuple", "entry", "tmatch_idx"]


def load(inputs, uboone_bdt=None, nproc=None):
    """(candidates with truth and cat, mcnu with cat, POT). With uboone_bdt,
    the candidates also get the MicroBooNE BDT scores (ub_*)."""
    frames, totals = dfio.load(dfio.expand_inputs(inputs), keys=["wc", "mcnu", "wcevt"], nproc=nproc)
    w, m, ev = frames["wc"], frames["mcnu"], frames["wcevt"]
    m.index = m.index.set_names(NU_KEY)
    if len(ev):
        m["wc_valid"] = ev.valid.reindex(m.index.droplevel("tmatch_idx")).fillna(0).astype(int).values
    else:
        print("Warning: no wcevt table (older inputs); taking every event as wc_valid.")
        m["wc_valid"] = 1
    m["cat"] = cuts.truth_category(m)
    w = w.join(m.drop(columns=["cat", "wc_valid"]), on=NU_KEY)
    w["wc_valid"] = 1  # a candidate implies Wire-Cell output
    w["cat"] = cuts.truth_category(w)
    if uboone_bdt:
        from analysis_village.SinglePhotonWC import microboone_bdt
        w = w.join(microboone_bdt.score(w, uboone_bdt, nproc=nproc))
    return w, m, totals.pot.sum()


WEIGHT_KEYS = ("wgt", "slimwgt")   # full / slim systematic weights, index as mcnu
SYST_GROUPS = {"flux": "flux", "genie": "cross section (GENIE)", "g4": "hadron reinteractions (Geant4)"}


def load_weights(inputs, keep=None, nproc=None):
    """The wgt table (else slimwgt) indexed like mcnu (NU_KEY), None if the
    inputs have neither. keep: an index of NU_KEY tuples to keep (rows are
    selected per file, so the full table is never in memory)."""
    files = dfio.expand_inputs(inputs)
    key = next((k for k in WEIGHT_KEYS if k in dfio.file_keys(files[0])), None)
    if key is None:
        return None

    def select(k, d):
        if k != key:
            return d
        d.index = d.index.set_names(NU_KEY)
        return d if keep is None else d[d.index.isin(keep)]
    frames, _ = dfio.load(files, keys=[key], nproc=nproc, select=select)
    print("systematic weights: '%s' table, %d neutrinos x %d weights" % (key, *frames[key].shape))
    return frames[key]


def syst_group(name):
    if name == "Flux" or name.endswith("_Flux"):
        return "flux"
    if name.startswith("GENIE"):
        return "genie"
    if "Geant4" in name:
        return "g4"
    return name


def syst_unc(sel, wgt, groups):
    """Fractional uncertainty per group on the rate of the candidates sel,
    summed within each label of groups (a Series over sel, e.g. the energy
    bin): DataFrame (label x group) plus "syst" (in quadrature) and "n" (MC
    candidates). Candidates without a weight row (cosmics) weigh 1."""
    W = wgt.reindex(pd.MultiIndex.from_frame(sel.reset_index()[NU_KEY].fillna(-1).astype(int))).fillna(1.0)
    W.index = sel.index
    N = W.groupby(groups.to_numpy()).sum()
    n = groups.value_counts()
    nom = n.reindex(N.index).astype(float)
    var = {}
    for s in N.columns.get_level_values(0).unique():
        v = N[s]
        if "cv" in v:
            d = np.maximum(abs(v.get("ps1", v.cv) - v.cv), abs(v.get("ms1", v.cv) - v.cv)) ** 2
        elif "morph" in v:
            d = (v.morph - nom) ** 2
        else:
            d = v.sub(nom, axis=0).pow(2).mean(axis=1)
        g = syst_group(s)
        var[g] = var.get(g, 0) + d
    out = pd.DataFrame({g: np.sqrt(var[g]) / nom for g in var})
    out["syst"] = np.sqrt(sum(var.values())) / nom
    out["n"] = nom
    return out


def n_signal_nu(p):
    """True signal neutrinos with at least one candidate in p."""
    s = p[cuts.is_signal_cat(p.cat)]
    return len(s.reset_index().drop_duplicates(NU_KEY))


def _binom(k, n):
    return 100 * np.sqrt(k / n * (1 - k / n) / n) if n else np.nan


def _row(name, p, nsel, nsig, pot):
    n = p.cat.value_counts().reindex(list(cuts.CATEGORY_LABELS)).fillna(0) / pot * 1e20
    nsp = cuts.is_signal_cat(p.cat).sum()
    return dict(cut=name, **n.rename(cuts.CATEGORY_LABELS), total=n.sum(),
                eff=100 * nsel / nsig if nsig else np.nan, eff_err=_binom(nsel, nsig),
                purity=100 * nsp / len(p) if len(p) else 0, purity_err=_binom(nsp, len(p)))


def cutflow(w, m, pot, chain):
    """Rates per 1e20 POT, eff and purity after each cut. Event-cut rows
    count true neutrinos (keeping them, before candidates); the rest count
    candidates. Every efficiency is over all true signal neutrinos."""
    nsig = cuts.is_signal_cat(m.cat).sum()
    rows, keep = [], pd.Series(True, index=m.index)
    events = [k for k in chain if k in cuts.EVENT_CUTS]
    for name in events:
        keep &= chain[name](m)
        rows.append(_row(name, m[keep], cuts.is_signal_cat(m.cat[keep]).sum(), nsig, pot))
    mask = pd.Series(True, index=w.index)
    for name, cut in [("all candidates", None)] + [(k, c) for k, c in chain.items() if k not in events]:
        if cut is not None:
            mask &= cut(w)
        p = w[mask]
        rows.append(_row(name, p, n_signal_nu(p), nsig, pot))
    return pd.DataFrame(rows).set_index("cut")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--inputs", nargs="+", required=True, help=".df files, directories or globs")
    ap.add_argument("-c", "--chain", nargs="+", default=["presel"], choices=list(cuts.CHAINS))
    ap.add_argument("--uboone-bdt", help="directory of the MicroBooNE BDT XMLs (needed by the uboone chains)")
    args = ap.parse_args()
    w, m, pot = load(args.inputs, args.uboone_bdt)
    pd.set_option("display.width", 250)
    for c in args.chain:
        print("\n[%s] candidates per 1e20 POT (sample %.3g POT); eff and purity in %%" % (c, pot))
        print(cutflow(w, m, pot, cuts.CHAINS[c]).T.round(1).to_string())


if __name__ == "__main__":
    main()
