#!/usr/bin/env python3
"""Cut flows (cuts.CHAINS: our preselection, the MicroBooNE reference
selection, or both) by truth category.

Input: configs/sp_mc.py outputs (wc + mcnu). Candidates join to the truth on
(file, __ntuple, entry, tmatch_idx); unmatched ones are cosmic. Event
cuts (cuts.EVENT_CUTS, e.g. wc_valid from wcevt) come first and are rows
of true neutrinos, before "all candidates".
  efficiency  true signal neutrinos with >= 1 passing candidate matched to
              them / all true signal neutrinos
  purity      signal candidates / passing candidates

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


def load(inputs, uboone_bdt=None):
    """(candidates with truth and cat, mcnu with cat, POT). With uboone_bdt,
    the candidates also get the MicroBooNE BDT scores (ub_*)."""
    frames, totals = dfio.load(dfio.expand_inputs(inputs), keys=["wc", "mcnu", "wcevt"])
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
        w = w.join(microboone_bdt.score(w, uboone_bdt))
    return w, m, totals.pot.sum()


def n_signal_nu(p):
    """True signal neutrinos with at least one candidate in p."""
    s = p[cuts.is_signal_cat(p.cat)]
    return len(s.reset_index().drop_duplicates(NU_KEY))


def _row(name, p, nsel, nsig, pot):
    n = p.cat.value_counts().reindex(list(cuts.CATEGORY_LABELS)).fillna(0) / pot * 1e20
    return dict(cut=name, **n.rename(cuts.CATEGORY_LABELS), total=n.sum(),
                eff=100 * nsel / nsig if nsig else np.nan,
                purity=100 * cuts.is_signal_cat(p.cat).mean() if len(p) else 0)


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
