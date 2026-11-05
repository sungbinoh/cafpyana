#!/usr/bin/env python3
"""Cut flows (cuts.CHAINS: our preselection, the MicroBooNE reference
selection, or both) by truth category.

Input: configs/sp_mc.py outputs (wc + mcnu). Candidates join to the truth on
(file, __ntuple, entry, tmatch_idx); unmatched ones are cosmic.
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
    frames, totals = dfio.load(dfio.expand_inputs(inputs), keys=["wc", "mcnu"])
    w, m = frames["wc"], frames["mcnu"]
    m.index = m.index.set_names(NU_KEY)
    m["cat"] = cuts.truth_category(m)
    w = w.join(m.drop(columns="cat"), on=NU_KEY)
    w["cat"] = cuts.truth_category(w)
    if uboone_bdt:
        from analysis_village.SinglePhotonWC import microboone_bdt
        w = w.join(microboone_bdt.score(w, uboone_bdt))
    return w, m, totals.pot.sum()


def n_signal_nu(p):
    """True signal neutrinos with at least one candidate in p."""
    s = p[cuts.is_signal_cat(p.cat)]
    return len(s.reset_index().drop_duplicates(NU_KEY))


def cutflow(w, m, pot, chain):
    nsig = cuts.is_signal_cat(m.cat).sum()
    rows, mask = [], pd.Series(True, index=w.index)
    for name, cut in [("all candidates", None)] + list(chain.items()):
        if cut is not None:
            mask &= cut(w)
        p = w[mask]
        n = p.cat.value_counts().reindex(list(cuts.CATEGORY_LABELS)).fillna(0) / pot * 1e20
        rows.append(dict(cut=name, **n.rename(cuts.CATEGORY_LABELS), total=n.sum(),
                         eff=100 * n_signal_nu(p) / nsig if nsig else np.nan,
                         purity=100 * cuts.is_signal_cat(p.cat).mean() if len(p) else 0))
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
