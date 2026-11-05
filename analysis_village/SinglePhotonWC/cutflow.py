#!/usr/bin/env python3
"""Cut flow of the reco preselection (cuts.PRESELECTION) by truth category.

Input: configs/sp_mc.py outputs (wc + mcnu). Candidates join to the truth on
(file, __ntuple, entry, tmatch_idx); unmatched ones are cosmic.
  efficiency  true signal neutrinos with >= 1 passing candidate matched to
              them / all true signal neutrinos
  purity      signal candidates / passing candidates

  python analysis_village/SinglePhotonWC/cutflow.py -i <dir or .df files>
"""

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from analysis_village.SinglePhotonWC import cuts  # noqa: E402
from analysis_village.NCDeltaWC import dfio  # noqa: E402

NU_KEY = ["file", "__ntuple", "entry"]


def load(inputs, extra=()):
    """(candidates with truth and cat, mcnu with cat, POT). extra: more
    candidate-indexed keys to join, e.g. "wctag"."""
    frames, totals = dfio.load(dfio.expand_inputs(inputs), keys=["wc", "mcnu"] + list(extra))
    w, m = frames["wc"], frames["mcnu"]
    for k in extra:
        w = w.join(frames[k])
    m.index = m.index.set_names(NU_KEY + ["tmatch_idx"])
    m["cat"] = cuts.truth_category(m)
    w = w.join(m.drop(columns="cat"), on=NU_KEY + ["tmatch_idx"])
    w["cat"] = cuts.truth_category(w)
    return w, m, totals.pot.sum()


def cutflow(w, m, pot):
    cats = list(cuts.CATEGORY_LABELS)
    sig_nu = m[m.cat.isin(cuts.SIGNAL_CATEGORIES)]
    rows, mask = [], pd.Series(True, index=w.index)
    steps = [("all candidates", None)] + list(cuts.PRESELECTION.items())
    for name, cut in steps:
        if cut is not None:
            mask &= cut(w)
        p = w[mask]
        n = p.cat.value_counts().reindex(cats).fillna(0) / pot * 1e20
        nsig = p.cat.isin(cuts.SIGNAL_CATEGORIES)
        matched = p[nsig].reset_index().set_index(NU_KEY + ["tmatch_idx"]).index.unique()
        rows.append(dict(cut=name, **{cuts.CATEGORY_LABELS[c]: n[c] for c in cats},
                         total=n.sum(), eff=100 * len(matched) / len(sig_nu),
                         purity=100 * nsig.mean() if len(p) else 0))
    return pd.DataFrame(rows).set_index("cut")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--inputs", nargs="+", required=True, help=".df files, directories or globs")
    args = ap.parse_args()
    w, m, pot = load(args.inputs)
    t = cutflow(w, m, pot)
    pd.set_option("display.width", 250)
    print("Candidates per 1e20 POT (sample %.3g POT); eff and purity in %%" % pot)
    print(t.T.round(1).to_string())


if __name__ == "__main__":
    main()
