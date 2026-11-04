#!/usr/bin/env python3
"""Train one BDT per channel, after the manual selection.

Pipeline: topology + pre-selection + energy window + manual_selection, then
the BDT. The BDT replaces the cut-based stage; it does not sit on top of it.

Files split 50/50 by parity of position in the sorted input list. Slices from
one file share an event, so splitting by file keeps the halves independent.

  TRAIN  even files, restricted to true NC Delta->N gamma (signal) and true
         NC 1pi0 (background). This is the discrimination that matters and it
         removes the ~45k easy slices that swamped the loss before.
  TEST   odd files, every category that survives the manual selection -- the
         real background mix the selection has to cope with.

  python train_bdt.py --cache <pickle> -o out/
"""

import argparse
import os
import pickle
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.utils.class_weight import compute_sample_weight

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import cuts                                                    # noqa: E402
import shower_energy                                           # noqa: E402
from dfio import expand_inputs, load                           # noqa: E402
from features import feature_columns, derived, sanitize        # noqa: E402
from manual_selection import stage                             # noqa: E402


def topology_only(df, topo, cuts):
    """Topology alone. The dataframes are built with pre-selection applied, so
    FV and the clear-cosmic veto are baked in and cannot be loosened here."""
    return (cuts.topology_1g1p(df) if topo == "1g1p"
            else cuts.topology_1g0p(df)).to_numpy()

SIGNAL = 0
NC_PI0 = (1, 2)            # coherent and non-coherent

PARAMS = {
    "objective": "binary:logistic",
    "eval_metric": ["auc", "logloss", "aucpr"],   # last drives early stopping
    "device": "cpu",
    "max_depth": 3,
    "eta": 0.05,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "min_child_weight": 5,
    "reg_lambda": 2.0,
}
NUM_ROUNDS = 300
EARLY_STOP = 30


def build(df, topo, pre="manual"):
    """(X, y, category, parity) after `pre`: "manual" or "topology"."""
    m = stage(df, topo, cuts) if pre == "manual" else topology_only(df, topo, cuts)
    d = df[m]
    cols = feature_columns(d)
    X = d[cols].astype(float).copy()
    for name, v in derived(d, cuts, shower_energy).items():
        X[name] = v
    X = sanitize(X)
    cat = cuts.truth_category(d).to_numpy()
    par = np.asarray(d.index.get_level_values("file")) % 2
    return X, (cat == SIGNAL).astype(int), cat, par, d


def learning_curve(hist, best, topo, out):
    metrics = list(hist["train"])
    fig, axes = plt.subplots(1, len(metrics), figsize=(5.2 * len(metrics), 4.2))
    for ax, m in zip(np.atleast_1d(axes), metrics):
        for split, c in (("train", "#4C78A8"), ("test", "#E45756")):
            ax.plot(hist[split][m], color=c, lw=1.8, label=split)
        ax.axvline(best, ls="--", c="k", lw=1.2)
        ax.set_xlabel("boosting round"); ax.set_ylabel(m)
        ax.set_title("%s  %s  (best round %d)" % (topo, m, best), fontsize=10)
        ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out, dpi=110, bbox_inches="tight"); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--inputs", nargs="+")
    ap.add_argument("--cache")
    ap.add_argument("-o", "--outdir", default="bdt_out")
    ap.add_argument("--weights", choices=("balanced", "none"), default="none")
    ap.add_argument("--pre", choices=("manual", "topology"), default="manual",
                    help="what to apply before the BDT. topology drops the "
                         "energy window and the manual selection.")
    ap.add_argument("--test-bkg", choices=("all", "ncpi0"), default="all",
                    help="ncpi0 restricts the TEST set the same way, which "
                         "measures the pure signal-vs-NC-1pi0 ceiling rather "
                         "than performance against the real background mix.")
    ap.add_argument("--train-bkg", choices=("all", "ncpi0"), default="all",
                    help="all: signal vs every surviving background (default). "
                         "ncpi0: signal vs true NC 1pi0 only -- tried and much "
                         "worse, the model stops rejecting CC nue because it "
                         "never sees an electron shower.")
    args = ap.parse_args()

    slc = (pickle.load(open(args.cache, "rb"))["slc"] if args.cache
           else load(expand_inputs(args.inputs), keys=("slc",))[0]["slc"])
    os.makedirs(args.outdir, exist_ok=True)

    for topo in cuts.TOPOLOGIES:
        X, y, cat, par, d = build(slc, topo, pre=args.pre)
        # train: even files, signal vs NC 1pi0 only
        tr = par == 0
        if args.train_bkg == "ncpi0":
            tr = tr & ((cat == SIGNAL) | np.isin(cat, NC_PI0))
        te = par == 1
        if args.test_bkg == "ncpi0":
            te = te & ((cat == SIGNAL) | np.isin(cat, NC_PI0))
        print("\n=== %s  %d features  (pre: %s, train bkg: %s, test bkg: %s)"
              % (topo, X.shape[1], args.pre, args.train_bkg, args.test_bkg))
        print("   train (even): %d signal, %d background"
              % ((y[tr] == 1).sum(), (y[tr] == 0).sum()))
        print("   test  (odd):  %d signal, %d background"
              % ((y[te] == 1).sum(), (y[te] == 0).sum()), flush=True)

        w = compute_sample_weight("balanced", y[tr]) if args.weights == "balanced" else None
        dtr = xgb.DMatrix(X[tr], label=y[tr], weight=w, feature_names=list(X.columns))
        dte = xgb.DMatrix(X[te], label=y[te], feature_names=list(X.columns))
        hist = {}
        bst = xgb.train(PARAMS, dtr, NUM_ROUNDS, evals=[(dtr, "train"), (dte, "test")],
                        evals_result=hist, early_stopping_rounds=EARLY_STOP,
                        verbose_eval=False)
        learning_curve(hist, bst.best_iteration, topo,
                       os.path.join(args.outdir, "learning_curve_%s.png" % topo))
        s_tr, s_te = bst.predict(dtr), bst.predict(dte)
        print("   AUC   train %.4f test %.4f  |  AUCPR train %.4f test %.4f  (round %d)"
              % (roc_auc_score(y[tr], s_tr), roc_auc_score(y[te], s_te),
                 average_precision_score(y[tr], s_tr),
                 average_precision_score(y[te], s_te), bst.best_iteration))

        # head to head with the cut-based selection, on the test half
        yt, st = y[te], s_te
        nbkg = int((yt == 0).sum())
        print("   %-28s %5s %7s %9s" % ("", "sig", "bkg", "purity"))
        for name in cuts.SELECTIONS:
            sel = cuts.selection(d, topo, cutbased=cuts.selection_config(name, topo)).to_numpy()[te]
            S, B = int((sel & (yt == 1)).sum()), int((sel & (yt == 0)).sum())
            print("   %-28s %5d %7d %8.2f%%" % ("cuts: " + name, S, B, 100 * S / max(S + B, 1)))
            if S:
                thr = np.sort(st[yt == 1])[::-1][S - 1]
                k = st >= thr
                Sb, Bb = int((k & (yt == 1)).sum()), int((k & (yt == 0)).sum())
                print("   %-28s %5d %7d %8.2f%%" % ("   BDT at the same yield", Sb, Bb,
                                                    100 * Sb / max(Sb + Bb, 1)))
        print("   signal efficiency vs background rejection on the test half:")
        for eff in (0.95, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40):
            thr = np.quantile(st[yt == 1], 1. - eff)
            rej = float(np.mean(st[yt == 0] < thr))
            kept = int(round((1 - rej) * nbkg))
            print("      eff %2.0f%% -> rejects %5.1f%% of background (%d of %d left)"
                  % (100 * eff, 100 * rej, kept, nbkg))
        imp = sorted(bst.get_score(importance_type="gain").items(), key=lambda kv: -kv[1])[:10]
        print("   top features: " + ", ".join("%s %.0f" % kv for kv in imp))
        bst.save_model(os.path.join(args.outdir, "bdt_%s.json" % topo))
    print("\nwrote models to", args.outdir)


if __name__ == "__main__":
    main()
