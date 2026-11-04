#!/usr/bin/env python3
"""Score the saved BDTs against the cut-based selection on the test half.

At the cut-based selection's own signal efficiency, which keeps less
background? That is the only comparison that matters; AUC over the whole
pre-selection sample is dominated by easy cosmics.
"""
import os, pickle, sys
import numpy as np, xgboost as xgb
from sklearn.metrics import roc_auc_score
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import cuts, shower_energy                                    # noqa: E402
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from features import feature_columns, derived, sanitize       # noqa: E402
from train_bdt import build                                   # noqa: E402

slc = pickle.load(open(sys.argv[1], "rb"))["slc"]
for topo in cuts.TOPOLOGIES:
    X, y, par, d = build(slc, topo)
    te = par == 1
    bst = xgb.Booster(); bst.load_model("bdt_out/bdt_%s.json" % topo)
    s = bst.predict(xgb.DMatrix(X[te], feature_names=list(X.columns)))
    yt = y[te]
    nsig, nbkg = int(yt.sum()), int((yt == 0).sum())
    print("=== %s  test half: %d signal, %d background  (AUC %.4f)"
          % (topo, nsig, nbkg, roc_auc_score(yt, s)))
    print("   %-14s %7s %8s %9s %9s" % ("", "sig", "bkg", "purity", "bkg rej"))
    for name in cuts.SELECTIONS:
        sel = cuts.selection(d, topo, cutbased=cuts.selection_config(name, topo)).values[te]
        S, B = int((sel & (yt == 1)).sum()), int((sel & (yt == 0)).sum())
        eff = S / max(nsig, 1)
        print("   %-14s %7d %8d %8.2f%% %8.2f%%   <- cuts, eff %.1f%%"
              % (name, S, B, 100 * S / max(S + B, 1), 100 * (1 - B / nbkg), 100 * eff))
        # BDT threshold that keeps the same number of signal events
        thr = np.sort(s[yt == 1])[::-1][S - 1] if S else 1.
        k = s >= thr
        Sb, Bb = int((k & (yt == 1)).sum()), int((k & (yt == 0)).sum())
        print("   %-14s %7d %8d %8.2f%% %8.2f%%   <- BDT at the same signal yield"
              % ("", Sb, Bb, 100 * Sb / max(Sb + Bb, 1), 100 * (1 - Bb / nbkg)))
    print()
