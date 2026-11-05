#!/usr/bin/env python3
"""BDT (XGBoost) separating the inclusive single-photon signal from every
background, on candidates passing cuts.PRESELECTION.

Inputs are reco only: every numeric Wire-Cell candidate column except
identifiers, MC-only truth-match columns, absolute position and time
(see feature_columns). With --tagger, also the wctag table (every scalar
T_tagger leaf, configs/sp_mc_bdt.py).

Two-fold by file: each half of the files is scored by a model trained on
the other half, so every candidate gets an unbiased score. A tenth of the
training files is held out for early stopping.

Efficiency is per true signal neutrino (as cutflow.py); purity per candidate.

  python analysis_village/SinglePhotonWC/train_bdt.py -i <.df files or dirs> -o bdt_out/ [--tagger]
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import xgboost as xgb
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from analysis_village.SinglePhotonWC import cuts  # noqa: E402
from analysis_village.SinglePhotonWC.cutflow import load, NU_KEY  # noqa: E402

# Not inputs: identifiers, MC-only truth match, detector position and time
# (the vertex enters as distance to the wall), sample bookkeeping.
EXCLUDE = re.compile(
    r"^(run|subrun|evt|n_nu|cluster_id|matched_flash_gid|kine_main_vertex_id|"
    r"tmatch_.*|pf_.*|nu_[xyz]|kine_nu_[xyz]_corr|flash_time_us|"
    r"tag_(run|subrun|event|nu_[xyz]|.*cluster_id|.*flash.*|.*time.*|match_.*))$")
PARAMS = dict(objective="binary:logistic", eval_metric="aucpr", tree_method="hist",
              max_depth=6, learning_rate=0.05, subsample=0.8, colsample_bytree=0.5,
              min_child_weight=5, nthread=os.cpu_count())
GROUP_COLORS = [("Other backgrounds", "#a3a29d"), ("CC numu (mu >= 100 MeV)", "#eda100"),
                ("Signal (inclusive 1g)", "#2a78d6")]


def feature_columns(df):
    num = df.select_dtypes(include=[np.number, bool]).columns
    cols = [c for c in num if not c.startswith("true_") and c != "cat" and not EXCLUDE.match(c)]
    return [c for c in cols if df[c].nunique(dropna=False) > 1]


def add_derived(df):
    df["vtx_dwall"] = cuts.dist_to_active_boundary(df.nu_x, df.nu_y, df.nu_z)
    df["E_em_frac"] = df.E_em_1 / df.kine_reco_Enu
    return df


def train_fold(X, y, train, seed):
    files = X.index.get_level_values("file")
    rng = np.random.default_rng(seed)
    tr_files = np.unique(files[train])
    val_files = rng.choice(tr_files, size=max(1, len(tr_files) // 10), replace=False)
    val = train & np.isin(files, val_files)
    fit = train & ~val
    dfit = xgb.DMatrix(X[fit], label=y[fit])
    dval = xgb.DMatrix(X[val], label=y[val])
    return xgb.train(PARAMS, dfit, num_boost_round=2000, evals=[(dval, "val")],
                     early_stopping_rounds=100, verbose_eval=False)


def working_points(p, nsig_nu, pot, effs=(0.40, 0.30, 0.20, 0.10, 0.05)):
    """Score cut for each target efficiency, with what survives."""
    sig = p.cat.isin(cuts.SIGNAL_CATEGORIES)
    s = p[sig].sort_values("bdt", ascending=False)
    first = ~s.reset_index().duplicated(NU_KEY + ["tmatch_idx"]).values  # one per true neutrino
    eff_curve = np.cumsum(first) / nsig_nu
    rows = []
    for e in effs:
        if eff_curve[-1] < e:
            continue
        cut = s.bdt.values[np.searchsorted(eff_curve, e)]
        q = p[p.bdt >= cut]
        qs = q.cat.isin(cuts.SIGNAL_CATEGORIES)
        sc = 1e20 / pot
        rows.append(dict(eff=100 * e, bdt_cut=cut, purity=100 * qs.mean(), sig=qs.sum() * sc,
                         bkg=(~qs).sum() * sc, numucc=(q.cat == 6).sum() * sc,
                         nc_2shw=(q.cat == 4).sum() * sc, zero_shw=(q.cat == 5).sum() * sc,
                         out_av=(q.cat == 8).sum() * sc, cosmic=(q.cat == 9).sum() * sc,
                         nue=(q.cat == 7).sum() * sc))
    return pd.DataFrame(rows)


def plots(p, imp, wp, nsig_nu, outdir):
    plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "font.size": 9,
                         "axes.grid": True, "grid.color": "#e6e5e1", "axes.axisbelow": True,
                         "legend.frameon": False})
    sig = p.cat.isin(cuts.SIGNAL_CATEGORIES)
    G = [p[~sig & (p.cat != 6)], p[p.cat == 6], p[sig]]
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.6))
    bins = np.linspace(0, 1, 51)
    axs[0].hist([g.bdt for g in G], bins=bins, stacked=True, color=[c for _, c in GROUP_COLORS],
                label=["%s (%d)" % (n, len(g)) for (n, _), g in zip(GROUP_COLORS, G)])
    axs[0].hist(G[2].bdt, bins=bins, histtype="step", lw=1.8, ls="--", color="#2a78d6",
                weights=np.full(len(G[2]), 20.), label="Signal x 20")
    axs[0].set_yscale("log")
    axs[0].set_xlabel("BDT score (out-of-fold)")
    axs[0].set_ylabel("preselected candidates")
    axs[0].set_title("BDT score", loc="left")
    axs[0].legend(fontsize=8)

    s = p[sig].sort_values("bdt", ascending=False)
    first = ~s.reset_index().duplicated(NU_KEY + ["tmatch_idx"]).values
    cuts_ = np.quantile(p.bdt, np.linspace(0.5, 0.99999, 300))
    eff = [first[s.bdt.values >= c].sum() / nsig_nu for c in cuts_]
    pur = [sig[p.bdt >= c].mean() for c in cuts_]
    axs[1].plot(100 * np.array(eff), 100 * np.array(pur), color="#2a78d6", lw=2)
    axs[1].set_xlabel("signal efficiency [% of all true signal]")
    axs[1].set_ylabel("purity [%]")
    axs[1].set_title("Efficiency vs purity", loc="left")

    top = imp.head(25)[::-1]
    axs[2].barh(top.index, top.values, color="#2a78d6", height=0.6)
    axs[2].set_xlabel("total gain (summed over folds), normalised")
    axs[2].set_title("Top 25 inputs", loc="left")
    axs[2].tick_params(axis="y", labelsize=7)
    axs[2].grid(axis="y", visible=False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "bdt_summary.png"), dpi=130)
    fig.savefig(os.path.join(outdir, "bdt_summary.pdf"))
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--inputs", nargs="+", required=True)
    ap.add_argument("-o", "--outdir", default="bdt_out")
    ap.add_argument("--tagger", action="store_true", help="also use the wctag (T_tagger) table")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    w, m, pot = load(args.inputs, extra=["wctag"] if args.tagger else [])
    nsig_nu = m.cat.isin(cuts.SIGNAL_CATEGORIES).sum()
    p = add_derived(w[cuts.preselection(w)].copy())
    feats = feature_columns(p)
    print("%d preselected candidates, %d signal; %d input variables" % (
        len(p), p.cat.isin(cuts.SIGNAL_CATEGORIES).sum(), len(feats)))
    X = p[feats].astype(np.float32)
    y = p.cat.isin(cuts.SIGNAL_CATEGORIES).astype(int).values
    fold = (p.index.get_level_values("file") % 2).values

    p["bdt"] = np.nan
    imp = pd.Series(0., index=feats)
    for k in (0, 1):
        bst = train_fold(X, y, fold != k, seed=k)
        p.loc[fold == k, "bdt"] = bst.predict(xgb.DMatrix(X[fold == k]), iteration_range=(0, bst.best_iteration + 1))
        imp = imp.add(pd.Series(bst.get_score(importance_type="total_gain")), fill_value=0)
        bst.save_model(os.path.join(args.outdir, "bdt_fold%d.json" % k))
        print("fold %d: %d trees, AUC %.4f" % (k, bst.best_iteration + 1, roc_auc_score(y[fold == k], p.bdt[fold == k])))
    imp = (imp / imp.sum()).sort_values(ascending=False)

    print("AUC (all, out-of-fold) %.4f" % roc_auc_score(y, p.bdt))
    pre_eff = p[y == 1].reset_index().drop_duplicates(NU_KEY + ["tmatch_idx"]).shape[0] / nsig_nu
    print("Preselection efficiency %.1f%%, purity %.2f%%" % (100 * pre_eff, 100 * y.mean()))
    wp = working_points(p, nsig_nu, pot)
    pd.set_option("display.width", 200)
    print("Working points (rates per 1e20 POT):")
    print(wp.round(3).to_string(index=False))
    print("Top inputs:")
    print((100 * imp.head(20)).round(1).to_string())

    p[["cat", "bdt"]].to_pickle(os.path.join(args.outdir, "scores.pkl"))
    imp.to_csv(os.path.join(args.outdir, "importance.csv"))
    wp.to_csv(os.path.join(args.outdir, "working_points.csv"), index=False)
    plots(p, imp, wp, nsig_nu, args.outdir)
    print("Wrote %s" % args.outdir)


if __name__ == "__main__":
    main()
