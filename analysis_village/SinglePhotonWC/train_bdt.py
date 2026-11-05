#!/usr/bin/env python3
"""BDT (XGBoost) separating the inclusive single-photon signal from every
background, on candidates passing cuts.PRESELECTION.

Inputs are reco only: every numeric Wire-Cell candidate column (summary
columns plus every scalar T_tagger leaf) except identifiers, MC-only
truth-match columns, absolute position and time (EXCLUDE). --uboone-bdt adds
the four MicroBooNE BDT scores; --features uboone keeps only the MicroBooNE
BDT inputs, i.e. retrains them on SBND.

Two-fold by file: each half of the files is scored by a model trained on
the other half, so every candidate gets an unbiased score. A tenth of the
training files is held out for early stopping.

Efficiency is per true signal neutrino (as cutflow.py); purity per candidate.

  python analysis_village/SinglePhotonWC/train_bdt.py -i <.df files or dirs> -o bdt_out/ \
      [--uboone-bdt <xml dir>] [--features uboone]
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
from analysis_village.SinglePhotonWC.cutflow import load, n_signal_nu  # noqa: E402

# Not inputs: identifiers, MC-only truth match, detector position and time
# (the vertex enters as distance to the wall).
EXCLUDE = re.compile(
    r"^(run|subrun|evt|n_nu|cluster_id|matched_flash_gid|kine_main_vertex_id|"
    r"tmatch_.*|nu_[xyz]|kine_nu_[xyz]_corr|flash_time_us|flash_group|vertex_moved_cluster)$")
PARAMS = dict(objective="binary:logistic", eval_metric="aucpr", tree_method="hist",
              max_depth=6, learning_rate=0.05, subsample=0.8, colsample_bytree=0.5,
              min_child_weight=5, nthread=os.cpu_count())
EFFS = [0.40, 0.30, 0.20, 0.10, 0.05]


def feature_columns(df, only=None):
    cols = [c for c in df.select_dtypes(include=[np.number, bool]).columns
            if not c.startswith("true_") and c != "cat" and not EXCLUDE.match(c)]
    if only is not None:
        cols = [c for c in cols if c in only]
    return [c for c in cols if df[c].nunique(dropna=False) > 1]


def train_fold(X, y, train, seed):
    files = X.index.get_level_values("file")
    tr_files = np.unique(files[train])
    val_files = np.random.default_rng(seed).choice(tr_files, size=max(1, len(tr_files) // 10), replace=False)
    val = train & np.isin(files, val_files)
    fit = train & ~val
    return xgb.train(PARAMS, xgb.DMatrix(X[fit], label=y[fit]), num_boost_round=2000,
                     evals=[(xgb.DMatrix(X[val], label=y[val]), "val")],
                     early_stopping_rounds=100, verbose_eval=False)


def eff_purity(p, nsig_nu, cut):
    q = p[p.bdt >= cut]
    return n_signal_nu(q) / nsig_nu, cuts.is_signal_cat(q.cat).mean()


def working_points(p, nsig_nu, pot):
    """For each target efficiency, the score cut and what survives (per 1e20 POT)."""
    sig = cuts.is_signal_cat(p.cat)
    grid = np.unique(np.quantile(p.bdt[sig], np.linspace(0, 1, 501)))
    effs = np.array([eff_purity(p, nsig_nu, c)[0] for c in grid])
    rows = []
    for e in EFFS:
        if effs.max() < e:
            continue
        cut = grid[effs >= e].max()
        q = p[p.bdt >= cut]
        n = q.cat.value_counts().reindex(list(cuts.CATEGORY_LABELS)).fillna(0) * 1e20 / pot
        rows.append(dict(eff=100 * eff_purity(p, nsig_nu, cut)[0], bdt_cut=cut,
                         purity=100 * cuts.is_signal_cat(q.cat).mean(), **n.rename(cuts.CATEGORY_LABELS)))
    return pd.DataFrame(rows)


def plots(p, imp, nsig_nu, outdir):
    sig = cuts.is_signal_cat(p.cat)
    groups = [("Other backgrounds", p[~sig & (p.cat != 6)], "#a3a29d"),
              ("CC numu (mu >= 100 MeV)", p[p.cat == 6], "#eda100"),
              ("Signal (inclusive 1g)", p[sig], "#2a78d6")]
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.6))
    bins = np.linspace(0, 1, 51)
    axs[0].hist([g.bdt for _, g, _ in groups], bins=bins, stacked=True, color=[c for *_, c in groups],
                label=["%s (%d)" % (n, len(g)) for n, g, _ in groups])
    axs[0].hist(p.bdt[sig], bins=bins, histtype="step", lw=1.8, ls="--", color="#2a78d6",
                weights=np.full(sig.sum(), 20.), label="Signal x 20")
    axs[0].set(yscale="log", xlabel="BDT score (out-of-fold)", ylabel="preselected candidates")
    axs[0].set_title("BDT score", loc="left")
    axs[0].legend(fontsize=8)

    ep = np.array([eff_purity(p, nsig_nu, c) for c in np.quantile(p.bdt, np.linspace(0.5, 0.9999, 150))])
    axs[1].plot(100 * ep[:, 0], 100 * ep[:, 1], color="#2a78d6", lw=2)
    axs[1].set(xlabel="signal efficiency [% of all true signal]", ylabel="purity [%]")
    axs[1].set_title("Efficiency vs purity", loc="left")

    top = imp.head(25)[::-1]
    axs[2].barh(top.index, top.values, color="#2a78d6", height=0.6)
    axs[2].set_xlabel("total gain (summed over folds), normalised")
    axs[2].set_title("Top 25 inputs", loc="left")
    axs[2].tick_params(axis="y", labelsize=7)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, "bdt_summary." + ext), dpi=130)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--inputs", nargs="+", required=True)
    ap.add_argument("-o", "--outdir", default="bdt_out")
    ap.add_argument("--uboone-bdt", help="MicroBooNE BDT XML directory: adds the four ub_* scores")
    ap.add_argument("--features", choices=["all", "uboone"], default="all",
                    help="uboone: only the inputs of the four MicroBooNE BDTs (needs --uboone-bdt)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    w, m, pot = load(args.inputs, args.uboone_bdt)
    nsig_nu = cuts.is_signal_cat(m.cat).sum()
    p = w[cuts.passes(w)].copy()
    p["vtx_dwall"] = cuts.dist_to_active_boundary(p.nu_x, p.nu_y, p.nu_z)
    p["E_em_frac"] = p.E_em_1 / p.kine_reco_Enu
    only = None
    if args.features == "uboone":
        from analysis_village.SinglePhotonWC import microboone_bdt
        only = {c for cols in microboone_bdt.input_columns(args.uboone_bdt).values() for c in cols}
    feats = feature_columns(p, only)
    X = p[feats].astype(np.float32)
    y = cuts.is_signal_cat(p.cat).astype(int).values
    fold = (p.index.get_level_values("file") % 2).values
    print("%d preselected candidates, %d signal; %d input variables" % (len(p), y.sum(), len(feats)))

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
    print("Preselection efficiency %.1f%%, purity %.2f%%" % (100 * n_signal_nu(p) / nsig_nu, 100 * y.mean()))
    wp = working_points(p, nsig_nu, pot)
    pd.set_option("display.width", 250)
    print("Working points (rates per 1e20 POT):")
    print(wp.round(3).to_string(index=False))
    print("Top inputs (% of gain):")
    print((100 * imp.head(20)).round(1).to_string())

    p[["cat", "bdt"]].to_pickle(os.path.join(args.outdir, "scores.pkl"))
    imp.to_csv(os.path.join(args.outdir, "importance.csv"))
    wp.to_csv(os.path.join(args.outdir, "working_points.csv"), index=False)
    plots(p, imp, nsig_nu, args.outdir)
    print("Wrote %s" % args.outdir)


if __name__ == "__main__":
    main()
