#!/usr/bin/env python3
"""Post-selection plots of the presel+uboone chain (cuts.py), stacked by
truth category, per 1e20 POT with MC-statistics errors on the total:

  final_shower_energy.png    reco energy of the counted shower
                             (shw_sp_energy) after the full selection
  final_scores_nminus1.png   the score cuts (WC numu/nue, uB numu/other/
                             ncpi0 BDTs), each with every other cut applied

Input: configs/sp_mc.py outputs (wc + mcnu + wcevt) and the MicroBooNE BDT
XMLs (the ub_* scores).

  python analysis_village/SinglePhotonWC/plot_selection.py -i <dir or .df files> \
      --uboone-bdt <xml dir> -o sp_selection_plots/
"""

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from analysis_village.SinglePhotonWC import cuts  # noqa: E402
from analysis_village.SinglePhotonWC.cutflow import load  # noqa: E402

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

# Stacked groups: label, truth categories, colour
GROUPS = [
    ("Signal: NC $\\pi^0$ 1$\\gamma$", [1], "#2a78d6"),
    ("Signal: NC $\\Delta\\to N\\gamma$", [0], "#e87ba4"),
    ("Signal: other ($\\nu_\\mu$CC, NC other, out-of-FV)", [2, 3, 10], "#6250d6"),
    ("$\\geq$2 visible photons", [4], "#eb6834"),
    ("CC $\\nu_e$", [7], "#1baf7a"),
    ("CC $\\nu_\\mu$ ($\\mu\\geq$100 MeV)", [6], "#eda100"),
    ("Cosmic / unmatched", [9], "#e34948"),
    ("Other bkg (0 showers, outside AV, dirt, soft 1$\\gamma$)", [5, 8, 11, 12], "#008300"),
]

# Score cuts: cut name in cuts.UB_SELECTION, column, label, bins, threshold, kept side
SCORES = [
    ("wc_numu", "numu_score", "Wire-Cell $\\nu_\\mu$ CC score", np.linspace(-3, 3, 31), cuts.WC_NUMU_SCORE_MAX, "<"),
    ("ub_numu", "ub_numu", "uB $\\nu_\\mu$ BDT", np.linspace(-2, 4, 31), cuts.UB_NUMU_MIN, ">"),
    ("ub_other", "ub_other", "uB 'other' BDT", np.linspace(-3, 3, 31), cuts.UB_OTHER_MIN, ">"),
    ("ub_ncpi0", "ub_ncpi0", "uB NC $\\pi^0$ BDT", np.linspace(-2, 3.5, 23), cuts.UB_NCPI0_MIN, ">"),
    ("wc_nue", "nue_score", "Wire-Cell $\\nu_e$ score", np.linspace(-15, 15, 31), cuts.WC_NUE_SCORE_MAX, "<"),
]


def stack(ax, q, var, bins, sc, logy=False):
    x = q[var].clip(bins[0], bins[-1] - 1e-6)   # under/overflow into the end bins
    vals = [x[q.cat.isin(cs)] for _, cs, _ in GROUPS]
    ax.hist(vals, bins=bins, stacked=True, color=[c for *_, c in GROUPS], weights=[np.full(len(v), sc) for v in vals],
            label=["%s (%.0f)" % (g, len(v) * sc) for (g, *_), v in zip(GROUPS, vals)], edgecolor="white", linewidth=0.6)
    n, _ = np.histogram(x, bins)
    ax.errorbar(0.5 * (bins[1:] + bins[:-1]), n * sc, np.sqrt(n) * sc, fmt="none", ecolor="#444", lw=1)
    if logy:
        ax.set_yscale("log")
    ax.grid(axis="y", alpha=0.25, lw=0.5)


def plot_energy(w, chain, sc, out):
    sel = w[cuts.passes(w, chain)]
    s = cuts.is_signal_cat(sel.cat)
    fig, ax = plt.subplots(figsize=(9, 5.6))
    stack(ax, sel, "shw_sp_energy", np.arange(0, 1201, 60), sc)
    ax.set_xlabel("reco energy of the counted shower, shw_sp_energy [MeV] (overflow in last bin)")
    ax.set_ylabel("candidates / 1e20 POT / 60 MeV")
    ax.set_title("After full selection: %.0f signal + %.0f background / 1e20 POT (%d MC candidates); purity %.0f%%"
                 % (s.sum() * sc, (~s).sum() * sc, len(sel), 100 * s.mean()), fontsize=10, loc="left")
    ax.legend(frameon=False, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "final_shower_energy.png"), dpi=150)
    plt.close(fig)


def plot_scores(w, chain, sc, pot, out):
    fig, axs = plt.subplots(2, 3, figsize=(18, 10))
    axs = axs.ravel()
    for a, (cut, var, lab, bins, thr, sense) in zip(axs, SCORES):
        q = w[cuts.passes(w, {k: v for k, v in chain.items() if k != cut})]
        stack(a, q, var, bins, sc, logy=True)
        a.axvline(thr, color="#222", lw=1.6, ls="--")
        a.text(thr, 0.95, "  keep $\\rightarrow$" if sense == ">" else "$\\leftarrow$ keep  ", transform=a.get_xaxis_transform(),
               ha="left" if sense == ">" else "right", va="center", fontsize=9.5, color="#222")
        sig, keep = cuts.is_signal_cat(q.cat), chain[cut](q)
        a.set_title("%s  (cut %s %g)\nall other cuts applied: this cut keeps %.0f%% of signal, %.0f%% of background"
                    % (lab, sense, thr, 100 * keep[sig].mean(), 100 * keep[~sig].mean()), fontsize=9.5, loc="left")
        a.set_xlabel(var)
        a.set_ylabel("candidates / 1e20 POT")
        a.set_ylim(bottom=3)
    h, l = axs[0].get_legend_handles_labels()
    axs[-1].axis("off")
    axs[-1].legend(h, [x.rsplit(" (", 1)[0] for x in l], frameon=False, fontsize=10, loc="center")
    fig.suptitle("Score distributions with every other selection cut applied (N$-$1), %.2e POT; error bars: MC statistics" % pot,
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(os.path.join(out, "final_scores_nminus1.png"), dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--inputs", nargs="+", required=True, help=".df files, directories or globs")
    ap.add_argument("--uboone-bdt", required=True, help="directory of the MicroBooNE BDT XMLs")
    ap.add_argument("-o", "--outdir", default="sp_selection_plots")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    w, m, pot = load(args.inputs, args.uboone_bdt)
    chain, sc = cuts.CHAINS["presel+uboone"], 1e20 / pot
    plot_energy(w, chain, sc, args.outdir)
    plot_scores(w, chain, sc, pot, args.outdir)
    print("wrote final_shower_energy.png, final_scores_nminus1.png to %s (%.3g POT)" % (args.outdir, pot))


if __name__ == "__main__":
    main()
