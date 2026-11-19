#!/usr/bin/env python3
"""Plots of the presel+uboone chain (cuts.py), per 1e20 POT:

  final_shower_energy.png    shw_sp_energy after the full selection
  final_scores_nminus1.png   each score cut with every other cut applied
  selection_sequential.png   each cut variable with the earlier cuts
                             applied, signal vs background
  presel_shower_multiplicity.png
                             true EM particles vs pf_n_em_primary after
                             preselection, row-normalised

The first two are stacked by truth category with the MC statistics on the
total as a hatched band and no under/overflow bins.

Input: configs/sp_mc.py outputs (wc + mcnu + wcevt) and the MicroBooNE BDT
XMLs (the ub_* scores).

  python analysis_village/SinglePhotonWC/plot_selection.py -i <dir or .df files> \
      --uboone-bdt <xml dir> -o sp_selection_plots/
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from analysis_village.SinglePhotonWC import cuts  # noqa: E402
from analysis_village.SinglePhotonWC.cutflow import load  # noqa: E402

plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), "dune.mplstyle"))

GROUPS = [   # label, truth categories, colour
    ("Signal: NC $\\pi^0$ 1$\\gamma$", [1], "#2a78d6"),
    ("Signal: NC $\\Delta\\to N\\gamma$", [0], "#e87ba4"),
    ("Signal: other 1$\\gamma$", [3], "#6250d6"),
    ("$\\geq$2 reconstructable photons (truth)", [4], "#eb6834"),
    ("1$\\gamma$ + escaping $\\gamma$ (not contained)", [13], "#a05a2c"),
    ("CC $\\nu_e$", [7], "#1baf7a"),
    ("Primary $\\mu$ (CC $\\nu_\\mu$)", [6], "#eda100"),
    ("Cosmic / unmatched", [9], "#e34948"),
    ("Other bkg (0 showers, outside AV, dirt, 1$\\gamma$ outside FV)", [5, 8, 11, 12], "#008300"),
]

SCORES = [   # cut name in cuts.UB_SELECTION, column, label, bins, threshold, kept side
    ("wc_numu", "numu_score", "Wire-Cell $\\nu_\\mu$ CC score", np.linspace(-3, 3, 31), cuts.WC_NUMU_SCORE_MAX, "<"),
    ("ub_numu", "ub_numu", "uB $\\nu_\\mu$ BDT", np.arange(-3.6, 4.01, 0.2), cuts.UB_NUMU_MIN, ">"),
    ("ub_other", "ub_other", "uB 'other' BDT", np.arange(-4.4, 3.01, 0.2), cuts.UB_OTHER_MIN, ">"),
    ("ub_ncpi0", "ub_ncpi0", "uB NC $\\pi^0$ BDT", np.arange(-3.0, 3.51, 0.25), cuts.UB_NCPI0_MIN, ">"),
    ("narrow_shower", "shw1_a68", "counted shower 68% opening angle [deg]", np.arange(0, 40.1, 1.5), cuts.SHOWER_A68_MAX, "<"),
    ("shower_e_max", "shw_sp_energy", "counted shower reco energy [MeV]", np.arange(0, 2501, 100), cuts.SHOWER_E_MAX, "<"),
    ("wc_nue", "nue_score", "Wire-Cell $\\nu_e$ score", np.arange(-15.5, 17.6, 1.0), cuts.WC_NUE_SCORE_MAX, "<"),
]


def keep_marker(ax, thr, sense, **kw):
    """Dashed line at the threshold, with an arrow to the kept side."""
    ax.axvline(thr, color="#222", lw=1.6, ls="--")
    ax.text(thr, 0.95, "  keep $\\rightarrow$" if sense == ">" else "$\\leftarrow$ keep  ", transform=ax.get_xaxis_transform(),
            ha="left" if sense == ">" else "right", va="center", fontsize=9.5, **kw)


def stack(ax, q, var, bins, sc):
    x = var(q) if callable(var) else q[var]
    out = ((x < bins[0]) | (x >= bins[-1])).sum()
    if out:
        print("warning: %d candidates outside the range, not drawn" % out)
    vals = [x[q.cat.isin(cs)] for _, cs, _ in GROUPS]
    ax.hist(vals, bins=bins, stacked=True, color=[c for *_, c in GROUPS], weights=[np.full(len(v), sc) for v in vals],
            label=["%s (%.0f)" % (g, len(v) * sc) for (g, *_), v in zip(GROUPS, vals)], edgecolor="white", linewidth=0.6)
    n, _ = np.histogram(x, bins)
    lo, hi = (n - np.sqrt(n)) * sc, (n + np.sqrt(n)) * sc     # MC statistics on the total, shaded
    ax.fill_between(bins, np.r_[lo, lo[-1]], np.r_[hi, hi[-1]], step="post", facecolor="none", edgecolor="#3d3c39",
                    hatch="////", linewidth=0, label="MC stat. uncertainty")
    ax.grid(axis="y", alpha=0.25, lw=0.5)


def plot_energy(w, chain, sc, out):
    sel = w[cuts.passes(w, chain)]
    s = cuts.is_signal_cat(sel.cat)
    fig, ax = plt.subplots(figsize=(9, 5.6))
    stack(ax, sel, "shw_sp_energy", np.arange(0, 2401, 100), sc)
    ax.set_xlabel("reco energy of the counted shower, shw_sp_energy [MeV]")
    ax.set_ylabel("candidates / 1e20 POT / 100 MeV")
    ax.set_title("After full selection: %.0f signal + %.0f background / 1e20 POT (%d MC candidates); purity %.0f%%"
                 % (s.sum() * sc, (~s).sum() * sc, len(sel), 100 * s.mean()), fontsize=10, loc="left")
    ax.legend(frameon=False, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "final_shower_energy.png"), dpi=150)
    plt.close(fig)


def plot_scores(w, chain, sc, pot, out):
    fig, axs = plt.subplots(2, 3, figsize=(18, 10))
    axs = axs.ravel()
    for a, (cut, var, lab, bins, thr, sense) in zip(axs, [x for x in SCORES if x[0] in chain]):
        q = w[cuts.passes(w, {k: v for k, v in chain.items() if k != cut})]
        stack(a, q, var, bins, sc)
        keep_marker(a, thr, sense, color="#222")
        sig, keep = cuts.is_signal_cat(q.cat), chain[cut](q)
        a.set_title("%s  (cut %s %g)\nall other cuts applied: this cut keeps %.0f%% of signal, %.0f%% of background"
                    % (lab, sense, thr, 100 * keep[sig].mean(), 100 * keep[~sig].mean()), fontsize=9.5, loc="left")
        a.set_xlabel(var if isinstance(var, str) else lab)
        a.set_ylabel("candidates / 1e20 POT")
    h, l = axs[0].get_legend_handles_labels()
    axs[-1].axis("off")
    axs[-1].legend(h, [x.rsplit(" (", 1)[0] for x in l], frameon=False, fontsize=10, loc="center")
    fig.suptitle("Score distributions with every other selection cut applied (N$-$1), %.2e POT; hatched band: MC statistics" % pot,
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(os.path.join(out, "final_scores_nminus1.png"), dpi=140)
    plt.close(fig)


def _fv_wall(df):
    return cuts.dist_to_active_boundary(df.nu_x, df.nu_y, df.nu_z)


def _fv_cathode(df):
    return df.nu_x.abs()


# cut name, column or function, label, bins, threshold, kept side. fv is two
# panels: the distance to the AV walls, then (after it) |x| to the cathode.
SEQUENTIAL = [
    ("fv", _fv_wall, "reco vertex distance to the AV wall [cm]", np.linspace(-20, 200, 45), cuts.FV_INSET, ">"),
    ("fv", _fv_cathode, "reco vertex |x|, distance to the cathode [cm]", np.linspace(0, 200, 41), cuts.CATHODE_INSET, ">"),
] + SCORES + [
    ("one_shower", "shw_sp_n_20br1_showers", "SP-tagger showers (shw_sp_n_20br1_showers)", np.arange(-0.5, 6.5), None, "== 1"),
    ("em_showers", lambda df: df.n_shower - df.n_em, "showers labelled $\\mu/\\pi/p$ (n_shower $-$ n_em)", np.arange(-0.5, 5.5), 0.5, "<"),
    ("no_second_em", lambda df: (df.E_em_all_2 if "E_em_all_2" in df else df.E_em_2).fillna(0), "2nd EM shower, anywhere [MeV]",
     np.arange(0, 201, 10), cuts.SECOND_EM_MAX, "<"),
    ("back_traced", "shw1_back_dist", "counted shower traced back to the AV wall [cm]", np.linspace(0, 500, 26), cuts.BACK_DIST_MIN, ">"),
]
SIG_BKG = [("signal", "#2a78d6"), ("background", "#a3a29d")]


def signal_scale(smax, bmax):
    """1, 2, 5, 10, ... so the signal peak reaches about half the background
    peak on a linear axis."""
    if smax <= 0 or smax >= 0.3 * bmax:
        return 1
    k = 0.5 * bmax / smax
    e = 10 ** np.floor(np.log10(k))
    return max(1, max(m * e for m in (1, 2, 5) if m * e <= k))


def plot_sequential(w, chain, sc, pot, out):
    """Each cut's variable for the candidates that passed every earlier cut,
    signal and background stacked; the title gives what the cut removes."""
    names = list(chain)
    seq = [x for x in SEQUENTIAL if x[0] in chain]
    nrow = -(-(len(seq) + 1) // 4)   # one extra panel for the legend
    fig, axs = plt.subplots(nrow, 4, figsize=(30, 5.7 * nrow))
    axs = axs.ravel()
    for a, (cut, var, lab, bins, thr, sense) in zip(axs, seq):
        before = names[names.index(cut) - 1]
        q = w[cuts.passes(w, chain, upto=before)]
        if var is _fv_cathode:
            q = q[_fv_wall(q) >= cuts.FV_INSET]
        x = (var(q) if callable(var) else q[var])
        x = np.clip(np.asarray(x, float), bins[0], bins[-1] - 1e-6)
        sig = cuts.is_signal_cat(q.cat).to_numpy()
        hb, _, _ = a.hist(x[~sig], bins=bins, color=SIG_BKG[1][1], edgecolor="white", linewidth=0.6,
                          weights=np.full((~sig).sum(), sc), label="background")
        hs, _ = np.histogram(x[sig], bins, weights=np.full(sig.sum(), sc))
        k = signal_scale(hs.max(), hb.max())
        a.stairs(hs * k, bins, color=SIG_BKG[0][1], lw=2.5, label="signal" + (" × %g" % k if k > 1 else ""))
        a.legend(fontsize=12, loc="upper right")
        a.grid(axis="y", alpha=0.25, lw=0.5)
        if thr is None:   # one_shower: keep the == 1 bin
            keep = x == 1
            a.axvspan(0.5, 1.5, color="#2a78d6", alpha=0.08)
            a.text(1, 0.95, "keep", transform=a.get_xaxis_transform(), ha="center", va="center", fontsize=9.5)
        else:
            keep = x > thr if sense == ">" else x < thr
            keep_marker(a, thr, sense)
        part = {_fv_wall: " (wall part)", _fv_cathode: " (cathode part)"}.get(var, "") if callable(var) else ""
        a.set_title("%s%s, after %s: %.0f sig + %.0f bkg / 1e20 POT\nremoves %.0f%% of signal, %.0f%% of background reaching it"
                    % (cut, part, "fv wall part" if part == " (cathode part)" else before, sig.sum() * sc, (~sig).sum() * sc,
                       100 * (~keep[sig]).mean(), 100 * (~keep[~sig]).mean()), fontsize=10, loc="left")
        a.set_xlabel(lab)
        a.set_ylabel("candidates / 1e20 POT")
    for a in axs[len(seq):]:
        a.axis("off")
    from matplotlib.patches import Patch
    axs[len(seq)].legend([Patch(color=SIG_BKG[1][1]), plt.Line2D([], [], color=SIG_BKG[0][1], lw=2.5)],
                                ["background", "signal (scaled up where the panel says × k)"], frameon=False, fontsize=13, loc="center")
    fig.suptitle("Selection variables in cut order: each panel has every earlier cut applied (%.2e POT). "
                 "Candidates; wc_valid and vertex are not shown." % pot, fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(os.path.join(out, "selection_sequential.png"), dpi=130)
    plt.close(fig)


MULT_TEXT = (
    "RECO  pf_n_em_primary\n"
    "Wire-Cell particle-flow tree (T_pf):\n"
    "e/$\\gamma$ showers (PDG 11) with no real\n"
    "particle above them in the tree:\n"
    "  - attached at the neutrino vertex, or\n"
    "  - displaced, linked to the vertex through a\n"
    "    pseudo photon / $\\pi^0$\n"
    "Showers below a $\\mu$, $\\pi$, p or a pseudo\n"
    "neutron are secondaries and not counted.\n"
    "No energy threshold.\n\n"
    "TRUTH  per matched true neutrino\n"
    "  - photons from $\\pi^0$/$\\eta$/hadron decays,\n"
    "    any energy\n"
    "  - GENIE-primary photons $\\geq$ 20 MeV\n"
    "  - primary e$^\\pm$\n"
    "  no containment, conversion or merging\n\n"
    "Candidates: after preselection (wc_valid,\n"
    "vertex, 10 cm FV, 5 cm cathode), matched to a\n"
    "true neutrino (%d; cosmics excluded).")


def plot_multiplicity(w, pot, out, var="pf_n_em_primary", L=5):
    p = w[cuts.passes(w, cuts.PRESELECTION)]
    nu = p[p.true_E.notna() & (p.tmatch_idx >= 0) & p[var].notna()]
    true_em = nu.true_ngamma_20MeV + nu.true_ndecayg_soft + nu.true_ne
    M = pd.crosstab(true_em.clip(upper=L), nu[var].clip(upper=L)).reindex(index=range(L + 1), columns=range(L + 1), fill_value=0)
    R = M.div(M.sum(1).replace(0, np.nan), axis=0)
    tk = [str(i) for i in range(L)] + ["$\\geq$%d" % L]
    fig = plt.figure(figsize=(14.5, 7.4))
    a = fig.add_axes([0.03, 0.11, 0.47, 0.74], anchor="E")
    cax = fig.add_axes([0.515, 0.11, 0.012, 0.74])
    im = a.imshow(R.values, origin="lower", cmap="Blues", vmin=0, vmax=1)
    for i in range(L + 1):
        for k in range(L + 1):
            if M.values[i, k]:
                a.text(k, i, "%.0f%%\n%d" % (100 * R.values[i, k], M.values[i, k]), ha="center", va="center", fontsize=8.5,
                       color="white" if R.values[i, k] > 0.55 else "#222")
    a.set_xticks(range(L + 1)); a.set_xticklabels(tk); a.set_yticks(range(L + 1))
    a.set_yticklabels(["%s  (%d)" % (t, n) for t, n in zip(tk, M.sum(1))])
    a.set_xlabel("RECO: %s" % var); a.set_ylabel("TRUE: prompt photons + primary e$^\\pm$  (candidates)")
    d0 = M.values[0, 0] / M.values[0].sum(); d1 = np.trace(M.values[1:L, 1:L]) / M.values[1:].sum()
    a.set_title("reco = truth: %.0f%% of candidates with true $\\geq$1, %.0f%% with true 0" % (100 * d1, 100 * d0), fontsize=10.5)
    fig.colorbar(im, cax=cax, label="fraction of the true row")
    fig.text(0.6, 0.85, MULT_TEXT % len(nu), va="top", fontsize=10)
    fig.suptitle("Wire-Cell primary e/$\\gamma$ shower count vs truth after preselection (%.2e POT)\n"
                 "row-normalised %%, MC candidates below" % pot, fontsize=11.5)
    fig.savefig(os.path.join(out, "presel_shower_multiplicity.png"), dpi=140)
    plt.close(fig)
    return R


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
    plot_sequential(w, chain, sc, pot, args.outdir)
    if "pf_n_em_primary" in w and w.pf_n_em_primary.notna().any():
        plot_multiplicity(w, pot, args.outdir)
    else:
        print("no pf_n_em_primary (dataframes made before T_pf): skipping presel_shower_multiplicity.png")
    print("wrote the plots to %s (%.3g POT)" % (args.outdir, pot))


if __name__ == "__main__":
    main()
