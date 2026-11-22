#!/usr/bin/env python3
"""Plots of the presel+uboone chain (cuts.py), per 1e20 POT:

  selection_eff_purity.png   efficiency and purity after each stage
  selection_composition.png  truth make-up of the candidates passing each stage
  selection_cutflow.csv      the numbers of both (cutflow.cutflow)
  final_shower_energy.png    shw_sp_energy after the full selection
  final_frac_uncertainty.png fractional uncertainty on the selected rate vs
                             shw_sp_energy: systematics by group (if the
                             inputs have wgt / slimwgt tables), MC statistics,
                             and the data statistics expected at --data-pot
  final_scores_nminus1.png   each score cut with every other cut applied
  selection_sequential.png   each cut variable with the earlier cuts
                             applied, signal vs background
  presel_shower_multiplicity.png
                             true EM particles vs pf_n_em_primary after
                             preselection, row-normalised

Stacked plots are by truth category with the uncertainty on the total as a
hatched band (MC statistics, plus systematics where available) and no
under/overflow bins. Efficiency and purity errors are binomial.

Input: configs/sp_mc.py outputs (wc + mcnu + wcevt + wgt; or slimwgt) and
the MicroBooNE BDT XMLs (the ub_* scores). --cache DIR keeps the loaded and
scored candidates (w_ub.pkl, m.pkl, pot.txt) for later runs; the weights
are always read from the inputs.

  python analysis_village/SinglePhotonWC/plot_selection.py -i <dir or .df files> \
      --uboone-bdt <xml dir> -o sp_selection_plots/ [--cache <dir>]
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
from analysis_village.SinglePhotonWC import cutflow  # noqa: E402

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


ENERGY_BINS = np.arange(0, cuts.SHOWER_E_MAX + 1, 100)   # shw_sp_energy after the full selection
UNC_BINS = np.r_[np.arange(0, 601, 100), 800, cuts.SHOWER_E_MAX]   # coarser where the sample is thin


def stack(ax, q, var, bins, sc, syst=None):
    """Stack by GROUPS; hatched band: MC statistics, in quadrature with the
    fractional systematic per bin syst if given."""
    x = var(q) if callable(var) else q[var]
    out = ((x < bins[0]) | (x >= bins[-1])).sum()
    if out:
        print("warning: %d candidates outside the range, not drawn" % out)
    vals = [x[q.cat.isin(cs)] for _, cs, _ in GROUPS]
    ax.hist(vals, bins=bins, stacked=True, color=[c for *_, c in GROUPS], weights=[np.full(len(v), sc) for v in vals],
            label=["%s (%.0f)" % (g, len(v) * sc) for (g, *_), v in zip(GROUPS, vals)], edgecolor="white", linewidth=0.6)
    n, _ = np.histogram(x, bins)
    err = np.sqrt(n + (0 if syst is None else (n * np.nan_to_num(syst)) ** 2))
    lo, hi = (n - err) * sc, (n + err) * sc
    ax.fill_between(bins, np.r_[lo, lo[-1]], np.r_[hi, hi[-1]], step="post", facecolor="none", edgecolor="#3d3c39",
                    hatch="////", linewidth=0, label="MC stat. uncertainty" if syst is None else "MC stat. + syst. uncertainty")
    ax.grid(axis="y", alpha=0.25, lw=0.5)


def plot_energy(w, chain, sc, out, syst=None):
    sel = w[cuts.passes(w, chain)]
    s = cuts.is_signal_cat(sel.cat)
    fig, ax = plt.subplots(figsize=(9, 5.6))
    stack(ax, sel, "shw_sp_energy", ENERGY_BINS, sc, syst)
    ax.set_xlabel("reco energy of the counted shower, shw_sp_energy [MeV]")
    ax.set_ylabel("candidates / 1e20 POT / 100 MeV")
    ax.set_ylim(0, 1.8 * ax.get_ylim()[1])   # room for the legend
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
        mid = 0.5 * (bins[1:] + bins[:-1])   # MC statistics
        a.errorbar(mid, hb, np.sqrt(hb * sc), fmt="none", ecolor="#52514e", elinewidth=1)
        a.errorbar(mid, hs * k, np.sqrt(hs * sc) * k, fmt="none", ecolor=SIG_BKG[0][1], elinewidth=1.2)
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
                         ["background", "signal (scaled up where the panel says × k)\nerror bars: MC statistics"],
                         frameon=False, fontsize=13, loc="center")
    fig.suptitle("Selection variables in cut order: each panel has every earlier cut applied (%.2e POT). "
                 "Candidates; wc_valid and vertex are not shown." % pot, fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(os.path.join(out, "selection_sequential.png"), dpi=130)
    plt.close(fig)

LAB = {   # stage labels of the eff / purity and composition plots
    "wc_valid": "WC output\n(true ν)", "all candidates": "all WC\ncandidates", "vertex": "vertex",
    "fv": "FV 10 cm\n= end of presel.", "wc_numu": "WC νμ\nscore < %g" % cuts.WC_NUMU_SCORE_MAX,
    "ub_numu": "uB νμ\nBDT > %g" % cuts.UB_NUMU_MIN, "ub_other": "uB other\nBDT > %g" % cuts.UB_OTHER_MIN,
    "ub_ncpi0": "uB NC π⁰\nBDT > %g" % cuts.UB_NCPI0_MIN, "narrow_shower": "shower angle\n< %g°" % cuts.SHOWER_A68_MAX,
    "shower_e_max": "shower E\n< %d MeV" % cuts.SHOWER_E_MAX, "wc_nue": "WC νe\nscore < %g" % cuts.WC_NUE_SCORE_MAX,
    "one_shower": "exactly 1 shower\n(none μ/π/p)", "no_second_em": "no 2nd EM\nshower > %g MeV" % cuts.SECOND_EM_MAX,
    "back_traced": "back-traced\n> %g cm" % cuts.BACK_DIST_MIN,
}
COMPOSITION = [   # label, truth categories, colour
    ("signal: 1 reconstructable γ (contained; NC π⁰ with the other γ < 20 MeV)", cuts.SIGNAL_CATEGORIES, "#2a78d6"),
    ("≥2 reconstructable γ in the TPC (truth; mostly NC π⁰)", [4], "#eb6834"),
    ("primary μ (CC νμ)", [6], "#eda100"),
    ("1γ + escaping γ (not contained)", [13], "#5c3317"),
    ("CC νe", [7], "#1baf7a"),
    ("outside AV, not 1γ", [8], "#e87ba4"),
    ("no reconstructable photon", [5], "#7f5bd5"),
    ("1γ from outside the AV (dirt)", [12], "#a05a2c"),
    ("1γ, not signal (vertex or conv. < 10 cm in)", [11], "#008300"),
    ("cosmic / unmatched", [9], "#a3a29d"),
]
INK = "#1f1f1e"


def display_chain(chain):
    """The chain for the per-stage plots: the μ/π/p-shower veto folded into
    "exactly one shower" (cuts are ANDed: only the per-stage numbers change)."""
    c = dict(chain)
    if "one_shower" in c and "em_showers" in c:
        one, em = c["one_shower"], c.pop("em_showers")
        c["one_shower"] = lambda df: one(df) & em(df)
    return c


def plot_eff_purity(t, pot, out):
    x = np.arange(len(t))
    fig, axs = plt.subplots(2, 1, figsize=(13, 7.5), sharex=True)
    for ax, col, ylab, top in [(axs[0], "eff", "signal efficiency [%]", 112), (axs[1], "purity", "purity [%]", 50)]:
        ax.errorbar(x, t[col], t[col + "_err"].fillna(0), fmt="o-", color="#2a78d6", lw=2, ms=6, capsize=3)
        for i, (v, e) in enumerate(zip(t[col], t[col + "_err"].fillna(0))):
            ax.annotate("%.1f" % v, (i, v + e), textcoords="offset points", xytext=(0, 5), ha="center", fontsize=8, color=INK)
        ax.set_ylabel(ylab)
        ax.set_ylim(0, top)
    axs[1].set_xticks(x, [LAB.get(n, n) for n in t.index], fontsize=8.5)
    axs[0].set_title("Selection efficiency and purity by stage (SBND MC, %.2g POT; signal: 1 contained photon, vertex and "
                     "conversion in the 10 cm FV); binomial MC-statistics errors" % pot, loc="left", color=INK, fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "selection_eff_purity.png"), dpi=140)
    plt.close(fig)


def plot_composition(t, pot, out):
    lab = cuts.CATEGORY_LABELS
    R = pd.DataFrame({g: t[[lab[c] for c in cs]].sum(1) for g, cs, _ in COMPOSITION})
    V = R.div(R.sum(1), axis=0) * 100
    fig, ax = plt.subplots(figsize=(11, 7))
    y = np.arange(len(R))[::-1]
    left = np.zeros(len(R))
    for g, _, c in COMPOSITION:
        ax.barh(y, V[g], left=left, color=c, height=0.66, edgecolor="#fcfcfb", linewidth=1.5, label=g)
        left += V[g].values
    ax.grid(axis="y", visible=False)
    ax.set_xlim(0, 100)
    ax.set_xlabel("share of candidates passing the stage [%]")
    ax.set_yticks(y, [LAB.get(n, n).replace("\n", " ") for n in R.index])
    for yi, (n, r) in zip(y, t.iterrows()):
        ax.text(101, yi, "%.1f ± %.1f%%" % (r.purity, r.purity_err) if r.purity_err > 0.05 else "%.1f%%" % r.purity,
                va="center", fontsize=8.5, color=INK)
    ax.text(101, y[0] + 0.7, "signal", fontsize=8.5, color=INK)
    h, l = ax.get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=9, loc="lower center", ncol=3)
    fig.suptitle("What survives each stage: composition of the passing candidates (%.2g POT)" % pot, x=0.01, ha="left",
                 color=INK, fontsize=12)
    fig.tight_layout(rect=(0, 0.11, 0.94, 1))
    fig.savefig(os.path.join(out, "selection_composition.png"), dpi=140)
    plt.close(fig)
    return R


UNC_STYLE = {   # component: colour, line style, width
    "total": ("#1f1f1e", "-", 3), "flux": ("#2a78d6", "-", 2), "genie": ("#eb6834", "-", 2), "g4": ("#1baf7a", "-", 2),
    "MC stat": ("#7f5bd5", "--", 2), "data stat": ("#a3a29d", ":", 2),
}


def frac_uncertainty(sel, wgt, pot, data_pot, bins=ENERGY_BINS):
    """Fractional uncertainty on the selected rate per shw_sp_energy bin and
    integrated ("all" row): systematic groups (wgt may be None), MC
    statistics, total (both in quadrature) and the data statistics at each
    of data_pot."""
    b = pd.cut(sel.shw_sp_energy, bins, right=False, labels=False)
    n = b.value_counts().reindex(range(len(bins) - 1), fill_value=0).astype(float)
    U = pd.DataFrame(index=list(range(len(bins) - 1)) + ["all"])
    U["n"] = pd.concat([n, pd.Series({"all": n.sum()})])
    if wgt is not None:
        ok = b.notna()
        bins_syst = cutflow.syst_unc(sel[ok], wgt, b[ok].astype(int))
        all_syst = cutflow.syst_unc(sel[ok], wgt, pd.Series("all", index=sel[ok].index))
        U = U.join(pd.concat([bins_syst, all_syst]).drop(columns="n"))
    U["MC stat"] = 1 / np.sqrt(U.n.replace(0, np.nan))
    U["total"] = np.sqrt(U.get("syst", 0) ** 2 + U["MC stat"] ** 2)
    for dp in data_pot:
        U["data stat %.0e POT" % dp] = 1 / np.sqrt(U.n.replace(0, np.nan) * dp / pot)
    return U


def plot_frac_uncertainty(w, chain, pot, out, wgt=None, data_pot=(1e20,), bins=UNC_BINS):
    sel = w[cuts.passes(w, chain)]
    U = frac_uncertainty(sel, wgt, pot, data_pot, bins)
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    comps = [c for c in ["total", "flux", "genie", "g4", "MC stat"] if c in U] + [c for c in U if c.startswith("data stat")]
    for i, c in enumerate(comps):
        col, ls, lw = UNC_STYLE.get(c, UNC_STYLE["data stat"])
        name = {"total": "total (syst. ⊕ MC stat.)" if "syst" in U else "total (MC stat. only)",
                "MC stat": "MC statistics (%.2g POT sample)" % pot}.get(c, cutflow.SYST_GROUPS.get(c, c))
        if c.startswith("data stat"):
            name, ls = "expected data statistics, %.0e POT" % float(c.split()[2]), [":", "-."][i % 2]
        v = U.loc[list(range(len(bins) - 1)), c].to_numpy(float)
        ax.stairs(100 * v, bins, color=col, ls=ls, lw=lw, label="%s: %.1f%% integrated" % (name, 100 * U.loc["all", c]))
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0, None)
    ax.set_xlabel("reco energy of the counted shower, shw_sp_energy [MeV]")
    ax.set_ylabel("fractional uncertainty [%]")
    ax.set_title("Uncertainty on the selected rate after the full selection (%d MC candidates)%s" % (len(sel), "" if "syst" in U else
                 ": no systematic weights in the inputs"), fontsize=10, loc="left")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "final_frac_uncertainty.png"), dpi=150)
    plt.close(fig)
    U.to_csv(os.path.join(out, "final_frac_uncertainty.csv"), float_format="%.4g")
    return U


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


def load_cached(inputs, uboone_bdt, cache=None, nproc=None):
    """cutflow.load, through the pickles in cache (written on the first run)."""
    if cache and os.path.exists(os.path.join(cache, "w_ub.pkl")):
        print("reading the candidates from %s" % cache)
        return (pd.read_pickle(os.path.join(cache, "w_ub.pkl")), pd.read_pickle(os.path.join(cache, "m.pkl")),
                float(open(os.path.join(cache, "pot.txt")).read()))
    w, m, pot = cutflow.load(inputs, uboone_bdt, nproc=nproc)
    if cache:
        os.makedirs(cache, exist_ok=True)
        w.to_pickle(os.path.join(cache, "w_ub.pkl"))
        m.to_pickle(os.path.join(cache, "m.pkl"))
        open(os.path.join(cache, "pot.txt"), "w").write(repr(pot))
    return w, m, pot


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--inputs", nargs="+", required=True, help=".df files, directories or globs")
    ap.add_argument("--uboone-bdt", required=True, help="directory of the MicroBooNE BDT XMLs")
    ap.add_argument("-o", "--outdir", default="sp_selection_plots")
    ap.add_argument("--cache", help="directory for the loaded candidates (read if present, else written)")
    ap.add_argument("--data-pot", type=float, nargs="*", default=[1e20, 1e21],
                    help="exposures for the expected data statistics in final_frac_uncertainty.png")
    ap.add_argument("-j", "--nproc", type=int, help="worker processes (default: up to 8)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    w, m, pot = load_cached(args.inputs, args.uboone_bdt, args.cache, args.nproc)
    chain, sc = cuts.CHAINS["presel+uboone"], 1e20 / pot
    sel = w[cuts.passes(w, chain)]
    wgt = cutflow.load_weights(args.inputs, keep=pd.MultiIndex.from_frame(sel.reset_index()[cutflow.NU_KEY]),
                               nproc=args.nproc)
    if wgt is None:
        print("no wgt / slimwgt table in the inputs: MC statistics only")
    t = cutflow.cutflow(w, m, pot, display_chain(chain))
    t.to_csv(os.path.join(args.outdir, "selection_cutflow.csv"), float_format="%.4g")
    plot_eff_purity(t, pot, args.outdir)
    plot_composition(t, pot, args.outdir)
    plot_frac_uncertainty(w, chain, pot, args.outdir, wgt, args.data_pot)
    U = frac_uncertainty(sel, wgt, pot, [], ENERGY_BINS)
    plot_energy(w, chain, sc, args.outdir, U.loc[list(range(len(ENERGY_BINS) - 1)), "syst"].to_numpy(float)
                if "syst" in U else None)
    plot_scores(w, chain, sc, pot, args.outdir)
    plot_sequential(w, chain, sc, pot, args.outdir)
    if "pf_n_em_primary" in w and w.pf_n_em_primary.notna().any():
        plot_multiplicity(w, pot, args.outdir)
    else:
        print("no pf_n_em_primary (dataframes made before T_pf): skipping presel_shower_multiplicity.png")
    pd.set_option("display.width", 250)
    print(t[["eff", "eff_err", "purity", "purity_err", "total"]].round(2).T.to_string())
    print(pd.read_csv(os.path.join(args.outdir, "final_frac_uncertainty.csv"), index_col=0).round(3).to_string())
    print("wrote the plots to %s (%.3g POT)" % (args.outdir, pot))


if __name__ == "__main__":
    main()
