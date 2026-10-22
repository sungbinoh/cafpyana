#!/usr/bin/env python3
"""Truth-composition plots for the NCDelta dataframes, stacked by the truth
categories in cuts.py. Per channel (1g1p, 1g0p):

  composition_trueE_<ch>.png  true neutrino energy
  composition_recoE_<ch>.png  reconstructed visible energy, E_gamma + T_p
  invmass_1g1p.png            proton + photon invariant mass (1g1p only)
  cutvars_<ch>.png            N-1 panel per cut: every other cut applied,
                              threshold drawn, rejected side shaded

Per-category yields are in each legend, totals and purity in the titles.

Visible energy, because the outgoing neutrino of an NC event is not measured.
Cosmic slices have no true neutrino energy, so in the true-energy plots they
appear only as a legend yield.

The energy and mass plots show two panels: before the cut-based selection
(topology + pre-selection) and after it.

Inputs (-i) are directories, .df files, glob patterns or .list files;
unreadable files are skipped with a warning and their POT is not counted.
"""

import argparse
import math
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Sibling modules are imported from this file's directory, so the script runs
# the same inside cafpyana and in a standalone copy of these files.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import cuts                            # noqa: E402
from shower_energy import visible_energy   # noqa: E402
from dfio import expand_inputs, load   # noqa: E402

DEFAULT_INPUT = "/pnfs/sbnd/scratch/users/nrowe/cafpyana_out/dfs/2026_10_17_164529__ncdelta.df/"

DELTA_MASS = 1.232  # GeV, PDG Delta(1232)
COSMIC = 7
SIGNAL = 0
COLORS = {
    0: "#d62728",  # NC Delta->Ngamma
    1: "#9467bd",  # NC 1pi0 coherent
    2: "#1f77b4",  # NC 1pi0 non-coherent
    3: "#17becf",  # CC numu 1pi0
    4: "#2ca02c",  # CC nue/nuebar
    5: "#7f7f7f",  # BNB other
    6: "#8c564b",  # Dirt
    7: "#bcbd22",  # Cosmic
}
STAGE_TITLES = {
    "preselection": "{topo} before selection",
    "selection": "{topo} after selection",
}
STAGE_ORDER = ("preselection", "selection")


# ---------------------------------------------------------------------------
# Derived quantities
# ---------------------------------------------------------------------------

# --- Cut variables: what to plot for each stage name in cuts.stage_cuts()
# One row per cut: stage name, the cuts.py config key, the quantity, the axis
# range, the label, and which side passes. A float range means (0, that x the
# threshold), for scores whose scale follows the cut. Thresholds come from
# cuts.py, so the plots follow it. Cuts absent from a config are not drawn.
_CUT_PLOTS = (
    ("nu_score",          "nu_score_min",        lambda d: d.nu_score,
     (0, 1),      "Pandora neutrino score",                     "above"),
    ("flash_match",       "bary_score_min",      lambda d: d.bary_score,
     (0, 1),      "Barycenter flash-match score",               "above"),
    ("opt0_score",        "opt0_score_min",      lambda d: d.opt0_score,
     6.,          "OpT0 flash-match score",                     "above"),
    ("shw_dedx",          "shw_dedx_min",        lambda d: d.shw_dedx_best,
     (0, 10),     "Shower dE/dx, best plane [MeV/cm]",          "above"),
    ("shw_dedx_max",      "shw_dedx_max",        lambda d: d.shw_dedx_best,
     (0, 20),     "Shower dE/dx, best plane [MeV/cm]",          "below"),
    ("shw_conv_dist",     "shw_conv_dist_min",   lambda d: d.shw_conv_dist,
     (0, 100),    "Shower conversion distance [cm]",            "above"),
    ("shw_conv_dist_max", "shw_conv_dist_max",   lambda d: d.shw_conv_dist,
     4.,          "Shower conversion distance [cm]",            "below"),
    ("shw_len",           "shw_len_min",         lambda d: d.shw_len,
     (0, 200),    "Shower length [cm]",                         "above"),
    ("shw_trackScore",    "shw_trackScore_max",  lambda d: d.shw_trackScore,
     (0, 0.5),    "Shower track/shower score",                  "below"),
    ("shw_open_angle",    "shw_open_angle_max",  lambda d: d.shw_open_angle,
     (0, 1),      "Shower opening angle [rad]",                 "below"),
    ("shw_nhits",         "shw_nhits_min",
     lambda d: d.shw_nhits_I0 + d.shw_nhits_I1 + d.shw_nhits_I2,
     (0, 1000),   "Shower hits, all planes",                    "above"),
    ("extrap_det_path",   "extrap_det_path_min", lambda d: cuts.extrap_det_path(d),
     (0, 600),    "Extrapolated detector path length [cm]",     "above"),
    ("shw_contained",     "shw_start_dist_min",  lambda d: d.shw_start_dist_to_wall,
     (-10, 100),  "Shower start distance to wall [cm]",         "above"),
    ("vtx_contained",     "vtx_dist_min",        lambda d: d.vtx_dist_to_wall,
     (-10, 100),  "Vertex distance to wall [cm]",               "above"),
    ("trk_contained",     "trk_dist_min",
     lambda d: np.minimum(d.trk_start_dist_to_wall, d.trk_end_dist_to_wall),
     (-10, 100),  "Track min distance to wall [cm]",            "above"),
    ("trk_ke",            "trk_ke_max",          lambda d: d.trk_ke_p,
     (0, 1.0),    "Proton KE from range [GeV]",                 "below"),
    ("trk_calo_E",        "trk_calo_E_max",      lambda d: d.trk_calo_E,
     (0, 1.0),    "Track calorimetric energy [GeV]",            "below"),
    ("trk_attached",      "trk_start_dist_to_vtx_max",
     lambda d: d.trk_start_dist_to_vtx,
     (0, 30),     "Track start distance to vertex [cm]",        "below"),
    ("trk_pida",          "trk_pida_max",        lambda d: d.trk_pida,
     (0, 40),     r"Track PIDA [MeV/cm$^{1.42}$]",              "below"),
)


def cut_var_specs(topo, C=None):
    """{stage name: plot spec} for the cuts enabled in config `C`."""
    C = cuts.CUTBASED[topo] if C is None else C
    s = {}
    for stage, key, var, rng, label, keep in _CUT_PLOTS:
        thr = C.get(key)
        if thr is None:
            continue
        s[stage] = dict(var=var, label=label, cuts=[thr], keep=keep,
                        range=(0, rng * thr) if isinstance(rng, float) else rng)

    # Cuts whose panel needs more than one threshold or a custom axis.
    if C.get("abs_cos_trk_shw_max") is not None:
        a = C["abs_cos_trk_shw_max"]
        s["trk_shw_angle"] = dict(var=lambda d: d.cos_trk_shw_implied, range=(-1, 1),
                                  label=r"cos$\theta$(track, vertex$\to$shower)",
                                  cuts=[-a, a], keep="between")
    if C.get("n_nonprim_shw_max") is not None:
        n = C["n_nonprim_shw_max"]
        s["no_second_shw"] = dict(var=lambda d: d.n_nonprim_shw, range=(-0.5, 5.5),
                                  nbins=6, label="Non-primary showers in slice",
                                  cuts=[n + 0.5], keep="below",
                                  cut_labels=[r"$\leq$ %d" % n])
    if C.get("delta_mass_range") is not None:
        s["delta_mass"] = dict(var=lambda d: d.delta_mass, range=(0.9, 2.0),
                               label="Proton + photon invariant mass [GeV]",
                               cuts=list(C["delta_mass_range"]), keep="between")
    if C.get("trk_proton_like"):
        s["trk_proton_like"] = dict(var=lambda d: d.trk_chi2_proton - d.trk_chi2_muon,
                                    range=(-100, 100), label=r"$\chi^2_p - \chi^2_\mu$",
                                    cuts=[0.], keep="below")

    lo_E, hi_E = cuts.energy_window(topo, C)
    s["energy_window"] = dict(var=lambda d: cuts.shw_energy_corrected(d), range=(0, 1.2),
                              label="Shower energy, scaled best plane (signal window) [GeV]",
                              cuts=[lo_E, hi_E], keep="between",
                              cut_labels=[r"$\geq$ %g GeV" % lo_E, "< %g GeV" % hi_E])
    return s


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def stack(ax, values, cat, weight, bins, include_cosmic=True, legend=True,
          legend_fontsize=7):
    """Stacked histogram of `values` by truth category. Returns total yield.

    Values outside the bins are put in the first/last bin; NaNs are dropped.
    """
    lo, hi = bins[0], bins[-1]
    eps = 1e-9 * (hi - lo)
    v = values.clip(lower=lo, upper=hi - eps)
    finite = v.notna().values

    hists, colors, labels = [], [], []
    for code, label in cuts.CATEGORY_LABELS.items():
        m = (cat == code).values & finite
        n = weight[m].sum()
        if n == 0 or (code == COSMIC and not include_cosmic):
            continue
        hists.append((v[m], weight[m]))
        colors.append(COLORS[code])
        labels.append("%s (%.1f)" % (label, n))

    if hists:
        ax.hist([h[0] for h in hists], bins=bins, weights=[h[1] for h in hists],
                stacked=True, histtype="stepfilled", color=colors, label=labels,
                edgecolor="black", linewidth=0.4)

    if not include_cosmic:
        n_cos = weight[(cat == COSMIC).values].sum()
        if n_cos > 0:
            ax.bar(np.nan, np.nan, color=COLORS[COSMIC], label="Cosmic (%.1f, not shown)" % n_cos)

    ax.set_xlim(lo, hi)
    if legend:
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=legend_fontsize, loc="upper right", frameon=False)
        else:
            ax.text(0.5, 0.5, "no slices", transform=ax.transAxes, ha="center", va="center")
    return weight[finite].sum()


def stage_title(name, topo, cat, w):
    total = w.sum()
    purity = w[(cat == SIGNAL).values].sum() / total if total > 0 else 0.
    return "%s\ntotal %.1f, purity %.1f%%" % (STAGE_TITLES[name].format(topo=topo), total, 100. * purity)


def stage_figure(stages, topo, values_fn, bins, xlabel, ylabel, out,
                 include_cosmic=True, logy=False, vline=None):
    names = [n for n in STAGE_ORDER if n in stages]
    fig, axes = plt.subplots(1, len(names), figsize=(5.4 * len(names), 5), squeeze=False)
    for ax, name in zip(axes[0], names):
        d, w, cat = stages[name]
        stack(ax, values_fn(d), cat, w, bins, include_cosmic=include_cosmic)
        ax.set_title(stage_title(name, topo, cat, w), fontsize=10)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if logy:
            ax.set_yscale("log")
        if vline is not None:
            x, text = vline
            ax.axvline(x, color="black", linestyle="--", linewidth=1.2)
            ax.text(x, 0.97, " " + text, transform=ax.get_xaxis_transform(),
                    ha="left", va="top", fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("Wrote", out)


def _unit(label):
    m = re.search(r"\[([^\]]+)\]\s*$", label)
    return m.group(1) if m else ""


def cut_labels(spec):
    """Text for each cut line: operator for the kept side, value, unit."""
    if "cut_labels" in spec:
        return spec["cut_labels"]
    unit = _unit(spec["label"])
    fmt = lambda v: ("%g %s" % (v, unit)).strip()
    c = spec["cuts"]
    if spec["keep"] == "above":
        return ["> " + fmt(c[0])]
    if spec["keep"] == "below":
        return ["< " + fmt(c[0])]
    return ["> " + fmt(c[0]), "< " + fmt(c[1])]


def shade_rejected(ax, spec):
    """Shade the rejected side, draw each cut line, and label its value."""
    lo, hi = spec["range"]
    c = spec["cuts"]
    style = dict(color="gray", alpha=0.25, zorder=0, linewidth=0)
    if spec["keep"] == "above":
        ax.axvspan(lo, c[0], **style)
        sides = ["left"]           # text sits on the kept (right) side
    elif spec["keep"] == "below":
        ax.axvspan(c[0], hi, **style)
        sides = ["right"]
    else:
        ax.axvspan(lo, c[0], **style)
        ax.axvspan(c[1], hi, **style)
        sides = ["left", "right"]  # both labels inside the kept window

    for x, text, side in zip(c, cut_labels(spec), sides):
        ax.axvline(x, color="black", linestyle="--", linewidth=1.0)
        ax.text(x, 0.97, text, transform=ax.get_xaxis_transform(), rotation=90,
                ha="left" if side == "left" else "right", va="top", fontsize=8,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=1.0))


def cutvars_figure(df, w, topo, nbins, ylabel, out, logy=False, cfg=None):
    specs = cut_var_specs(topo, cfg)
    nm1 = cuts.nminus1_masks(df, topo, cutbased=cfg if cfg is not None else True)
    names = [n for n in nm1 if n in specs]
    if not names:
        return

    ncol = 4
    nrow = math.ceil(len(names) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.6 * nrow), squeeze=False)
    for ax, name in zip(axes.flat, names):
        spec = specs[name]
        m = nm1[name]
        d, ww = df[m], w[m]
        cat = cuts.truth_category(d)
        bins = np.linspace(*spec["range"], spec.get("nbins", nbins) + 1)
        stack(ax, spec["var"](d), cat, ww, bins, legend=False)
        shade_rejected(ax, spec)
        ax.set_xlabel(spec["label"], fontsize=9)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.set_title("%s: keep %s\n(N-1: %.1f slices)"
                     % (name, " and ".join(cut_labels(spec)), ww.sum()), fontsize=9)
        ax.tick_params(labelsize=8)
        if logy:
            ax.set_yscale("log")
    for ax in list(axes.flat)[len(names):]:
        ax.axis("off")

    handles = [Patch(facecolor=COLORS[c], edgecolor="black", linewidth=0.4, label=l)
               for c, l in cuts.CATEGORY_LABELS.items()]
    handles.append(Patch(facecolor="gray", alpha=0.25, label="Rejected by this cut"))
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=9, frameon=False)
    fig.suptitle("%s cut variables, all other cuts applied" % topo, fontsize=12)
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("Wrote", out)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--inputs", nargs="+", default=[DEFAULT_INPUT],
                    help="grid output directories, .df files, glob patterns (quote them), "
                         "or a .list/.txt of paths (default: %(default)s)")
    ap.add_argument("-o", "--outdir", default="ncdelta_plots")
    ap.add_argument("--selection", default="balanced", choices=sorted(cuts.SELECTIONS),
                    help="which working point in cuts.SELECTIONS to plot (default: %(default)s)")
    ap.add_argument("--pot", type=float, default=None,
                    help="scale MC to this POT (default: unscaled MC counts)")
    ap.add_argument("--emax", type=float, default=3.0, help="upper edge of the true-energy axis [GeV]")
    ap.add_argument("--recoemax", type=float, default=1.5, help="upper edge of the visible-energy axis [GeV]")
    ap.add_argument("--nbins", type=int, default=30)
    ap.add_argument("--logy", action="store_true")
    args = ap.parse_args()

    files = expand_inputs(args.inputs)
    frames, pot_by_file = load(files, keys=("slc",))
    df = frames["slc"]
    mc_pot = pot_by_file.sum()
    print("Read %d slices, MC POT = %.3g" % (len(df), mc_pot))

    if args.pot is not None:
        if mc_pot <= 0:
            sys.exit("MC POT is zero; cannot scale to --pot.")
        scale = args.pot / mc_pot
        ylabel = "Slices / bin (%.2g POT)" % args.pot
    else:
        scale = 1.
        ylabel = "Slices / bin (MC, %.2g POT)" % mc_pot
    w = pd.Series(scale, index=df.index)

    os.makedirs(args.outdir, exist_ok=True)
    true_bins = np.linspace(0., args.emax, args.nbins + 1)
    reco_bins = np.linspace(0., args.recoemax, args.nbins + 1)
    mass_bins = np.arange(0.900, 2.000 + 1e-9, 0.025)   # GeV

    for topo in cuts.TOPOLOGIES:
        cfg = cuts.selection_config(args.selection, topo)
        pre = cuts.PRESELECTION[topo](df)
        sel = cuts.selection(df, topo, cutbased=cfg)
        stages = {}
        for name, m in (("preselection", pre), ("selection", sel)):
            d = df[m]
            stages[name] = (d, w[m], cuts.truth_category(d))

        stage_figure(stages, topo, lambda d: d.true_E, true_bins,
                           "True neutrino energy [GeV]", ylabel,
                           os.path.join(args.outdir, "composition_trueE_%s.png" % topo),
                           include_cosmic=False, logy=args.logy)

        reco_label = (r"Reconstructed visible energy $E_\gamma + T_p$ [GeV]" if topo == "1g1p"
                      else r"Reconstructed visible energy $E_\gamma$ [GeV]")
        stage_figure(stages, topo, lambda d: visible_energy(d), reco_bins,
                           reco_label, ylabel,
                           os.path.join(args.outdir, "composition_recoE_%s.png" % topo),
                           logy=args.logy)

        if topo == "1g1p":            # needs a proton candidate
            stage_figure(stages, topo, lambda d: d.delta_mass, mass_bins,
                               r"Proton + photon invariant mass $M_{p\gamma}$ [GeV]", ylabel,
                               os.path.join(args.outdir, "invmass_%s.png" % topo),
                               logy=args.logy, vline=(DELTA_MASS, r"$\Delta(1232)$"))

        cutvars_figure(df, w, topo, args.nbins, ylabel,
                       os.path.join(args.outdir, "cutvars_%s.png" % topo), logy=args.logy,
                       cfg=cfg)


if __name__ == "__main__":
    main()
