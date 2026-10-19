#!/usr/bin/env python3
"""Truth-composition plots for the NCDelta dataframes.

For each channel (1g1p, 1g0p), stacked by the truth categories in cuts.py:

  composition_trueE_<ch>.png  true neutrino energy
  composition_recoE_<ch>.png  reconstructed visible energy, E_gamma + T_p
  invmass_1g1p.png            proton + photon invariant mass, with the
                              Delta(1232) mass marked (1g1p only)
  cutvars_<ch>.png            N-1 distribution of each cut variable: every
                              other cut applied, threshold drawn, rejected
                              side shaded
  yields_<ch>.csv             yields per category at each stage

The energy and mass plots have three stages:
  dataframe     every slice in the .df (already topology + pre-selection
                for ncdelta_mc.py; all non-clear-cosmic slices for
                ncdelta_mc_nopresel.py)
  preselection  topology + the fiducial volume, which is all the
                pre-selection now is (see cuts.py)
  selection     full selection

Visible energy is used because the outgoing neutrino in an NC event is not
measured, so the true neutrino energy cannot be reconstructed.

Cosmic slices have no true neutrino energy. In the true-energy plots they
appear only as a yield in the legend; elsewhere they are stacked normally.

Pre-selection is now only the fiducial volume, which depends on all three
vertex coordinates and is not drawn; every cut that does get an N-1 panel
lives in CUTBASED_*. On dataframes made with ncdelta_mc.py the FV cut is
already applied, so use the ncdelta_mc_nopresel.py ones to see what it
removes.

Inputs (-i) can be grid output directories, .df files, glob patterns, or
.list/.txt files of paths. A directory contributes every .df file in it.
Files that cannot be read (e.g. from crashed grid jobs) are skipped with a
warning, and their POT is not counted.

Example (from the repository root):
  python analysis_village/NCDelta/plot_composition.py \\
      -i /pnfs/sbnd/scratch/users/nrowe/cafpyana_out/dfs/2026_10_17_164529__ncdelta.df/ \\
      --pot 1e20 -o plots/
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
from dfio import expand_inputs, load   # noqa: E402

DEFAULT_INPUT = "/pnfs/sbnd/scratch/users/nrowe/cafpyana_out/dfs/2026_10_17_164529__ncdelta.df/"

DELTA_MASS = 1232.0  # MeV, PDG Delta(1232)
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
    "dataframe": "All slices in dataframe",
    "preselection": "{topo} topology + pre-selection",
    "selection": "{topo} full selection",
}


# ---------------------------------------------------------------------------
# Derived quantities
# ---------------------------------------------------------------------------

def visible_energy(df):
    """Shower energy plus proton KE from range [MeV]; shower only without a track."""
    return df.shw_energy_best + df.trk_ke_p.where(df.n_trk > 0, 0.).fillna(0.)


# ---------------------------------------------------------------------------
# Cut variables: what to plot for each stage name in cuts.stage_cuts()
# ---------------------------------------------------------------------------
# keep: "above" (pass > cut), "below" (pass < cut), "between" (lo < x < hi).
# Thresholds are read from cuts.py so the plots follow any changes there.

def cut_var_specs(topo, C=None):
    C = cuts.CUTBASED[topo] if C is None else C
    lo_E, hi_E = cuts.energy_window(topo, C)
    s = {}

    # These used to be pre-selection; they now live in CUTBASED_* too, and
    # each is drawn only when it is enabled there. The pre-selection proper is
    # the fiducial volume, which depends on all three coordinates and is not
    # drawn (see the module docstring).
    if topo == "1g1p":
        if C.get("vtx_dist_min") is not None:
            s["vtx_contained"] = dict(var=lambda d: d.vtx_dist_to_wall, range=(-10, 100),
                                      label="Vertex distance to wall [cm]",
                                      cuts=[C["vtx_dist_min"]], keep="above")
        if C.get("trk_dist_min") is not None:
            s["trk_contained"] = dict(var=lambda d: np.minimum(d.trk_start_dist_to_wall, d.trk_end_dist_to_wall),
                                      range=(-10, 100), label="Track min distance to wall [cm]",
                                      cuts=[C["trk_dist_min"]], keep="above")
        if C.get("trk_dedx_min") is not None:
            s["trk_dedx"] = dict(var=lambda d: d.trk_dedx_trunc, range=(0, 15),
                                 label="Track truncated-mean dE/dx [MeV/cm]",
                                 cuts=[C["trk_dedx_min"]], keep="above")
        if C.get("abs_cos_trk_shw_max") is not None:
            s["trk_shw_angle"] = dict(var=lambda d: d.cos_trk_shw_implied, range=(-1, 1),
                                      label=r"cos$\theta$(track, vertex$\to$shower)",
                                      cuts=[-C["abs_cos_trk_shw_max"], C["abs_cos_trk_shw_max"]],
                                      keep="between")
        if C.get("shw_start_dist_min") is not None:
            s["shw_contained"] = dict(var=lambda d: d.shw_start_dist_to_wall, range=(-10, 100),
                                      label="Shower start distance to wall [cm]",
                                      cuts=[C["shw_start_dist_min"]], keep="above")
        if C.get("trk_ke_max") is not None:
            s["trk_ke"] = dict(var=lambda d: d.trk_ke_p, range=(0, 1000),
                               label="Proton KE from range [MeV]",
                               cuts=[C["trk_ke_max"]], keep="below")
        if C.get("trk_calo_E_max") is not None:
            s["trk_calo_E"] = dict(var=lambda d: d.trk_calo_E, range=(0, 1000),
                                   label="Track calorimetric energy [MeV]",
                                   cuts=[C["trk_calo_E_max"]], keep="below")
    else:
        if C.get("shw_start_dist_min") is not None:
            s["shw_contained"] = dict(var=lambda d: d.shw_start_dist_to_wall, range=(-10, 100),
                                      label="Shower start distance to wall [cm]",
                                      cuts=[C["shw_start_dist_min"]], keep="above")

    # cut-based stage (only stages enabled in CUTBASED appear in the masks)
    if C.get("nu_score_min") is not None:
        s["nu_score"] = dict(var=lambda d: d.nu_score, range=(0, 1),
                             label="Pandora neutrino score", cuts=[C["nu_score_min"]], keep="above")
    if C.get("bary_score_min") is not None:
        s["flash_match"] = dict(var=lambda d: d.bary_score, range=(0, 1),
                                label="Barycenter flash-match score",
                                cuts=[C["bary_score_min"]], keep="above")
    if C.get("shw_dedx_min") is not None:
        s["shw_dedx"] = dict(var=lambda d: d.shw_dedx_best, range=(0, 10),
                             label="Shower dE/dx, best plane [MeV/cm]",
                             cuts=[C["shw_dedx_min"]], keep="above")
    if C.get("shw_conv_dist_min") is not None:
        s["shw_conv_dist"] = dict(var=lambda d: d.shw_conv_dist, range=(0, 100),
                                  label="Shower conversion distance [cm]",
                                  cuts=[C["shw_conv_dist_min"]], keep="above")
    if C.get("shw_conv_dist_max") is not None:
        s["shw_conv_dist_max"] = dict(var=lambda d: d.shw_conv_dist,
                                      range=(0, 4 * C["shw_conv_dist_max"]),
                                      label="Shower conversion distance [cm]",
                                      cuts=[C["shw_conv_dist_max"]], keep="below")
    if C.get("shw_len_min") is not None:
        s["shw_len"] = dict(var=lambda d: d.shw_len, range=(0, 200),
                            label="Shower length [cm]",
                            cuts=[C["shw_len_min"]], keep="above")
    if C.get("shw_nhits_min") is not None:
        s["shw_nhits"] = dict(var=lambda d: d.shw_nhits_I0 + d.shw_nhits_I1 + d.shw_nhits_I2,
                              range=(0, 1000), label="Shower hits, all planes",
                              cuts=[C["shw_nhits_min"]], keep="above")
    if C.get("shw_open_angle_max") is not None:
        s["shw_open_angle"] = dict(var=lambda d: d.shw_open_angle, range=(0, 1),
                                   label="Shower opening angle [rad]",
                                   cuts=[C["shw_open_angle_max"]], keep="below")
    if C.get("opt0_score_min") is not None:
        # The score runs to ~5000, so the axis is scaled to the cut instead:
        # what matters here is the rejected low tail, not the bulk.
        s["opt0_score"] = dict(var=lambda d: d.opt0_score,
                               range=(0, 6 * C["opt0_score_min"]),
                               label="OpT0 flash-match score",
                               cuts=[C["opt0_score_min"]], keep="above")
    if C.get("n_nonprim_shw_max") is not None:
        s["no_second_shw"] = dict(var=lambda d: d.n_nonprim_shw, range=(-0.5, 5.5), nbins=6,
                                  label="Non-primary showers in slice",
                                  cuts=[C["n_nonprim_shw_max"] + 0.5], keep="below",
                                  cut_labels=[r"$\leq$ %d" % C["n_nonprim_shw_max"]])
    if C.get("delta_mass_range") is not None:
        s["delta_mass"] = dict(var=lambda d: d.delta_mass, range=(900, 2000),
                               label="Proton + photon invariant mass [MeV]",
                               cuts=list(C["delta_mass_range"]), keep="between")
    if C.get("trk_proton_like"):
        s["trk_proton_like"] = dict(var=lambda d: d.trk_chi2_proton - d.trk_chi2_muon, range=(-100, 100),
                                    label=r"$\chi^2_p - \chi^2_\mu$", cuts=[0.], keep="below")
    if C.get("trk_start_dist_to_vtx_max") is not None:
        s["trk_attached"] = dict(var=lambda d: d.trk_start_dist_to_vtx, range=(0, 30),
                                 label="Track start distance to vertex [cm]",
                                 cuts=[C["trk_start_dist_to_vtx_max"]], keep="below")
    if C.get("trk_pida_max") is not None:
        s["trk_pida"] = dict(var=lambda d: d.trk_pida, range=(0, 40),
                             label="Track PIDA [MeV/cm$^{1.42}$]",
                             cuts=[C["trk_pida_max"]], keep="below")

    s["energy_window"] = dict(var=lambda d: d.shw_energy_best, range=(0, 1000),
                              label="Shower energy, best plane (signal window) [MeV]",
                              cuts=[lo_E, hi_E], keep="between",
                              cut_labels=[r"$\geq$ %g MeV" % lo_E, "< %g MeV" % hi_E])
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


def three_stage_figure(stages, topo, values_fn, bins, xlabel, ylabel, out,
                       include_cosmic=True, logy=False, vline=None):
    names = [n for n in ("dataframe", "preselection", "selection") if n in stages]
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


def yield_table(stages):
    rows = {}
    for name, (d, w, cat) in stages.items():
        rows[name] = {label: w[(cat == code).values].sum()
                      for code, label in cuts.CATEGORY_LABELS.items()}
        rows[name]["total"] = w.sum()
    return pd.DataFrame(rows).T


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
    mass_bins = np.arange(900., 2000. + 1e-6, 25.)

    for topo in cuts.TOPOLOGIES:
        cfg = cuts.selection_config(args.selection, topo)
        pre = cuts.PRESELECTION[topo](df)
        sel = cuts.selection(df, topo, cutbased=cfg)
        stages = {}
        for name, m in (("dataframe", slice(None)), ("preselection", pre), ("selection", sel)):
            d = df[m]
            stages[name] = (d, w[m], cuts.truth_category(d))

        three_stage_figure(stages, topo, lambda d: d.true_E / 1000., true_bins,
                           "True neutrino energy [GeV]", ylabel,
                           os.path.join(args.outdir, "composition_trueE_%s.png" % topo),
                           include_cosmic=False, logy=args.logy)

        reco_label = (r"Reconstructed visible energy $E_\gamma + T_p$ [GeV]" if topo == "1g1p"
                      else r"Reconstructed visible energy $E_\gamma$ [GeV]")
        three_stage_figure(stages, topo, lambda d: visible_energy(d) / 1000., reco_bins,
                           reco_label, ylabel,
                           os.path.join(args.outdir, "composition_recoE_%s.png" % topo),
                           logy=args.logy)

        if topo == "1g1p":
            # Needs a proton candidate, so only the 1g1p stages are drawn.
            mass_stages = {k: v for k, v in stages.items() if k != "dataframe"}
            three_stage_figure(mass_stages, topo, lambda d: d.delta_mass, mass_bins,
                               r"Proton + photon invariant mass $M_{p\gamma}$ [MeV]", ylabel,
                               os.path.join(args.outdir, "invmass_%s.png" % topo),
                               logy=args.logy, vline=(DELTA_MASS, r"$\Delta(1232)$"))

        cutvars_figure(df, w, topo, args.nbins, ylabel,
                       os.path.join(args.outdir, "cutvars_%s.png" % topo), logy=args.logy,
                       cfg=cfg)

        table = yield_table(stages)
        print("\n%s yields by truth category:\n%s\n" % (topo, table.round(2).to_string()))
        table.to_csv(os.path.join(args.outdir, "yields_%s.csv" % topo))


if __name__ == "__main__":
    main()
