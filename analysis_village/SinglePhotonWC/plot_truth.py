#!/usr/bin/env python3
"""Truth-level, pre-selection look at the inclusive single-photon signal and
its backgrounds (signal definition and categories: cuts.py).

Input: .df files (or directories) with the "mcnu" table from
configs/sp_mc.py. Every true neutrino counts, so these are rates before any
reconstruction; cosmics (no true neutrino) are not in it.

  python analysis_village/SinglePhotonWC/plot_truth.py -i <dir or .df files> -o sp_truth_plots/
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from analysis_village.SinglePhotonWC import cuts  # noqa: E402
from analysis_village.NCDeltaWC import dfio  # noqa: E402

INK2 = "#52514e"
plt.rcParams.update({
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": "#e6e5e1", "grid.linewidth": 0.6,
    "axes.axisbelow": True, "font.size": 9, "legend.frameon": False,
})
SIG_COLOR, BKG_COLOR = "#2a78d6", "#a3a29d"

# Compared samples: signal solid, background dashed
GROUPS = {
    "signal: NC Delta":      dict(cats=[0], color="#2a78d6", ls="-"),
    "signal: pi0 (NC+CC)":   dict(cats=[1, 2], color="#eb6834", ls="-"),
    "bkg: >=2 showers":      dict(cats=[4], color="#1baf7a", ls="--"),
    "bkg: CC numu mu>=100":  dict(cats=[6], color="#eda100", ls="--"),
    "bkg: CC nue":           dict(cats=[7], color="#7f5bd5", ls="--"),
}
PI0_REASONS = ["photons merged (< 20 deg)", "other photon escapes / < 3 cm from wall",
               "other photon < 20 MeV", "other"]


def groups(m):
    return {k: m[m.cat.isin(v["cats"])] for k, v in GROUPS.items()}


def hist(ax, G, col, bins, xlabel, scale=1., discrete=False):
    """Area-normalised steps; out-of-range values go in the first / last bin."""
    for name, d in G.items():
        v = d[col].dropna() * scale
        if not len(v):
            continue
        ax.hist(v.clip(bins[0], bins[-1] - (0.5 if discrete else 0)), bins=bins, histtype="step", lw=2,
                density=True, label="%s (%d)" % (name, len(v)),
                color=GROUPS[name]["color"], ls=GROUPS[name]["ls"])
    ax.set_xlabel(xlabel + " (last bin: overflow)")
    ax.set_ylabel("area-normalised")


def pi0_reason(d):
    r = np.select([d.true_nmerged > 0, d.true_ngamma_out > 0, d.true_ndecayg_soft > 0],
                  PI0_REASONS[:3], PI0_REASONS[3])
    return pd.Series(r, index=d.index)


def page_rates(pdf, m, pot):
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw=dict(width_ratios=[1.4, 1]))
    ax = axs[0]
    rate = m.cat.value_counts().reindex(cuts.CATEGORY_LABELS.keys()).fillna(0).drop(9) / pot * 1e20
    y = np.arange(len(rate))
    ax.barh(y, rate.values, height=0.65,
            color=[SIG_COLOR if c in cuts.SIGNAL_CATEGORIES else BKG_COLOR for c in rate.index])
    for yi, v in zip(y, rate.values):
        ax.text(v, yi, "  %.3g" % v, va="center", color=INK2, fontsize=8)
    ax.set_yticks(y, [cuts.CATEGORY_LABELS[c] for c in rate.index])
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("true neutrinos per 1e20 POT (log)")
    ax.set_title("All true neutrinos by category (blue: signal)", loc="left")
    ax.grid(axis="y", visible=False)

    ax = axs[1]
    sp = m[m.cat.isin([1, 2])]
    v = pi0_reason(sp).value_counts(normalize=True).reindex(PI0_REASONS).fillna(0)
    ax.barh(np.arange(len(v)), v.values, color=SIG_COLOR, height=0.6)
    for yi, w in enumerate(v.values):
        ax.text(w + 0.01, yi, "%.0f%%" % (100 * w), va="center", color=INK2)
    ax.set_yticks(np.arange(len(v)), v.index)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of pi0 signal (%d)" % len(sp))
    ax.set_title("pi0 signal: why only one shower", loc="left")
    ax.grid(False)
    fig.subplots_adjust(left=0.15, right=0.97, bottom=0.12, wspace=0.75)
    pdf.savefig(fig)
    return fig


def page_kinematics(pdf, G):
    fig, axs = plt.subplots(2, 3, figsize=(15, 8.4))
    hist(axs[0, 0], G, "true_shw1_E", np.linspace(0, 1000, 41), "leading shower E [MeV]", scale=1000.)
    axs[0, 0].set_title("Leading reconstructable shower energy", loc="left")
    hist(axs[0, 1], G, "true_shw1_gap", np.linspace(0, 100, 41), "vertex to shower start [cm]")
    axs[0, 1].set_title("Leading shower conversion distance", loc="left")
    hist(axs[0, 2], G, "true_shw1_costh", np.linspace(-1, 1, 41), "leading shower cos(theta_beam)")
    axs[0, 2].set_title("Leading shower angle", loc="left")
    hist(axs[1, 0], G, "true_E", np.linspace(0, 3, 31), "true E$_\\nu$ [GeV]")
    axs[1, 0].set_title("Neutrino energy", loc="left")
    hist(axs[1, 1], G, "true_np_35MeV", np.arange(0, 6) - 0.5, "protons above 35 MeV KE", discrete=True)
    axs[1, 1].set_title("Proton multiplicity", loc="left")
    hist(axs[1, 2], G, "true_nshower", np.arange(0, 6) - 0.5, "reconstructable showers", discrete=True)
    axs[1, 2].set_title("Shower multiplicity", loc="left")
    axs[0, 0].legend(fontsize=8)
    fig.tight_layout()
    pdf.savefig(fig)
    return fig


def page_extra(pdf, m):
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.2))
    G = groups(m)
    hist(axs[0], {k: G[k] for k in ["signal: pi0 (NC+CC)", "bkg: >=2 showers"]}, "true_shw2_E",
         np.linspace(0, 400, 41), "sub-leading shower E [MeV]", scale=1000.)
    axs[0].set_title("Sub-leading shower (signal: none)", loc="left")
    axs[0].legend(fontsize=8)
    cc1 = m[(m.true_iscc == 1) & (m.true_nshower == 1) & (m.true_ne == 0) & (m.cat != 8)]
    axs[1].hist((cc1.true_mu_KE * 1000).clip(0, 999), bins=np.linspace(0, 1000, 41), histtype="step", lw=2,
                color=INK2)
    axs[1].axvline(cuts.MUON_KE_MAX * 1000, color="#eb6834", lw=1.5)
    axs[1].set_xlabel("leading muon KE [MeV] (last bin: overflow)")
    axs[1].set_ylabel("true neutrinos")
    axs[1].set_title("CC numu with one shower: muon KE (line: 100 MeV)", loc="left")
    hist(axs[2], G, "true_npipm", np.arange(0, 5) - 0.5, "charged pions (any KE)", discrete=True)
    axs[2].set_title("Charged-pion multiplicity", loc="left")
    fig.tight_layout()
    pdf.savefig(fig)
    return fig


def summary(m, pot):
    print("POT %.4g, %d true neutrinos" % (pot, len(m)))
    n = m.cat.value_counts()
    sig = n.reindex(cuts.SIGNAL_CATEGORIES).fillna(0).sum()
    for c, lab in cuts.CATEGORY_LABELS.items():
        if c == 9:
            continue
        k = n.get(c, 0)
        frac = "  %4.1f%% of signal" % (100 * k / sig) if c in cuts.SIGNAL_CATEGORIES and sig else ""
        print("  %-28s %8d  %10.4g / 1e20 POT%s" % (lab, k, k / pot * 1e20, frac))
    s = m[m.cat.isin(cuts.SIGNAL_CATEGORIES)]
    if len(s):
        print("  signal: shower E median %.0f MeV, gap median %.1f cm, %d%% have a merged photon pair"
              % (1000 * s.true_shw1_E.median(), s.true_shw1_gap.median(), 100 * (s.true_shw1_nphot > 1).mean()))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--inputs", nargs="+", required=True, help=".df files, directories or globs")
    ap.add_argument("-o", "--outdir", default="sp_truth_plots")
    args = ap.parse_args()

    frames, totals = dfio.load(dfio.expand_inputs(args.inputs), keys=["mcnu"])
    m = frames["mcnu"]
    if "true_nshower" not in m:
        sys.exit("mcnu has no single-photon truth columns; make it with SinglePhotonWC/configs/sp_mc.py.")
    m["cat"] = cuts.truth_category(m)
    pot = totals.pot.sum()
    summary(m, pot)

    os.makedirs(args.outdir, exist_ok=True)
    out = os.path.join(args.outdir, "sp_truth.pdf")
    inav = m[m.cat != 8]
    with PdfPages(out) as pdf:
        figs = [page_rates(pdf, m, pot), page_kinematics(pdf, groups(inav)), page_extra(pdf, inav)]
        for i, fig in enumerate(figs):
            fig.savefig(os.path.join(args.outdir, "sp_truth_%02d.png" % i), dpi=130)
            plt.close(fig)
    print("Wrote %s (%d pages)" % (out, len(figs)))


if __name__ == "__main__":
    main()
