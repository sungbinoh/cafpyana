#!/usr/bin/env python3
"""Truth-level look at the NC Delta -> N gamma signal and the NC pi0 background.

Input: .df files (or directories of them) with an "mcnu" table: one row per
true neutrino with the make_ncdeltawc_truthdf columns. Tables from
configs/wc_candidates_mc.py (make_truth.py) also have photon gaps, the pi0
decay photons and the leading proton; their pages are added when those
columns exist. Only neutrinos with a true vertex in the active volume count.

  signal   cuts.is_signal, and a clean final state: exactly one photon above
           20 MeV, no charged pion above 30 MeV KE (0 pi0 is in is_signal)
  NC 1pi0  NC, exactly one pi0, no charged pion above 30 MeV KE

Both are split by the number of protons above PROTON_COL's threshold.
Distributions are area-normalised; the legend gives each sample's count.

  python analysis_village/NCDeltaWC/plot_truth_study.py -i <dir or .df files> -o truth_plots/
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
from analysis_village.NCDeltaWC import cuts, dfio  # noqa: E402

INK2 = "#52514e"
plt.rcParams.update({
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": "#e6e5e1", "grid.linewidth": 0.6,
    "axes.axisbelow": True, "font.size": 9, "legend.frameon": False,
})

# Samples: signal solid, background dashed, so identity is not colour alone
STYLE = {
    "1g0p":       dict(color="#2a78d6", ls="-"),
    "1g1p":       dict(color="#eb6834", ls="-"),
    "NC1pi0 0p":  dict(color="#1baf7a", ls="--"),
    "NC1pi0 1p":  dict(color="#eda100", ls="--"),
}
PROTON_THRESHOLDS = [("true_np", "any"), ("true_np_20MeV", "> 20 MeV"),
                     ("true_np_35MeV", "> 35 MeV"), ("true_np_50MeV", "> 50 MeV")]
GENIE_MODES = {0: "QE", 1: "Res", 2: "DIS", 3: "Coh", 10: "MEC"}


def samples(m, pcol):
    clean = (m.true_npi_30MeV == 0)
    sig = cuts.is_signal(m) & (m.true_ngamma_20MeV == 1) & clean
    pi0 = (m.true_isnc == 1) & (m.true_npi0 == 1) & clean & cuts.true_in_active(m)
    np_ = m[pcol]
    return {"1g0p": m[sig & (np_ == 0)], "1g1p": m[sig & (np_ == 1)],
            "NC1pi0 0p": m[pi0 & (np_ == 0)], "NC1pi0 1p": m[pi0 & (np_ == 1)]}


def hist(ax, S, col, bins, xlabel, scale=1., names=None):
    """Step histograms; col may be a list (values pooled). Out-of-range values
    go in the first / last bin."""
    cols = col if isinstance(col, list) else [col]
    for name in names or S:
        v = pd.concat([S[name][c] for c in cols]).dropna() * scale
        if not len(v):
            continue
        ax.hist(v.clip(bins[0], bins[-1]), bins=bins, histtype="step", lw=2, density=True,
                label="%s (%d)" % (name, len(v)), **STYLE[name])
    ax.set_xlabel(xlabel + " (last bin: overflow)")
    ax.set_ylabel("area-normalised")
    ax.legend(fontsize=8)


def page_topology(pdf, m, sig_all):
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axs[0]
    rows = [(c, lab) for c, lab in PROTON_THRESHOLDS if c in m]
    fr = pd.DataFrame({lab: sig_all[c].clip(upper=2).value_counts(normalize=True).reindex([0, 1, 2]).fillna(0)
                       for c, lab in rows}).T
    y = np.arange(len(fr))
    left = np.zeros(len(fr))
    for k, (lab, color) in enumerate(zip(["0p", "1p", "2+p"], ["#2a78d6", "#eb6834", "#1baf7a"])):
        ax.barh(y, fr[k], left=left, color=color, edgecolor="white", lw=2, label=lab, height=0.6)
        for yi, (l, w) in enumerate(zip(left, fr[k])):
            if w > 0.06:
                ax.text(l + w / 2, yi, "%.0f%%" % (100 * w), ha="center", va="center", color="white", fontsize=9)
        left += fr[k]
    ax.set_yticks(y, ["protons " + i for i in fr.index])
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of clean NC Delta -> N gamma (%d)" % len(sig_all))
    ax.grid(False)
    ax.legend(ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    ax.set_title("Signal proton multiplicity vs. threshold", loc="left", pad=22)

    ax = axs[1]
    pi = m[(m.true_isnc == 1) & (m.true_npi0 == 1) & cuts.true_in_active(m)]
    v = pi.true_genie_mode.map(GENIE_MODES).fillna("other").value_counts(normalize=True)
    ax.barh(v.index, v.values, color="#1baf7a", height=0.6)
    for yi, w in enumerate(v.values):
        ax.text(w + 0.01, yi, "%.1f%%" % (100 * w), va="center", color=INK2)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of NC 1pi0 (%d)" % len(pi))
    ax.set_title("NC 1pi0 by GENIE mode", loc="left")
    ax.grid(False)
    fig.tight_layout()
    pdf.savefig(fig)
    return fig


def page_kinematics(pdf, S, pcol_label):
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))
    hist(axs[0], S, "true_leading_gamma_E", np.linspace(0, 800, 41), "leading GENIE photon E [MeV]",
         scale=1000., names=["1g0p", "1g1p"])
    axs[0].set_title("Signal photon energy", loc="left")
    hist(axs[1], S, "true_E", np.linspace(0, 3, 31), "true E$_\\nu$ [GeV]")
    axs[1].set_title("Neutrino energy", loc="left")
    hist(axs[2], S, "true_w", np.linspace(1.0, 2.0, 41), "true W [GeV]")
    axs[2].set_title("Hadronic invariant mass", loc="left")
    fig.suptitle("protons counted %s" % pcol_label, x=0.01, ha="left", color=INK2)
    fig.tight_layout()
    pdf.savefig(fig)
    return fig


def page_photons(pdf, S):
    """make_truth.py columns: gaps, pi0 photons, leading proton."""
    fig, axs = plt.subplots(2, 3, figsize=(14, 8))
    hist(axs[0, 0], S, "g1_gap", np.linspace(0, 100, 41), "signal photon conversion gap [cm]", names=["1g0p", "1g1p"])
    axs[0, 0].set_title("Signal photon gap (converts in TPC)", loc="left")
    hist(axs[0, 1], S, ["pi0_g1_gap", "pi0_g2_gap"], np.linspace(0, 100, 41), "pi0 photon gap [cm]",
         names=["NC1pi0 0p", "NC1pi0 1p"])
    axs[0, 1].set_title("pi0 photon gap (both photons, converts in TPC)", loc="left")
    hist(axs[0, 2], S, "pi0_g2_E", np.linspace(0, 400, 41), "pi0 sub-leading photon E [MeV]", scale=1000.,
         names=["NC1pi0 0p", "NC1pi0 1p"])
    axs[0, 2].set_title("pi0 sub-leading photon energy", loc="left")
    hist(axs[1, 0], S, "pi0_open", np.linspace(0, 180, 37), "pi0 opening angle [deg]", names=["NC1pi0 0p", "NC1pi0 1p"])
    axs[1, 0].set_title("pi0 opening angle", loc="left")
    hist(axs[1, 1], S, "pi0_asym", np.linspace(0, 1, 26), "pi0 (E1 - E2) / (E1 + E2)", names=["NC1pi0 0p", "NC1pi0 1p"])
    axs[1, 1].set_title("pi0 energy asymmetry", loc="left")
    hist(axs[1, 2], S, "p1_KE", np.linspace(0, 600, 31), "leading proton KE [MeV]", scale=1000.,
         names=["1g1p", "NC1pi0 1p"])
    axs[1, 2].set_title("Leading proton KE (1p samples)", loc="left")
    fig.tight_layout()
    pdf.savefig(fig)
    return fig


def summary(S, pot):
    print("POT %.3g" % pot)
    for name, d in S.items():
        print("  %-10s %7d  (%.1f per 1e20 POT)" % (name, len(d), len(d) / pot * 1e20))
    if "g1_gap" in S["1g0p"]:
        for name in ["1g0p", "1g1p"]:
            g = S[name].g1_gap
            print("  %s gap: median %.1f cm, < 3 cm %.0f%%, lost/out of TPC %.0f%%" % (
                name, g.median(), 100 * (g < 3).mean(), 100 * g.isna().mean()))
        for name in ["NC1pi0 0p", "NC1pi0 1p"]:
            d = S[name]
            print("  %s: E2 < 30 MeV %.0f%%, open < 20 deg %.0f%%, a photon leaves TPC %.0f%%" % (
                name, 100 * (d.pi0_g2_E < 0.03).mean(), 100 * (d.pi0_open < 20).mean(),
                100 * ((d.pi0_g1_lost > 0) | (d.pi0_g2_lost > 0)).mean()))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--inputs", nargs="+", required=True, help=".df files, directories or globs")
    ap.add_argument("-o", "--outdir", default="truth_plots")
    ap.add_argument("--proton-col", default=None,
                    help="proton count for 0p/1p (default true_np_35MeV if present, else true_np_20MeV)")
    args = ap.parse_args()

    frames, totals = dfio.load(dfio.expand_inputs(args.inputs), keys=["mcnu"])
    m = frames["mcnu"]
    m = m[cuts.true_in_active(m)]
    pcol = args.proton_col or ("true_np_35MeV" if "true_np_35MeV" in m else "true_np_20MeV")
    plab = dict(PROTON_THRESHOLDS).get(pcol, pcol)
    S = samples(m, pcol)
    sig_all = m[cuts.is_signal(m) & (m.true_ngamma_20MeV == 1) & (m.true_npi_30MeV == 0)]
    summary(S, totals.pot.sum())

    os.makedirs(args.outdir, exist_ok=True)
    out = os.path.join(args.outdir, "truth_study.pdf")
    with PdfPages(out) as pdf:
        figs = [page_topology(pdf, m, sig_all), page_kinematics(pdf, S, plab)]
        if "g1_gap" in m:
            figs.append(page_photons(pdf, S))
        for i, fig in enumerate(figs):
            fig.savefig(os.path.join(args.outdir, "truth_%02d.png" % i), dpi=130)
            plt.close(fig)
    print("Wrote %s (%d pages)" % (out, len(figs)))


if __name__ == "__main__":
    main()
