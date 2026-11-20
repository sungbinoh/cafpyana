"""True neutrino energy vs (MC-scaled) flash PE for the selected 1uNp MC CV.

One 2D histogram per run set -- SBND, ICARUS Run 2, ICARUS Run 4 -- of the
selected 1u1p + 1uMp events (the signal-box selection, see
signalbox_systematics_gumple.add_selection) in the MC CV sample.

flash_maxpe is the per-run scaled MC PE: loaddf.load_one multiplies the MC flash
PE by the all-trigger data/MC scale (loaddf.MC_PE_SCALE: SBND 0.647, ICARUS Run2 0.574, Run4 0.369)
when the frame is loaded with the matching detector=. The flash cut is NOT
applied (PE is the plotted axis); its threshold is drawn as a dashed line.
x = scaled flash PE (linear), y = true neutrino energy.

Outputs, in <plotdir>/{png,pdf}/:
  <det>_trueE_vs_flashpe          weighted counts / 1e20 POT, log colour scale
and <plotdir>/summary.txt.

    python flashpe_vs_trueE_gumple.py [--df-dir DIR] [--plotdir DIR] [--nsbnd N]
"""

import argparse
import os
import sys
import warnings

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
# repo root first (pyanalib, makedf package), then ../gumple (gumple_cuts),
# and the gump dir LAST: gump/makedf.py would shadow the makedf package
sys.path.insert(0, _ROOT)
sys.path.insert(1, os.path.abspath(os.path.join(_HERE, "..", "gumple")))
sys.path.append(_HERE)

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

import loaddf
import signalbox_systematics_gumple as sb
import gumple_cuts as gc

DEFAULT_DF_DIR = "/Users/gputnam/Work/osc/sbn-rewgted-25/"
DEFAULT_PLOTDIR = "/Users/gputnam/Work/osc/cafpyana/plots-gumple-2026-10-24-rewgted25/flashpe_vs_trueE/"

POT = 1e20
DETECTORS = ["SBND", "ICARUS Run2", "ICARUS Run4"]
# gumple_cuts.flash_cut thresholds on the scaled flash_maxpe
FLASH_THRESHOLD = gc.FLASH_PE_THRESHOLD
DET_LABEL = {"SBND": "SBND", "ICARUS Run2": "ICARUS Run 2", "ICARUS Run4": "ICARUS Run 4"}

E_BINS = np.arange(0., 3.0 + 1e-9, 0.05)
FONTSIZE = 14


# NB: keep at module scope -- loaddf's cache key depends on the preselection's qualname
def FV(df):
    return sb.FV(df)


def mc_files(df_dir, detector, nsbnd):
    if detector == "SBND":
        return [df_dir + "SBNDMCCV_%i.df" % i for i in range(nsbnd)]
    if detector == "ICARUS Run2":
        return [df_dir + "ICARUSRun2_SpringMCOverlay_rewgt_%i.df" % i for i in range(3)]
    if detector == "ICARUS Run4":
        return [df_dir + "ICARUSRun4_SpringMCOverlay_rewgt_%i.df" % i for i in range(4)]
    raise ValueError(detector)


def load_selected(df_dir, detector, nsbnd, log):
    files = mc_files(df_dir, detector, nsbnd)
    for f in files:
        assert os.path.exists(f), f
    df, _, pot = loaddf.loadl(files, njob=len(files), detector=detector, preselection=FV,
                              reweight_aFF=True, include_syst=False)
    loaddf.scale_pot(df, pot, POT)
    sb.add_selection(df)
    n1p, nmp = int(df.selected_noflash_1p.sum()), int(df.selected_noflash_mp.sum())
    sel = df[df.selected_noflash_1p | df.selected_noflash_mp]
    nnan = int(sel.true_E.isna().sum())
    sel = sel[sel.true_E.notna()]
    thr = FLASH_THRESHOLD[detector]
    w = sel.glob_scale
    frac_above = float(w[sel.flash_maxpe > thr].sum() / w.sum())
    stats = dict(pot=pot, n1p=n1p, nmp=nmp, nnan=nnan, nsel=len(sel),
                 wsel=float(w.sum()), frac_above=frac_above,
                 pe_median=float(sel.flash_maxpe.median()))
    log("%s: %d files, MC POT %.3e; selected (no flash cut) 1p=%d Mp=%d, dropped %d with no true nu; "
        "weighted %.1f / 1e20 POT; above %g PE: %.1f%%; median PE %.0f" % (
            detector, len(files), pot, n1p, nmp, nnan, stats["wsel"], thr,
            100 * frac_above, stats["pe_median"]))
    return sel[["true_E", "flash_maxpe", "glob_scale"]].copy(), stats


def pe_bins(sel):
    # linear bins from 0 to the 99.5% PE quantile, rounded up to 5000 PE
    hi = 5000. * np.ceil(np.quantile(sel.flash_maxpe, 0.995) / 5000.)
    return np.linspace(0., hi, 61)


def draw(ax, sel, detector):
    xbins = pe_bins(sel)
    x, y, w = sel.flash_maxpe.to_numpy(), sel.true_E.to_numpy(), sel.glob_scale.to_numpy()
    x = np.clip(x, xbins[0], xbins[-1] * 0.999)  # keep overflow PE on the plot edge
    H, _, _ = np.histogram2d(x, y, bins=[xbins, E_BINS], weights=w)
    clabel = "Events / $10^{20}$ POT"
    Hm = np.ma.masked_where(H <= 0, H)
    pos = H[H > 0]
    norm = mcolors.LogNorm(vmin=max(pos.min(), pos.max() * 1e-4), vmax=pos.max())
    mesh = ax.pcolormesh(xbins, E_BINS, Hm.T, norm=norm, cmap="viridis")
    cb = ax.figure.colorbar(mesh, ax=ax)
    cb.set_label(clabel, fontsize=FONTSIZE - 2)

    ax.axvline(FLASH_THRESHOLD[detector], color="red", ls="--", lw=2)
    ax.set_xlim(xbins[0], xbins[-1])
    ax.set_ylim(E_BINS[0], E_BINS[-1])
    ax.set_xlabel("Scaled MC flash PE", fontsize=FONTSIZE)
    ax.set_ylabel("True $E_\\nu$ [GeV]", fontsize=FONTSIZE)
    ax.set_title("%s MC, $1\\mu Np$ Selection" % DET_LABEL[detector], fontsize=FONTSIZE)


def save(fig, plotdir, name):
    for ext in ["png", "pdf"]:
        fig.savefig(os.path.join(plotdir, ext, name + "." + ext), bbox_inches="tight", dpi=150)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--df-dir", default=DEFAULT_DF_DIR)
    ap.add_argument("--plotdir", default=DEFAULT_PLOTDIR)
    ap.add_argument("--nsbnd", type=int, default=4, help="number of SBNDMCCV_<i>.df files to load")
    args = ap.parse_args(argv)
    df_dir = args.df_dir if args.df_dir.endswith("/") else args.df_dir + "/"
    for ext in ["png", "pdf"]:
        os.makedirs(os.path.join(args.plotdir, ext), exist_ok=True)

    lines = []

    def log(msg):
        print(msg, flush=True)
        lines.append(msg)

    log("DF_DIR: %s" % df_dir)
    sels, stats = {}, {}
    for det in DETECTORS:
        sels[det], stats[det] = load_selected(df_dir, det, args.nsbnd, log)

    for det in DETECTORS:
        fig, ax = plt.subplots(figsize=(8, 6))
        draw(ax, sels[det], det)
        save(fig, args.plotdir, "%s_trueE_vs_flashpe" % det.replace(" ", "-"))

    with open(os.path.join(args.plotdir, "summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
