"""Impact of the cross-section central-value weights on the GUMPLE selection.

For SBND and ICARUS, 1u1p and 1uMp, draws the selected reco-neutrino-energy
spectrum with and without the CV weights (loaddf.xsec_cv_rwgt: MvA ZExp,
CCQE correction, PiProFix) and their ratio. The prediction is the selected
neutrino MC + dirt, as the signal-box CV (beam-off is not included); only the
neutrino MC carries CV weights (cvwgt == 1 for dirt), and "without" is
glob_scale / cvwgt.

Inputs, selection, binning and POT match signalbox_systematics_gumple.py.

    python cvweight_impact_gumple.py [--df-dir DIR] [--plotbase DIR]

FV stays at module scope of __main__ so the loaddf cache keys match the
signalbox_systematics_gumple.py runs.
"""

import argparse
import glob
import os
import re
import sys
import warnings

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, _ROOT)
sys.path.insert(1, os.path.abspath(os.path.join(_HERE, "..", "gumple")))
sys.path.append(_HERE)

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

import gumple_cuts as gc
import loaddf
import signalbox_systematics_gumple as sb

DEFAULT_DF_DIR = "/Users/gputnam/Work/osc/sbn-rewgted-24/"
DEFAULT_PLOTBASE = "/Users/gputnam/Work/osc/cafpyana/plots-gumple-2026-10-13-cvwgt-impact"

FONTSIZE = sb.FONTSIZE
RATIO_RANGE = 0.15


# ---------------------------------------------------------------------------
# Preselection (module scope: part of the loaddf cache key)
# ---------------------------------------------------------------------------
def FV(df):
    return gc.slcfv_cut(df) & df.cut_contained

def FVSBND(df):
    return FV(df)

def FVICARUSRun2(df):
    return FV(df)

def FVICARUSRun4(df):
    return FV(df)


# ---------------------------------------------------------------------------
def load_sbnd(DF_DIR):
    files = [f for f in glob.glob(DF_DIR + "SBNDMCCV_*.df") if re.search(r"SBNDMCCV_\d+\.df$", f)]
    files = sorted(files, key=lambda f: int(re.search(r"SBNDMCCV_(\d+)\.df$", f).group(1)))
    assert files, "no SBNDMCCV_<i>.df in " + DF_DIR
    print("SBND CV: %i files" % len(files), flush=True)
    mc, _, pot = loaddf.loadl(files, njob=min(len(files), 6), detector="SBND",
                              preselection=FVSBND, reweight_aFF=True, pot_univ=True)
    dirt, _, dirtpot = loaddf.load(DF_DIR + "SBND_SpringLowEMC.df", detector="SBND",
                                   preselection=FVSBND, include_syst=False)
    loaddf.scale_pot(mc, pot, sb.SGOAL_POT)
    loaddf.scale_pot(dirt, dirtpot, sb.SGOAL_POT)
    return mc, dirt


def load_icarus(DF_DIR):
    parts = {"mc": [], "dirt": []}
    for run, presel, pot_goal, ncv in [(2, FVICARUSRun2, sb.IRUN2_POT, 3), (4, FVICARUSRun4, sb.IRUN4_POT, 4)]:
        det = "ICARUS Run%i" % run
        cv = [DF_DIR + "ICARUSRun%i_SpringMCOverlay_rewgt_%i.df" % (run, i) for i in range(ncv)]
        mc, _, pot = loaddf.loadl(cv, njob=len(cv), detector=det,
                                  preselection=presel, reweight_aFF=True, pot_univ=True)
        dirt, _, dirtpot = loaddf.loadl([DF_DIR + "ICARUSRun%i_Spring_Overlay_Dirt.df" % run], njob=1,
                                        detector=det, preselection=presel, include_syst=False)
        loaddf.scale_pot(mc, pot, pot_goal)
        loaddf.scale_pot(dirt, dirtpot, pot_goal)
        parts["mc"].append(mc)
        parts["dirt"].append(dirt)
    return tuple(pd.concat(parts[k]).reset_index(drop=True) for k in ["mc", "dirt"])


def combine(mc, dirt):
    """Selected prediction: neutrino MC (minus its dirt rows) + dirt."""
    if "dirt" not in mc.columns:
        mc["dirt"] = False
    dirt["dirt"] = True
    for d in [mc, dirt]:
        sb.add_selection(d)
    cols = ["selected_" + s for s in sb.PSELS] + [sb.var, "glob_scale", "cvwgt"]
    df = pd.concat([mc.loc[~mc.dirt, cols], dirt[cols]]).reset_index(drop=True)
    assert not (df.selected_1p & df.selected_mp).any()
    assert np.allclose(dirt.cvwgt, 1)
    return df


def plot(df, det, sel, outdir):
    xmin = sb.MP_XMIN if sel == "mp" else None
    b = sb.bins if xmin is None else sb.bins[sb.bins >= xmin - 1e-9]
    c = df["selected_" + sel]
    w_cv = df.loc[c, "glob_scale"]
    w_nocv = w_cv / df.loc[c, "cvwgt"]
    n_cv = np.histogram(df.loc[c, sb.var], bins=b, weights=w_cv)[0]
    n_nocv = np.histogram(df.loc[c, sb.var], bins=b, weights=w_nocv)[0]
    ratio = np.divide(n_cv, n_nocv, out=np.ones_like(n_cv), where=n_nocv > 0)
    print("%-6s %-2s: with CV %.1f, without CV %.1f, ratio %.4f" %
          (det, sel, w_cv.sum(), w_nocv.sum(), w_cv.sum() / w_nocv.sum()), flush=True)
    for lo, hi, a, n, r in zip(b[:-1], b[1:], n_cv, n_nocv, ratio):
        print("    [%.2f, %.2f] %10.1f %10.1f %8.4f" % (lo, hi, a, n, r))

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, sharex=True, figsize=(6.4, 6.4),
                                         gridspec_kw={"height_ratios": [3, 1], "hspace": 0})
    centers = (b[1:] + b[:-1]) / 2
    ax_top.hist(centers, bins=b, weights=n_nocv, histtype="step", linewidth=2, color="tab:gray",
                label="Without CV weights")
    ax_top.hist(centers, bins=b, weights=n_cv, histtype="step", linewidth=2, color="tab:blue",
                label="With CV weights")
    ax_top.set_ylim(0, 1.3 * max(n_cv.max(), n_nocv.max()))
    pot = sb.SGOAL_POT if det == "SBND" else sb.IGOAL_POT
    sb.TAG = sb.PSEL_LABEL[sel]
    sb.add_style(ax_top, "", det=det, ylabel="Events / $%g \\times 10^{20}$ POT" % (pot / 1e20),
                 legend_loc="upper right")
    ax_bot.stairs(ratio, b, baseline=None, linewidth=2, color="tab:blue")
    ax_bot.axhline(1, color="gray", linestyle="--", linewidth=1)
    # common ratio range, so the four panels compare directly
    assert np.all(np.abs(ratio - 1) < RATIO_RANGE)
    ax_bot.set_ylim(1 - RATIO_RANGE, 1 + RATIO_RANGE)
    ax_bot.set_yticks([0.9, 1, 1.1])
    ax_top.set_yticks([t for t in ax_top.get_yticks() if 0 < t <= ax_top.get_ylim()[1]])
    ax_bot.set_xlim(b[0], b[-1])
    ax_bot.set_xlabel("Reco. Neutrino Energy [GeV]", fontsize=FONTSIZE, fontweight="bold")
    ax_bot.set_ylabel("With / Without", fontsize=FONTSIZE, fontweight="bold")
    ax_bot.tick_params(axis="both", which="both", direction="in", length=6, width=1.5,
                       labelsize=FONTSIZE, top=True, right=True)
    for spine in ax_bot.spines.values():
        spine.set_linewidth(1.5)
    fig.align_ylabels([ax_top, ax_bot])

    name = "%s_%s_cvwgt_impact" % (det, sel)
    fig.savefig(os.path.join(outdir, "png", name + ".png"), bbox_inches="tight", dpi=150)
    fig.savefig(os.path.join(outdir, "pdf", name + ".pdf"), bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--df-dir", default=DEFAULT_DF_DIR)
    ap.add_argument("--plotbase", default=DEFAULT_PLOTBASE)
    args = ap.parse_args()

    DF_DIR = args.df_dir if args.df_dir.endswith("/") else args.df_dir + "/"
    for d in ["png", "pdf"]:
        os.makedirs(os.path.join(args.plotbase, d), exist_ok=True)
    print("DF_DIR:", DF_DIR)
    print("PLOTBASE:", args.plotbase)

    dfs = {"SBND": combine(*load_sbnd(DF_DIR)), "ICARUS": combine(*load_icarus(DF_DIR))}
    for det, df in dfs.items():
        for sel in sb.PSELS:
            plot(df, det, sel, args.plotbase)


if __name__ == "__main__":
    main()
