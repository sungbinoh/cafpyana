"""Signal-box systematics for the GUMPLE dfs, 1p and Mp, plus the joint 1p+Mp
SBND -> ICARUS constraint.

Script port of nb/SignalBoxSystematics-GUMPLE.ipynb (left untouched). The
notebook evaluates one proton multiplicity per run (PROTON_SEL); this script
loads the production once and then

  1. makes every notebook plot for 1p (n_pfp == 2) and for Mp (n_pfp > 2),
     into <plotbase>/signalbox/{1p,mp}/{png,pdf}/, and
  2. treats the two multiplicities together (<plotbase>/signalbox/joint/):
     the joint covariance of [SBND 1u1p, SBND 1uMp, ICARUS 1u1p, ICARUS 1uMp]
     and the ICARUS 1p / Mp uncertainties conditioned on SBND 1p AND Mp.

The joint covariance reuses syst.ConcatBins with a per-block cut list, so each
systematic's universes are shared between the 1p and Mp blocks and every
1p <-> Mp correlation comes out of the same machinery as the per-sample one.

Changes relative to the notebook (besides the multiplicity handling):
  - default production sbn-rewgted-24;
  - SBND CV is every SBNDMCCV_<i>.df in --df-dir (13 shards on -24; --nsbnd N
    keeps the first N);
  - SBND DENT is DENT2 vs Nom2 (was DENT vs Nom). The SBND detector set holds
    no SCE sample (0xSCE / 2xSCE are not used);
  - the ICARUS Run 4 trigger ratio plot uses the r4 down universe (the
    notebook paired r4_up with r2_dn).

    python signalbox_systematics_gumple.py [--df-dir DIR] [--plotbase DIR] [--nsbnd N]

Run from anywhere; FV stays at module scope of __main__ so the loaddf cache
keys match between runs.
"""

import argparse
import glob
import os
import re
import sys
import time
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
import matplotlib.cm as cm
import matplotlib.colors as mcolors

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

import gumple_cuts as gc
import loaddf
import syst

DEFAULT_DF_DIR = "/Users/gputnam/Work/osc/sbn-rewgted-24/"
DEFAULT_PLOTBASE = "/Users/gputnam/Work/osc/cafpyana/plots-gumple-2026-10-05-rewgted24-allSBNDCV"

RECO = "PANDORA"

# Proton-multiplicity selections. gc.trk_cut already requires n_pfp >= 2.
PSEL_CUTS = {
    "1p": lambda df: df.n_pfp == 2,
    "mp": lambda df: df.n_pfp > 2,
}
PSELS = list(PSEL_CUTS)
PSEL_LABEL = {"1p": "1$\\mu$1p", "mp": "1$\\mu$Mp"}

IRUN2_POT = 2e20
IRUN4_POT = 3e20
IGOAL_POT = IRUN2_POT + IRUN4_POT
SGOAL_POT = 1e20

FONTSIZE = 14

var = "nu_E_calo"
wgt = "glob_scale"
# fine binning (as the notebook, restored 2026-09-02)
bins = np.array([0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.25, 1.5])
centers = (bins[1:] + bins[:-1]) / 2
nbin = len(bins) - 1
DROP_EDGE_BINS = bins[0] == 0.0
EDGE = slice(1, -1) if DROP_EDGE_BINS else slice(None)

labels = [
    "Flux",
    "G4",
    "XSec",
    "POT Norm.",
    "Detector",
    "Dirt",
    "Beam Off",
    "Stat",
]

# plotting state set per output directory
PLOTDIR = None
TAG = ""
# lower edge of the drawn neutrino-energy axis (None: the full binning); the
# Mp plots start at 0.5 GeV, below which there are almost no Mp events
XMIN = None
MP_XMIN = 0.5


# ---------------------------------------------------------------------------
# Preselection (module scope: part of the loaddf cache key)
# ---------------------------------------------------------------------------
def FV(df):
    is_spine = "SPINE" in RECO
    ret = gc.slcfv_cut(df) & df.cut_contained
    if is_spine:
        ret = ret & (df.is_time_contained)
    return ret

def FVSBND(df):
    return FV(df)

def FVICARUSRun2(df):
    return FV(df)

def FVICARUSRun4(df):
    return FV(df)


def TrueAV(df):
    vtx = pd.DataFrame({'detector': df.detector, 'Run': df.Run,
                        'x': df.true_vtx_x, 'y': df.true_vtx_y, 'z': df.true_vtx_z}, index=df.index)
    return gc._fv_cut(vtx, 0, 0, 0, 0)

def OOAVSBND(df):
    return ~np.isnan(df.true_vtx_x) & ~TrueAV(df)

def OOAVICARUS(df):
    return ~np.isnan(df.true_vtx_x) & ~TrueAV(df)


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------
def add_style(ax, xlabel, title=None, det="ICARUS", ylabel='Events / $10^{20}$ POT', legend_loc=None,
              legend_ncol=1, legend_title=None, legend_fontsize=FONTSIZE-1):
    if title is None:
        title = TAG
    ax.tick_params(axis='both', which='both', direction='in', length=6, width=1.5, labelsize=FONTSIZE, top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    ax.set_xlabel(xlabel, fontsize=FONTSIZE, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=FONTSIZE, fontweight='bold')
    ax.set_title(f"$\\bf{{{det}}}$  {title}", fontsize=FONTSIZE+2)
    handles, labs = ax.get_legend_handles_labels()
    ax.legend(handles, [abbrev(l) for l in labs], fontsize=legend_fontsize, loc=legend_loc, ncol=legend_ncol,
              title=legend_title, title_fontsize=FONTSIZE)


# legend-only abbreviations (the full names still build the file names), so
# the wide 3-column ICARUS breakdown legend fits inside the axes
LEGEND_ABBREV = [("Constant", "Const."), ("East", "E"), ("West", "W"), ("Cathode", "Cath.")]


def abbrev(label):
    for a, b in LEGEND_ABBREV:
        label = label.replace(a, b)
    return label


def B(b):
    """The drawn part of an energy binning: the bins above XMIN (the Mp plots
    exclude the low bins from the histogram rather than cutting the x axis)."""
    b = np.asarray(b)
    return b if XMIN is None else b[b >= XMIN - 1e-9]


LEGEND_PAD = 0.03  # clearance between the highest line under the legend and the legend (axes fraction)


def _step_segments(ax):
    """(x0, x1, y) of every horizontal segment drawn in `ax` (step histograms, lines)."""
    xys = [p.get_xy() for p in ax.patches if hasattr(p, "get_xy")]
    xys += [l.get_xydata() for l in ax.lines]
    segs = []
    for xy in xys:
        xy = np.asarray(xy, dtype=float)
        if len(xy) < 2:
            continue
        x0, x1, y0, y1 = xy[:-1, 0], xy[1:, 0], xy[:-1, 1], xy[1:, 1]
        m = np.isfinite(x0) & np.isfinite(x1) & np.isfinite(y0) & np.isfinite(y1)
        # horizontal segments carry the bin contents; vertical ones join them
        h = m & (x1 != x0)
        segs.append(np.c_[np.minimum(x0, x1)[h], np.maximum(x0, x1)[h], np.maximum(y0, y1)[h]])
    return np.concatenate(segs) if segs else np.zeros((0, 3))


def _legend_ylim(ax):
    """(bottom, top) that puts every drawn line below the legend of `ax`.
    The bottom is kept (0 on linear axes); None if nothing to fit."""
    leg = ax.get_legend()
    segs = _step_segments(ax)
    log = ax.get_yscale() == "log"
    if log:
        segs = segs[segs[:, 2] > 0]
    if segs.size == 0:
        return None
    bot = ax.get_ylim()[0]
    if not log:
        bot = min(bot, 0) if np.all(segs[:, 2] >= 0) else bot
    f = (np.log10 if log else (lambda v: v))
    fb = f(bot)
    ymax_all = f(segs[:, 2].max())
    top = fb + 1.05*(ymax_all - fb)
    if leg is not None:
        fig = ax.figure
        fig.canvas.draw()
        bb = leg.get_window_extent().transformed(ax.transAxes.inverted())
        xl = ax.transData.inverted().transform(ax.transAxes.transform([[bb.x0, 0], [bb.x1, 0]]))[:, 0]
        under = segs[(segs[:, 1] > xl.min()) & (segs[:, 0] < xl.max())]
        if under.size and bb.y0 > LEGEND_PAD + 0.1:
            top = max(top, fb + (f(under[:, 2].max()) - fb)/(bb.y0 - LEGEND_PAD))
    return bot, (10**top if log else top)


def _fit_legend_ylim(fig):
    for ax in fig.axes:
        if ax.get_legend() is None:
            continue
        # twice: a loc="best" legend can move once the range changes
        for _ in range(2):
            r = _legend_ylim(ax)
            if r is not None:
                ax.set_ylim(*r)


# figures held back so that equivalent SBND / ICARUS plots share a y range
_PENDING = {}


def savefig(name, dirs=None, group=None, **kw):
    fig = plt.gcf()
    if group is not None:
        _PENDING.setdefault(group, []).append((fig, name, dirs, kw))
        return
    _fit_legend_ylim(fig)
    _write(fig, name, dirs, kw)


def flush(group):
    """Save the figures of `group` with the widest y range any of them needs."""
    items = _PENDING.pop(group, [])
    for fig, *_ in items:
        _fit_legend_ylim(fig)
    lims = [fig.axes[0].get_ylim() for fig, *_ in items]
    if lims:
        lo, hi = min(l[0] for l in lims), max(l[1] for l in lims)
        for fig, *_ in items:
            fig.axes[0].set_ylim(lo, hi)
    for fig, name, dirs, kw in items:
        _write(fig, name, dirs, kw)


def _write(fig, name, dirs, kw):
    for d in (dirs or [PLOTDIR]):
        fig.savefig(os.path.join(d, "png", name + ".png"), **kw)
        fig.savefig(os.path.join(d, "pdf", name + ".pdf"), **kw)
    plt.close(fig)


def clean(name):
    return name.translate(dict.fromkeys(map(ord, "$\\_{} ="), None))


def hist_dropedges(centers, bins, weights, **kw):
    return plt.hist(np.asarray(centers)[EDGE], bins=B(np.asarray(bins)[EDGE]),
                    weights=np.asarray(weights)[EDGE], **kw)


def ratioplot(detvars, names, title, cut, presel=None, scales=None, ratio_range=0.1):
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, sharex=True,
                                         gridspec_kw={"height_ratios": [3, 1], "hspace": 0})
    cv_vals = None
    for i, d in enumerate(detvars):
        c = d[cut]
        if presel is not None:
            c = c & presel(d)
        n, _ = np.histogram(d.loc[c, var], bins=bins, weights=d.loc[c, wgt])
        n = n * (scales[i] if scales is not None else 1.)
        ax_top.hist(centers[EDGE], bins=B(bins[EDGE]), weights=n[EDGE], color="tab:gray" if i == 0 else None,
                    histtype="step", linewidth=2, label=names[i])
        if i == 0:
            cv_vals = n
        else:
            ratio = np.divide(n, cv_vals, out=np.ones_like(n, dtype=float), where=cv_vals > 0)
            ax_bot.hist(centers[EDGE], bins=B(bins[EDGE]), weights=ratio[EDGE], histtype="step", linewidth=2)

    ax_bot.axhline(1, color="gray", linestyle="--", linewidth=1)
    add_style(ax_top, "", det=title)
    ax_top.set_xlabel("")
    ax_bot.set_xlabel("Reco. Neutrino Energy [GeV]", fontsize=FONTSIZE, fontweight="bold")
    ax_bot.set_ylabel("Var. / Nom.", fontsize=FONTSIZE, fontweight="bold")
    ax_bot.tick_params(axis="both", which="both", direction="in", length=6, width=1.5, labelsize=FONTSIZE, top=True, right=True)
    ax_bot.set_ylim([1-ratio_range/2., 1+ratio_range/2.])
    for spine in ax_bot.spines.values():
        spine.set_linewidth(1.5)
    return ax_top, ax_bot


def ratioplot_varyselect(df, selects, names, title, cut):
    dfs = []
    for s in selects:
        d = df[[var, wgt]].copy()
        d[cut] = df[s]
        dfs.append(d)
    return ratioplot(dfs, names, title, cut)


def corr_f(cov):
    err = np.sqrt(np.diag(cov))
    err[err == 0] = 1
    err_inv = np.diag(1/err)
    return err_inv@cov@err_inv


def plot_corr_blocks(C, names, title, fname, dirs=None):
    """Block correlation matrix, one energy-binned panel per block pair (the
    notebook's 2x2 SBND/ICARUS layout, generalized to any number of blocks).
    Block 0 sits bottom-left; rows go up, columns go right."""
    nb = len(names)
    size = (6.4, 4.8) if nb == 2 else (10, 8.5)
    fig, axs = plt.subplots(nb, nb, layout='constrained', figsize=size,
                            gridspec_kw={'wspace': 0.0, 'hspace': 0.0})
    fig.get_layout_engine().set(w_pad=0 / 72, h_pad=0 / 72, hspace=0, wspace=0)
    ticks = np.linspace(bins[0], bins[-1], 7)

    def trim(k):
        # drop the tick that would collide with the neighbouring panel
        if k == 0:
            return ticks[:-1]
        if k == nb - 1:
            return ticks[1:]
        return ticks[1:-1]

    for i in range(nb):
        for j in range(nb):
            a = axs[i][j]
            ic = nb - 1 - i
            a.pcolormesh(bins, bins, C[nbin*ic:nbin*(ic+1), nbin*j:nbin*(j+1)], cmap="bwr",
                         vmin=-1, vmax=1, linewidth=0, rasterized=True)
            a.tick_params(axis='both', which='major', labelsize=FONTSIZE-2 if nb == 2 else FONTSIZE-4)
            if nb > 2:
                a.axhline(bins[0], color="black", linewidth=0.5)
                a.axvline(bins[0], color="black", linewidth=0.5)
            if j == 0:
                a.set_yticks(trim(ic))
                a.set_ylabel(names[ic], size=FONTSIZE-2)
            else:
                a.set_yticks([])
            if i == nb - 1:
                a.set_xticks(trim(j))
                a.set_xlabel(names[j], size=FONTSIZE-2)
            else:
                a.set_xticks([])

    fig.suptitle(title, size=FONTSIZE+2)
    fig.supxlabel("Reco. Neutrino Energy [GeV]", size=FONTSIZE)
    fig.supylabel("Reco. Neutrino Energy [GeV]", size=FONTSIZE)
    sm = cm.ScalarMappable(cmap="bwr", norm=mcolors.Normalize(vmin=-1, vmax=1))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=[a for ax in axs for a in ax])
    cbar.set_label("Correlation", rotation=270, size=FONTSIZE-2, labelpad=12)
    savefig(fname, dirs=dirs)


def ratio_cov_full(x, y, cov):
    """Covariance of r = x / y given the full covariance of (x, y)."""
    n = len(x)
    assert cov.shape == (2*n, 2*n)
    y_safe = np.where(np.abs(y) < 1e-12, 1e-12, y)
    J = np.hstack([np.diag(1.0 / y_safe), np.diag(-x / y_safe**2)])
    return J @ cov @ J.T


def cond(cov, nx, what):
    c, _, k = syst.conditional_constraint(cov, nx)
    print("  conditional constraint %-28s cond(Cxx) = %.3g" % (what, k))
    return c


def makedirs(d):
    os.makedirs(os.path.join(d, "png"), exist_ok=True)
    os.makedirs(os.path.join(d, "pdf"), exist_ok=True)


# ---------------------------------------------------------------------------
# Selection columns
# ---------------------------------------------------------------------------
def signalbox_base(df):
    # the notebook's signalbox() minus the proton-multiplicity cut; the flash
    # cut is applied separately so the trigger systematic can move it
    return gc.cosmic_cut(df) & gc.pid_cut(df) & gc.trk_cut(df) & df.cut_cathode


def add_selection(df):
    """selected_noflash_<sel> / selected_<sel> for every multiplicity."""
    if len(df) == 0:  # gumple cuts are not safe on zero-row frames
        for s in PSELS:
            df["selected_noflash_" + s] = False
            df["selected_" + s] = False
        return df
    base = signalbox_base(df)
    fl = gc.flash_cut(df)
    for s, f in PSEL_CUTS.items():
        df["selected_noflash_" + s] = base & f(df)
        df["selected_" + s] = df["selected_noflash_" + s] & fl
    return df


# TRIGGER SYSTEMATIC: flash-PE scale-factor uncertainty (notebook cell 46). The
# uncertainty is |nominal s - GUMPLE Contained-stage fit s|, from
# plots-gumple-2026-10-20-rewgted25/flashpe_fit/gumple (as in mcdata_comparison_gumple.py)
TRIG_PE_SCALE     = {1: 0.642, 2: 0.632, 4: 0.358}
TRIG_PE_SCALE_UNC = {1: 0.024, 2: 0.155, 4: 0.065}
TRIG_THRESHOLD    = {1: 2000, 2: 5000., 4: 1000.}

def trig_flash_cut_var(df, updn, runs=None):
    sel = pd.Series(False, index=df.index)
    for run, T in TRIG_THRESHOLD.items():
        f = TRIG_PE_SCALE_UNC[run] / TRIG_PE_SCALE[run]
        Tvar = T / (1 - updn*f) if (runs is None or run in runs) else T
        sel = sel | ((df.Run == run) & (df.flash_maxpe > Tvar))
    return sel


def tcol(name, sel):
    """Trigger-universe cut column for `sel`; a per-block list for the joint fit."""
    if sel == "joint":
        return [name + "_" + p for p in PSELS]
    return name + "_" + sel


# ---------------------------------------------------------------------------
# Calorimetric variations
# ---------------------------------------------------------------------------
chi2vars = ["mu_chi2_of_mu_cand", "prot_chi2_of_mu_cand", "mu_chi2_of_prot_cand", "prot_chi2_of_prot_cand"]
chi2labels = ["$\\chi^2_\\mu$ of Muon Candidate", "$\\chi^2_p$ of Muon Candidate",
              "$\\chi^2_\\mu$ of Proton Candidate", "$\\chi^2_p$ of Proton Candidate"]
chi2bins = [np.linspace(0, 60, 21), np.linspace(0, 350, 61)[:30],
            np.linspace(0, 60, 21), np.linspace(0, 350, 61)[:30]]

SCHI2_VARIATIONS      = ["smear13", "hi", "alpha_p", "beta_p", "R_p"]
SCHI2_VARIATION_NAMES = ["dE/dx Constant Res.", "Gain Hi", "EMB $\\alpha$", "EMB $\\beta$", "EMB $R$"]
ICHI2_VARIATIONS      = ["smear13", "hi", "dedxbias"]
ICHI2_VARIATION_NAMES = ["dE/dx Constant Res.", "Gain Hi", "dE/dx Bias"]
chi2_variations = list(dict.fromkeys(SCHI2_VARIATIONS + ICHI2_VARIATIONS))

SCHI2_DETVAR_NAMES = ["dE/dx Constant Res.", "Calo. Gain", "EMB $\\alpha$", "EMB $\\beta$", "EMB $R$"]
ICHI2_DETVAR_NAMES = ["dE/dx Constant Res.", "Calo. Gain", "dE/dx Bias"]
S_TRIGGER_NAMES = ["Trigger"]
I_TRIGGER_NAMES = ["Trigger Run 2", "Trigger Run 4"]

def v_calovariation(df, variation):
    df = df[[c for c in df.columns if "univ" not in c]].copy()
    for c in chi2vars:
        df[c] = df[c.replace("chi2", "chi2%s" % variation)]
    return df


def select_loose(df):
    # loose selection for the area-normalized chi2 shape plots
    return gc.cosmic_cut(df) & gc.pid_cut(df)


# ---------------------------------------------------------------------------
def main():
    global PLOTDIR, TAG, XMIN

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--df-dir", default=DEFAULT_DF_DIR)
    ap.add_argument("--plotbase", default=DEFAULT_PLOTBASE)
    ap.add_argument("--nsbnd", type=int, default=None,
                    help="use only the first N SBNDMCCV_<i>.df shards (default: all)")
    args = ap.parse_args()

    DF_DIR = args.df_dir if args.df_dir.endswith("/") else args.df_dir + "/"
    OUTBASE = os.path.join(args.plotbase, "signalbox")
    SELDIRS = {s: os.path.join(OUTBASE, s) for s in PSELS}
    JOINTDIR = os.path.join(OUTBASE, "joint")
    for d in list(SELDIRS.values()) + [JOINTDIR]:
        makedirs(d)
    print("DF_DIR:", DF_DIR)
    print("PLOTBASE:", args.plotbase)
    t0 = time.time()

    # --- inputs ---------------------------------------------------------------
    SCV_FILES = [f for f in glob.glob(DF_DIR + "SBNDMCCV_*.df") if re.search(r"SBNDMCCV_\d+\.df$", f)]
    SCV_FILES = sorted(SCV_FILES, key=lambda f: int(re.search(r"SBNDMCCV_(\d+)\.df$", f).group(1)))[:args.nsbnd]
    assert SCV_FILES, "no SBNDMCCV_<i>.df in " + DF_DIR
    print("SBND CV: %i files" % len(SCV_FILES), flush=True)
    SDIRT_FILE = DF_DIR + "SBND_SpringLowEMC.df"
    SBEAMOFF_FILE = DF_DIR + "SBND_SpringBNBOffData.df"
    # NB no SBND SCE sample in the detector set (0xSCE / 2xSCE deliberately unused)
    SDETVARS = [
        ([DF_DIR + "SBND_SpringMC_Nom2.df"], [DF_DIR + "SBND_SpringMC_DENT2.df"]),
        ([DF_DIR + "SBND_SpringMC_WMNom.df"], [DF_DIR + "SBND_SpringMC_WMXThetaXW.df"]),
        ([DF_DIR + "SBND_SpringMC_WMNom.df"], [DF_DIR + "SBND_SpringMC_WMYZ.df"]),
    ]
    SDETVAR_NAMES = ["DENT", "WM $X\\theta_{xw}$", "WM $YZ$"]

    IRUN2_CV_FILES = [DF_DIR + "ICARUSRun2_SpringMCOverlay_rewgt_%i.df" % i for i in range(3)]
    IRUN4_CV_FILES = [DF_DIR + "ICARUSRun4_SpringMCOverlay_rewgt_%i.df" % i for i in range(4)]
    IRUN2_DIRT_FILES = [DF_DIR + "ICARUSRun2_Spring_Overlay_Dirt.df"]
    IRUN4_DIRT_FILES = [DF_DIR + "ICARUSRun4_Spring_Overlay_Dirt.df"]
    IRUN2_BEAMOFF_FILE = DF_DIR + "ICARUS_SpringRun2BNBOff_unblind.df"
    IRUN4_BEAMOFF_FILE = DF_DIR + "ICARUS_SpringRun4BNBOff_ReCAF2026.df"
    # no ICARUS *_Overlay_Nom.df: the detvar CV partner is the rewgt CV sample
    IRUN2_DETVARS = [(IRUN2_CV_FILES, [DF_DIR + "ICARUSRun2_Spring_Overlay_%s.df" % v]) for v in ["SCE", "WMXThXW", "WMYZ"]]
    IRUN4_DETVARS = [(IRUN4_CV_FILES, [DF_DIR + "ICARUSRun4_Spring_Overlay_%s.df" % v]) for v in ["SCE", "WMXThXW", "WMYZ"]]
    IDETVAR_NAMES = ["SCE", "WM $X\\theta_{xw}$", "WM $YZ$"]

    for f in (SCV_FILES + [SDIRT_FILE, SBEAMOFF_FILE] + [x for p in SDETVARS for l in p for x in l] +
              IRUN2_CV_FILES + IRUN4_CV_FILES + IRUN2_DIRT_FILES + IRUN4_DIRT_FILES +
              [IRUN2_BEAMOFF_FILE, IRUN4_BEAMOFF_FILE] + [x for p in IRUN2_DETVARS + IRUN4_DETVARS for x in p[1]]):
        assert os.path.exists(f), f

    # --- load -----------------------------------------------------------------
    Idf_r2, _, Ipot_r2 = loaddf.loadl(IRUN2_CV_FILES, njob=len(IRUN2_CV_FILES), detector="ICARUS Run2",
                                      preselection=FVICARUSRun2, reweight_aFF=True, pot_univ=True)
    Idf_r4, _, Ipot_r4 = loaddf.loadl(IRUN4_CV_FILES, njob=len(IRUN4_CV_FILES), detector="ICARUS Run4",
                                      preselection=FVICARUSRun4, reweight_aFF=True, pot_univ=True)
    loaddf.scale_pot(Idf_r2, Ipot_r2, IRUN2_POT)
    loaddf.scale_pot(Idf_r4, Ipot_r4, IRUN4_POT)
    Idf = pd.concat([Idf_r2, Idf_r4]).reset_index(drop=True)
    del Idf_r2, Idf_r4

    Idirt_r2, _, Idirtpot_r2 = loaddf.loadl(IRUN2_DIRT_FILES, njob=1, detector="ICARUS Run2",
                                            preselection=FVICARUSRun2, include_syst=False)
    Idirt_r4, _, Idirtpot_r4 = loaddf.loadl(IRUN4_DIRT_FILES, njob=1, detector="ICARUS Run4",
                                            preselection=FVICARUSRun4, include_syst=False)
    loaddf.scale_pot(Idirt_r2, Idirtpot_r2, IRUN2_POT)
    loaddf.scale_pot(Idirt_r4, Idirtpot_r4, IRUN4_POT)
    Idirt = pd.concat([Idirt_r2, Idirt_r4]).reset_index(drop=True)

    Ioffbeam_r2, _, Ioffbeampot_r2 = loaddf.load(IRUN2_BEAMOFF_FILE, detector="ICARUS Run2", offbeampot=True,
                                                 preselection=FVICARUSRun2, include_syst=False, load_truth=False, match_Enu=False)
    Ioffbeam_r4, _, Ioffbeampot_r4 = loaddf.load(IRUN4_BEAMOFF_FILE, detector="ICARUS Run4", offbeampot=True,
                                                 preselection=FVICARUSRun4, include_syst=False, load_truth=False, match_Enu=False)
    loaddf.scale_pot(Ioffbeam_r2, Ioffbeampot_r2, IRUN2_POT)
    loaddf.scale_pot(Ioffbeam_r4, Ioffbeampot_r4, IRUN4_POT)
    Ioffbeam = pd.concat([Ioffbeam_r2, Ioffbeam_r4]).reset_index(drop=True)
    print("ICARUS loaded (%.0f s)" % (time.time() - t0), flush=True)

    Sdf, _, Spot = loaddf.loadl(SCV_FILES, njob=min(len(SCV_FILES), 6), detector="SBND",
                                preselection=FVSBND, reweight_aFF=True, pot_univ=True)
    Sdirt, _, Sdirtpot = loaddf.load(SDIRT_FILE, detector="SBND", preselection=FVSBND, include_syst=False)
    Soffbeam, _, Soffbeampot = loaddf.load(SBEAMOFF_FILE, detector="SBND", preselection=FVSBND, offbeampot=True,
                                           include_syst=False, load_truth=False, match_Enu=False)
    print("SBND loaded (%.0f s)" % (time.time() - t0), flush=True)

    def load_detvar(detvar, detector):
        presel = {"SBND": FVSBND, "ICARUS Run2": FVICARUSRun2, "ICARUS Run4": FVICARUSRun4}[detector]
        df, match, pots = zip(*[loaddf.loadl(flist, preselection=presel, include_syst=False, detector=detector)
                                for flist in detvar])
        df, pots = loaddf.match_common_evts(match, df, pots)
        return df, pots

    Sdetvars, Sdetvar_pots = zip(*[load_detvar(d, "SBND") for d in SDETVARS])
    Ir2detvars, Ir2detvar_pots = zip(*[load_detvar(d, "ICARUS Run2") for d in IRUN2_DETVARS])
    Ir4detvars, Ir4detvar_pots = zip(*[load_detvar(d, "ICARUS Run4") for d in IRUN4_DETVARS])
    print("detvars loaded (%.0f s)" % (time.time() - t0), flush=True)

    for c in Sdf.columns:
        if "_univ" in c:
            Sdirt[c] = 1
    for c in Idf.columns:
        if "_univ" in c:
            Idirt[c] = 1
    if "dirt" not in Sdf.columns:
        Sdf["dirt"] = False
        Sdirt["dirt"] = True
    if "dirt" not in Idf.columns:
        Idf["dirt"] = False
        Idirt["dirt"] = True

    for i in range(len(Ir2detvars)):
        for j in range(len(Ir2detvars[i])):
            loaddf.scale_pot(Ir2detvars[i][j], Ir2detvar_pots[i][j], IRUN2_POT)
    for i in range(len(Ir4detvars)):
        for j in range(len(Ir4detvars[i])):
            loaddf.scale_pot(Ir4detvars[i][j], Ir4detvar_pots[i][j], IRUN4_POT)
    Idetvars = [[pd.concat([Ir2detvars[i][j], Ir4detvars[i][j]]).reset_index(drop=True)
                 for j in range(len(Ir4detvars[i]))] for i in range(len(Ir4detvars))]
    del Ir2detvars, Ir4detvars

    loaddf.scale_pot(Sdf, Spot, SGOAL_POT)
    loaddf.scale_pot(Sdirt, Sdirtpot, SGOAL_POT)
    loaddf.scale_pot(Soffbeam, Soffbeampot, SGOAL_POT)
    for i in range(len(Sdetvars)):
        for j in range(len(Sdetvars[i])):
            loaddf.scale_pot(Sdetvars[i][j], Sdetvar_pots[i][j], SGOAL_POT)

    Sdf = pd.concat([Sdf[~Sdf.dirt], Sdirt])
    Idf = pd.concat([Idf[~Idf.dirt], Idirt])

    # dirt carries no chi2-variation columns: fall back to the CV chi2
    for c in chi2vars:
        for chi2_var in chi2_variations:
            vc = c.replace("chi2", "chi2" + chi2_var)
            Sdirt[vc] = Sdirt[c]
            Idirt[vc] = Idirt[c]
            Sdf[vc] = Sdf[vc].fillna(Sdf[c]) if vc in Sdf.columns else Sdf[c]
            Idf[vc] = Idf[vc].fillna(Idf[c]) if vc in Idf.columns else Idf[c]

    COMMON_DIRS = list(SELDIRS.values())

    # --- selection-independent chi2 shape plots (written to both sel dirs) -----
    TAG = ""
    Ssel = select_loose(Sdf)
    for c, l, b in zip(chi2vars, chi2labels, chi2bins):
        plt.figure()
        plt.hist(Sdf.loc[Ssel, c], bins=b, density=True, histtype="step", linewidth=2, label="CV")
        for chi2_var, label in zip(SCHI2_VARIATIONS, SCHI2_VARIATION_NAMES):
            plt.hist(Sdf.loc[Ssel, c.replace("chi2", "chi2" + chi2_var)], bins=b, density=True, histtype="step", linewidth=2, label=label)
        add_style(plt.gca(), l, det="SBND", ylabel="Area Normalized")
        savefig("SBND_%s_variations" % c, dirs=COMMON_DIRS, bbox_inches="tight")

        plt.figure()
        plt.hist(Idf[c], bins=b, density=True, histtype="step", linewidth=2, label="CV")
        for chi2_var, label in zip(ICHI2_VARIATIONS, ICHI2_VARIATION_NAMES):
            plt.hist(Idf[c.replace("chi2", "chi2" + chi2_var)], bins=b, density=True, histtype="step", linewidth=2, label=label)
        add_style(plt.gca(), l, det="ICARUS", ylabel="Area Normalized")
        savefig("ICARUS_%s_variations" % c, dirs=COMMON_DIRS, bbox_inches="tight")

        for j, varname in enumerate(SDETVAR_NAMES):
            dv = Sdetvars[j]
            plt.figure()
            plt.hist(dv[0].loc[select_loose(dv[0]), c], bins=b, density=True, histtype="step", linewidth=2, label="CV")
            for d in dv[1:]:
                plt.hist(d.loc[select_loose(d), c], bins=b, density=True, histtype="step", linewidth=2, label=varname)
            add_style(plt.gca(), l, det="SBND", ylabel="Area Normalized")
            savefig("SBND_detvar_%s_variation_%s" % (c, clean(varname)), dirs=COMMON_DIRS, bbox_inches="tight")

            plt.figure()
            plt.hist(dv[0][c], bins=b, histtype="step", linewidth=2, label="CV", weights=dv[0].glob_scale)
            for d in dv[1:]:
                plt.hist(d[c], bins=b, histtype="step", linewidth=2, label=varname, weights=d.glob_scale)
            add_style(plt.gca(), l, det="SBND", ylabel="Events / 1e20 POT")
            savefig("SBND_detvar_absnorm_%s_variation_%s" % (c, clean(varname)), dirs=COMMON_DIRS, bbox_inches="tight")

        for j, varname in enumerate(IDETVAR_NAMES):
            dv = Idetvars[j]
            plt.figure()
            plt.hist(dv[0][c], bins=b, density=True, histtype="step", linewidth=2, label="CV")
            for d in dv[1:]:
                plt.hist(d[c], bins=b, density=True, histtype="step", linewidth=2, label=varname)
            add_style(plt.gca(), l, det="ICARUS", ylabel="Area Normalized")
            savefig("ICARUS_detvar_%s_variation_%s" % (c, clean(varname)), dirs=COMMON_DIRS, bbox_inches="tight")
    del Ssel

    # --- universes / selection columns ----------------------------------------
    Schi2_detvars = [v_calovariation(Sdf, V) for V in SCHI2_VARIATIONS]
    Ichi2_detvars = [v_calovariation(Idf, V) for V in ICHI2_VARIATIONS]

    assert np.all(np.isnan(Sdirt[~OOAVSBND(Sdirt)].true_vtx_x))
    assert len(Idirt) == 0 or np.all(np.isnan(Idirt[~OOAVICARUS(Idirt)].true_vtx_x))

    for d in ([Sdf, Sdirt, Soffbeam, Idf, Idirt, Ioffbeam] + Schi2_detvars + Ichi2_detvars +
              [x for dv in Sdetvars for x in dv] + [x for dv in Idetvars for x in dv]):
        add_selection(d)

    for d in [Sdf, Sdirt, Idf, Idirt]:
        d["true_E"] = d.true_E.fillna(-1)
    for d in [Soffbeam, Ioffbeam]:
        d["true_E"] = -1
    for d in Schi2_detvars + Ichi2_detvars + [x for dv in Sdetvars for x in dv] + [x for dv in Idetvars for x in dv]:
        d["true_E"] = d.true_E.fillna(-1)

    # 1p and Mp partition the multiplicity-inclusive selection
    for d in [Sdf, Idf]:
        assert not (d.selected_1p & d.selected_mp).any()

    # trigger (flash-PE scale) universes, per multiplicity
    for s in PSELS:
        nf = "selected_noflash_" + s
        Sdf["selected_trig_up_" + s] = Sdf[nf] & trig_flash_cut_var(Sdf, +1)
        Sdf["selected_trig_dn_" + s] = Sdf[nf] & trig_flash_cut_var(Sdf, -1)
        Idf["selected_trig_r2_up_" + s] = Idf[nf] & trig_flash_cut_var(Idf, +1, runs=[2])
        Idf["selected_trig_r2_dn_" + s] = Idf[nf] & trig_flash_cut_var(Idf, -1, runs=[2])
        Idf["selected_trig_r4_up_" + s] = Idf[nf] & trig_flash_cut_var(Idf, +1, runs=[4])
        Idf["selected_trig_r4_dn_" + s] = Idf[nf] & trig_flash_cut_var(Idf, -1, runs=[4])
        sel = "selected_" + s
        assert not (Sdf["selected_trig_up_" + s] & ~Sdf[sel]).any()
        assert not (Sdf[sel] & ~Sdf["selected_trig_dn_" + s]).any()
        assert not ((Idf["selected_trig_r2_up_" + s] | Idf["selected_trig_r4_up_" + s]) & ~Idf[sel]).any()
        assert not (Idf[sel] & ~(Idf["selected_trig_r2_dn_" + s] & Idf["selected_trig_r4_dn_" + s])).any()

    # TRACK-SPLITTING SYSTEMATIC (ICARUS), notebook cell 48
    fv = gc.ICARUSRun4FVCuts
    TRACKSPLIT_PLANES = [
        ("Z=0",          "z", 0.0),
        ("East Cathode", "x", 0.5*(fv["C0"]["x"]["min"] + fv["C0"]["x"]["max"])),
        ("West Cathode", "x", 0.5*(fv["C1"]["x"]["min"] + fv["C1"]["x"]["max"])),
    ]
    TRACKSPLIT_FRAC = {"Z=0": 0.1143, "East Cathode": 0.0446, "West Cathode": 0.0314}
    Isplit_systs = []
    I_TRACKSPLIT_NAMES = []
    for name, dim, coord in TRACKSPLIT_PLANES:
        sdf, crosses = syst.split_tracks(Idf, dim, coord, runs=[2, 4])
        add_selection(sdf)
        Isplit_systs.append(syst.TrackSplittingSystematic(Idf, sdf, crosses, TRACKSPLIT_FRAC[name]))
        I_TRACKSPLIT_NAMES.append("Trk. Split %s" % name)

    # BINDING-ENERGY SYSTEMATIC, notebook cell 50
    BE_SHIFT = 0.025
    BE_COLS = ["selected_" + s for s in PSELS] + [
        "glob_scale", "genie_mode", "mu_E", "p_E",
        "mu_dir_x", "mu_dir_y", "mu_dir_z", "p_dir_x", "p_dir_y", "p_dir_z",
        "psum_E", "psum_p", "psum_dir_x", "psum_dir_y", "psum_dir_z", "n_pfp",
        "nu_E_calo", "del_p", "del_Tp", "del_phi"]
    Sbedf = syst.shift_binding_energy(Sdf[BE_COLS], BE_SHIFT, fraction=0.5)
    Ibedf = syst.shift_binding_energy(Idf[BE_COLS], BE_SHIFT, fraction=0.5)
    S_BE_syst = syst.SampleSystematic(Sbedf)
    I_BE_syst = syst.SampleSystematic(Ibedf)
    SXSecWgtSyst = loaddf.XSecSystematic(Sdf)
    IXSecWgtSyst = loaddf.XSecSystematic(Idf)

    S_cosmic = Sdf[Sdf.true_iscosmic]
    S_ooav = Sdf[OOAVSBND(Sdf)]
    I_ooav = Idf[OOAVICARUS(Idf)]
    print("selection + universes ready (%.0f s)" % (time.time() - t0), flush=True)

    def build_systematics(sel):
        """Notebook cell 51, with the trigger universes for `sel` ("1p", "mp"
        or "joint" -> per-block cut lists)."""
        Ssystematics = [
            loaddf.FluxSystematic(Sdf),
            loaddf.G4Systematic(Sdf),
            syst.SystematicList([SXSecWgtSyst, S_BE_syst]),
            syst.NormalizationSystematic(0.005),
            syst.SystematicList(
                [syst.SampleSystematic(dv[1:], cvdf=dv[0]) for dv in Sdetvars] +
                [syst.SampleSystematic(d, norm=(2. if n == "smear13" else 1.)) for d, n in zip(Schi2_detvars, SCHI2_VARIATIONS)] +
                [syst.SelectionSystematic(Sdf, [tcol("selected_trig_up", sel), tcol("selected_trig_dn", sel)])] +
                [syst.SystSampleSystematic(S_cosmic, norm=0.107)] +
                [syst.StatSampleSystematic(Sdf)]
            ),
            syst.SystSampleSystematic(S_ooav),
            syst.StatSampleSystematic(Soffbeam),
        ]
        Isystematics = [
            loaddf.FluxSystematic(Idf),
            loaddf.G4Systematic(Idf),
            syst.SystematicList([IXSecWgtSyst, I_BE_syst]),
            syst.NormalizationSystematic(0.005),
            syst.SystematicList(
                [syst.SampleSystematic(dv[1:], cvdf=dv[0]) for dv in Idetvars] +
                [syst.SampleSystematic(d) for d in Ichi2_detvars] +
                [syst.SelectionSystematic(Idf, [tcol("selected_trig_r2_up", sel), tcol("selected_trig_r2_dn", sel)]),
                 syst.SelectionSystematic(Idf, [tcol("selected_trig_r4_up", sel), tcol("selected_trig_r4_dn", sel)])] +
                Isplit_systs +
                [syst.StatSampleSystematic(Idf)]
            ),
            syst.SystSampleSystematic(I_ooav),
            syst.StatSampleSystematic(Ioffbeam, norm=0.5),  # TODO: change after unblinding (2x stats)
        ]
        systematics = [
            syst.CorrelatedSystematic(Ssystematics[0], Isystematics[0]),
            syst.CorrelatedSystematic(Ssystematics[1], Isystematics[1]),
            syst.SystematicList([
                syst.CorrelatedSystematic(SXSecWgtSyst, IXSecWgtSyst),
                syst.CorrelatedSystematic(S_BE_syst, I_BE_syst),
            ]),
            # POT norm: SBND Run1 correlated with ICARUS Run4, ICARUS Run2 uncorrelated
            syst.SystematicList([
                syst.UnCorrelatedSystematic(syst.NormalizationSystematic(0), syst.NormalizationSystematic(0.005*(2./5.))),
                syst.CorrelatedSystematic(syst.NormalizationSystematic(0.005), syst.NormalizationSystematic(0.005*(3./5.))),
            ]),
            syst.UnCorrelatedSystematic(Ssystematics[4], Isystematics[4]),
            syst.UnCorrelatedSystematic(Ssystematics[5], Isystematics[5]),
            syst.UnCorrelatedSystematic(Ssystematics[6], Isystematics[6]),
        ]
        return Ssystematics, Isystematics, systematics

    # --- per-multiplicity plots -------------------------------------------------
    results = {}
    for s in PSELS:
        PLOTDIR = SELDIRS[s]
        TAG = PSEL_LABEL[s]
        XMIN = MP_XMIN if s == "mp" else None
        cut = "selected_" + s
        print("\n==== %s -> %s" % (s, PLOTDIR), flush=True)
        results[s] = run_selection(s, cut, locals_=dict(
            Sdf=Sdf, Idf=Idf, Soffbeam=Soffbeam, Ioffbeam=Ioffbeam, Sdetvars=Sdetvars, Idetvars=Idetvars,
            Schi2_detvars=Schi2_detvars, Ichi2_detvars=Ichi2_detvars, Isplit_systs=Isplit_systs,
            I_TRACKSPLIT_NAMES=I_TRACKSPLIT_NAMES, SDETVAR_NAMES=SDETVAR_NAMES, IDETVAR_NAMES=IDETVAR_NAMES,
            BE_COLS=BE_COLS, S_BE_syst=S_BE_syst, I_BE_syst=I_BE_syst, SXSecWgtSyst=SXSecWgtSyst,
            IXSecWgtSyst=IXSecWgtSyst, build_systematics=build_systematics))
        print("%s done (%.0f s)" % (s, time.time() - t0), flush=True)

    # --- joint 1p + Mp ------------------------------------------------------------
    PLOTDIR = JOINTDIR
    TAG = ""
    XMIN = None
    print("\n==== joint -> %s" % JOINTDIR, flush=True)
    var2 = [var, var]
    bins2 = syst.ConcatBins([bins, bins])
    cut2 = ["selected_" + p for p in PSELS]
    NCVj = np.concatenate([results["1p"]["SCV"], results["mp"]["SCV"],
                           results["1p"]["ICV"], results["mp"]["ICV"]])
    _, _, systematicsJ = build_systematics("joint")
    covsJ = [sj.cov(var2, cut2, bins2, NCVj) for sj in systematicsJ]
    covsJ.append(np.diag(NCVj))
    n = nbin
    total = np.sum(covsJ, axis=0)

    # the joint diagonal blocks must reproduce the per-multiplicity covariances
    for s, idx in [("1p", np.r_[0:n, 2*n:3*n]), ("mp", np.r_[n:2*n, 3*n:4*n])]:
        for l, cj, cs in zip(labels, covsJ, results[s]["covs"]):
            sub = cj[np.ix_(idx, idx)]
            err = np.max(np.abs(sub - cs)) / max(np.max(np.abs(cs)), 1e-30)
            print("  joint vs %s %-10s max rel. diff %.2e" % (s, l, err))
            if err > 1e-6:
                print("  WARNING: joint %s block differs from the per-sample covariance" % l)

    np.savez(os.path.join(JOINTDIR, "joint_covariances.npz"), bins=bins, NCV=NCVj,
             labels=np.array(labels), covs=np.array(covsJ),
             blocks=np.array(["SBND 1u1p", "SBND 1uMp", "ICARUS 1u1p", "ICARUS 1uMp"]))

    # (1) 4-block correlation matrices
    block_names = ["SBND " + PSEL_LABEL["1p"], "SBND " + PSEL_LABEL["mp"],
                   "ICARUS " + PSEL_LABEL["1p"], "ICARUS " + PSEL_LABEL["mp"]]
    for k, (l, f) in enumerate([("Flux", "flux"), ("Geant4", "g4"), ("XSec", "xsec"),
                                ("POT Norm.", "norm"), ("Detector", "detector")]):
        plot_corr_blocks(corr_f(covsJ[k]), block_names, "%s Systematic Correlation" % l,
                         "joint_%s_correlation" % f)
    plot_corr_blocks(corr_f(total), block_names, "Total Uncertainty Correlation", "joint_total_correlation")
    plot_corr_blocks(corr_f(total - covsJ[7]), block_names, "Total Systematic Correlation",
                     "joint_totalsyst_correlation")

    # (2) ICARUS 1p and Mp conditioned on SBND 1p + Mp
    omit = [("No Flux", covsJ[0]), ("No G4", covsJ[1]), ("No XSec", covsJ[2]),
            ("No Detector", covsJ[4] + covsJ[5] + covsJ[6]), ("No POT", covsJ[3]), ("No Stat.", covsJ[7])]
    condJ = [(name, cond(total - c, 2*n, "joint " + name)) for name, c in omit]
    condJ.append(("All", cond(total, 2*n, "joint All")))
    ICVj = {"1p": results["1p"]["ICV"], "mp": results["mp"]["ICV"]}

    summary = ["ICARUS fractional uncertainty per reco-energy bin (%s)" % os.path.basename(args.plotbase), ""]
    for bi, s in enumerate(PSELS):
        XMIN = MP_XMIN if s == "mp" else None
        blk = slice(bi*n, (bi + 1)*n)
        ICV = ICVj[s]
        plt.figure()
        for name, c in condJ:
            style = dict(color="black") if name == "All" else {}
            hist_dropedges(centers, bins, np.sqrt(np.diag(c)[blk])/ICV, histtype="step", linewidth=2, label=name, **style)
        own = results[s]["cond_cov_all"]
        hist_dropedges(centers, bins, np.sqrt(np.diag(own))/ICV, histtype="step", linewidth=2, color="black",
                       linestyle=":", label="All, SBND %s only" % PSEL_LABEL[s])
        add_style(plt.gca(), "Reco. Neutrino Energy [GeV]", det="ICARUS", title=PSEL_LABEL[s],
                  ylabel="Frac. Unc. on ICARUS %s\nConditioned by SBND 1p+Mp" % PSEL_LABEL[s],
                  legend_loc="upper left", legend_ncol=2, legend_fontsize=11)
        savefig("condconstraint_joint_ICARUS%s" % ("1p" if s == "1p" else "Mp"), bbox_inches="tight")

        unc = np.sqrt(np.diag(total)[2*n + bi*n: 2*n + (bi + 1)*n])/ICV
        summary.append("ICARUS %s" % s)
        summary.append("%-12s %10s %10s %14s %14s" % ("bin [GeV]", "ICV", "unconstr.", "SBND same-mult", "SBND 1p+Mp"))
        allJ = dict(condJ)["All"]
        for b in range(n):
            summary.append("%5.2f-%-5.2f %10.1f %10.4f %14.4f %14.4f" % (
                bins[b], bins[b+1], ICV[b], unc[b], np.sqrt(own[b, b])/ICV[b], np.sqrt(allJ[blk, blk][b, b])/ICV[b]))
        summary.append("")
    # (2b) the same, unfolded to true neutrino energy: each ICARUS multiplicity
    # is unfolded with its own reco -> true migration (as the per-sample
    # condconstraint_signalbox_systematics_trueE), i.e. a block-diagonal
    # migration over the conditioned [ICARUS 1p, ICARUS Mp] covariance
    true_bins = np.array([-2] + list(bins) + [np.inf])
    tbc = (true_bins[:-1] + true_bins[1:])/2
    plot_bins = true_bins[1:-1]
    Imig1, ImigM = results["1p"]["Imig"], results["mp"]["Imig"]
    Imigj = np.block([[Imig1, np.zeros((Imig1.shape[0], ImigM.shape[1]))],
                      [np.zeros((ImigM.shape[0], Imig1.shape[1])), ImigM]])
    nt = Imig1.shape[0]
    summary.append("ICARUS unfolded (true neutrino energy)")
    summary.append("%-12s %8s %14s %14s   %8s %14s %14s" % ("bin [GeV]", "1p ICVT", "1p same-mult", "1p joint",
                                                           "Mp ICVT", "Mp same-mult", "Mp joint"))
    unf = {}
    for bi, s in enumerate(PSELS):
        XMIN = MP_XMIN if s == "mp" else None
        blk = slice(bi*nt, (bi + 1)*nt)
        ICVT = results[s]["ICVT"]
        Imig = results[s]["Imig"]
        plt.figure()
        for name, c in condJ:
            u = np.sqrt(np.diag(Imigj@c@Imigj.T)[blk])/ICVT
            if name == "All":
                unf[s] = u
            plt.hist(tbc[1:-1], bins=B(plot_bins), weights=u[1:-1], histtype="step", linewidth=2, label=name,
                     **(dict(color="black") if name == "All" else {}))
        own = np.sqrt(np.diag(Imig@results[s]["cond_cov_all"]@Imig.T))/ICVT
        unf[s + "_own"] = own
        plt.hist(tbc[1:-1], bins=B(plot_bins), weights=own[1:-1], histtype="step", linewidth=2, color="black",
                 linestyle=":", label="All, SBND %s only" % PSEL_LABEL[s])
        add_style(plt.gca(), "True Neutrino Energy [GeV]", det="ICARUS", title=PSEL_LABEL[s],
                  ylabel="Unfolded Frac. Unc. on ICARUS %s\nConditioned by SBND 1p+Mp" % PSEL_LABEL[s],
                  legend_loc="upper left", legend_ncol=2, legend_fontsize=11)
        savefig("condconstraint_joint_ICARUS%s_trueE" % ("1p" if s == "1p" else "Mp"), bbox_inches="tight")
    for b in range(1, nt - 1):
        summary.append("%5.2f-%-5.2f %8.1f %14.4f %14.4f   %8.1f %14.4f %14.4f" % (
            true_bins[b], true_bins[b+1],
            results["1p"]["ICVT"][b], unf["1p_own"][b], unf["1p"][b],
            results["mp"]["ICVT"][b], unf["mp_own"][b], unf["mp"][b]))
    summary.append("")

    with open(os.path.join(JOINTDIR, "summary.txt"), "w") as f:
        f.write("\n".join(summary) + "\n")
    print("\n".join(summary))
    print("DONE (%.0f s)" % (time.time() - t0), flush=True)


def run_selection(s, cut, locals_):
    """Every saved plot of the notebook for one proton multiplicity."""
    L = locals_
    Sdf, Idf, Soffbeam, Ioffbeam = L["Sdf"], L["Idf"], L["Soffbeam"], L["Ioffbeam"]
    Sdetvars, Idetvars = L["Sdetvars"], L["Idetvars"]
    Schi2_detvars, Ichi2_detvars = L["Schi2_detvars"], L["Ichi2_detvars"]
    Isplit_systs, I_TRACKSPLIT_NAMES = L["Isplit_systs"], L["I_TRACKSPLIT_NAMES"]
    SDETVAR_NAMES, IDETVAR_NAMES = L["SDETVAR_NAMES"], L["IDETVAR_NAMES"]
    BE_COLS = L["BE_COLS"]
    xE = "Reco. Neutrino Energy [GeV]"

    # --- diagnostics ---
    for det, df, cols in [("SBND", Sdf, {1: ("selected_trig_up_", "selected_trig_dn_")}),
                          ("ICARUS", Idf, {2: ("selected_trig_r2_up_", "selected_trig_r2_dn_"),
                                           4: ("selected_trig_r4_up_", "selected_trig_r4_dn_")})]:
        for run, (up, dn) in cols.items():
            inrun = df.Run == run
            wcv = df.loc[df[cut] & inrun, wgt].sum()
            wup = df.loc[df[up + s] & inrun, wgt].sum()
            wdn = df.loc[df[dn + s] & inrun, wgt].sum()
            print("  %s Run %i trigger: selected %.2f -> %+.2f%% up / %+.2f%% dn" %
                  (det, run, wcv, 100*(wup/wcv - 1), 100*(wdn/wcv - 1)))
    for name, sp in zip(I_TRACKSPLIT_NAMES, Isplit_systs):
        wcv = Idf.loc[Idf[cut], wgt].sum()
        w = Idf[wgt].to_numpy()*(1 - sp.frac*sp.crosses)
        wuniv = w[Idf[cut].to_numpy()].sum() + sp.frac*sp.splitdf.loc[sp.splitdf[cut], wgt].sum()
        print("  %-26s: selected %.2f -> %.2f (%+.3f%%)" % (name, wcv, wuniv, 100*(wuniv/wcv - 1)))

    # --- alt samples ---
    for det, df, off, ooav in [("ICARUS", Idf, Ioffbeam, OOAVICARUS), ("SBND", Sdf, Soffbeam, OOAVSBND)]:
        plt.figure()
        m = df[cut] & ~ooav(df)
        plt.hist(df.loc[m, var], bins=B(bins[EDGE]), weights=df.loc[m, wgt], histtype="step", linewidth=2, label="CV")
        m = df[cut] & ooav(df)
        plt.hist(df.loc[m, var], bins=B(bins[EDGE]), weights=df.loc[m, wgt], histtype="step", linewidth=2, label="OOAV")
        plt.hist(off.loc[off[cut], var], bins=B(bins[EDGE]), weights=off.loc[off[cut], wgt], histtype="step", linewidth=2, label="Beam Off")
        add_style(plt.gca(), xE, det=det)
        plt.yscale("log")
        savefig("%s_altsamples" % det, bbox_inches="tight")

    # --- weight universes ---
    for det, df in [("SBND", Sdf), ("ICARUS", Idf)]:
        c = df[cut]
        for kind, cols, lab in [("flux", ["flux_univ%i" % i for i in range(20)], "Flux Univ's"),
                                ("g4", ["g4_univ%i" % i for i in range(20)], "G4 Univ's"),
                                ("xsec", ["%s_univ" % x for x in loaddf.xsec_syst], "XSec Syst's")]:
            plt.figure()
            plt.hist(df.loc[c, var], bins=B(bins), weights=df.loc[c, wgt], histtype="step", linewidth=2, label="CV")
            for i, col in enumerate(cols):
                w = df[wgt]*df[col]
                plt.hist(df.loc[c, var], bins=B(bins), weights=w[c], histtype="step", linewidth=1, color="gray",
                         label=lab if i == 0 else None)
            add_style(plt.gca(), xE, det=det)
            savefig("%s_%s_universes" % (det, kind), bbox_inches="tight")

    # --- binding energy ---
    shifts = [0.005, 0.01, 0.025]
    for det, df in [("SBND", Sdf), ("ICARUS", Idf)]:
        be_shifted = [syst.shift_binding_energy(df[BE_COLS], sh, fraction=1) for sh in shifts]
        be_shifted += [syst.shift_binding_energy(df[BE_COLS], shifts[-1], fraction=0.5)]
        top, bot = ratioplot([df] + be_shifted,
                             ["CV"] + ["$\\Delta E_\\mathrm{miss}$: %i MeV" % (sh*1000) for sh in shifts] +
                             ["50%% $\\Delta E_\\mathrm{miss}$: %i MeV" % (shifts[-1]*1000)],
                             det, cut, ratio_range=1)
        top.set_title("%s %s Shifted $E_\\mathrm{miss}$" % (det, TAG))
        savefig("%s_be_shift_variations" % det, bbox_inches="tight")
        del be_shifted

    # --- detector-variation / trigger / split ratio plots ---
    varnames = ["CV", "$+1\\sigma$", "$-1\\sigma$"]
    for name, dv in zip(SDETVAR_NAMES, Sdetvars):
        top, _ = ratioplot(dv, varnames, "SBND", cut)
        top.set_title("SBND %s %s" % (TAG, name), fontsize=14)
        savefig("SBND_detvar_ratio_%s" % clean(name), bbox_inches="tight")
    for name, dv in zip(IDETVAR_NAMES, Idetvars):
        for runtag, f in [("", None), ("r2", 2), ("r4", 4)]:
            d = dv if f is None else [x[x.Run == f] for x in dv]
            det = "ICARUS" if f is None else "ICARUS Run%i" % f
            top, _ = ratioplot(d, varnames, det, cut)
            top.set_title("%s %s %s" % (det, TAG, name), fontsize=14)
            savefig("ICARUS%s_detvar_ratio_%s" % (runtag, clean(name)), bbox_inches="tight")

    top, _ = ratioplot_varyselect(Sdf, [cut, "selected_trig_up_" + s, "selected_trig_dn_" + s], varnames, "SBND", cut)
    top.set_title("SBND %s Trigger" % TAG, fontsize=14)
    savefig("SBND_Trig_variations", bbox_inches="tight")
    for r in (2, 4):
        top, _ = ratioplot_varyselect(Idf[Idf.Run == r], [cut, "selected_trig_r%i_up_%s" % (r, s), "selected_trig_r%i_dn_%s" % (r, s)],
                                      varnames, "ICARUS Run%i" % r, cut)
        top.set_title("ICARUS Run %i %s Trigger" % (r, TAG), fontsize=14)
        savefig("ICARUSr%i_Trig_variations" % r, bbox_inches="tight")

    for name, sp in zip(I_TRACKSPLIT_NAMES, Isplit_systs):
        base = Idf[[cut, var, wgt]].copy()
        base[wgt] = Idf[wgt].to_numpy()*(1 - sp.frac*sp.crosses)
        spl = sp.splitdf[[cut, var, wgt]].copy()
        spl[wgt] = spl[wgt]*sp.frac
        top, _ = ratioplot([Idf, pd.concat([base, spl])], ["CV", "Split Universe"], "ICARUS", cut)
        top.set_title("ICARUS %s %s" % (TAG, name), fontsize=14)
        savefig("ICARUS_tracksplit_%s" % clean(name), bbox_inches="tight")

    for det, df, dvs, names in [("SBND", Sdf, Schi2_detvars, SCHI2_DETVAR_NAMES), ("ICARUS", Idf, Ichi2_detvars, ICHI2_DETVAR_NAMES)]:
        plt.figure()
        plt.hist(df.loc[df[cut], var], bins=B(bins), weights=df.loc[df[cut], wgt], histtype="step", linewidth=2, label="CV")
        for d, l in zip(dvs, names):
            plt.hist(d.loc[d[cut], var], bins=B(bins), weights=d.loc[d[cut], wgt], histtype="step", linewidth=2, label=l)
        add_style(plt.gca(), xE, det=det)
        savefig("%s_calo_variations" % det, bbox_inches="tight")

    # --- covariances ---
    Ssystematics, Isystematics, systematics = L["build_systematics"](s)
    SCV = np.histogram(Sdf.loc[Sdf[cut], var], bins=bins, weights=Sdf.loc[Sdf[cut], wgt])[0]
    ICV = np.histogram(Idf.loc[Idf[cut], var], bins=bins, weights=Idf.loc[Idf[cut], wgt])[0]
    print("  SBND CV %.1f  ICARUS CV %.1f events" % (SCV.sum(), ICV.sum()))
    Scovs = [x.cov(var, cut, bins, SCV) for x in Ssystematics] + [np.diag(SCV)]
    Icovs = [x.cov(var, cut, bins, ICV) for x in Isystematics] + [np.diag(ICV)]
    S_detsyst_covs = [x.cov(var, cut, bins, SCV) for x in Ssystematics[labels.index("Detector")].systs]
    I_detsyst_covs = [x.cov(var, cut, bins, ICV) for x in Isystematics[labels.index("Detector")].systs]
    S_DET_NAMES = SDETVAR_NAMES + SCHI2_DETVAR_NAMES + S_TRIGGER_NAMES + ["Corsika", "MC Stat"]
    I_DET_NAMES = IDETVAR_NAMES + ICHI2_DETVAR_NAMES + I_TRIGGER_NAMES + I_TRACKSPLIT_NAMES + ["MC Stat"]

    # --- budgets ---
    for det, covs_, CV in [("SBND", Scovs, SCV), ("ICARUS", Icovs, ICV)]:
        plt.figure()
        for c, l in zip(covs_, labels):
            hist_dropedges(centers, bins, np.sqrt(np.diag(c))/CV, label=l, histtype="step", linewidth=2)
        add_style(plt.gca(), xE, det=det, ylabel="Fractional Uncertainty", legend_loc="upper center", legend_ncol=2)
        plt.ylim(bottom=0)
        savefig("%s_signalbox_systematics" % det, group="budget", bbox_inches="tight")

        plt.figure()
        m = CV > 0
        for c, l in zip(covs_, labels):
            plt.hist(centers[m], bins=B(bins), weights=(np.sqrt(np.diag(c))/np.where(m, CV, 1))[m], label=l, histtype="step", linewidth=2)
        add_style(plt.gca(), xE, det=det, ylabel="Fractional Uncertainty", legend_loc="upper center", legend_ncol=2)
        plt.yscale("log")
        plt.ylim(bottom=1e-3)
        savefig("%s_signalbox_systematics_log" % det, group="budget_log", bbox_inches="tight")

    flush("budget")
    flush("budget_log")

    def breakdown(det, detcovs, names, covs_, CV, keep, legend_title, fname, with_tail):
        plt.figure()
        combined = np.zeros(covs_[0].shape)
        for c, l in zip(detcovs, names):
            if keep is None or l in keep:
                hist_dropedges(centers, bins, np.sqrt(np.diag(c))/CV, label=l, histtype="step", linewidth=2,
                               linestyle="--" if l == "MC Stat" else "-")
                combined += c
        if with_tail:
            for c, l in zip(covs_[5:], labels[5:]):
                hist_dropedges(centers, bins, np.sqrt(np.diag(c))/CV, label=l, histtype="step", linewidth=2, linestyle="--")
                combined += c
        hist_dropedges(centers, bins, np.sqrt(np.diag(combined))/CV, label="Total", histtype="step", color="black", linewidth=2, linestyle=":")
        add_style(plt.gca(), xE, det=det, ylabel="Fractional Uncertainty", legend_loc="upper center", legend_ncol=3,
                  legend_title=legend_title, legend_fontsize=10)
        plt.ylim(bottom=0)
        savefig(fname, group=fname.split("_", 1)[1], bbox_inches="tight")

    breakdown("SBND", S_detsyst_covs, S_DET_NAMES, Scovs, SCV, None, "Uncorrelated Uncertainties", "SBND_signalbox_systematics_uncorr", True)
    breakdown("SBND", S_detsyst_covs, S_DET_NAMES, Scovs, SCV, SDETVAR_NAMES + S_TRIGGER_NAMES + ["Corsika"], "Det. Var. Uncertainties", "SBND_signalbox_systematics_detvar", False)
    breakdown("SBND", S_detsyst_covs, S_DET_NAMES, Scovs, SCV, SCHI2_DETVAR_NAMES, "Calorimetric Uncertainties", "SBND_signalbox_systematics_calo", False)
    breakdown("ICARUS", I_detsyst_covs, I_DET_NAMES, Icovs, ICV, None, "Uncorrelated Uncertainties", "ICARUS_signalbox_systematics_uncorr", True)
    breakdown("ICARUS", I_detsyst_covs, I_DET_NAMES, Icovs, ICV, IDETVAR_NAMES + I_TRIGGER_NAMES + I_TRACKSPLIT_NAMES, "Det. Var. Uncertainties", "ICARUS_signalbox_systematics_detvar", False)
    breakdown("ICARUS", I_detsyst_covs, I_DET_NAMES, Icovs, ICV, ICHI2_DETVAR_NAMES, "Calorimetric Uncertainties", "ICARUS_signalbox_systematics_calo", False)
    for g in ("signalbox_systematics_uncorr", "signalbox_systematics_detvar", "signalbox_systematics_calo"):
        flush(g)

    # --- joint SBND/ICARUS covariance ---
    CVc = np.concatenate((SCV, ICV))
    covs = [x.cov(var, cut, bins, CVc) for x in systematics] + [np.diag(CVc)]
    for k, (l, f) in enumerate([("Flux", "flux"), ("Geant4", "g4"), ("XSec", "xsec"),
                                ("POT Norm.", "norm"), ("Detector", "detector")]):
        plot_corr_blocks(corr_f(covs[k]), ["SBND", "ICARUS"], "%s %s Systematic Correlation" % (TAG, l), "%s_correlation" % f)

    ratio = SCV / ICV
    tot = np.sum(covs, axis=0)
    parts = [("Flux", covs[0]), ("G4", covs[1]), ("XSec", covs[2]), ("Detector", covs[4] + covs[5] + covs[6]),
             ("POT Norm.", covs[3]), ("Stat.", covs[7]), ("All", tot)]
    plt.figure()
    for l, c in parts:
        hist_dropedges(centers, bins, np.sqrt(np.diag(ratio_cov_full(SCV, ICV, c)))/ratio, histtype="step", linewidth=2, label=l)
    add_style(plt.gca(), xE, det="SBND / ICARUS", ylabel="Uncertainty on SBND/ICARUS Ratio", legend_loc="upper left", legend_ncol=2)
    plt.ylim(bottom=0)
    savefig("ratio_signalbox_systematics", bbox_inches="tight")

    conds = [("No " + ("POT" if l == "POT Norm." else l), cond(tot - c, nbin, "%s No %s" % (s, l))) for l, c in parts[:-1]]
    cond_cov_all = cond(tot, nbin, "%s All" % s)
    conds.append(("All", cond_cov_all))
    plt.figure()
    for l, c in conds:
        hist_dropedges(centers, bins, np.sqrt(np.diag(c))/ICV, histtype="step", linewidth=2, label=l)
    add_style(plt.gca(), xE, det="", ylabel="Frac. Unc. on ICARUS\nConditioned by SBND", legend_loc="upper left", legend_ncol=2)
    savefig("condconstraint_signalbox_systematics", bbox_inches="tight")

    # --- unfolded (true energy) ---
    true_bins = np.array([-2] + list(bins) + [np.inf])
    tbc = (true_bins[:-1] + true_bins[1:])/2
    SCV2 = np.histogramdd([Sdf.loc[Sdf[cut], var], Sdf.loc[Sdf[cut], "true_E"]], bins=[bins, true_bins], weights=Sdf.loc[Sdf[cut], wgt])[0]
    ICV2 = np.histogramdd([Idf.loc[Idf[cut], var], Idf.loc[Idf[cut], "true_E"]], bins=[bins, true_bins], weights=Idf.loc[Idf[cut], wgt])[0]
    SCVT = np.histogram(Sdf.loc[Sdf[cut], "true_E"], bins=true_bins, weights=Sdf.loc[Sdf[cut], wgt])[0]
    ICVT = np.histogram(Idf.loc[Idf[cut], "true_E"], bins=true_bins, weights=Idf.loc[Idf[cut], wgt])[0]
    Smig = SCV2.T / SCV2.sum(axis=1)
    Imig = ICV2.T / ICV2.sum(axis=1)
    migration = np.block([[Smig, np.zeros(Smig.shape)], [np.zeros(Imig.shape), Imig]])
    unfold = lambda C: migration@C@migration.T
    uS = lambda C: Smig@C@Smig.T
    uI = lambda C: Imig@C@Imig.T

    ratioT = SCVT / ICVT
    ratioT = ratioT / ratioT.mean()
    plot_bins = true_bins[1:-1]
    xT = "True Neutrino Energy [GeV]"
    plt.figure()
    for l, c in parts:
        plt.hist(tbc[1:-1], bins=B(plot_bins), weights=(np.sqrt(np.diag(ratio_cov_full(SCVT, ICVT, unfold(c))))/ratioT)[1:-1],
                 histtype="step", linewidth=2, label=l)
    add_style(plt.gca(), xT, det="SBND / ICARUS", ylabel="Unfolded Unc. on SBND/ICARUS Ratio", legend_loc="upper center", legend_ncol=2)
    plt.ylim(bottom=0)
    savefig("ratio_signalbox_systematics_trueE", bbox_inches="tight")

    for det, detcovs, names, covs_, CVT, u in [("ICARUS", I_detsyst_covs, I_DET_NAMES, Icovs, ICVT, uI),
                                              ("SBND", S_detsyst_covs, S_DET_NAMES, Scovs, SCVT, uS)]:
        plt.figure()
        combined = np.zeros(covs_[0].shape)
        for c, l in zip(detcovs, names):
            plt.hist(tbc[1:-1], bins=B(plot_bins), weights=(np.sqrt(np.diag(u(c)))/CVT)[1:-1], label=l, histtype="step",
                     linewidth=2, linestyle="--" if l == "MC Stat" else "-")
            combined += c
        for c, l in zip(covs_[5:], labels[5:]):
            plt.hist(tbc[1:-1], bins=B(plot_bins), weights=(np.sqrt(np.diag(u(c)))/CVT)[1:-1], label=l, histtype="step", linewidth=2, linestyle="--")
            combined += c
        plt.hist(tbc[1:-1], bins=B(plot_bins), weights=(np.sqrt(np.diag(u(combined)))/CVT)[1:-1], label="Total", histtype="step", color="black", linewidth=2, linestyle=":")
        add_style(plt.gca(), xT, det=det, ylabel="Unfolded Fractional Uncertainty", legend_loc="upper center", legend_ncol=3,
                  legend_title="Uncorrelated Uncertainties", legend_fontsize=10)
        plt.ylim(bottom=0)
        savefig("%s_signalbox_systematics_trueE_uncorr" % det, group="trueE_uncorr", bbox_inches="tight")
    flush("trueE_uncorr")

    plt.figure()
    for l, c in conds:
        plt.hist(tbc[1:-1], bins=B(plot_bins), weights=(np.sqrt(np.diag(uI(c)))/ICVT)[1:-1], histtype="step", linewidth=2, label=l)
    add_style(plt.gca(), "True. Neutrino Energy [GeV]", det="", ylabel="Unfolded Frac. Unc. on ICARUS\nConditioned by SBND",
              legend_loc="upper right", legend_ncol=2)
    plt.ylim(bottom=0)
    savefig("condconstraint_signalbox_systematics_trueE", bbox_inches="tight")

    return dict(SCV=SCV, ICV=ICV, covs=covs, cond_cov_all=cond_cov_all, Imig=Imig, ICVT=ICVT)


if __name__ == "__main__":
    main()
