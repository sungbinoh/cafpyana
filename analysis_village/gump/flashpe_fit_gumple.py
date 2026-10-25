"""MC flash-PE scale fit, on the rewgted-25 productions.

Script version of nb/FlashMCDataComparison.ipynb (the study behind the
loaddf.MC_PE_SCALE factors). For each of SBND, ICARUS Run 2 and ICARUS Run 4 the
raw MC max in-time flash PE is scaled by s and fit to the on - off data flash PE
distribution. The fit, plots and bootstrap follow the notebook, with a few changes:

  * Input: event-level flash PE instead of the old per-flash tables (see the
    schemas below).
  * A normalization nuisance alpha on the neutrino MC (CV + dirt) with a Gaussian
    prior of --norm-unc (default 20%), in lieu of the full systematic model:
        chi2(s, a) = sum_bins (N_data - (1+a) P(s))^2 / V + (a / norm_unc)^2
    with V = data stat (ON + OFF_w^2 OFF) + MC stat (sum w^2). alpha is profiled
    analytically. The notebook's fit (no norm uncertainty, data stat only) is
    also run as a reference.
  * Uncertainties, as in the notebook: delta-chi2 = 1 from the profiled chi2 scan
    (the notebook's Hessian), bootstrap over the fit bins, and the
    flash-selection systematic. The notebook took |s_max - s_all|, with s_all
    fit to every in-time flash; no schema here stores individual flashes,
    so s_all is fit to flash_sumpe (sum of in-time flash PE) instead.
    Total = bootstrap (+) flash-selection.
  * The notebook's flash-time plots are dropped: no schema stores flash times.

Two input schemas (--schema):
  gumple     full GUMPLE evt DFs (sbn-rewgted-25/), slices at the Contained stage
             (mcdata_comparison_gumple.FV, no flash cut). flash_maxpe is the max
             in-time flash PE of the event (ICARUS: of the vertex cryostat, see
             loaddf). MC CV weights from the wgt table.
  flashonly  flash-only DFs (sbn-rewgted-25-flashonly/): one flashevt row per
             trigger, no slices, so every triggered event enters ("All Triggers").
             ICARUS flash_maxpe is the max over both cryostats (no vertex to pick
             one). Data cuts (good runs, on-beam spill beam quality) and the
             loaddf dedup are applied here; no wgt table -> cvwgt = 1. The Run 4
             on-beam file is a 10% subsample of the full ReCAF2026 stream
             (make_unblind_subsample.py), as for gumple.
Input files come from mcdata_comparison_gumple.detector_files, filtered to those
present in --df-dir.

Outputs, in <plotdir>/{png,pdf}/, prefixed by detector:
  _flash_totalpe            max flash PE, MC at the fitted s (+-norm_unc band)
  _flash_totalpe_current    the same at the loaddf.MC_PE_SCALE currently applied
  _flash_sumpe              sum of in-time flash PE, MC at its own fitted s
  _flash_scale_bootstrap    bootstrap distribution of s
  _flash_scale_chi2scan     profiled delta-chi2(s), with and without norm unc.
and <plotdir>/flashpe_fit.tsv, <plotdir>/summary.txt.

Trigger systematic: the CV MC PE scale (loaddf.MC_PE_SCALE) is the flashonly
(all-trigger) fit, and its uncertainty (loaddf.MC_PE_SCALE_UNC) is |all-trigger -
Contained (gumple) fit|. Once both schemas have run into sibling <plotdir>s, the
parent directory gets trigger_scale_summary.tsv and trigger_scale_table.tex (rows
for the technote's tab:trigger-systematic); a warning is logged if loaddf's
constants do not match. Keep loaddf and the technote table in sync with them.

    python flashpe_fit_gumple.py --schema {gumple,flashonly} [--df-dir DIR] [--plotdir DIR]
                                 [--detectors ...] [--nboot N] [--norm-unc F] [--full-onbeam]
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

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

import h5py

import mcdata_comparison_gumple as mc
import dataquality
import loaddf

FONTSIZE = 14
DETECTORS = ["SBND", "ICARUS Run2", "ICARUS Run4"]
DET_LABEL = {"SBND": "SBND", "ICARUS Run2": "ICARUS Run 2", "ICARUS Run4": "ICARUS Run 4"}

PLOTBASE = "/Users/gputnam/Work/osc/cafpyana/plots-gumple-2026-10-24-rewgted25/flashpe_fit/"
SCHEMA = {
    "gumple": dict(df_dir="/Users/gputnam/Work/osc/sbn-rewgted-25/", reweight_aFF=True, stage="Contained"),
    "flashonly": dict(df_dir="/Users/gputnam/Work/osc/sbn-rewgted-25-flashonly/", reweight_aFF=False,
                      stage="All Triggers"),
}

# Fit configuration from the notebook
# The fit (and the shaded band on the plots) starts at the trigger-emulation
# flash-PE threshold, gumple_cuts.FLASH_PE_THRESHOLD (SBND 2000, ICARUS Run 2 4000,
# Run 4 2000 PE = the lowest 4 bins the notebook skipped).
FIT_MIN_PE = mc.gc.FLASH_PE_THRESHOLD
S_RANGE = (0.1, 0.8)
S_GRID = np.linspace(*S_RANGE, 3501)  # step 2e-4, well below the fit uncertainty
SEED = 42
BINS_PE = {
    "SBND": np.linspace(0, 30000, 61),
    "ICARUS Run2": np.linspace(0, 40000, 41),
    "ICARUS Run4": np.linspace(0, 30000, 61),
}
VARS = {"flash_maxpe": "Max In-Time Flash Total PE",
        "flash_sumpe": "Sum of In-Time Flash PE"}

MC_COLORS = ["#315031", "#8B4513"]
MC_LABELS = ["MC Simulation", "Add. Dirt MC"]


# NB: keep at module scope, named FV: loaddf's cache key is the preselection's
# __module__.__qualname__, and other scripts register different bodies under the
# same "__main__.FV" key, so the mask is re-applied after every load (refv).
def FV(df):
    return mc.FV(df)


def refv(d):
    return d[FV(d)].copy() if len(d) else d


# ============================================================
# Loading
# ============================================================
def existing(files, log):
    keep = [f for f in files if os.path.exists(f)]
    for f in files:
        if f not in keep:
            log("  missing, skipped: %s" % f)
    return keep


def load_flashonly(fname, det, kind):
    """flashevt frame of one flash-only .df, all splits -> (df[VARS], POT).

    kind: "on" (good runs + spill beam quality), "off" (good runs) or "mc" (POT
    from hdr, times the SBND MC POT factor). Duplicated events are dropped as in
    loaddf.load_one (keep=False; (run, evt), plus nu_E0 for MC), and the MC POT
    scaled by the surviving fraction. The PE is raw (no MC_PE_SCALE)."""
    assert dataquality.external_beam_file(fname) is None or kind != "on", \
        "%s has no spill information of its own" % fname
    tag = os.path.basename(fname)
    out, pot = [], 0.
    for hk in dataquality._keys(fname, "hdr"):
        idf = int(hk.split("_")[-1])
        hdr = pd.read_hdf(fname, hk)
        df = pd.read_hdf(fname, "flashevt_%i" % idf)[list(VARS)].fillna(0.).astype(float)
        df = df.reindex(hdr.index).fillna(0.)

        keep = np.ones(len(df), dtype=bool)
        if kind != "mc":
            keep &= np.asarray(dataquality.data_quality_cut(hdr.run, det))
        if kind == "on":
            spill = dataquality.event_spill_fom(dataquality.spill_table(fname, idf, det), hdr,
                                                cols=("TOR860", "TOR875", "THCURR", "FOM"))
            keep &= dataquality.beam_quality_cut(spill, det, prefix="").to_numpy()
            print("[%s idf=%i] %d/%d events mapped to a spill" % (tag, idf, spill.mapped.sum(), len(spill)))
        if kind != "mc":
            print("[%s idf=%i] data cuts (BQ=%s, DQ=True): kept %d/%d events" % (
                tag, idf, kind == "on", keep.sum(), len(df)))

        match = hdr[["run", "evt"]].copy()
        dedup_cols = ["run", "evt"]
        if kind == "mc":
            with h5py.File(fname, "r") as f:
                has_mcnu = ("mcnu_%i" % idf) in f
            match["nu_E0"] = (pd.read_hdf(fname, "mcnu_%i" % idf).nu_E.groupby(level=[0, 1]).max()
                              .reindex(hdr.index) if has_mcnu else np.nan)
            dedup_cols.append("nu_E0")
        dup = match.duplicated(subset=dedup_cols, keep=False).to_numpy()
        print("[%s idf=%i] dedup: dropped %d duplicated %s keys (%d hdr rows)" % (
            tag, idf, match[dup][dedup_cols].drop_duplicates().shape[0], tuple(dedup_cols), dup.sum()))

        if kind == "mc":
            pot += hdr.pot.sum() * (1. - dup.sum()/len(hdr)) * loaddf._mc_pot_scale(det, True)
        out.append(df[keep & ~dup])
    return pd.concat(out), pot


def load_detector(det, df_dir, schema, full_onbeam, log):
    """-> dict(norm, ON, OFF, MC=[cv, dirt]) with raw (unscaled) MC flash PE."""
    files = mc.detector_files(det, df_dir, full_onbeam=full_onbeam)
    log("on-beam: %s" % os.path.basename(files["ONBEAM"]))
    offbeam = existing(files["OFFBEAM_FILES"], log)
    mc_files = existing(files["MC_FILES"], log)
    dirt_files = existing(files["DIRT_FILES"], log)
    log("off-beam: %s" % ", ".join(map(os.path.basename, offbeam)))
    log("MC CV: %s" % ", ".join(map(os.path.basename, mc_files)))
    log("dirt: %s" % ", ".join(map(os.path.basename, dirt_files)))

    norm = dataquality.data_norm(det, files["ONBEAM"], offbeam, log=lambda m: log("    " + str(m)))
    POT = norm["POT"]

    cols = list(VARS)
    if schema == "flashonly":
        ON = load_flashonly(files["ONBEAM"], det, "on")[0]
        OFF = pd.concat([load_flashonly(f, det, "off")[0] for f in offbeam])
        mcs = []
        for flist in (mc_files, dirt_files):
            if not flist:
                mcs.append(None)
                continue
            parts = [load_flashonly(f, det, "mc") for f in flist]
            d = pd.concat([p[0] for p in parts])
            pot = sum(p[1] for p in parts)
            d["cvwgt"] = 1.
            loaddf.scale_pot(d, pot, POT)
            log("  MC %s: POT %.4e (scale %.4g), %d events" % (
                os.path.basename(flist[0]), pot, POT/pot, len(d)))
            mcs.append(d[cols + ["glob_scale"]])
        log("data POT %.4e, OFF_w %.5g; on-beam %d, off-beam %d events" % (POT, norm["OFF_w"], len(ON), len(OFF)))
        return dict(norm=norm, ON=ON, OFF=OFF, MC=mcs)

    ON, _, _ = loaddf.load(files["ONBEAM"], load_truth=False, include_syst=False,
                           detector=det, preselection=FV, match_Enu=False)
    OFF = pd.concat([loaddf.load(f, load_truth=False, include_syst=False, detector=det,
                                 preselection=FV, offbeampot=True, match_Enu=False)[0]
                     for f in offbeam])
    ON, OFF = refv(ON)[cols], refv(OFF)[cols]

    mcs = []
    for flist, kw in ((mc_files, dict(reweight_aFF=SCHEMA[schema]["reweight_aFF"])), (dirt_files, {})):
        if not flist:
            mcs.append(None)
            continue
        d, _, pot = loaddf.loadl(flist, njob=min(len(flist), 4), detector=det,
                                 preselection=FV, include_syst=False, **kw)
        loaddf.scale_pot(d, pot, POT)
        d = refv(d)
        # undo the MC PE scale loaddf applied to flash_maxpe (flash_sumpe is raw)
        d["flash_maxpe"] = d["flash_maxpe"] / loaddf.MC_PE_SCALE[det]
        log("  MC %s: POT %.4e (scale %.4g), %d slices, sum cvwgt/N = %.4f" % (
            os.path.basename(flist[0]), pot, POT/pot, len(d), d.cvwgt.mean()))
        mcs.append(d[cols + ["glob_scale"]])

    log("data POT %.4e, OFF_w %.5g; on-beam %d, off-beam %d slices" % (POT, norm["OFF_w"], len(ON), len(OFF)))
    return dict(norm=norm, ON=ON, OFF=OFF, MC=mcs)


# ============================================================
# Fit
# ============================================================
class SortedMC:
    """Sorted PE + cumulative weights (and weights^2) for fast scaled histograms
    via searchsorted, as in the notebook's _mc_pred_fast."""
    def __init__(self, pe, w):
        o = np.argsort(pe)
        self.pe = np.asarray(pe)[o]
        w = np.asarray(w)[o]
        self.cw = np.concatenate([[0], np.cumsum(w)])
        self.cw2 = np.concatenate([[0], np.cumsum(w**2)])

    def hist(self, s, bins):
        i = np.searchsorted(self.pe, bins / s)
        return np.diff(self.cw[i]), np.diff(self.cw2[i])


def data_hist(d, var, bins):
    NON, _ = np.histogram(d["ON"][var], bins=bins)
    NOFF, _ = np.histogram(d["OFF"][var], bins=bins)
    w = d["norm"]["OFF_w"]
    return NON - NOFF*w, NON + NOFF*w**2


def fit_mask(var_data, bins, pe_min):
    """Bins entering the fit: those starting at or above pe_min, and not empty."""
    return (var_data > 0) & (bins[:-1] >= pe_min)


class Fit:
    """Profiled chi2(s) on S_GRID for one detector / variable.

    P[i] and W2[i] (MC prediction and sum w^2) are tabulated on S_GRID once, so
    every refit (bootstrap, with or without norm unc.) is array arithmetic."""
    def __init__(self, d, var, bins, norm_unc, pe_min, mc_stat=True):
        self.bins, self.norm_unc = bins, norm_unc
        self.N, self.Vd = data_hist(d, var, bins)
        self.mask = fit_mask(self.Vd, bins, pe_min)
        self.mcs = [SortedMC(m[var].values, m["glob_scale"].values) for m in d["MC"] if m is not None]
        P = np.zeros((len(S_GRID), len(bins) - 1))
        W2 = np.zeros_like(P)
        for m in self.mcs:
            for i, s in enumerate(S_GRID):
                p, w2 = m.hist(s, bins)
                P[i] += p
                W2[i] += w2
        self.P = P
        self.V = self.Vd[None, :] + (W2 if mc_stat else 0.)
        self.ndof = int(self.mask.sum()) - 1

    def chi2(self, counts=None):
        """-> (chi2[s], alpha[s]) with alpha profiled; counts = bootstrap bin multiplicities."""
        c = (np.ones(len(self.N)) if counts is None else counts) * self.mask
        P, V = self.P, self.V
        with np.errstate(divide="ignore", invalid="ignore"):
            iv = np.where(self.mask[None, :], c[None, :] / V, 0.)
        r = self.N[None, :] - P
        if self.norm_unc > 0:
            a = np.sum(iv*P*r, axis=1) / (np.sum(iv*P*P, axis=1) + 1./self.norm_unc**2)
            pen = (a/self.norm_unc)**2
        else:
            a = np.zeros(len(S_GRID))
            pen = 0.
        chi2 = np.sum(iv*(r - a[:, None]*P)**2, axis=1) + pen
        return chi2, a

    @staticmethod
    def minimum(chi2):
        """Grid minimum, refined by a parabola through the neighbouring points."""
        i = int(np.clip(np.argmin(chi2), 1, len(S_GRID) - 2))
        y0, y1, y2 = chi2[i-1:i+2]
        den = y0 - 2*y1 + y2
        h = S_GRID[1] - S_GRID[0]
        off = 0.5*h*(y0 - y2)/den if den > 0 else 0.
        return S_GRID[i] + np.clip(off, -h, h), i

    def result(self):
        chi2, a = self.chi2()
        s, i = self.minimum(chi2)
        # The scan is jagged (MC events hop across bin edges as s varies), so the
        # delta-chi2 = 1 width is taken from a parabola fit to the delta-chi2 < 9 region
        sel = chi2 - chi2[i] < 9.
        c2 = np.polyfit(S_GRID[sel] - s, chi2[sel], 2)[0]
        err = 1./np.sqrt(c2) if c2 > 0 else np.nan
        return dict(s=s, err=err, chi2=chi2[i], ndof=self.ndof, alpha=a[i], scan=chi2)

    def bootstrap(self, nboot, rng):
        """Resample which bins enter the fit (multinomial bin multiplicities), refit."""
        nb = len(self.N)
        out = np.empty(nboot)
        for k in range(nboot):
            counts = np.bincount(rng.integers(0, nb, size=nb), minlength=nb)
            out[k] = self.minimum(self.chi2(counts)[0])[0]
        return out


# ============================================================
# Plots
# ============================================================
def style(ax):
    ax.tick_params(axis="both", which="both", direction="in", length=6,
                   width=1.5, labelsize=FONTSIZE, top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)


def save(fig, plotdir, stem):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(plotdir, ext, stem + "." + ext), bbox_inches="tight", dpi=150)
    plt.close(fig)


def ratio_plot(d, det, var, bins, s, norm_unc, stage, alpha=None, text=None, title=""):
    """Stacked MC (CV + dirt, PE scaled by s) vs data (ON - OFF_w*OFF) with a
    data/MC ratio panel and a +-norm_unc band on the MC. alpha, if given, draws
    the post-fit (1+alpha) MC total as a dashed line."""
    N, Vd = data_hist(d, var, bins)
    Nerr = np.sqrt(Vd)
    hists, colors, labels = [], [], []
    for m, c, l in zip(d["MC"], MC_COLORS, MC_LABELS):
        if m is None:
            continue
        hists.append(SortedMC(m[var].values, m["glob_scale"].values).hist(s, bins)[0])
        colors.append(c)
        labels.append(l)
    tot = np.sum(hists, axis=0)
    centers = (bins[:-1] + bins[1:])/2

    fig, (ax0, ax1) = plt.subplots(2, 1, height_ratios=[3, 1], sharex=True, figsize=(7.5, 6.5))
    ax0.hist([centers]*len(hists), bins=bins, weights=hists, histtype="stepfilled",
             stacked=True, color=colors, alpha=0.7, label=labels)
    band = lambda ax, lo, hi, **kw: ax.fill_between(bins, np.append(lo, lo[-1]), np.append(hi, hi[-1]),
                                                    step="post", **kw)
    band(ax0, tot*(1 - norm_unc), tot*(1 + norm_unc), facecolor="none", edgecolor="gray",
         hatch="///", linewidth=0, label="$\\pm$%d%% $\\nu$ norm." % round(100*norm_unc))
    if alpha is not None:
        ax0.stairs(tot*(1 + alpha), bins, color="black", linestyle="--", linewidth=1.5,
                   label="Post-fit norm. (%+.1f%%)" % (100*alpha))
    ax0.errorbar(centers, N, Nerr, color="black", linestyle="none", marker=".", label="Data (ON $-$ OFF)")

    with np.errstate(invalid="ignore", divide="ignore"):
        den = np.where(tot > 0, tot, np.nan)
        ax1.errorbar(centers, N/den, Nerr/den, color="black", linestyle="none", marker=".")
    band(ax1, np.full(len(tot), 1 - norm_unc), np.full(len(tot), 1 + norm_unc),
         facecolor="none", edgecolor="gray", hatch="///", linewidth=0)
    if alpha is not None:
        ax1.axhline(1 + alpha, color="black", linestyle="--", linewidth=1.5)
    ax1.axhline(1, color="red", linestyle="--")
    ax1.set_ylim([0.5, 1.5])
    for ax in (ax0, ax1):
        ax.axvspan(bins[0], FIT_MIN_PE[det], color="gray", alpha=0.15, linewidth=0)
        style(ax)

    ax0.set_ylabel("Events / %.2f$\\times10^{19}$ POT" % (d["norm"]["POT"]/1e19),
                   fontsize=FONTSIZE, fontweight="bold")
    ax1.set_ylabel("Data / MC", fontsize=FONTSIZE, fontweight="bold")
    ax1.set_xlabel(VARS[var], fontsize=FONTSIZE, fontweight="bold")
    ax1.set_xlim(bins[0], bins[-1])
    ax0.set_title("$\\bf{%s}$ %s %s" % (DET_LABEL[det].replace(" ", "\\ "), stage, title), fontsize=FONTSIZE + 2)
    ax0.legend(fontsize=FONTSIZE - 3, loc="upper right")
    if text:
        ax0.text(0.97, 0.45, text, transform=ax0.transAxes, ha="right", va="top",
                 fontsize=FONTSIZE - 3, bbox=dict(boxstyle="round", fc="wheat", alpha=0.8))
    plt.subplots_adjust(hspace=0.05)
    return fig


def plot_bootstrap(det, s_boot, s, err):
    fig, ax = plt.subplots()
    ax.hist(s_boot, bins=30, color="#315031", alpha=0.7, edgecolor="black")
    ax.axvline(s, color="red", linestyle="--", label="Best fit: %.4f" % s)
    ax.axvline(s - err, color="red", linestyle=":", alpha=0.5)
    ax.axvline(s + err, color="red", linestyle=":", alpha=0.5)
    ax.set_xlabel("PE Scale Factor $s$", fontsize=FONTSIZE)
    ax.set_ylabel("Bootstrap Samples", fontsize=FONTSIZE)
    ax.set_title("%s Bootstrap Distribution (N=%d), $\\sigma$ = %.4f" % (DET_LABEL[det], len(s_boot), err),
                 fontsize=FONTSIZE)
    ax.legend(fontsize=FONTSIZE)
    style(ax)
    return fig


def plot_scan(det, fits, current):
    fig, ax = plt.subplots()
    for (label, r), ls in zip(fits, ["-", "--", ":"]):
        ax.plot(S_GRID, r["scan"] - r["chi2"], linestyle=ls, linewidth=2,
                label="%s: %.4f$\\pm$%.4f" % (label, r["s"], r["err"]))
    ax.axhline(1, color="gray", linewidth=1)
    ax.axvline(current, color="red", linestyle="--", linewidth=1.5, label="loaddf scale %.3f" % current)
    s0 = fits[0][1]["s"]
    w = max(0.05, 6*max(r["err"] for _, r in fits))
    ax.set_xlim(s0 - w, s0 + w)
    ax.set_ylim(0, 10)
    ax.set_xlabel("PE Scale Factor $s$", fontsize=FONTSIZE)
    ax.set_ylabel("$\\Delta\\chi^2$ (profiled)", fontsize=FONTSIZE)
    ax.set_title("%s flash PE scale fit" % DET_LABEL[det], fontsize=FONTSIZE)
    ax.legend(fontsize=FONTSIZE - 4)
    style(ax)
    return fig


# ============================================================
# Main
# ============================================================
def run_detector(det, d, args, log):
    bins = BINS_PE[det]
    current = loaddf.MC_PE_SCALE[det]
    stage = SCHEMA[args.schema]["stage"]
    tag = det.replace(" ", "-")

    pe_min = FIT_MIN_PE[det]
    fit_max = Fit(d, "flash_maxpe", bins, args.norm_unc, pe_min)
    r_max = fit_max.result()
    r_ref = Fit(d, "flash_maxpe", bins, 0., pe_min, mc_stat=False).result()   # the notebook's fit
    r_sum = Fit(d, "flash_sumpe", bins, args.norm_unc, pe_min).result()

    s_boot = fit_max.bootstrap(args.nboot, np.random.default_rng(SEED))
    boot_err = np.std(s_boot)
    flash_syst = abs(r_max["s"] - r_sum["s"])
    total = np.hypot(boot_err, flash_syst)

    log("max-flash fit:  s = %.4f +- %.4f (parabolic dchi2=1, stat+norm)  chi2/ndof = %.1f/%d  alpha = %+.3f" % (
        r_max["s"], r_max["err"], r_max["chi2"], r_max["ndof"], r_max["alpha"]))
    log("  notebook fit (no norm unc, data stat): s = %.4f +- %.4f  chi2/ndof = %.1f/%d" % (
        r_ref["s"], r_ref["err"], r_ref["chi2"], r_ref["ndof"]))
    log("sum-PE fit:     s = %.4f +- %.4f  chi2/ndof = %.1f/%d  alpha = %+.3f" % (
        r_sum["s"], r_sum["err"], r_sum["chi2"], r_sum["ndof"], r_sum["alpha"]))
    log("bootstrap (N=%d): sigma_s = %.4f" % (args.nboot, boot_err))
    log("flash-selection syst |s_max - s_sum| = %.4f" % flash_syst)
    log("FINAL %s: s = %.4f +- %.4f (boot) +- %.4f (flash sel.) = +- %.4f (total);  loaddf has %.3f (%+.1f sigma)" % (
        det, r_max["s"], boot_err, flash_syst, total, current, (current - r_max["s"])/total))

    fmt = "PE scale = %.3f $\\pm$ %.3f\n$\\chi^2$/ndof = %.1f/%d\nfit norm. $\\alpha$ = %+.3f"
    save(ratio_plot(d, det, "flash_maxpe", bins, r_max["s"], args.norm_unc, stage, alpha=r_max["alpha"],
                    text=fmt % (r_max["s"], total, r_max["chi2"], r_max["ndof"], r_max["alpha"])),
         args.plotdir, tag + "_flash_totalpe")
    cur_chi2 = fit_max.chi2()[0][np.argmin(np.abs(S_GRID - current))]
    save(ratio_plot(d, det, "flash_maxpe", bins, current, args.norm_unc, stage,
                    text="PE scale = %.3f (loaddf)\nprofiled $\\chi^2$/ndof = %.1f/%d" % (current, cur_chi2, r_max["ndof"]),
                    title="(current scale)"),
         args.plotdir, tag + "_flash_totalpe_current")
    save(ratio_plot(d, det, "flash_sumpe", bins, r_sum["s"], args.norm_unc, stage, alpha=r_sum["alpha"],
                    text=fmt % (r_sum["s"], r_sum["err"], r_sum["chi2"], r_sum["ndof"], r_sum["alpha"]),
                    title="(sum PE)"),
         args.plotdir, tag + "_flash_sumpe")
    save(plot_bootstrap(det, s_boot, r_max["s"], boot_err), args.plotdir, tag + "_flash_scale_bootstrap")
    save(plot_scan(det, [("max PE, %d%% norm" % round(100*args.norm_unc), r_max),
                         ("max PE, no norm (notebook)", r_ref),
                         ("sum PE, %d%% norm" % round(100*args.norm_unc), r_sum)], current),
         args.plotdir, tag + "_flash_scale_chi2scan")

    return dict(detector=det, schema=args.schema, s=r_max["s"], err_dchi2=r_max["err"],
                err_boot=boot_err, err_flashsel=flash_syst, err_total=total,
                alpha=r_max["alpha"], chi2=r_max["chi2"], ndof=r_max["ndof"],
                s_sumpe=r_sum["s"], alpha_sumpe=r_sum["alpha"],
                s_nonorm=r_ref["s"], err_nonorm=r_ref["err"], chi2_nonorm=r_ref["chi2"],
                s_loaddf=current, chi2_loaddf=cur_chi2,
                POT=d["norm"]["POT"], OFF_w=d["norm"]["OFF_w"])


def write_scale_summary(base, log):
    """trigger_scale_summary.tsv / trigger_scale_table.tex in `base` from the
    flashonly (all-trigger, CV) and gumple (Contained) fits, if both have run."""
    tsvs = {sch: os.path.join(base, sch, "flashpe_fit.tsv") for sch in ("flashonly", "gumple")}
    if not all(os.path.exists(f) for f in tsvs.values()):
        return
    a = pd.read_csv(tsvs["flashonly"], sep="\t").set_index("detector")
    c = pd.read_csv(tsvs["gumple"], sep="\t").set_index("detector")
    rows, tex = [], []
    for det in DETECTORS:
        if det not in a.index or det not in c.index:
            continue
        run = {v: k for k, v in loaddf.RUN_DETECTOR.items()}[det]
        unc = abs(a.s[det] - c.s[det])
        rows.append(dict(detector=det, s_alltrig=a.s[det], err_alltrig=a.err_total[det],
                         s_contained=c.s[det], err_contained=c.err_total[det], trig_unc=unc,
                         threshold=mc.TRIG_THRESHOLD[run],
                         loaddf_scale=loaddf.MC_PE_SCALE[det], loaddf_unc=loaddf.MC_PE_SCALE_UNC[det]))
        tex.append("        %s & $%.3f \\pm %.3f$ & $%.3f \\pm %.3f$ & %.3f & %d \\\\" % (
            DET_LABEL[det] + (" Run 1" if det == "SBND" else ""), a.s[det], a.err_total[det],
            c.s[det], c.err_total[det], unc, mc.TRIG_THRESHOLD[run]))
        if (round(a.s[det], 3), round(unc, 3)) != (loaddf.MC_PE_SCALE[det], loaddf.MC_PE_SCALE_UNC[det]):
            log("WARNING %s: loaddf MC_PE_SCALE/_UNC = %.3f/%.3f, fits give %.3f/%.3f -- update loaddf "
                "and the technote's tab:trigger-systematic" % (
                    det, loaddf.MC_PE_SCALE[det], loaddf.MC_PE_SCALE_UNC[det], a.s[det], unc))
    pd.DataFrame(rows).to_csv(os.path.join(base, "trigger_scale_summary.tsv"), sep="\t",
                              index=False, float_format="%.6g")
    with open(os.path.join(base, "trigger_scale_table.tex"), "w") as f:
        f.write("        %% Sample & All-trigger scale (CV) & Contained scale & Trigger syst. & Threshold [PE]\n")
        f.write("\n".join(tex) + "\n")
    log("trigger systematic summary: %s" % os.path.join(base, "trigger_scale_summary.tsv"))
    log("\n".join(tex))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--schema", required=True, choices=sorted(SCHEMA))
    ap.add_argument("--df-dir", default=None)
    ap.add_argument("--plotdir", default=None)
    ap.add_argument("--detectors", nargs="+", default=DETECTORS, choices=DETECTORS)
    ap.add_argument("--nboot", type=int, default=1000)
    ap.add_argument("--norm-unc", type=float, default=0.20)
    ap.add_argument("--full-onbeam", action="store_true")
    args = ap.parse_args(argv)
    df_dir = args.df_dir or SCHEMA[args.schema]["df_dir"]
    df_dir = df_dir if df_dir.endswith("/") else df_dir + "/"
    args.plotdir = args.plotdir or os.path.join(PLOTBASE, args.schema)

    for sub in ("pdf", "png"):
        os.makedirs(os.path.join(args.plotdir, sub), exist_ok=True)
    summary = []

    def log(msg):
        print(msg, flush=True)
        summary.append(str(msg))

    log("Flash PE scale fit, %s stage, schema=%s" % (SCHEMA[args.schema]["stage"], args.schema))
    log("DF_DIR: %s" % df_dir)
    log("norm. unc. on nu MC: %.0f%%; N_BOOT = %d; fit from the flash-PE threshold %s; s range %s" % (
        100*args.norm_unc, args.nboot, FIT_MIN_PE, S_RANGE))
    log("MC CV weights: %s" % ("reweight_aFF (wgt table)" if SCHEMA[args.schema]["reweight_aFF"]
                               else "none (no wgt table; cvwgt = 1)"))
    log("flash-selection syst uses flash_sumpe in place of the notebook's all-flash fit; "
        "the notebook's flash-time plots are not made (no flash times in any schema)")

    rows = []
    for det in args.detectors:
        log("\n==== %s" % det)
        d = load_detector(det, df_dir, args.schema, args.full_onbeam, log)
        rows.append(run_detector(det, d, args, log))
        del d

    tsv = os.path.join(args.plotdir, "flashpe_fit.tsv")
    pd.DataFrame(rows).to_csv(tsv, sep="\t", index=False, float_format="%.6g")
    log("\nfigures (pdf+png) in %s; table %s" % (args.plotdir, tsv))
    write_scale_summary(os.path.dirname(os.path.normpath(args.plotdir)), log)
    with open(os.path.join(args.plotdir, "summary.txt"), "w") as f:
        f.write("\n".join(summary) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
