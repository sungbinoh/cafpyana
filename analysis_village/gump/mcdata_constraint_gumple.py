"""MC/data comparison with a conditioned (Gaussian) constraint (GUMPLE).

GUMPLE (sbn-rewgted-19+) port of mcdata_constraint.py, built on
mcdata_comparison_gumple (gumple_cuts selection, -19 file layout, renamed
columns). Two additions over the reCAF version:

  * This script now also owns the plain MC/data comparisons of the variables
    REMOVED from mcdata_comparison_gumple -- mu_p, nu_E_calo, nu_E_ccqe and
    nu_E_frac_diff -- in BOTH area and absolute normalization, at the chosen
    selection stage (their plot specs come from mdc.CONSTRAINT_VARS).
  * The proton multiplicity is selected with PROTON_SEL / --proton-sel
    ("1p" -> the "PID 1p" stage (n_pfp == 2), ">1p" -> "PID Np" (n_pfp > 2),
    ">=1p" -> multiplicity-inclusive), as in SignalBoxSystematics-GUMPLE.
    BLINDING: only "1p" may currently be run -- anything else raises
    BlindedError before any data is touched.

Runs the selection for one detector and makes three constraint figures, in the
same stacked-MC + ratio style as mcdata_comparison_gumple.py:

  1. MC/data of the PLOT variable, unconstrained covariance
  2. MC/data of the CONSTRAINT variable
  3. MC/data of the PLOT variable, covariance conditioned on the constraint
     variable, with the constraint-shifted MC central value overlaid

The two variables are histogrammed into one concatenated bin vector
(syst.ConcatBins), so a single Systematic.cov call returns the joint block
covariance

    [[Cov(x,x), Cov(x,y)],
     [Cov(y,x), Cov(y,y)]]        x = constraint variable, y = plot variable

from which syst.conditional_constraint gives the Schur complement
Cyy - Cyx Cxx^-1 Cxy and the gain matrix K = Cyx Cxx^-1 that shifts the central
value by K (data_x - MC_x).

The Cxx block used for the constraint includes the DATA statistical covariance
of the constraint variable (with its cross terms to the plot variable): a finite
data sample cannot constrain systematics beyond its own statistical power.

Normalization is absolute (POT), not area: a constraint that moves the central
value is only meaningful on an absolute scale, and area normalization is not
defined over a concatenated two-variable vector. Absolute normalization exposes
the data rate, so this script is SBND-only for as long as ICARUS is blinded
(mcdata_comparison.ABSNORM_DETECTORS).

Example:

    python mcdata_constraint_gumple.py -d SBND --plot-var nu_E_calo --constraint-var nu_E_ccqe

Band conventions differ between the figures, deliberately, and are restated on
each one: figures 1 and 2 keep the mcdata_comparison convention (the band is the
MC systematic only; the data statistical uncertainty enters the chi2 alone),
while figure 3's band is the POST-CONSTRAINT TOTAL, which necessarily contains
the data statistics that limit the constraint. The summary file prints the
per-bin fractional uncertainties for all three so the shrinkage can be read
apples-to-apples.
"""

import argparse
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (_HERE, os.path.abspath(os.path.join(_HERE, "..", ".."))):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import syst
import mcdata_comparison_gumple as mdc


DEFAULT_PLOTDIR = "/Users/gputnam/Work/osc/cafpyana/plots-gumple-2026-08-26/constraint/"

# Proton-multiplicity selection (as SignalBoxSystematics-GUMPLE's PROTON_SEL):
# maps onto the mcdata_comparison_gumple PID stages. BLINDING: only "1p" may
# currently be run; every other value raises BlindedError.
PROTON_SEL = "1p"
PROTON_SEL_STAGE = {
    "1p":   "PID 1p",   # n_pfp == 2, the GUMP selection
    ">1p":  "PID Np",   # n_pfp > 2, the MAPLE selection
    ">=1p": None,       # multiplicity-inclusive; no stage column yet (blinded anyway)
}


class BlindedError(RuntimeError):
    """A blinded proton selection was requested."""


# The GUMP signal-box binning (historically nb/SignalBoxSystematics cell 53 /
# enrico_gump.ENU_BINS): fine through the flux peak, coarse in the tail. NB the
# GUMPLE signal-box notebook has since reverted to its coarse 8-bin array; this
# stays available here via --gump-bins.
GUMP_BINS = np.array([0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8,
                      0.85, 0.9, 0.95, 1.0, 1.25, 1.5])


def density_scale(bins, per_mev):
    """Per-bin factor turning raw bin contents into events per `per_mev` MeV.

    With variable-width bins the raw contents are not comparable between bins,
    so everything drawn is divided by the bin width and multiplied by a common
    reference width.
    """
    return (per_mev/1000.)/(np.asarray(bins)[1:] - np.asarray(bins)[:-1])


def ylabel_for(per_mev, POT):
    """Y-axis wording for fixed-width-normalized contents."""
    return "Events / %g MeV / %.1f$\\times 10^{19}$ POT" % (per_mev, POT/1e19)


def scale_plotdata(p, s, ylabel=None):
    """Rescale a plotdata dict by a per-bin factor, in place.

    The covariances transform as D C D, so the chi2 is invariant -- this is a
    display normalization and must not move any number it reports.
    """
    D = np.diag(s)
    p["NMC_total"] = p["NMC_total"]*s
    p["NMC_breakdown"] = p["NMC_breakdown"]*s          # (ncategory, nbin)
    p["NData"] = p["NData"]*s
    p["NDataErr"] = p["NDataErr"]*s
    for k in ("cov", "cov_w_stat"):
        if k in p and p[k] is not None:
            p[k] = D @ p[k] @ D
    if ylabel:
        p["ylabel"] = ylabel
    return p


def _strip_unit(label):
    """Axis label without its trailing unit, for use inside a title."""
    return re.sub(r"\s*\[[^\]]*\]\s*$", "", label)


def _plotvar_spec(detector, var):
    """(bins, label) for a variable, from the mcdata_comparison_gumple plot
    lists -- including the variables removed from that script's own
    comparisons (mdc.CONSTRAINT_VARS)."""
    if var in mdc.CONSTRAINT_VARS:
        return mdc.CONSTRAINT_VARS[var]
    if var not in mdc.PLOTVARS:
        raise ValueError("unknown variable %r; choose from %s"
                         % (var, ", ".join(list(mdc.CONSTRAINT_VARS) + mdc.PLOTVARS)))
    i = mdc.PLOTVARS.index(var)
    return mdc.plot_bins(detector)[i], mdc.LABELS[i]


def removed_var_comparisons(detector, cut, plotdir, no_save, log):
    """The plain MC/data comparisons removed from mcdata_comparison_gumple:
    mu_p and the neutrino-energy estimators, area- AND absolute-normalized,
    at the chosen selection stage. Uses the same make_plot_data/_save_plots
    machinery (and file naming) mcdata_comparison used for them."""
    vars_ = list(mdc.CONSTRAINT_VARS)
    n = 0
    for areanorm, sfx in ((True, ""), (False, "_absnorm")):
        store = {cut: {}}
        for v in vars_:
            bins, label = mdc.CONSTRAINT_VARS[v]
            store[cut][v] = mdc.make_plot_data(
                v, bins, cut, "glob_scale", mdc.breakdown_mode, areanorm,
                mdc.mode_labels, mdc.mode_colors, label, cut, fillna=-1,
                det=detector)
        if not no_save:
            n += mdc._save_plots(store, [cut], vars_, plotdir, suffix=sfx)
        for v in vars_:
            c2, nd, pv = mdc.chi2_ndof_pvalue(store[cut][v])
            log("  %-16s %-9s chi2 %7.2f / %2d   p = %s"
                % (v, "shape" if areanorm else "absolute", c2, nd,
                   "%.3g" % pv if np.isfinite(pv) else "--"))
    log("Removed-variable comparisons: %d figures (area + absolute) at %r"
        % (n, cut))
    return n


def data_stat_cov_blocks(ONdf, OFFdf, OFF_w, blocks, fillna=-1.):
    """Joint statistical covariance of the (ON - OFF_w x OFF) data vector over
    an arbitrary list of concatenated blocks.

    `blocks` is a list of (var, bins, cut-column). The same data events can fill
    several blocks, so element (a, b) of block (i, j) counts the events in bin a
    of block i AND bin b of block j -- the 2-D histogram over the rows passing
    both blocks' cuts. Disjoint del_p slices share no events and so come out
    uncorrelated, while the two variables within one slice do not.
    """
    sizes = [len(b) - 1 for _, b, _ in blocks]
    edges = np.cumsum([0] + sizes)
    C = np.zeros((edges[-1], edges[-1]))

    for frame, w in ((ONdf, 1.0), (OFFdf, OFF_w**2)):
        vals = {}
        for v, _, _ in blocks:
            if v not in vals:
                vals[v] = np.asarray(frame[v].fillna(fillna))
        masks = [np.asarray(frame[c], dtype=bool) for _, _, c in blocks]
        for i, (vi, bi, _) in enumerate(blocks):
            for j, (vj, bj, _) in enumerate(blocks):
                m = masks[i] & masks[j]
                if not m.any():
                    continue
                h = np.histogramdd([vals[vi][m], vals[vj][m]], bins=[bi, bj])[0]
                C[edges[i]:edges[i+1], edges[j]:edges[j+1]] += w*h

    return 0.5*(C + C.T)


def data_stat_cov(ONdf, OFFdf, OFF_w, cut, varx, binsx, vary, binsy, fillna=-1.):
    """Joint statistical covariance of the (ON - OFF_w x OFF) data vector over
    the concatenated [x-bins, y-bins].

    The same data events fill both blocks, so the blocks are correlated: element
    (a, b) of block (i, j) is the number of events in bin a of var i AND bin b of
    var j -- the 2-D histogram -- with the off-beam sample entering at OFF_w^2
    (it is subtracted with weight OFF_w, so its Poisson variance scales as the
    square). For i == j this reduces to the familiar diagonal
    diag(N_ON + OFF_w^2 N_OFF) that mcdata_comparison uses for NDataErr.
    """
    nx, ny = len(binsx) - 1, len(binsy) - 1
    vs = [(varx, binsx), (vary, binsy)]
    C = np.zeros((nx + ny, nx + ny))
    edges = [0, nx, nx + ny]

    for frame, w in ((ONdf, 1.0), (OFFdf, OFF_w**2)):
        m = frame[cut]
        for i, (vi, bi) in enumerate(vs):
            for j, (vj, bj) in enumerate(vs):
                h = np.histogramdd([frame.loc[m, vi].fillna(fillna),
                                    frame.loc[m, vj].fillna(fillna)],
                                   bins=[bi, bj])[0]
                C[edges[i]:edges[i+1], edges[j]:edges[j+1]] += w*h

    return 0.5*(C + C.T)


def draw_figures(plotdir, detector, cut, var_y, var_x, label_x,
                 pd_y, pd_x, pd_c, NMC_y_shift, no_save, log,
                 stem_suffix="", ax0_note=None):
    """The three figures for one (plot var, constraint var) pair.

    stem_suffix tags the output files for a del_p slice ("_lodp", ...) and
    ax0_note is the slice annotation drawn on the right, below the legend.
    Shared by run(), run_dp() and replot().
    """
    detname = detector.replace(" ", "-")
    cutstem = cut.replace(" ", "").replace(".", "").lower()
    sfx = stem_suffix

    # the constraint is called out on the plot rather than in the title
    constraint_note = "Uncertainty constrained by %s" % _strip_unit(label_x)

    figs = [
        ("%s_%s_%s_unconstrained%s" % (detname, cutstem, var_y, sfx), pd_y, None),
        ("%s_%s_%s_constraintvar%s" % (detname, cutstem, var_x, sfx), pd_x, None),
        ("%s_%s_%s_constrained_by_%s%s" % (detname, cutstem, var_y, var_x, sfx), pd_c,
         {"N": NMC_y_shift, "label": "Const. MC", "color": "black",
          "linestyle": "--"}),
    ]

    n_saved = 0
    for stem, pdata, overlay in figs:
        fig, ax0, ax1 = mdc.ratio_plot(plt, pdata, overlay=overlay)
        # right-hand side of the main panel, below the stack legend; the
        # constraint note drops a line when a del_p annotation takes the top slot
        # Anchor the notes just under the stack legend by measuring it, rather
        # than to a fixed fraction: the legend's height depends on how many
        # breakdown categories are populated. zorder beats the legend, and the
        # white backing keeps them legible where they cross the top of the stack.
        ynote = 0.54
        leg = ax0.get_legend()
        if leg is not None:
            fig.canvas.draw()   # legend has no extent until the figure is laid out
            bb = leg.get_window_extent().transformed(ax0.transAxes.inverted())
            ynote = bb.y0 - 0.03
        for txt, size, color in ((ax0_note, mdc.FONTSIZE-2, "black"),
                                 (constraint_note if overlay is not None else None,
                                  mdc.FONTSIZE-4, "gray")):
            if not txt:
                continue
            ax0.text(0.95, ynote, txt, verticalalignment="top",
                     horizontalalignment="right", fontsize=size, color=color,
                     transform=ax0.transAxes, zorder=6,
                     bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.5))
            ynote -= 0.085
        if not no_save:
            fig.savefig(os.path.join(plotdir, "pdf", stem + ".pdf"), bbox_inches="tight")
            fig.savefig(os.path.join(plotdir, "png", stem + ".png"), bbox_inches="tight")
            n_saved += 1
            log("Wrote %s.{pdf,png}" % stem)
        plt.close(fig)

    return n_saved


def replot(npz_path, args):
    """Rebuild the three figures from a saved .npz, with no reload.

    Everything the figures need is in the archive, so cosmetic changes cost
    seconds instead of a full detector load.
    """
    z = np.load(npz_path, allow_pickle=False)
    det, cut = str(z["detector"]), str(z["cut"])
    var_x, var_y = str(z["var_x"]), str(z["var_y"])
    label_x, label_y = str(z["label_x"]), str(z["label_y"])
    nx = len(z["bins_x"]) - 1
    POT = float(z["POT"])

    # archives store already-normalized contents, so the y-axis wording rides along
    ylab = str(z["ylabel"]) if "ylabel" in z.files else None

    def mk(bins, label, bkd, NMC, ND, NDErr, cov, cov_stat, chi2ndof):
        return {"det": det, "title": cut, "xlabel": label, "bins": bins,
                "areanorm": False, "breakdown_labels": mdc.mode_labels,
                "breakdown_colors": mdc.mode_colors, "NMC_breakdown": bkd,
                "NMC_total": NMC, "NData": ND, "NDataErr": NDErr,
                "cov": cov, "cov_w_stat": cov_stat, "ylabel": ylab,
                "chi2": float(chi2ndof[0]), "ndof": int(chi2ndof[1]), "POT": POT}

    Csyst, Cstat = z["Csyst"], z["Cstat"]
    # The del_p-sliced archives store the whole 6-block matrix, so index the
    # blocks for THIS slice rather than assuming a 2-block layout.
    ny = len(z["bins_y"]) - 1
    if Csyst.shape[0] == nx + ny:
        bx_, by_ = slice(0, nx), slice(nx, nx + ny)
    else:
        ns = Csyst.shape[0]//(nx + ny)          # slices per variable
        # the tag is the last component of e.g. "_gump_lodp"
        i = [t for t, _, _ in mdc.DP_SLICES].index(
            str(z["stem_suffix"]).rsplit("_", 1)[-1])
        bx_, by_ = slice(i*nx, (i+1)*nx), slice(ns*nx + i*ny, ns*nx + (i+1)*ny)

    cov_x, cov_y = Csyst[bx_, bx_], Csyst[by_, by_]
    covw_x, covw_y = cov_x + Cstat[bx_, bx_], cov_y + Cstat[by_, by_]
    cov_c = z["cov_cond"]

    # Recompute every chi2 from the stored blocks. Each figure then gets its
    # own -- the constraint variable used to inherit the plot variable's -- and
    # the numbers cannot drift from the covariance actually drawn.
    NMC_x, NMC_y, NMC_shift = z["NMC_x"], z["NMC_y"], z["NMC_y_shift"]
    ND_x, ND_y = z["NData_x"], z["NData_y"]
    c_x = mdc.f_chi2(NMC_x, ND_x, covw_x)
    c_y = mdc.f_chi2(NMC_y, ND_y, covw_y)
    c_c = mdc.f_chi2(NMC_shift, ND_y, cov_c)

    pd_y = mk(z["bins_y"], label_y, z["NMC_breakdown_y"], NMC_y,
              ND_y, z["NDataErr_y"], cov_y, covw_y, c_y)
    pd_x = mk(z["bins_x"], label_x, z["NMC_breakdown_x"], NMC_x,
              ND_x, z["NDataErr_x"], cov_x, covw_x, c_x)
    pd_c = mk(z["bins_y"], label_y, z["NMC_breakdown_y"], NMC_y,
              ND_y, z["NDataErr_y"], cov_c, cov_c, c_c)

    plotdir = args.plotdir
    os.makedirs(os.path.join(plotdir, "png"), exist_ok=True)
    os.makedirs(os.path.join(plotdir, "pdf"), exist_ok=True)

    def log(m):
        print(m, flush=True)

    # del_p-sliced archives carry their file tag and slice annotation
    sfx = str(z["stem_suffix"]) if "stem_suffix" in z.files else ""
    note = str(z["ax0_note"]) if "ax0_note" in z.files else None

    log("Replotting %s (%s, %s, plot=%s constrained by %s)%s"
        % (npz_path, det, cut, var_y, var_x, " " + note if note else ""))
    log("  chi2  %-12s %7.2f / %d" % (var_y, c_y[0], c_y[1]))
    log("  chi2  %-12s %7.2f / %d" % (var_x, c_x[0], c_x[1]))
    log("  chi2  %-12s %7.2f / %d" % ("constrained", c_c[0], c_c[1]))
    n = draw_figures(plotdir, det, cut, var_y, var_x, label_x,
                     pd_y, pd_x, pd_c, z["NMC_y_shift"], args.no_save, log,
                     stem_suffix=sfx, ax0_note=note)
    log("Saved %d figures (pdf+png) to %s" % (n, plotdir))
    return 0


def run(detector, args, ctx=None):
    plotdir = args.plotdir
    os.makedirs(plotdir, exist_ok=True)
    os.makedirs(os.path.join(plotdir, "png"), exist_ok=True)
    os.makedirs(os.path.join(plotdir, "pdf"), exist_ok=True)

    var_y, var_x = args.plot_var, args.constraint_var
    cut = args.cut
    if cut not in mdc.CUTNAMES:
        raise ValueError("unknown cut %r; choose from %s" % (cut, mdc.CUTNAMES))
    if var_x == var_y:
        raise ValueError("the constraint variable must differ from the plot variable")

    bins_y, label_y = _plotvar_spec(detector, var_y)
    bins_x, label_x = _plotvar_spec(detector, var_x)
    if args.gump_bins:
        bins_x = bins_y = GUMP_BINS
    nx, ny = len(bins_x) - 1, len(bins_y) - 1

    # ---- load + build systematics (absolute normalization) ----
    # a caller running several rounds passes the shared ctx to load only once
    if ctx is None:
        ctx = mdc.prepare_detector(detector, args)
    log = ctx["log"]
    df, ONdf, OFFdf = ctx["df"], ctx["ONdf"], ctx["OFFdf"]
    OFF_w, POT = ctx["OFF_w"], ctx["POT"]

    # ---- the plain comparisons removed from mcdata_comparison_gumple ----
    log("")
    log("Removed-variable MC/data comparisons (proton sel %r, stage %r):"
        % (args.proton_sel, cut))
    removed_var_comparisons(detector, cut, plotdir, args.no_save, log)

    log("")
    log("=" * 70)
    log("CONDITIONED CONSTRAINT")
    log("  proton selection    : %s" % args.proton_sel)
    log("  cut                 : %s" % cut)
    log("  plot variable   (y) : %s  (%d bins)" % (var_y, ny))
    log("  constraint var  (x) : %s  (%d bins)" % (var_x, nx))
    log("=" * 70)

    # The CCQE estimator (and any other variable with a restricted domain) is
    # NaN where its denominator is non-positive; with fillna=-1 those events
    # land in underflow and drop out of that histogram. They still count in the
    # other block, so the constraint sees fewer events than the selection does.
    sel = df[cut]
    for name, v in (("plot", var_y), ("constraint", var_x)):
        nnan = int(df.loc[sel, v].isna().sum())
        log("  %-10s variable %-12s: %d / %d selected MC rows are NaN (dropped "
            "from that histogram)" % (name, v, nnan, int(sel.sum())))

    # ---- figures 1 and 2: the two single-variable comparisons ----
    pd_y = mdc.make_plot_data(var_y, bins_y, cut, "glob_scale", mdc.breakdown_mode,
                              False, mdc.mode_labels, mdc.mode_colors, label_y,
                              cut, fillna=-1, det=detector)
    pd_x = mdc.make_plot_data(var_x, bins_x, cut, "glob_scale", mdc.breakdown_mode,
                              False, mdc.mode_labels, mdc.mode_colors, label_x,
                              cut, fillna=-1, det=detector)

    NMC_y, NMC_x = pd_y["NMC_total"], pd_x["NMC_total"]
    ND_y, ND_x = pd_y["NData"], pd_x["NData"]

    # ---- joint covariance over the concatenated [x-bins, y-bins] ----
    jbins = syst.ConcatBins([bins_x, bins_y])
    NCV = np.concatenate((NMC_x, NMC_y))
    log("\nBuilding the joint (%d + %d) systematic covariance..." % (nx, ny))
    Csyst = ctx["systematics_by_cut"][cut].cov([var_x, var_y], cut, jbins, NCV, fillna=-1)

    # cross-check: the joint blocks must reproduce the single-variable covariances
    for name, blk, ref in (("xx", Csyst[:nx, :nx], pd_x["cov"]),
                           ("yy", Csyst[nx:, nx:], pd_y["cov"])):
        d = np.abs(blk - ref).max()
        scale = max(np.abs(ref).max(), 1e-30)
        log("  block %s vs the single-variable covariance: max |diff| = %.3g (rel %.2g)"
            % (name, d, d/scale))
        if d > 1e-6*scale:
            log("  !! WARNING: joint block does not reproduce the 1-D covariance")

    # ---- data statistical covariance (with cross terms) ----
    Cstat = data_stat_cov(ONdf, OFFdf, OFF_w, cut, var_x, bins_x, var_y, bins_y)
    # consistency with the error bars drawn on figures 1 and 2
    for name, got, ref in (("x", np.sqrt(np.diag(Cstat[:nx, :nx])), pd_x["NDataErr"]),
                           ("y", np.sqrt(np.diag(Cstat[nx:, nx:])), pd_y["NDataErr"])):
        d = np.abs(got - ref).max()
        log("  data stat %s block vs NDataErr: max |diff| = %.3g" % (name, d))

    # ---- fixed-width normalization (variable-width binning only) ----
    # Applied here, before the constraint: the covariances transform as D C D
    # and the central values as N s, so the constraint and every chi2 are
    # unchanged -- doing it up front just keeps one consistent set of arrays.
    ylabel = None
    if args.per_mev:
        sx, sy = density_scale(bins_x, args.per_mev), density_scale(bins_y, args.per_mev)
        ylabel = ylabel_for(args.per_mev, POT)
        dscale = np.concatenate((sx, sy))
        D = np.diag(dscale)
        Csyst, Cstat = D @ Csyst @ D, D @ Cstat @ D
        scale_plotdata(pd_x, sx, ylabel)
        scale_plotdata(pd_y, sy, ylabel)
        NMC_y, NMC_x = pd_y["NMC_total"], pd_x["NMC_total"]
        ND_y, ND_x = pd_y["NData"], pd_x["NData"]
        log("\nNormalizing variable-width bins to events / %g MeV" % args.per_mev)

    Ctot = Csyst + Cstat

    # ---- condition the plot variable on the constraint variable ----
    cov_cond, K, cond_number = syst.conditional_constraint(Ctot, nx)
    log("\ncond(Cxx) = %.4g%s" % (cond_number,
        "  (ill-conditioned: pseudo-inverse used)" if not np.isfinite(cond_number)
        or cond_number > 1e12 else ""))

    NMC_y_shift = NMC_y + K @ (ND_x - NMC_x)
    # A Gaussian constraint can push a near-empty bin below zero; f_chi2 drops
    # non-positive bins, so this silently costs a degree of freedom.
    _nneg = int((NMC_y_shift <= 0).sum())
    if _nneg:
        log("!! %d bin(s) shifted to <= 0 and are dropped from the chi2 (bins %s)"
            % (_nneg, np.flatnonzero(NMC_y_shift <= 0).tolist()))

    # ---- chi2 before and after ----
    chi2_0, ndof_0 = mdc.f_chi2(NMC_y, ND_y, pd_y["cov_w_stat"])
    chi2_c, ndof_c = mdc.f_chi2(NMC_y_shift, ND_y, cov_cond)
    p0 = dict(pd_y, chi2=chi2_0, ndof=ndof_0)
    pc = dict(pd_y, chi2=chi2_c, ndof=ndof_c)
    _, _, pval_0 = mdc.chi2_ndof_pvalue(p0)
    _, _, pval_c = mdc.chi2_ndof_pvalue(pc)

    log("")
    log("chi2 (unconstrained, syst + data stat) : %.2f / %d   p = %.3g"
        % (chi2_0, ndof_0, pval_0))
    log("chi2 (constrained, shifted CV)         : %.2f / %d   p = %.3g"
        % (chi2_c, ndof_c, pval_c))

    # ---- fractional uncertainties, all three on the same footing ----
    with np.errstate(divide="ignore", invalid="ignore"):
        f_syst = np.sqrt(np.diag(pd_y["cov"]))/NMC_y
        f_tot = np.sqrt(np.diag(pd_y["cov_w_stat"]))/NMC_y
        f_cond = np.sqrt(np.diag(cov_cond))/NMC_y
    log("")
    log("Fractional uncertainty on %s per bin:" % var_y)
    log("  %-28s %s" % ("bin [low, high)", " ".join("%7.4f-%-7.4f" % (a, b)
                                                    for a, b in zip(bins_y[:-1], bins_y[1:]))))
    for nm, arr in (("unconstrained syst only", f_syst),
                    ("unconstrained syst+stat", f_tot),
                    ("constrained  syst+stat", f_cond)):
        log("  %-28s %s" % (nm, " ".join("%15.4f" % a for a in arr)))
    log("  %-28s %s" % ("shift / MC", " ".join("%15.4f" % a for a in
                                               (NMC_y_shift - NMC_y)/NMC_y)))

    if np.any(f_cond > f_tot + 1e-9):
        log("  !! WARNING: the constraint increased the uncertainty in some bin")

    # ---- figure 3 ----
    pd_c = dict(pd_y)
    pd_c["cov"] = cov_cond
    pd_c["cov_w_stat"] = cov_cond          # data stat already inside
    pd_c["chi2"], pd_c["ndof"] = chi2_c, ndof_c
    pd_c["title"] = cut

    tag = "_gump" if args.gump_bins else ""
    n_saved = draw_figures(plotdir, detector, cut, var_y, var_x, label_x,
                           pd_y, pd_x, pd_c, NMC_y_shift, args.no_save, log,
                           stem_suffix=tag)

    detname = detector.replace(" ", "-")
    if not args.no_save:
        npz = os.path.join(plotdir, "%s_constraint_%s_by_%s%s.npz"
                           % (detname, var_y, var_x, tag))
        np.savez(npz, Csyst=Csyst, Cstat=Cstat, cov_cond=cov_cond, K=K,
                 NMC_x=NMC_x, NMC_y=NMC_y, NMC_y_shift=NMC_y_shift,
                 NData_x=ND_x, NData_y=ND_y,
                 NDataErr_x=pd_x["NDataErr"], NDataErr_y=pd_y["NDataErr"],
                 bins_x=bins_x, bins_y=bins_y, POT=POT,
                 NMC_breakdown_x=pd_x["NMC_breakdown"],
                 NMC_breakdown_y=pd_y["NMC_breakdown"],
                 # metadata, so --replot-from can rebuild the figures alone
                 detector=np.array(detector), cut=np.array(cut),
                 var_x=np.array(var_x), var_y=np.array(var_y),
                 label_x=np.array(label_x), label_y=np.array(label_y),
                 stem_suffix=np.array(tag), ylabel=np.array(ylabel or ""),
                 chi2_unconstrained=np.array([chi2_0, ndof_0]),
                 chi2_constrained=np.array([chi2_c, ndof_c]))
        log("Wrote %s" % npz)

        summary = os.path.join(plotdir, "%s_constraint%s_summary.txt" % (detname, tag))
        with open(summary, "w") as f:
            f.write("\n".join(ctx["summary_lines"]) + "\n")
        print("Wrote %s" % summary, flush=True)

    log("\nSaved %d figures (pdf+png) to %s" % (n_saved, plotdir))
    return 0


def run_dp(detector, args, ctx=None):
    """As run(), but sliced into the del_p bins of mcdata_comparison.DP_SLICES.

    Builds ONE joint covariance over the six blocks

        [x_lo, x_mid, x_hi, y_lo, y_mid, y_hi]

    (x = constraint variable, y = plot variable), so the cross-slice
    correlations are carried explicitly, then applies the conditional
    constraint within each del_p slice: y_i constrained by x_i, using that
    slice's 2-block sub-matrix. Nine figures, three per slice.
    """
    plotdir = args.plotdir
    for sub in ("png", "pdf"):
        os.makedirs(os.path.join(plotdir, sub), exist_ok=True)

    var_y, var_x = args.plot_var, args.constraint_var
    cut = args.cut
    if cut not in mdc.CUTNAMES:
        raise ValueError("unknown cut %r; choose from %s" % (cut, mdc.CUTNAMES))
    if var_x == var_y:
        raise ValueError("the constraint variable must differ from the plot variable")

    bins_y, label_y = _plotvar_spec(detector, var_y)
    bins_x, label_x = _plotvar_spec(detector, var_x)
    if args.gump_bins:
        bins_x = bins_y = GUMP_BINS
    nx, ny = len(bins_x) - 1, len(bins_y) - 1
    slices = mdc.DP_SLICES
    ns = len(slices)

    # a caller running several rounds passes the shared ctx to load only once
    if ctx is None:
        ctx = mdc.prepare_detector(detector, args)
    log = ctx["log"]
    df, ONdf, OFFdf = ctx["df"], ctx["ONdf"], ctx["OFFdf"]
    OFF_w, POT = ctx["OFF_w"], ctx["POT"]

    log("")
    log("=" * 70)
    log("CONDITIONED CONSTRAINT, SLICED IN del_p")
    log("  cut                 : %s" % cut)
    log("  plot variable   (y) : %s  (%d bins)" % (var_y, ny))
    log("  constraint var  (x) : %s  (%d bins)" % (var_x, nx))
    log("  del_p slices        : %s" % ", ".join(
        "%s [%.1f, %.1f)" % (t, lo, hi) for t, lo, hi in slices))
    log("=" * 70)

    # block layout: the three constraint slices, then the three plot slices
    cuts6 = [mdc.dp_cutname(cut, t) for t, _, _ in slices]*2
    vars6 = [var_x]*ns + [var_y]*ns
    bins6 = [bins_x]*ns + [bins_y]*ns
    sizes = [len(b) - 1 for b in bins6]
    edges = np.cumsum([0] + sizes)
    blk = [slice(edges[i], edges[i+1]) for i in range(2*ns)]

    for c in cuts6[:ns]:
        log("  %-24s: %7d selected MC rows, %6d on-beam" %
            (c, int(df[c].sum()), int(ONdf[c].sum())))

    # ---- per-block plot data (no systematics: the covariance comes from the
    # single joint call below, so the universe loop runs once, not six times) ----
    nosyst = syst.SystematicList([])
    pds = []
    for i in range(2*ns):
        lbl = label_x if i < ns else label_y
        pds.append(mdc.make_plot_data(vars6[i], bins6[i], cuts6[i], "glob_scale",
                                      mdc.breakdown_mode, False, mdc.mode_labels,
                                      mdc.mode_colors, lbl, cut, fillna=-1,
                                      det=detector, syst=nosyst))

    NCV = np.concatenate([p["NMC_total"] for p in pds])

    log("\nBuilding the joint %d-block (%d bin) systematic covariance..."
        % (2*ns, edges[-1]))
    # systematics_by_cut is keyed by a single cut column; for the joint vector
    # the trigger SelectionSystematic needs a per-block cut list per universe.
    base = ctx["systematics_by_cut"][cut].systs[0]
    trig = syst.SelectionSystematic(df, [[c + updn for c in cuts6]
                                         for updn in ("_trig_up", "_trig_dn")])
    joint = syst.SystematicList([base, trig])
    Csyst = joint.cov(vars6, cuts6, syst.ConcatBins(bins6), NCV, fillna=-1)

    # cross-check one block against an independent single-variable covariance
    ref = ctx["systematics_by_cut"][cuts6[0]].cov(vars6[0], cuts6[0], bins6[0],
                                                  NCV[blk[0]], fillna=-1)
    d = np.abs(Csyst[blk[0], blk[0]] - ref).max()
    log("  block 0 vs the single-variable covariance: max |diff| = %.3g (rel %.2g)"
        % (d, d/max(np.abs(ref).max(), 1e-30)))
    if d > 1e-6*max(np.abs(ref).max(), 1e-30):
        log("  !! WARNING: joint block does not reproduce the 1-D covariance")

    Cstat = data_stat_cov_blocks(ONdf, OFFdf, OFF_w,
                                 list(zip(vars6, bins6, cuts6)))

    # ---- fixed-width normalization (see run(); chi2 and constraint invariant) ----
    ylabel = None
    if args.per_mev:
        bscale = [density_scale(b, args.per_mev) for b in bins6]
        ylabel = ylabel_for(args.per_mev, POT)
        D = np.diag(np.concatenate(bscale))
        Csyst, Cstat = D @ Csyst @ D, D @ Cstat @ D
        for p, s_ in zip(pds, bscale):
            scale_plotdata(p, s_, ylabel)
        NCV = np.concatenate([p["NMC_total"] for p in pds])
        log("Normalizing variable-width bins to events / %g MeV" % args.per_mev)

    Ctot = Csyst + Cstat

    # ---- constrain within each del_p slice ----
    n_saved = 0
    per_slice = {}
    for i, (tag, lo, hi) in enumerate(slices):
        pd_x, pd_y = pds[i], pds[ns + i]
        idx = np.concatenate([np.arange(edges[i], edges[i+1]),
                              np.arange(edges[ns+i], edges[ns+i+1])])
        Ci = Ctot[np.ix_(idx, idx)]
        cov_cond, K, cn = syst.conditional_constraint(Ci, nx)

        NMC_x, NMC_y = pd_x["NMC_total"], pd_y["NMC_total"]
        ND_x, ND_y = pd_x["NData"], pd_y["NData"]
        NMC_y_shift = NMC_y + K @ (ND_x - NMC_x)
        # A Gaussian constraint can push a near-empty bin below zero; f_chi2
        # drops non-positive bins, so this silently costs a degree of freedom.
        nneg = int((NMC_y_shift <= 0).sum())
        if nneg:
            log("  !! %d bin(s) shifted to <= 0 and are dropped from the chi2 "
                "(bins %s)" % (nneg, np.flatnonzero(NMC_y_shift <= 0).tolist()))

        # the unconstrained figures keep the mcdata_comparison band convention
        # (systematic only); the constrained one shows the post-constraint total
        for p, s_ in ((pd_x, blk[i]), (pd_y, blk[ns+i])):
            p["cov"] = Csyst[s_, s_]
            p["cov_w_stat"] = Csyst[s_, s_] + Cstat[s_, s_]
            p["chi2"], p["ndof"] = mdc.f_chi2(p["NMC_total"], p["NData"], p["cov_w_stat"])

        chi2_c, ndof_c = mdc.f_chi2(NMC_y_shift, ND_y, cov_cond)
        pd_c = dict(pd_y)
        pd_c["cov"] = cov_cond
        pd_c["cov_w_stat"] = cov_cond
        pd_c["chi2"], pd_c["ndof"] = chi2_c, ndof_c

        _, _, pv0 = mdc.chi2_ndof_pvalue(pd_y)
        _, _, pvc = mdc.chi2_ndof_pvalue(pd_c)
        with np.errstate(divide="ignore", invalid="ignore"):
            f_tot = np.sqrt(np.diag(pd_y["cov_w_stat"]))/NMC_y
            f_cond = np.sqrt(np.diag(cov_cond))/NMC_y

        log("")
        log("--- %s  (%.1f <= del_p < %.1f) ---" % (tag, lo, hi))
        log("  cond(Cxx) = %.4g" % cn)
        log("  chi2 unconstrained : %7.2f / %d   p = %.3g" % (pd_y["chi2"], pd_y["ndof"], pv0))
        log("  chi2 constrained   : %7.2f / %d   p = %.3g" % (chi2_c, ndof_c, pvc))
        log("  frac unc syst+stat : %s" % " ".join("%6.4f" % a for a in f_tot))
        log("  frac unc constrained: %s" % " ".join("%6.4f" % a for a in f_cond))
        if np.any(f_cond > f_tot + 1e-9):
            log("  !! WARNING: the constraint increased the uncertainty in some bin")

        note = mdc.dp_annotation(lo, hi)
        sfx = ("_gump" if args.gump_bins else "") + "_" + tag
        n_saved += draw_figures(plotdir, detector, cut, var_y, var_x, label_x,
                                pd_y, pd_x, pd_c, NMC_y_shift, args.no_save, log,
                                stem_suffix=sfx, ax0_note=note)

        per_slice[tag] = dict(pd_x=pd_x, pd_y=pd_y, cov_cond=cov_cond, K=K,
                              NMC_y_shift=NMC_y_shift, chi2_c=(chi2_c, ndof_c),
                              note=note, sfx=sfx)

    detname = detector.replace(" ", "-")
    if not args.no_save:
        for tag, d_ in per_slice.items():
            px, py = d_["pd_x"], d_["pd_y"]
            npz = os.path.join(plotdir, "%s_constraint_%s_by_%s%s.npz"
                               % (detname, var_y, var_x, d_["sfx"]))
            np.savez(npz, Csyst=Csyst, Cstat=Cstat, cov_cond=d_["cov_cond"], K=d_["K"],
                     NMC_x=px["NMC_total"], NMC_y=py["NMC_total"],
                     NMC_y_shift=d_["NMC_y_shift"],
                     NData_x=px["NData"], NData_y=py["NData"],
                     NDataErr_x=px["NDataErr"], NDataErr_y=py["NDataErr"],
                     bins_x=bins_x, bins_y=bins_y, POT=POT,
                     NMC_breakdown_x=px["NMC_breakdown"],
                     NMC_breakdown_y=py["NMC_breakdown"],
                     detector=np.array(detector), cut=np.array(cut),
                     var_x=np.array(var_x), var_y=np.array(var_y),
                     label_x=np.array(label_x), label_y=np.array(label_y),
                     stem_suffix=np.array(d_["sfx"]), ax0_note=np.array(d_["note"]),
                     ylabel=np.array(ylabel or ""),
                     chi2_unconstrained=np.array([py["chi2"], py["ndof"]]),
                     chi2_constrained=np.array(d_["chi2_c"]))
            log("Wrote %s" % npz)

        summary = os.path.join(plotdir, "%s_constraint_dp%s_summary.txt"
                               % (detname, "_gump" if args.gump_bins else ""))
        with open(summary, "w") as f:
            f.write("\n".join(ctx["summary_lines"]) + "\n")
        print("Wrote %s" % summary, flush=True)

    log("\nSaved %d figures (pdf+png) to %s" % (n_saved, plotdir))
    return 0


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-d", "--detector", default="SBND",
                   help="Detector to run (default SBND). Absolute normalization "
                        "restricts this to %s." % mdc.ABSNORM_DETECTORS)
    p.add_argument("--plot-var", default="nu_E_calo",
                   help="Variable to plot and constrain (default nu_E_calo)")
    p.add_argument("--constraint-var", default="nu_E_ccqe",
                   help="Variable to constrain WITH (default nu_E_ccqe)")
    p.add_argument("--proton-sel", default=PROTON_SEL, choices=("1p", ">1p", ">=1p"),
                   help="Proton-multiplicity selection (default %r). BLINDED: "
                        "only '1p' may currently be run; anything else raises "
                        "BlindedError." % PROTON_SEL)
    p.add_argument("--cut", default=None,
                   help="Selection stage (default: the stage for --proton-sel, "
                        "i.e. %r for '1p')" % PROTON_SEL_STAGE["1p"])
    p.add_argument("--df-dir", default=mdc.DEFAULT_DF_DIR, help="Directory of input .df files")
    p.add_argument("--plotdir", default=DEFAULT_PLOTDIR, help="Output directory")
    p.add_argument("--nproc", type=int, default=1, help="Unused; kept for prepare_detector")
    p.add_argument("--no-dirt", action="store_true", help="Skip the dirt sample")
    p.add_argument("--no-save", action="store_true", help="Build plots but write nothing")
    p.add_argument("--dp-bins", action="store_true",
                   help="Run ONLY the del_p-sliced round: slice the comparison "
                        "into the del_p bins of mcdata_comparison.DP_SLICES (%s), "
                        "build one joint 6-block covariance and constrain within "
                        "each slice. Produces 9 figures tagged %s. By default "
                        "both the inclusive and the sliced rounds run (one "
                        "detector load); see also --no-dp."
                        % (", ".join("%s [%.1f,%.1f)" % s for s in mdc.DP_SLICES),
                           "/".join(t for t, _, _ in mdc.DP_SLICES)))
    p.add_argument("--no-dp", action="store_true",
                   help="Run only the inclusive round (skip the del_p slices)")
    p.add_argument("--gump-bins", action="store_true",
                   help="Use the GUMP signal-box binning (%s) for both variables "
                        "instead of the mcdata_comparison binning. Implies "
                        "--per-mev 50, and tags the outputs '_gump'."
                        % " ".join("%g" % b for b in GUMP_BINS))
    p.add_argument("--per-mev", type=float, default=None, metavar="MEV",
                   help="Normalize bin contents to events per MEV MeV. Needed "
                        "with variable-width bins, where raw contents are not "
                        "comparable between bins. Default 50 with --gump-bins, "
                        "off otherwise (equal-width bins need no rescaling).")
    p.add_argument("--replot-from", default=None, metavar="NPZ",
                   help="Rebuild the three figures from a previously written .npz "
                        "instead of reloading the detector (seconds, not minutes). "
                        "All other options except --plotdir/--no-save are ignored.")
    args = p.parse_args(argv)

    # BLINDING: everything except the 1p selection is blinded. Gate before any
    # other action (including --replot-from) so a blinded selection can never
    # produce output.
    if args.proton_sel != "1p":
        raise BlindedError(
            "Blinded: proton selection %r is blinded; only '1p' may currently "
            "be run" % args.proton_sel)

    # resolve the selection stage from the proton selection
    if args.cut is None:
        args.cut = PROTON_SEL_STAGE[args.proton_sel]

    # the GUMP binning is variable-width, so it needs a reference width
    if args.gump_bins and args.per_mev is None:
        args.per_mev = 50.

    if args.replot_from:
        return args

    if args.detector not in mdc.DETECTORS:
        p.error("unknown detector %r (choose from %s)" % (args.detector, mdc.DETECTORS))

    # Absolute normalization exposes the data rate; ICARUS is still blinded.
    if args.detector not in mdc.ABSNORM_DETECTORS:
        p.error("this comparison is absolute-normalized and so is only allowed for "
                "%s (requested %r)." % (mdc.ABSNORM_DETECTORS, args.detector))

    # consumed by mcdata_comparison.prepare_detector
    args.absnorm = True
    args.blinded = False
    args.skip_pdg = True

    return args


def main(argv=None):
    args = parse_args(argv)
    if args.replot_from:
        return replot(args.replot_from, args)
    if args.dp_bins:
        return run_dp(args.detector, args)
    # default: the inclusive round AND the del_p-sliced round (with its
    # per-slice constraint), sharing one detector load
    ctx = mdc.prepare_detector(args.detector, args)
    rc = run(args.detector, args, ctx)
    if not args.no_dp:
        rc = run_dp(args.detector, args, ctx) or rc
    return rc


if __name__ == "__main__":
    sys.exit(main())
