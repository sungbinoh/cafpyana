"""Track-splitting correction at the SBND APA gap.

SBND counterpart of TrackSplittingCorrection_GUMPLE.py (ICARUS Z=0 gap and
cathodes). Same measurement, same estimator, same two selection stages; the
counting, angular and stage-ratio machinery is imported from that module rather
than copied. See its docstring for the definition of

    f = [N_data(W) - a*N_MC(W)] / (a * D)

(W the window around the plane, a the sideband MC->data normalization, D the
at-risk MC population: plane-crossers + muons already ending in W).

The plane
---------
Each SBND anode plane is two APA frames side by side in z. In sbndcode
(sbndcode/Geometry/gdml/sbnd_v02_06_base.gdml) the frames are
dimAPAFrame_W_mm = 2565.677 wide, placed at z = whereTPC_Z +/- coordAPA_DZ with
coordAPA_DZ = dimAPAFrame_W_mm/2 + dimAPASpacing_mm/2 (dimAPASpacing_mm = 38),
i.e. symmetric about the TPC z-centre with a 3.8 cm gap between them. The TPC is
dimDetector_Z = 5010 mm long and starts at z ~ 0, so the gap is a vertical plane
(normal along z) at z = 250 cm, in both the East and the West TPC. This is the
same z = 250 at which gumple_cuts.SBNDFVCuts splits into lowYZ / highYZ.

The FV notch
------------
For z > 250 the SBND FV stops at y < 100 (highYZ; with the 10 cm inset, 90),
while for z < 250 it runs to y < 200. Above y = 90 the FV edge therefore sits
exactly on the APA gap and a muon there cannot cross it. The measurement is
restricted to the region where the FV is continuous across the gap: slc_vtx_y
and mu_end_y both < Y_MAX = 90 cm. The restriction is ANDed into every stage
mask, so it applies to the window, the at-risk denominator and the
normalization sideband alike.

Regions
-------
    both   both TPCs in one fit (the combined comparison)
    east   East TPC, slc_vtx_x < 0
    west   West TPC, slc_vtx_x > 0

TPC assignment follows mcdata_comparison_gumple.cryo_mask (slice vertex x). The
production cathode cut (cut_cathode, in FV) already removes cathode crossers, so
the vertex side is the end side. alpha is refitted per region; the window
half-width is chosen once, from the both-TPC simple-stage scan, and reused for
every region and stage so the per-TPC numbers and the stage ratios integrate
over the same window. As a cross-check the East and West results are also
combined as two independent samples (sum excess / sum a*D, as the ICARUS
Run2+Run4 combination) and compared with the single both-TPC fit.

Angular dependence: per region and stage, 4 equal-population quantile bins of
theta_z = arccos|mu_dir_z| over the at-risk MC, alpha refitted per bin
(ts.run_angular). The per-TPC overlay keeps each region's own quantile edges.

Data: the full SBND_SpringBNBData_FullOnBeam.df on-beam stream alone (never
summed with FixedDev/RollingDev -- they overlap). Beam- and data-quality cuts
are passed explicitly on the ON/OFF loads so the event selection matches the
POT from mdc.compute_norm (loaddf.load_one defaults them off). MC: all
SBNDMCCV_* shards.

    python analysis_village/gump/TrackSplittingCorrection_SBND_GUMPLE.py

NB: `FV` must stay a module-level function in `__main__` (loaddf cache key,
see TrackSplittingCorrection_GUMPLE.py).
"""

import argparse
import glob
import json
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (_HERE, os.path.abspath(os.path.join(_HERE, "..", ".."))):
    if _p not in sys.path:
        sys.path.insert(0, _p)

sys.path.insert(0, os.path.join(_HERE, "..", "gumple"))
import gumple_cuts as gc
import loaddf
import mcdata_comparison_gumple as mdc
import TrackSplittingCorrection_GUMPLE as ts

DETECTOR = "SBND"

DEFAULT_DF_DIR = "/Users/gputnam/Work/osc/sbn-rewgted-25/"
DEFAULT_PLOTDIR = ("/Users/gputnam/Work/osc/cafpyana/"
                   "plots-gumple-2026-10-24-rewgted25-sbnd-tracksplit/")

FONTSIZE = ts.FONTSIZE

# ---- geometry (see the module docstring) ----
APA_GAP_Z = 250.0
DIM = "z"
POINT_NAME = "APA Gap"
POINT_KEY = "apa"
# upper y edge of the region where the FV is continuous across the APA gap
Y_MAX = 90.0

# (key, label, side): side is the sign of slc_vtx_x, None for both TPCs. The
# first region/stage pair is the one the common window half-width comes from.
REGIONS = [
    ("both", "Both TPCs", None),
    ("east", "East TPC",  -1),
    ("west", "West TPC",  +1),
]
REGION_LABEL = dict((k, l) for (k, l, _) in REGIONS)
REGION_COLORS = {"both": "#315031", "east": "#1e3f54", "west": "#d54c28",
                 "east+west": "#7a5c9e"}

STAGES = ts.STAGES
STAGE_LABEL = ts.STAGE_LABEL
RATIO_STAGES = ts.RATIO_STAGES


# ============================================================
# Preselection and region masks
# ============================================================
# NB: keep this at module scope, named FV, with a body identical to
# mcdata_comparison_gumple.FV -- loaddf's cache key depends on its __qualname__.
def FV(df):
    # gumple_cuts equivalent of the old slcfv & mufv & pfv & cathode chain
    # (identical body to mcdata_comparison_gumple.FV -- cache-key contract).
    return gc.sanity_cut(df) & gc.slcfv_cut(df) & df.cut_contained & df.cut_cathode


def continuous_fv(df):
    """Below the FV notch: the FV is continuous across z = 250 here."""
    return ((df.slc_vtx_y < Y_MAX) & (df.mu_end_y < Y_MAX)).to_numpy()


def region_mask(df, side):
    if side is None:
        return np.ones(len(df), dtype=bool)
    return (np.sign(df.slc_vtx_x) == side).to_numpy()


# ============================================================
# Loading
# ============================================================
def mc_shards(df_dir):
    files = glob.glob(os.path.join(df_dir, "SBNDMCCV_*.df"))
    return sorted(files, key=lambda f: int(re.findall(r"_(\d+)\.df$", f)[0]))


def load_all(files, pot, include_dirt, log):
    """ts.load_all with this module's FV, and BQ/DQ on the data loads."""
    mc_files = files["MC_FILES"]
    log("loading MC (%d shards): %s"
        % (len(mc_files), ", ".join(os.path.basename(f) for f in mc_files)))
    mcdf, _, mcpot = loaddf.loadl(mc_files, njob=min(len(mc_files), 10),
                                  detector=DETECTOR, preselection=FV,
                                  reweight_aFF=True, include_syst=False)
    mcdf = ts._slim(mcdf)
    loaddf.scale_pot(mcdf, mcpot, pot)
    mcdf["dirt"] = False
    log("  MC: %d rows, POT %.4g (scale %.4g)" % (len(mcdf), mcpot, pot/mcpot))

    if include_dirt and files["DIRT_FILES"]:
        dirt_files = files["DIRT_FILES"]
        log("loading dirt: %s" % ", ".join(os.path.basename(f) for f in dirt_files))
        dirt, _, dirtpot = loaddf.loadl(dirt_files, njob=min(len(dirt_files), 10),
                                        detector=DETECTOR, preselection=FV,
                                        include_syst=False)
        dirt = ts._slim(dirt)
        loaddf.scale_pot(dirt, dirtpot, pot)
        dirt["dirt"] = True
        log("  dirt: %d rows, POT %.4g (scale %.4g)" % (len(dirt), dirtpot, pot/dirtpot))
        mcdf = pd.concat([mcdf, dirt], ignore_index=True)

    # beam_quality/data_quality: load_one defaults them off, while
    # mdc.compute_norm's POT applies both -- pass them so events and POT match
    log("loading ON: %s" % os.path.basename(files["ONBEAM"]))
    ondf, _, _ = loaddf.load(files["ONBEAM"], load_truth=False, include_syst=False,
                             detector=DETECTOR, preselection=FV, match_Enu=False,
                             beam_quality=True, data_quality=True)
    ondf = ts._slim(ondf)

    log("loading OFF: %s" % ", ".join(os.path.basename(f) for f in files["OFFBEAM_FILES"]))
    offs = [loaddf.load(f, load_truth=False, include_syst=False, detector=DETECTOR,
                        preselection=FV, offbeampot=True, match_Enu=False,
                        beam_quality=True, data_quality=True)
            for f in files["OFFBEAM_FILES"]]
    offdf = pd.concat([ts._slim(o[0]) for o in offs], ignore_index=True)

    return (mcdf.reset_index(drop=True), ondf.reset_index(drop=True),
            offdf.reset_index(drop=True))


# ============================================================
# Per region x stage
# ============================================================
def run_one(region, stage, frames, masks, norm_masks, off_w, pot, w_star, args, log):
    """Normalization, counting, plots and angular scan for one region and stage.

    `masks` already carry stage & y-restriction & region. Returns
    (result row, scan frame, angular rows, chosen w).
    """
    rkey, rlabel, _ = region
    skey, slabel, _ = stage
    mcdf, ondf, offdf = frames
    plotdir, dosave = args.plotdir, not args.no_save

    log("")
    log("=" * 70)
    log("%s  --  %s" % (rlabel, slabel))
    log("=" * 70)
    for nm, f in (("MC", mcdf), ("ON", ondf), ("OFF", offdf)):
        log("  %s selected: %d / %d" % (nm, int(masks[nm].sum()), len(f)))

    smc = mcdf[masks["MC"]].reset_index(drop=True)
    son = ondf[masks["ON"]].reset_index(drop=True)
    soff = offdf[masks["OFF"]].reset_index(drop=True)

    mc_r = ts._wsum(mcdf, masks["MC"] & norm_masks["MC"])
    on_r = float((masks["ON"] & norm_masks["ON"]).sum())
    off_r = float((masks["OFF"] & norm_masks["OFF"]).sum())
    alpha = (on_r - off_r*off_w)/mc_r
    alpha_var = (on_r + off_r*off_w**2)/mc_r**2
    log("normalization region (|mu_end_z - %.0f| >= %.0f cm):" % (APA_GAP_Z, ts.NORM_EXCLUDE_W))
    log("  N_ON = %.0f  N_OFF = %.0f  N_MC = %.1f  ->  alpha = %.4f +/- %.4f"
        % (on_r, off_r, mc_r, alpha, np.sqrt(alpha_var)))

    counter = ts.PointCounter(smc, son, soff, off_w, alpha, alpha_var, DIM, APA_GAP_Z)
    scan = pd.DataFrame([counter.result(w) for w in ts.SCAN_W])
    s, s_err = counter.shells(ts.SCAN_W)
    scan["shell_excess"] = np.append(s, np.nan)
    scan["shell_excess_err"] = np.append(s_err, np.nan)
    scan.insert(0, "point", POINT_NAME)
    scan.insert(0, "stage", skey)
    scan.insert(0, "region", rkey)

    if w_star is None:
        w_star, how = ts.choose_halfwidth(counter, ts.SCAN_W, log)
        log("  half-width chosen here and reused everywhere: %.1f cm (%s)" % (w_star, how))
    else:
        how = "common"
        log("  half-width: %.1f cm (common)" % w_star)

    res = counter.result(w_star)
    res.update(detector=rlabel, region=rkey, region_label=rlabel,
               stage=skey, stage_label=slabel,
               point=POINT_NAME, key=POINT_KEY, dim=DIM, coord=APA_GAP_Z,
               alpha=alpha, alpha_err=float(np.sqrt(alpha_var)),
               norm_on=on_r, norm_off=off_r, norm_mc=mc_r,
               halfwidth_how=how, POT=pot, OFF_w=off_w, y_max=Y_MAX,
               onbeam=args.onbeam_name)

    log("  N_ON = %.0f  N_OFF = %.0f  ->  N_data = %.1f +/- %.1f"
        % (res["N_on"], res["N_off"], res["N_data"], res["N_data_err"]))
    log("  N_MC (normalized) = %.1f" % res["N_mc_scaled"])
    log("  excess = %.1f +/- %.1f  (%.1f sigma)"
        % (res["excess"], res["excess_err"],
           res["excess"]/res["excess_err"] if res["excess_err"] > 0 else 0))
    log("  MC muons crossing the gap = %.1f ; at-risk total = %.1f"
        % (alpha*res["n_cross_mc"], res["denom_scaled"]))
    log("  SPLIT FRACTION = %.3f%% +/- %.3f%%"
        % (100*res["f_split"], 100*res["f_split_err"]))
    neigh = ["w=%.1f: %.2f%%" % (w, 100*counter.result(w)["f_split"])
             for w in (w_star - 1.0, w_star + 1.0) if w >= ts.SCAN_W[0]]
    log("    (stability: %s)" % "; ".join(neigh))

    tag = "SBND_%s_%s_%s" % (skey, rkey, POINT_KEY)
    name = "%s (%s)" % (POINT_NAME, rlabel)
    ts._save(ts.zoom_plot(smc, son, soff, off_w, alpha, DIM, APA_GAP_Z, w_star,
                          name, DETECTOR, skey),
             plotdir, tag + "_muend_zoom", dosave)
    ts._save(ts.scan_plot(scan, w_star, name, DETECTOR, skey),
             plotdir, tag + "_scan", dosave)

    arows = []
    if not args.no_angular:
        log("")
        log("  angular scan (%d bins of theta_z):" % args.nquantile)
        norm_sel = {nm: norm_masks[nm][masks[nm]] for nm in ("MC", "ON", "OFF")}
        arows = ts.run_angular(smc, son, soff, off_w, norm_sel, DIM, APA_GAP_Z,
                               w_star, args.nquantile, log)
        for r in arows:
            r.update(detector=rlabel, region=rkey, region_label=rlabel,
                     stage=skey, stage_label=slabel, point=POINT_NAME,
                     key=POINT_KEY, dim=DIM, coord=APA_GAP_Z,
                     nquantile=args.nquantile, POT=pot, OFF_w=off_w)
        if arows:
            ts._save(ts.angular_plot(arows, res, name, DETECTOR, skey, DIM),
                     plotdir, tag + "_f_vs_theta", dosave)
            for r in arows:
                lo, hi, bi = r["theta_lo"], r["theta_hi"], r["theta_bin"]
                last = bi == args.nquantile - 1

                def _sub(fr, lo=lo, hi=hi, last=last):
                    t = np.asarray(ts.theta_normal(fr, DIM))
                    return fr[(t >= lo) & ((t <= hi) if last else (t < hi))]

                ts._save(ts.zoom_plot(_sub(smc), _sub(son), _sub(soff), off_w,
                                      r["alpha"], DIM, APA_GAP_Z, w_star, name,
                                      DETECTOR, skey,
                                      annotate=r"$\theta_z$ = %.1f$-$%.1f$^\circ$"
                                               % (lo, hi)),
                         plotdir, "%s_muend_zoom_th%d" % (tag, bi), dosave)

    return res, scan, arows, w_star


# ============================================================
# Comparison plots
# ============================================================
def tpc_summary_plot(rows, combo, stages=(ts.WIDTH_STAGE,)):
    """f +/- sigma for both TPCs, East, West (and the E+W combination), at the
    simple-cosmic-rejection stage only by default."""
    stages = [s for s in stages if any(r["stage"] == s for r in rows)]
    groups = [(k, REGION_LABEL[k]) for (k, _, _) in REGIONS] + [("east+west", "East+West")]
    fig, ax = plt.subplots(figsize=(8, 5))
    for gi, (g, glabel) in enumerate(groups):
        src = combo if g == "east+west" else rows
        xs, ys, es = [], [], []
        for si, s in enumerate(stages):
            m = [r for r in src if r["region"] == g and r["stage"] == s]
            if not m:
                continue
            xs.append(si + (gi - 0.5*(len(groups) - 1))*0.16)
            ys.append(100*m[0]["f_split"])
            es.append(100*m[0]["f_split_err"])
        if xs:
            ax.errorbar(xs, ys, es, linestyle="none", capsize=4, markersize=7,
                        marker="s" if g == "east+west" else "o",
                        mfc="none" if g == "east+west" else REGION_COLORS[g],
                        color=REGION_COLORS[g], label=glabel)
    ax.set_xticks(range(len(stages)))
    ax.set_xticklabels([STAGE_LABEL[s] for s in stages], fontsize=FONTSIZE)
    ax.set_xlim([-0.5, len(stages) - 0.5])
    ax.axhline(0, color="gray", linewidth=0.8)
    ax.set_ylabel("Muons to split [%]", fontsize=FONTSIZE, fontweight='bold')
    ax.set_title("$\\bf{SBND}$ APA gap (z = %.0f cm, y < %.0f cm)" % (APA_GAP_Z, Y_MAX),
                 fontsize=FONTSIZE+2)
    ax.legend(fontsize=11)
    lo_y, hi_y = ax.get_ylim()
    ax.set_ylim([lo_y, hi_y + 0.3*(hi_y - lo_y)])
    ts._style(ax)
    return fig


def tpc_ratio_plot(ratios):
    """f(full GUMP) / f(simple cosmic rejection) per region."""
    fig, ax = plt.subplots(figsize=(7, 5))
    keys = [k for (k, _, _) in REGIONS]
    for i, k in enumerate(keys):
        m = [r for r in ratios if r["region"] == k]
        if not m or not np.isfinite(m[0]["ratio"]):
            continue
        ax.errorbar([i], [m[0]["ratio"]], [m[0]["ratio_err"]], linestyle="none",
                    marker="o", markersize=7, capsize=4, color=REGION_COLORS[k])
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels([REGION_LABEL[k] for k in keys], fontsize=FONTSIZE)
    ax.set_xlim([-0.5, len(keys) - 0.5])
    ax.axhline(1, color="red", linestyle="--", linewidth=1.2)
    ax.axhline(0, color="gray", linewidth=0.8)
    ax.set_ylabel("$f$(Full GUMP) / $f$(Simple Cos. Rej.)",
                  fontsize=FONTSIZE, fontweight='bold')
    ax.set_title("$\\bf{SBND}$ APA gap: stage ratio", fontsize=FONTSIZE+2)
    ts._style(ax)
    return fig


def tpc_angular_compare_plot(per, stage):
    """f vs theta_z overlaid for the regions at one stage (adapted from
    ts.angular_compare_plot). Each region keeps its own quantile edges."""
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for i, (rkey, (rows, incl)) in enumerate(per.items()):
        c = REGION_COLORS[rkey]
        lo = np.array([r["theta_lo"] for r in rows])
        hi = np.array([r["theta_hi"] for r in rows])
        f = np.array([100*r["f_split"] for r in rows])
        fe = np.array([100*r["f_split_err"] for r in rows])
        x = 0.5*(lo + hi) + (i - 0.5*(len(per) - 1))*1.0
        chi2, ndof, pval = ts.flat_line_chi2([r["f_split"] for r in rows],
                                             [r["f_split_err"] for r in rows])
        lab = REGION_LABEL[rkey]
        if ndof > 0:
            lab += ": flat $\\chi^2$/ndof = %.1f/%d (p = %.3f)" % (chi2, ndof, pval)
        ax.errorbar(x, f, yerr=fe, xerr=[x - lo, hi - x], color=c, linestyle="none",
                    marker="o", markersize=5, capsize=0, zorder=3, label=lab)
        if incl is not None:
            fi, fie = 100*incl["f_split"], 100*incl["f_split_err"]
            ax.axhspan(fi - fie, fi + fie, color=c, alpha=0.12, zorder=0)
            ax.axhline(fi, color=c, linewidth=1.0, linestyle="--", zorder=1)
    ax.axhline(0, color="gray", linewidth=0.8, zorder=0)
    ts._style(ax)
    ax.set_xlim([0, 90])
    ax.set_xlabel("$\\theta_{z}$ to plane normal [deg]", fontsize=FONTSIZE, fontweight='bold')
    ax.set_ylabel("Split fraction [%]", fontsize=FONTSIZE, fontweight='bold')
    ax.set_title("$\\bf{SBND}$ APA gap -- %s" % STAGE_LABEL[stage], fontsize=FONTSIZE+2)
    ax.legend(fontsize=10, loc="upper right")
    lo_y, hi_y = ax.get_ylim()
    ax.set_ylim([lo_y, hi_y + 0.45*(hi_y - lo_y)])
    return fig


# ============================================================
# Driver
# ============================================================
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--df-dir", default=DEFAULT_DF_DIR, help="Directory of input .df files")
    p.add_argument("--plotdir", default=DEFAULT_PLOTDIR, help="Output directory")
    p.add_argument("--halfwidth", type=float, default=None,
                   help="Override the automatic window half-width [cm]")
    p.add_argument("--stage", action="append", default=None,
                   help="Selection stage: %s, or 'all' (default). Repeatable."
                        % ", ".join(repr(k) for (k, _, _) in STAGES))
    p.add_argument("--nquantile", type=int, default=ts.N_QUANTILE,
                   help="Number of equal-population angular bins [default %d]"
                        % ts.N_QUANTILE)
    p.add_argument("--no-angular", action="store_true", help="Skip the angular scan")
    p.add_argument("--no-dirt", action="store_true", help="Skip the dirt sample")
    p.add_argument("--no-save", action="store_true", help="Compute but write nothing")
    args = p.parse_args(argv)

    keys = [k for (k, _, _) in STAGES]
    sts = keys if not args.stage or "all" in args.stage else args.stage
    bad = [s for s in sts if s not in keys]
    if bad:
        p.error("unknown stage(s) %s (choose from %s, or 'all')" % (bad, keys))
    stages = [s for s in STAGES if s[0] in sts]
    if stages[0][0] != ts.WIDTH_STAGE and args.halfwidth is None:
        p.error("the %r stage sets the window half-width; include it or pass --halfwidth"
                % ts.WIDTH_STAGE)
    if args.nquantile < 2:
        p.error("--nquantile must be at least 2")

    plotdir, dosave = args.plotdir, not args.no_save
    if dosave:
        for sub in ("", "png", "pdf"):
            os.makedirs(os.path.join(plotdir, sub), exist_ok=True)

    lines = []

    def log(msg):
        print(msg, flush=True)
        lines.append(str(msg))

    log("=" * 70)
    log("TrackSplittingCorrection: SBND APA gap (z = %.1f cm, y < %.0f cm)" % (APA_GAP_Z, Y_MAX))
    log("=" * 70)

    files = mdc.detector_files(DETECTOR, args.df_dir, full_onbeam=True)
    files["MC_FILES"] = mc_shards(args.df_dir)
    args.onbeam_name = os.path.basename(files["ONBEAM"])
    log("on-beam stream: %s  [FULL STATS]" % args.onbeam_name)

    norm = mdc.compute_norm(DETECTOR, files, log)
    off_w, pot = float(norm["OFF_w"]), float(norm["POT"])
    log("POT = %.4g  OFF_w = %.4f" % (pot, off_w))

    mcdf, ondf, offdf = load_all(files, pot, not args.no_dirt, log)
    frames = (mcdf, ondf, offdf)
    named = (("MC", mcdf), ("ON", ondf), ("OFF", offdf))

    norm_masks = {nm: ((f["mu_end_" + DIM] - APA_GAP_Z).abs() >= ts.NORM_EXCLUDE_W).to_numpy()
                  for nm, f in named}
    ycut = {nm: continuous_fv(f) for nm, f in named}
    for nm, f in named:
        log("  %s: %d / %d FV rows below y = %.0f" % (nm, int(ycut[nm].sum()), len(f), Y_MAX))
    stage_base = {skey: {nm: (sfn(f) & gc.flash_cut(f)).to_numpy() & ycut[nm]
                         for nm, f in named}
                  for (skey, _, sfn) in stages}

    masks = {}
    results, scans, angular = [], [], []
    w_star = args.halfwidth
    for region in REGIONS:
        rkey, _, side = region
        rm = {nm: region_mask(f, side) for nm, f in named}
        for stage in stages:
            skey = stage[0]
            masks[(rkey, skey)] = {nm: stage_base[skey][nm] & rm[nm] for nm in rm}
            res, scan, arows, w_star = run_one(region, stage, frames, masks[(rkey, skey)],
                                               norm_masks, off_w, pot, w_star, args, log)
            results.append(res)
            scans.append(scan)
            angular += arows

    # ---- East+West combined as two independent samples (cross-check) ----
    combo = []
    for (skey, slabel, _) in stages:
        parts = [r for r in results if r["stage"] == skey and r["region"] in ("east", "west")]
        if len(parts) != 2:
            continue
        excess = sum(r["excess"] for r in parts)
        denom = sum(r["denom_scaled"] for r in parts)
        var = sum(r["excess_err"]**2 for r in parts)
        var_a = sum((r["alpha_err"]*r["N_data"]/r["alpha"])**2 for r in parts)
        combo.append(dict(region="east+west", stage=skey, stage_label=slabel,
                          excess=excess, excess_err=np.sqrt(var), denom_scaled=denom,
                          f_split=excess/denom, f_split_err=np.sqrt((var + var_a)/denom**2)))

    # ---- stage ratio per region ----
    ratios = []
    a, b = RATIO_STAGES
    if all(s in [x[0] for x in stages] for s in (a, b)):
        log("")
        log("-" * 70)
        log("stage overlap (%s vs %s):" % (STAGE_LABEL[a], STAGE_LABEL[b]))
        for (rkey, rlabel, _) in REGIONS:
            ra = [r for r in results if r["region"] == rkey and r["stage"] == a][0]
            rb = [r for r in results if r["region"] == rkey and r["stage"] == b][0]
            cov_data, cov_norm = ts.stage_overlap(ondf, offdf, masks[(rkey, a)],
                                                  masks[(rkey, b)], norm_masks, off_w,
                                                  DIM, APA_GAP_Z, ra["halfwidth"])
            x = dict(detector=rlabel, cov_data=cov_data,
                     cov_alpha=cov_norm/(ra["norm_mc"]*rb["norm_mc"]))
            ratio, rerr, cov = ts.stage_ratio([ra], [rb], [x])
            ratios.append(dict(region=rkey, region_label=rlabel, f_a=ra["f_split"],
                               f_b=rb["f_split"], ratio=ratio, ratio_err=rerr, cov=cov,
                               cov_data=cov_data, cov_alpha=x["cov_alpha"]))
            log("  %-10s cov(N_data) = %.2f  cov(alpha) = %.3e  ratio = %.2f +/- %.2f"
                % (rlabel, cov_data, x["cov_alpha"], ratio, rerr))

    # ---- tables ----
    W = 92
    out = ["", "=" * W,
           "SBND APA-GAP TRACK SPLITTING  (z = %.0f cm, y < %.0f cm, w = %.1f cm)"
           % (APA_GAP_Z, Y_MAX, w_star),
           "fraction of gap-crossing MC muons to split; data-stat errors only",
           "=" * W]
    for (skey, slabel, _) in stages:
        out += ["", "STAGE: %s" % slabel, "-" * W,
                "%-18s %8s %10s %10s %12s %20s" % ("Region", "alpha", "N_ON", "excess",
                                                   "at-risk", "split fraction"),
                "-" * W]
        for r in [r for r in results if r["stage"] == skey]:
            out.append("%-18s %8.4f %10.0f %10.1f %12.1f   (%5.2f +/- %4.2f) %%"
                       % (r["region_label"], r["alpha"], r["N_on"], r["excess"],
                          r["denom_scaled"], 100*r["f_split"], 100*r["f_split_err"]))
        for c in [c for c in combo if c["stage"] == skey]:
            out.append("%-18s %8s %10s %10.1f %12.1f   (%5.2f +/- %4.2f) %%"
                       % ("East+West (sum)", "", "", c["excess"], c["denom_scaled"],
                          100*c["f_split"], 100*c["f_split_err"]))
        e = [r for r in results if r["stage"] == skey and r["region"] == "east"]
        w_ = [r for r in results if r["stage"] == skey and r["region"] == "west"]
        if e and w_:
            d = e[0]["f_split"] - w_[0]["f_split"]
            de = np.hypot(e[0]["f_split_err"], w_[0]["f_split_err"])
            out.append("East - West: (%.2f +/- %.2f) %%  = %.1f sigma"
                       % (100*d, 100*de, d/de if de > 0 else 0))
        out.append("-" * W)

    if ratios:
        out += ["", "STAGE RATIO  f(%s) / f(%s)  (overlap propagated)"
                % (STAGE_LABEL[b], STAGE_LABEL[a]), "-" * W]
        for r in ratios:
            out.append("%-18s %8.2f%% %8.2f%%   %6.2f +/- %.2f"
                       % (r["region_label"], 100*r["f_a"], 100*r["f_b"],
                          r["ratio"], r["ratio_err"]))
        out.append("-" * W)

    if angular:
        out += ["", "=" * W, "SPLIT FRACTION vs theta_z  (equal-population bins of "
                "the at-risk MC; alpha refitted per bin)", "=" * W]
        for (skey, slabel, _) in stages:
            for (rkey, rlabel, _) in REGIONS:
                rows_ = sorted([r for r in angular if r["stage"] == skey
                                and r["region"] == rkey], key=lambda r: r["theta_bin"])
                if not rows_:
                    continue
                chi2, ndof, pval = ts.flat_line_chi2([r["f_split"] for r in rows_],
                                                     [r["f_split_err"] for r in rows_])
                out += ["-" * W, "%s  [%s]" % (rlabel, slabel),
                        "%18s %8s %10s %12s %20s" % ("theta_z [deg]", "N_ON", "alpha",
                                                     "D (MC)", "split fraction")]
                for r in rows_:
                    out.append("%8.1f - %7.1f %8.0f %10.4f %12.1f   (%6.2f +/- %5.2f) %%"
                               % (r["theta_lo"], r["theta_hi"], r["N_on"], r["alpha"],
                                  r["denom_mc"], 100*r["f_split"], 100*r["f_split_err"]))
                if ndof > 0:
                    out.append("%18s  chi2/ndof = %.2f/%d,  p = %.4f"
                               % ("flat line:", chi2, ndof, pval))
        out.append("-" * W)

    for l in out:
        log(l)

    if not dosave:
        return 0

    ts._save(tpc_summary_plot(results, combo), plotdir, "SBND_tracksplit_summary_tpc", True)
    if ratios:
        ts._save(tpc_ratio_plot(ratios), plotdir, "SBND_tracksplit_stage_ratio", True)
    for (skey, _, _) in stages:
        per = {}
        for (rkey, _, _) in REGIONS:
            rows_ = sorted([r for r in angular if r["stage"] == skey and r["region"] == rkey],
                           key=lambda r: r["theta_bin"])
            incl = [r for r in results if r["stage"] == skey and r["region"] == rkey]
            if rows_:
                per[rkey] = (rows_, incl[0] if incl else None)
        if per:
            ts._save(tpc_angular_compare_plot(per, skey), plotdir,
                     "SBND_tracksplit_%s_f_vs_theta_tpc" % skey, True)

    pd.DataFrame(results + combo).to_csv(os.path.join(plotdir, "SBND_tracksplit.csv"), index=False)
    with open(os.path.join(plotdir, "SBND_tracksplit.json"), "w") as f:
        json.dump(results + combo, f, indent=2, default=ts._jsonable)
    pd.concat(scans, ignore_index=True).to_csv(os.path.join(plotdir, "SBND_scan.csv"),
                                               index=False)
    if ratios:
        pd.DataFrame(ratios).to_csv(os.path.join(plotdir, "SBND_stage_ratio.csv"), index=False)
    if angular:
        pd.DataFrame(angular).to_csv(os.path.join(plotdir, "SBND_tracksplit_angular.csv"),
                                     index=False)
        with open(os.path.join(plotdir, "SBND_tracksplit_angular.json"), "w") as f:
            json.dump(angular, f, indent=2, default=ts._jsonable)
    with open(os.path.join(plotdir, "tracksplit_summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\nwrote results to %s" % plotdir, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
