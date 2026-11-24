"""N-1 MC/data comparisons for the GUMPLE selections.

For every cut variable of the final GUMP 1mu1p and MAPLE 1muMp selections
(gumple_cuts.all_gump_cuts / all_maplemp_cuts), plot the MC/data comparison of
that variable with every OTHER cut applied. The machinery -- loading,
normalization, the full systematic covariance, the chi2 and the figure style --
is mcdata_comparison_gumple.py's: the N-1 stages are registered as extra stage
cuts before its prepare_detector() runs, so they get every systematic (including
the trigger universes and the ICARUS track-splitting re-evaluation) for free.

Cuts that CANNOT be N-1'd, because the production preselection (maple/makedf.py
selection="presel") already removed the failing slices, are listed in
NOT_N1_CUTS and written to the summary. In particular the muon candidate is the
longest track passing gc.get_base_muon_mask(level="trk"), so its Pandora
track-score / length / containment requirements are baked into the candidate
choice. The analysis track-score cut is the slice-level coltopo score
(gc.coltopo_cut, sbn-rewgted-26+), which IS N-1'd ("trkscore").

Mp follows gc.all_maplemp_cuts and so INCLUDES the far-shower veto, unlike the
"PID Mp" stage of mcdata_comparison_gumple.py.

Deliverables per detector, in --plotdir: pdf/ + png/ figures
<det>_n-1-<sel>_<var>[_bkdwnpdg|_bkdwntrkshw]
(_bkdwntrkshw: track-score plots only, track-like vs shower-like muon cand.), <det>_n-1_chi2.tsv and
<det>_n-1_summary.txt.

Options:
  --trkscore {coltopo,pandora}  muon-candidate track-score cut used in EVERY
      stage (and the one N-1'd as "trkscore"): the analysis coltopo cut
      (default) or the legacy Pandora mu_trackScore >= musel_track_score_min.
      The stored candidate is the legacy trackScore > 0 track, so both N-1.
  --full-onbeam  use the full-statistics on-beam stream (mcd.FULL_ONBEAM)
      instead of the prescaled/unblind default.
  --g4-scale S  Geant4 reinteraction uncertainty inflation, applied at LOAD time
      (loaddf.G4_UNC_SCALE: each universe weight w -> 1 + S*(w - 1)); default:
      the loaddf default (2). NB before 2026-11-13 this flag instead scaled the
      covariance by S^2 on top of the load-time scale.
  --sbnd-no-vtx-zgap  put the SBND slice-vertex APA gap (gc.SBND_VTX_ZGAP,
      250 < z < 254) back into the FV (cross-check).
  --exclude-tpc T / --only-tpc T  (ICARUS only) drop / keep only slices with
      their reco vertex in TPC T (EE, EW, WE, WW), via mcd.cryo_mask.
  --pass-only  compare only the part of each N-1 distribution that PASSES the
      removed cut (i.e. the full selection, histogrammed over the kept side of
      the cut). The kept range [lo, cut] or [cut, hi] of the nominal binning is
      re-binned with edges anchored exactly on the cut value and about the
      nominal bin width (see pass_bins); the shape-only normalization, the
      covariance and the chi2 then cover the pass region only. Figures get a
      "_passonly" suffix; the chi2 .tsv gains an n_pass_bins column.

    python analysis_village/gump/nminus1_gumple.py -d all --plotdir <dir>
"""

import argparse
import os
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from multiprocess import Pool
from tqdm.auto import tqdm

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import mcdata_comparison_gumple as mcd
gc = mcd.gc

# loaddf keys its load cache on the preselection's __module__.__qualname__;
# mcdata_comparison_gumple runs (and so caches) as "__main__.FV". FV's body is
# unchanged by importing it, so share that cache rather than reloading.
mcd.FV.__module__ = "__main__"

DEFAULT_PLOTDIR = "/Users/gputnam/Work/osc/cafpyana/plots-gumple-2026-11-03-rewgted26/evtsel_n-1/"
DEFAULT_DF_DIR = "/Users/gputnam/Work/osc/sbn-rewgted-26/"

# Cuts of the selection chain that the stored frames cannot relax.
# The muon track-score requirement (coltopo, slice level) is NOT among them: it
# is not part of the production's candidate choice.
NOT_N1_CUTS = [
    "sanity (slice vertex defined)",
    "slice-vertex FV (incl. ICARUS WW cable region)",
    "all-candidate containment (cut_contained)",
    "cathode crossing (cut_cathode)",
    "has_muon / cut_np (trk_cut)",
    "muon cand. length 40-400 cm, dist_start <= 10 cm, primary, contained, "
    "same-side (muon-candidate definition; only >400 cm and non-primary have "
    "any failing slices stored, ~1%)",
]


# ============================================================
# Cut pieces
# ============================================================
def c_angle(d):
    return d["mu_p_opening_angle_deg"] < gc._det_cut_th(d.detector, "max_opening_angle")


def c_crt(d):
    return d["crthit"] == 0


def c_mu_chi2u(d):
    return d["mu_chi2_of_mu_cand"] < gc._det_cut_th(d.detector, "musel_muscore_th")


def c_mu_chi2p(d):
    return d["prot_chi2_of_mu_cand"] > gc._det_cut_th(d.detector, "musel_pscore_th")


def c_p_chi2p(d):
    return gc.get_proton_mask(d)


def c_1p(d):
    return d.n_pfp == 2


def c_mp(d):
    return d.n_pfp > 2


def c_farshw(d):
    return np.isnan(d.max_far_shw_len)


# "coltopo" (analysis) or "pandora" (legacy mu_trackScore >= 0.5); set by
# parse_args before register_stages.
TRKSCORE = "coltopo"
# --only-keys: restrict the plotted N-1 cuts (None: all)
ONLY_KEYS = None
# --pass-only: histogram only the kept side of each removed cut
PASS_ONLY = False


def pass_bins(bins, line, side):
    """Bins covering only the kept side of a single cut at `line`.

    The nominal bins do not always have the cut on an edge (e.g. coltopo 0.89,
    the chi2 thresholds), so the kept range is re-binned: [lo, line] for "<",
    [line, hi] for ">", in n equal bins with n = round(range / nominal width)
    (at least 1). The bins then end exactly on the cut."""
    bins = np.asarray(bins, dtype=float)
    w = np.median(np.diff(bins))
    lo, hi = bins[0], bins[-1]
    if side == "<":
        a, b = lo, min(line, hi)
    elif side == ">":
        a, b = max(line, lo), hi
    else:
        raise ValueError("pass_bins: unsupported cut side %r" % side)
    n = max(1, int(round((b - a)/w)))
    return np.linspace(a, b, n + 1)


def c_trkscore(d):
    if TRKSCORE == "pandora":
        # NaN fails
        return d["mu_trackScore"] >= gc._det_cut_th(d.detector, "musel_track_score_min")
    return gc.coltopo_cut(d)


def c_base(d):
    # always applied: FV (== the load preselection), trk_cut and the
    # chi2-free muon-candidate mask (gc.get_base_muon_mask), minus its
    # track-score requirement, which is the removable c_trkscore
    return (mcd.FV(d) & gc.trk_cut(d)
            & d["mu_end_x"].notna() & d["mu_len"].notna()
            & (d["mu_dist_start"] <= 10.0)
            & d["mu_len"].between(gc._det_cut_th(d.detector, "musel_len_th_min"),
                                  gc._det_cut_th(d.detector, "musel_len_th_max"))
            & d["mu_prim_pfp"] & d["mu_contained10"]
            & (d["mu_end_x"] * d["slc_vtx_x"] > 0))


# The removable cuts of each selection, in gc.all_gump_cuts / all_maplemp_cuts
# order. "flash" is handled by mcd's stage builder (NO_FLASH_STAGES).
SELECTIONS = {
    "1p": dict(angle=c_angle, crt=c_crt, flash=None, trkscore=c_trkscore,
               mu_chi2u=c_mu_chi2u, mu_chi2p=c_mu_chi2p, p_chi2p=c_p_chi2p, npfp=c_1p),
    "mp": dict(flash=None, trkscore=c_trkscore, mu_chi2u=c_mu_chi2u, mu_chi2p=c_mu_chi2p,
               p_chi2p=c_p_chi2p, npfp=c_mp, farshw=c_farshw),
}
SEL_TITLES = {"1p": "GUMP 1$\\mu$1p", "mp": "MAPLE 1$\\mu$Mp"}


# Zoomed on the bulk of the selected events, with the nominal threshold on a
# bin edge.
FLASH_BINS = {
    "SBND":        np.linspace(0, 20000, 21),
    "ICARUS Run2": np.linspace(0, 40000, 17),
    "ICARUS Run4": np.linspace(0, 20000, 21),
}


def flash_bins(detector):
    return FLASH_BINS[detector]


# cut key -> (var, bins(detector), label, cut lines(detector, sel), kept side)
# kept side: "<" keeps the region left of the (single) line, ">" right of it,
# "in" between two lines.
CUTVARS = {
    "angle":    ("mu_p_opening_angle_deg", lambda det: np.linspace(0, 180, 19),
                 "$\\theta_{\\mu p}$ [$^\\circ$]",
                 lambda det, sel: [gc._det_cut_th(det, "max_opening_angle")], "<"),
    "crt":      ("crthit_f", lambda det: np.array([-0.5, 0.5, 1.5]),
                 "CRT Hit in Beam Window",
                 lambda det, sel: [0.5], "<"),
    "flash":    ("flash_maxpe", flash_bins, "Max. Flash PE",
                 lambda det, sel: [gc.FLASH_PE_THRESHOLD[det]], ">"),
    # coltopo is shower-like (muons near 0): the cut keeps the region below the line.
    # NaN scores (~0.3%) land in the fillna=-1 underflow, outside the bins.
    "trkscore": None,  # per TRKSCORE, see TRKSCORE_VARS
    "mu_chi2u": ("mu_chi2_of_mu_cand", lambda det: np.linspace(0, 60, 13),
                 "Muon Cand. $\\chi^2_\\mu$",
                 lambda det, sel: [gc._det_cut_th(det, "musel_muscore_th")], "<"),
    # 0-255 in 15-wide bins (the 255-270 bin dropped 2026-11-24)
    "mu_chi2p": ("prot_chi2_of_mu_cand", lambda det: np.linspace(0, 255, 18),
                 "Muon Cand. $\\chi^2_p$",
                 lambda det, sel: [gc._det_cut_th(det, "musel_pscore_th")], ">"),
    # 10-wide bins: SBND 0-200, ICARUS 0-150
    "p_chi2p":  ("prot_chi2_of_prot_cand", lambda det: (np.linspace(0, 200, 21) if det == "SBND"
                                                        else np.linspace(0, 150, 16)),
                 "Proton Cand. $\\chi^2_p$",
                 lambda det, sel: [gc._det_cut_th(det, "psel_pscore_th")], "<"),
    # n_pfp >= 2 always (cut_np, production preselection), so 1p == "n_pfp <= 2"
    "npfp":     ("n_pfp", lambda det: np.linspace(1.5, 4.5, 4),
                 "Number of Candidate PFPs",
                 lambda det, sel: [2.5],
                 None),  # side depends on sel, see cut_side()
    "farshw":   ("max_far_shw_len", lambda det: np.array([-10.] + list(np.linspace(0, 50, 11))),
                 "Far Shower Length [cm] (none: underflow)",
                 lambda det, sel: [0.], "<"),
}
TRKSCORE_VARS = {
    "coltopo": ("coltopo_trkscore", lambda det: np.linspace(0, 1, 21),
                "Muon Cand. Track/Shower Score (coltopo)",
                lambda det, sel: [gc._det_cut_th(det, "musel_coltopo_max")], "<"),
    # Pandora is track-like (muons near 1): the cut keeps the region above the line.
    "pandora": ("mu_trackScore", lambda det: np.linspace(0, 1, 21),
                "Muon Cand. Track Score (Pandora)",
                lambda det, sel: [gc._det_cut_th(det, "musel_track_score_min")], ">"),
}


def cutvar(key):
    return TRKSCORE_VARS[TRKSCORE] if key == "trkscore" else CUTVARS[key]


# Track-score breakdown of the muon candidate's true particle: track-like
# (mu / pi+- / p) vs shower-like (e+- / gamma). |pdg| is used, so both charges
# land in their class. "Other": anything else or no truth match (incl. dirt).
TRKSHW_PDGS = [(13, 211, 2212), (11, 22)]
trkshw_labels = ["Track-like ($\\mu$/$\\pi^\\pm$/$p$)", "Shower-like ($e^\\pm$/$\\gamma$)", "Other"]
trkshw_colors = ["#1e3f54", "#d54c28", "#c89648"]


def breakdown_trkshw(var, d):
    apdg = d["mu_true_pdg"].abs()
    masks = [apdg.isin(p) for p in TRKSHW_PDGS]
    return [var[m] for m in masks] + [var[~(masks[0] | masks[1])]]


PDG_CUTS = {"mu_chi2u", "mu_chi2p", "p_chi2p"}


def cut_side(key, sel):
    if key == "npfp":
        return "<" if sel == "1p" else ">"
    return cutvar(key)[4]


def n1_keys(detector, sel):
    keys = list(SELECTIONS[sel])
    # the ICARUS chi2u threshold is inf (cut disabled): nothing to N-1
    if "ICARUS" in detector:
        keys.remove("mu_chi2u")
    # SBND's production crthit is False for every slice (no CRT veto
    # information): the CRT N-1 stage would equal the full selection
    elif "crt" in keys:
        keys.remove("crt")
    if ONLY_KEYS is not None:
        keys = [k for k in keys if k in ONLY_KEYS]
    return keys


def stage_name(sel, key):
    return "N-1 %s %s" % (sel, key)


def make_stage(sel, removed):
    pieces = [f for k, f in SELECTIONS[sel].items() if k != removed and f is not None]

    def stage(d):
        m = c_base(d)
        for f in pieces:
            m = m & f(d)
        return m
    stage.__name__ = "n1_%s_%s" % (sel, removed)
    return stage


def register_stages(detector):
    """Append the N-1 stages to mcd's stage list. Returns the stage specs."""
    specs = []
    for sel in SELECTIONS:
        mcd.STAGE_CUTS.append(make_stage(sel, None))
        mcd.STAGE_CUTNAMES.append(stage_name(sel, "full"))
        for key in n1_keys(detector, sel):
            name = stage_name(sel, key)
            mcd.STAGE_CUTS.append(make_stage(sel, key))
            mcd.STAGE_CUTNAMES.append(name)
            if key == "flash":
                mcd.NO_FLASH_STAGES.add(name)
            var, binf, label, linef, _ = cutvar(key)
            bins, lines, side = binf(detector), linef(detector, sel), cut_side(key, sel)
            if PASS_ONLY:
                # the N-1 stage minus its failing side == the full selection
                bins = pass_bins(bins, lines[0], side)
            specs.append(dict(sel=sel, key=key, stage=name, var=var, bins=bins,
                              label=label, lines=lines, side=side))
    mcd.extravarnames.append("crthit_f")
    mcd.extravars.append(lambda d: d.crthit.astype(float))
    mcd.EXTRA_UNIV_COLS.extend(["crthit_f", "flash_maxpe", "n_pfp", "max_far_shw_len"])
    return specs


# ============================================================
# Pool workers (module level so fork'd children resolve them)
# ============================================================
def _worker(dat):
    (i, spec, kind) = dat
    pdg = kind.startswith("pdg")
    if kind == "trkshw":
        bkdwn, labels, colors = breakdown_trkshw, trkshw_labels, trkshw_colors
    elif pdg:
        bkdwn = mcd.breakdown_pdg_p if "prot_cand" in spec["var"] else mcd.breakdown_pdg_mu
        labels, colors = mcd.pdg_labels, mcd.pdg_colors
    else:
        bkdwn, labels, colors = mcd.breakdown_mode, mcd.mode_labels, mcd.mode_colors
    title = "%s, %s" % (SEL_TITLES[spec["sel"]], "N-1, pass cut" if PASS_ONLY else "N-1")
    # PID (chi2) comparisons: G4 / calorimetric covariance lines on the ratio panel
    comps = mcd.COMPONENT_SYSTS if (pdg and spec["key"] in PDG_CUTS) else None
    return i, kind, mcd.make_plot_data(spec["var"], spec["bins"], spec["stage"], "glob_scale",
                                       bkdwn, True, labels, colors, spec["label"], title,
                                       fillna=-1, det=mcd.DETECTOR, components=comps)


def draw_cut(ax0, ax1, spec):
    for x in spec["lines"]:
        for ax in (ax0, ax1):
            ax.axvline(x, color="k", linestyle="--", linewidth=1.5, zorder=10)
    # arrow(s) on the main panel pointing into the kept region
    lo, hi = spec["bins"][0], spec["bins"][-1]
    width = 0.08*(hi - lo)
    y = 0.45  # below the legend and chi2 text
    trans = ax0.get_xaxis_transform()
    lines = spec["lines"]
    arrows = []
    if spec["side"] == "<":
        arrows = [(lines[0], -width)]
    elif spec["side"] == ">":
        arrows = [(lines[0], +width)]
    elif spec["side"] == "in":
        arrows = [(lines[0], +min(width, 0.4*(lines[1] - lines[0]))),
                  (lines[1], -min(width, 0.4*(lines[1] - lines[0])))]
    for x, dx in arrows:
        ax0.annotate("", xy=(x + dx, y), xytext=(x, y), xycoords=trans, textcoords=trans,
                     arrowprops=dict(arrowstyle="->", color="k", lw=1.5), zorder=10)


def save_plot(pd_, spec, plotdir, suffix):
    fig, ax0, ax1 = mcd.ratio_plot(plt, pd_)
    if PASS_ONLY:
        suffix = suffix + "_passonly"  # the cut is a plot edge: no line/arrow
    else:
        draw_cut(ax0, ax1, spec)
    stem = "%s_n-1-%s_%s%s" % (pd_["det"].replace(" ", "-"), spec["sel"], spec["var"], suffix)
    fig.savefig(os.path.join(plotdir, "pdf", stem + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(plotdir, "png", stem + ".png"), bbox_inches="tight")
    plt.close(fig)


def run_detector(detector, args):
    specs = register_stages(detector)
    ctx = mcd.prepare_detector(detector, args)
    plotdir, detname, log = ctx["plotdir"], ctx["detname"], ctx["log"]
    log("\nG4 uncertainty: load-time scale loaddf.G4_UNC_SCALE = %g" % mcd.loaddf.G4_UNC_SCALE)
    log("SBND slice-vertex z gap excluded from FV: %s" % (gc.SBND_VTX_ZGAP,))
    df, ONdf, OFFdf, OFF_w = ctx["df"], ctx["ONdf"], ctx["OFFdf"], ctx["OFF_w"]

    log("\nN-1 track-score cut: %s   on-beam: %s" % (TRKSCORE, os.path.basename(mcd.detector_files(detector, args.df_dir, full_onbeam=args.full_onbeam)["ONBEAM"])))

    # ---- yields ----
    log("\nN-1 stage yields (MC: POT-scaled glob_scale sum; data: on-beam count, "
        "off-beam subtracted)")
    log("%-22s %12s %10s %10s %12s" % ("stage", "MC", "ON", "OFF", "ON-OFF*w"))
    for sel in SELECTIONS:
        for name in [stage_name(sel, "full")] + [s["stage"] for s in specs if s["sel"] == sel]:
            mc = df.loc[df[name], "glob_scale"].sum()
            non, noff = int(ONdf[name].sum()), int(OFFdf[name].sum())
            log("%-22s %12.1f %10d %10d %12.1f" % (name, mc, non, noff, non - noff*OFF_w))
    if TRKSCORE == "coltopo":  # the evtsel stages are built with coltopo
        log("cross-check, evtsel 'PID 1p' stage: MC %.1f  ON %d  OFF %d"
            % (df.loc[df["PID 1p"], "glob_scale"].sum(), ONdf["PID 1p"].sum(), OFFdf["PID 1p"].sum()))
        _mp = lambda f: f["PID Mp"] & c_farshw(f)
        log("cross-check, evtsel 'PID Mp' & far-shower veto: MC %.1f  ON %d  OFF %d"
            % (df.loc[_mp(df), "glob_scale"].sum(), _mp(ONdf).sum(), _mp(OFFdf).sum()))
    log("\nNot N-1'd (applied in the production preselection / muon-candidate definition):")
    for c in NOT_N1_CUTS:
        log("  - " + c)
    if "ICARUS" in detector:
        log("  - muon cand. chi2_mu (threshold is inf for ICARUS: cut disabled)")
    else:
        log("  - CRT hit veto (crthit is False for every SBND slice in the production)")

    # ---- plots ----
    inputs = []
    for i, s in enumerate(specs):
        inputs.append((i, s, "mode"))
        if s["key"] in PDG_CUTS:
            inputs.append((i, s, "pdg"))
        if s["key"] == "trkscore":
            inputs.append((i, s, "trkshw"))

    results = {}
    with Pool(min(args.nproc, len(inputs))) as p:
        for i, kind, pd_ in tqdm(p.imap_unordered(_worker, inputs), total=len(inputs)):
            results[(i, kind)] = pd_

    suffixes = {"mode": "", "pdg": "_bkdwnpdg", "trkshw": "_bkdwntrkshw"}
    rows = []
    for (i, kind) in sorted(results, key=lambda t: (t[0], t[1])):
        s, pd_ = specs[i], results[(i, kind)]
        if not args.no_save:
            save_plot(pd_, s, plotdir, suffixes[kind])
        chi2, ndof, pval = mcd.chi2_ndof_pvalue(pd_)
        rows.append(dict(detector=detector, selection=s["sel"], removed_cut=s["key"],
                         var=s["var"], breakdown=kind, norm="shape",
                         chi2=chi2, ndof=ndof, pvalue=pval,
                         pass_only=int(PASS_ONLY), n_pass_bins=len(s["bins"]) - 1))
    chi2df = pd.DataFrame(rows)
    log("\n" + chi2df[chi2df.breakdown == "mode"].to_string(index=False, float_format="%.3g"))
    log("\nSaved %d figures (pdf+png) to %s" % (0 if args.no_save else len(results), plotdir))

    if args.save_plotdata and not args.no_save:
        import pickle
        out = {(specs[i]["sel"], specs[i]["key"], kind): pd_ for (i, kind), pd_ in results.items()}
        fn = os.path.join(plotdir, "%s_n-1_plotdata.pkl" % detname)
        with open(fn, "wb") as f:
            pickle.dump(out, f)
        log("Saved plotdata (NMC_total, NData, NON, cov_w_stat, ...) to %s" % fn)

    if not args.no_save:
        chi2df.to_csv(os.path.join(plotdir, "%s_n-1_chi2.tsv" % detname), sep="\t",
                      index=False, float_format="%.6g")
        with open(os.path.join(plotdir, "%s_n-1_summary.txt" % detname), "w") as f:
            f.write("\n".join(ctx["summary_lines"]) + "\n")
    return 0


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-d", "--detector", action="append", default=None,
                   help="'SBND', 'ICARUS Run2', 'ICARUS Run4', or 'all'. Repeatable.")
    p.add_argument("--df-dir", default=DEFAULT_DF_DIR)
    p.add_argument("--plotdir", default=DEFAULT_PLOTDIR)
    p.add_argument("--nproc", type=int, default=8)
    p.add_argument("--no-dirt", action="store_true")
    p.add_argument("--no-save", action="store_true")
    p.add_argument("--trkscore", choices=sorted(TRKSCORE_VARS), default="coltopo",
                   help="muon-candidate track-score cut (default: coltopo)")
    p.add_argument("--full-onbeam", action="store_true",
                   help="use the full-statistics on-beam stream (mcd.FULL_ONBEAM)")
    p.add_argument("--g4-scale", type=float, default=None,
                   help="load-time G4 uncertainty scale (loaddf.G4_UNC_SCALE; default: loaddf's)")
    p.add_argument("--sbnd-no-vtx-zgap", action="store_true",
                   help="keep the SBND slice-vertex APA gap (gc.SBND_VTX_ZGAP) in the FV")
    p.add_argument("--save-plotdata", action="store_true",
                   help="pickle every comparison's histograms + covariance to <det>_n-1_plotdata.pkl")
    tpcgrp = p.add_mutually_exclusive_group()
    tpcgrp.add_argument("--only-tpc", choices=mcd.ICARUS_TPCS, default=None,
                   help="ICARUS only: keep only slices whose reco vertex is in this TPC")
    tpcgrp.add_argument("--exclude-tpc", choices=mcd.ICARUS_TPCS, default=None,
                   help="ICARUS only: drop slices whose reco vertex is in this TPC "
                        "(every MC/detvar/dirt/on/off frame, via mcd.cryo_mask)")
    p.add_argument("--run2-fv-include-ew", action="store_true",
                   help="ICARUS Run2: put TPC EW back into the FV (pre-2026-11-10 Run 2 FV, "
                        "gumple_cuts.RUN2_FV_INCLUDE_EW); combine with --only-tpc EW for EW only")
    p.add_argument("--pass-only", action="store_true",
                   help="histogram only the side of each removed cut that passes it")
    p.add_argument("--onbeam-file", default=None, metavar="NAME",
                   help="with --full-onbeam: use this on-beam .df (basename in --df-dir) instead of "
                        "mcd.FULL_ONBEAM[detector], e.g. SBND_Run1BNBData_Full2p5e20.df")
    p.add_argument("--offbeam-file", default=None, metavar="NAME",
                   help="use this off-beam .df (basename in --df-dir) instead of the detector's default "
                        "off-beam file(s), e.g. SBND_Run1BNBOffData_Full2p5e20.df")
    p.add_argument("--only-keys", default=None,
                   help="comma-separated N-1 cut keys to plot (default: all)")
    args = p.parse_args(argv)
    global TRKSCORE, ONLY_KEYS, PASS_ONLY
    TRKSCORE = args.trkscore
    PASS_ONLY = args.pass_only
    ONLY_KEYS = None if args.only_keys is None else args.only_keys.split(",")
    if args.onbeam_file is not None:
        if not args.full_onbeam:
            p.error("--onbeam-file needs --full-onbeam")
    if args.g4_scale is not None:
        # env so spawned loader workers see it; enters loaddf's MC cache key
        os.environ["GUMPLE_G4_UNC_SCALE"] = repr(args.g4_scale)
        mcd.loaddf.G4_UNC_SCALE = args.g4_scale
    if args.sbnd_no_vtx_zgap:
        # env for loader workers; the gap enters loaddf's SBND cache key
        os.environ["GUMPLE_SBND_VTX_ZGAP"] = "0"
        gc.SBND_VTX_ZGAP = None
    if args.run2_fv_include_ew:
        # Since f140cbc (2026-11-13) the default Run 2 FV already includes TPC EW
        # and gumple_cuts no longer reads RUN2_FV_INCLUDE_EW, so this flag would
        # silently do nothing. Use --exclude-tpc EW for the no-EW variant.
        p.error("--run2-fv-include-ew is obsolete: the Run 2 FV includes TPC EW by default "
                "(use --exclude-tpc EW to drop it)")
    # fields mcd.prepare_detector reads
    args.absnorm = False  # shape-normalized only
    args.cryostat = ("no-" + args.exclude_tpc if args.exclude_tpc else
                     "only-" + args.only_tpc if args.only_tpc else None)
    args.blinded = False
    args.offbeam_scale = 1.0

    dets = []
    for d in (args.detector or ["all"]):
        if d == "all":
            dets += mcd.DETECTORS
        elif d in mcd.DETECTORS:
            dets.append(d)
        else:
            p.error("unknown detector %r" % d)
    args.detector = dets
    if args.onbeam_file is not None:
        for d in dets:
            mcd.FULL_ONBEAM[d] = args.onbeam_file
    if args.offbeam_file is not None:
        for d in dets:
            mcd.OFFBEAM_OVERRIDE[d] = [args.offbeam_file]
    return args


def main(argv=None):
    args = parse_args(argv)

    if len(args.detector) == 1:
        return run_detector(args.detector[0], args)

    passthrough = ["--df-dir", args.df_dir, "--plotdir", args.plotdir, "--nproc", str(args.nproc),
                   "--trkscore", args.trkscore]
    if args.g4_scale is not None:
        passthrough += ["--g4-scale", repr(args.g4_scale)]
    if args.exclude_tpc:
        passthrough += ["--exclude-tpc", args.exclude_tpc]
    if args.only_tpc:
        passthrough += ["--only-tpc", args.only_tpc]
    if args.only_keys:
        passthrough += ["--only-keys", args.only_keys]
    if args.onbeam_file:
        passthrough += ["--onbeam-file", args.onbeam_file]
    if args.offbeam_file:
        passthrough += ["--offbeam-file", args.offbeam_file]
    for flag in ("no_dirt", "no_save", "full_onbeam", "save_plotdata", "sbnd_no_vtx_zgap", "pass_only"):
        if getattr(args, flag):
            passthrough.append("--" + flag.replace("_", "-"))
    rc = 0
    for det in args.detector:
        cmd = [sys.executable, os.path.abspath(__file__), "--detector", det] + passthrough
        print("\n>>> %s" % " ".join(repr(c) if " " in c else c for c in cmd), flush=True)
        r = subprocess.run(cmd)
        if r.returncode != 0:
            print("!!! %s FAILED (exit %d)" % (det, r.returncode), flush=True)
            rc = r.returncode
    return rc


if __name__ == "__main__":
    sys.exit(main())
