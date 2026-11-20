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
longest track passing gc.get_base_muon_mask, so its track-score / length /
containment requirements are baked into the candidate choice.

Mp follows gc.all_maplemp_cuts and so INCLUDES the far-shower veto, unlike the
"PID Mp" stage of mcdata_comparison_gumple.py.

Deliverables per detector, in --plotdir: pdf/ + png/ figures
<det>_n-1-<sel>_<var>[_bkdwnpdg], <det>_n-1_chi2.tsv and
<det>_n-1_summary.txt.

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

DEFAULT_PLOTDIR = "/Users/gputnam/Work/osc/cafpyana/plots-gumple-2026-10-20-rewgted25/evtsel_n-1/"
DEFAULT_DF_DIR = "/Users/gputnam/Work/osc/sbn-rewgted-25/"

# Cuts of the selection chain that the stored frames cannot relax.
# The muon track-score requirement is NOT among them: the stored candidate is
# chosen with a looser score (~9% of stored slices have mu_trackScore < 0.5).
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


def c_trkscore(d):
    return d["mu_trackScore"] >= gc._det_cut_th(d.detector, "musel_track_score_min")


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
    "trkscore": ("mu_trackScore", lambda det: np.linspace(0.4, 0.9, 11),
                 "Muon Cand. Track Score",
                 lambda det, sel: [gc._det_cut_th(det, "musel_track_score_min")], ">"),
    "mu_chi2u": ("mu_chi2_of_mu_cand", lambda det: np.linspace(0, 60, 13),
                 "Muon Cand. $\\chi^2_\\mu$",
                 lambda det, sel: [gc._det_cut_th(det, "musel_muscore_th")], "<"),
    "mu_chi2p": ("prot_chi2_of_mu_cand", lambda det: np.linspace(0, 270, 19),
                 "Muon Cand. $\\chi^2_p$",
                 lambda det, sel: [gc._det_cut_th(det, "musel_pscore_th")], ">"),
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
PDG_CUTS = {"mu_chi2u", "mu_chi2p", "p_chi2p"}


def cut_side(key, sel):
    if key == "npfp":
        return "<" if sel == "1p" else ">"
    return CUTVARS[key][4]


def n1_keys(detector, sel):
    keys = list(SELECTIONS[sel])
    # the ICARUS chi2u threshold is inf (cut disabled): nothing to N-1
    if "ICARUS" in detector:
        keys.remove("mu_chi2u")
    # SBND's production crthit is False for every slice (no CRT veto
    # information): the CRT N-1 stage would equal the full selection
    elif "crt" in keys:
        keys.remove("crt")
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
            var, binf, label, linef, _ = CUTVARS[key]
            specs.append(dict(sel=sel, key=key, stage=name, var=var, bins=binf(detector),
                              label=label, lines=linef(detector, sel),
                              side=cut_side(key, sel)))
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
    if pdg:
        bkdwn = mcd.breakdown_pdg_p if "prot_cand" in spec["var"] else mcd.breakdown_pdg_mu
        labels, colors = mcd.pdg_labels, mcd.pdg_colors
    else:
        bkdwn, labels, colors = mcd.breakdown_mode, mcd.mode_labels, mcd.mode_colors
    title = "%s, N-1" % SEL_TITLES[spec["sel"]]
    return i, kind, mcd.make_plot_data(spec["var"], spec["bins"], spec["stage"], "glob_scale",
                                       bkdwn, True, labels, colors, spec["label"], title,
                                       fillna=-1, det=mcd.DETECTOR)


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
    draw_cut(ax0, ax1, spec)
    stem = "%s_n-1-%s_%s%s" % (pd_["det"].replace(" ", "-"), spec["sel"], spec["var"], suffix)
    fig.savefig(os.path.join(plotdir, "pdf", stem + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(plotdir, "png", stem + ".png"), bbox_inches="tight")
    plt.close(fig)


def run_detector(detector, args):
    specs = register_stages(detector)
    ctx = mcd.prepare_detector(detector, args)
    plotdir, detname, log = ctx["plotdir"], ctx["detname"], ctx["log"]
    df, ONdf, OFFdf, OFF_w = ctx["df"], ctx["ONdf"], ctx["OFFdf"], ctx["OFF_w"]

    # ---- yields ----
    log("\nN-1 stage yields (MC: POT-scaled glob_scale sum; data: on-beam count, "
        "off-beam subtracted)")
    log("%-22s %12s %10s %10s %12s" % ("stage", "MC", "ON", "OFF", "ON-OFF*w"))
    for sel in SELECTIONS:
        for name in [stage_name(sel, "full")] + [s["stage"] for s in specs if s["sel"] == sel]:
            mc = df.loc[df[name], "glob_scale"].sum()
            non, noff = int(ONdf[name].sum()), int(OFFdf[name].sum())
            log("%-22s %12.1f %10d %10d %12.1f" % (name, mc, non, noff, non - noff*OFF_w))
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

    results = {}
    with Pool(min(args.nproc, len(inputs))) as p:
        for i, kind, pd_ in tqdm(p.imap_unordered(_worker, inputs), total=len(inputs)):
            results[(i, kind)] = pd_

    suffixes = {"mode": "", "pdg": "_bkdwnpdg"}
    rows = []
    for (i, kind) in sorted(results, key=lambda t: (t[0], t[1])):
        s, pd_ = specs[i], results[(i, kind)]
        if not args.no_save:
            save_plot(pd_, s, plotdir, suffixes[kind])
        chi2, ndof, pval = mcd.chi2_ndof_pvalue(pd_)
        rows.append(dict(detector=detector, selection=s["sel"], removed_cut=s["key"],
                         var=s["var"], breakdown=kind, norm="shape",
                         chi2=chi2, ndof=ndof, pvalue=pval))
    chi2df = pd.DataFrame(rows)
    log("\n" + chi2df[chi2df.breakdown == "mode"].to_string(index=False, float_format="%.3g"))
    log("\nSaved %d figures (pdf+png) to %s" % (0 if args.no_save else len(results), plotdir))

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
    args = p.parse_args(argv)
    # fields mcd.prepare_detector reads
    args.absnorm = False  # shape-normalized only
    args.cryostat = None
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
    return args


def main(argv=None):
    args = parse_args(argv)

    if len(args.detector) == 1:
        return run_detector(args.detector[0], args)

    passthrough = ["--df-dir", args.df_dir, "--plotdir", args.plotdir, "--nproc", str(args.nproc)]
    for flag in ("no_dirt", "no_save"):
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
