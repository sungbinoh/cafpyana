#!/usr/bin/env python3
"""Data / MC dQ/dx comparisons for calibration checks: the median dQ/dx in
bins of x, y, z, track angle and residual range, on-beam data against the
prediction (MC + off-beam, weighted as plot_data_mc.py), with a data /
prediction ratio panel. What a pattern points to:
  flat offset               gain / overall charge scale
  slope in drift distance   electron lifetime (|x| = 0 is the cathode)
  structure in y or z       YZ response (dead / noisy wires, field distortions)
  trend with the angle      recombination angle dependence, wire response
  offset growing with dQ/dx recombination (protons vs MIPs)
  one TPC only              that TPC's conditions (channels, response)

Candidates: --stage, default "muon" (cuts.CHAINS["presel+muon"]:
preselection and a muon-labelled track > 30 cm, to isolate muons);
"muon_nofv" without the FV, "cosmic" with the Wire-Cell cosmic tagger fired
too, or a plot_data_mc.py stage. summary.txt has the MC truth make-up of the
candidates.

Samples (SAMPLES):
  from the trkseg table (make_wc_candidates.make_wc_trkseg_df: one row per
  4 cm of each track, position and direction of that piece):
    mu_mip       muon-labelled tracks, residual range > RR_MIP_MIN (MIP region)
    proton       proton-labelled tracks
    all_tracks   every track
  from the candidate table (any dataframes):
    trk1_vtx     longest track median dQ/dx at the reco vertex (no direction)
    mu1_angle    longest muon-labelled track median dQ/dx, with the direction
                 of the cosmic tagger's main track (cosmict_2_theta / _phi:
                 polar angle to z, azimuth from x; the same track as the
                 trkseg direction to |cos| > 0.99 for 3/4 of the candidates),
                 position at the reco vertex

Plots, per sample: dqdx_<sample>_dist.png (dQ/dx distribution, prediction
scaled to the data), dqdx_<sample>_<coordinate>.png (median profile, with
the MC-only and off-beam-only medians as thin lines) and
dqdx_<sample>_<coordinate>_bytpc.png (the profile and ratio for each TPC:
east x < 0, west x > 0; by the piece's x for trkseg samples, by the reco
vertex x for the candidate-level ones, so a track crossing the cathode
counts in its vertex TPC). Errors on a median: 1.253 x (16-84% half width) /
sqrt(effective entries); track pieces of one track are correlated, so these
are underestimates. dQ/dx in MIPs: q / nq / DQDX_MIP of the Wire-Cell track
points. Nothing is calibrated (see wc_calo.py).

  python analysis_village/SinglePhotonWC/plot_dqdx_data_mc.py \\
      --mc <MC dir or .df> --onbeam <on-beam> --offbeam <off-beam> -o dqdx_plots/ [--area-norm]
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from analysis_village.SinglePhotonWC import cuts, dfio  # noqa: E402
from analysis_village.SinglePhotonWC.cutflow import load as load_mc  # noqa: E402
from analysis_village.SinglePhotonWC.make_wc_candidates import DQDX_MIP, PID  # noqa: E402
from analysis_village.SinglePhotonWC.plot_data_mc import (STAGES, UBOONE_BDT_DIR, OFFBEAM_COLOR,  # noqa: E402
                                                          load_data, normalisation)

plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), "dune.mplstyle"))

_mu = "μ track > %g cm" % cuts.MUON_TRACK_LEN_MIN
STAGES = {"muon": (cuts.CHAINS["presel+muon"], "preselection + " + _mu),
          "muon_nofv": (cuts.CHAINS["presel+muon_nofv"], "preselection without the FV + " + _mu),
          "cosmic": (cuts.CHAINS["presel+cosmic"], "preselection + cosmic tag + " + _mu),
          **STAGES}
SEG_KEY = "trkseg"
CAND_LEVELS = ["file", "__ntuple", "entry", "nu_index"]
RR_MIP_MIN = 30.0   # cm, residual range above which a muon is a MIP (no Bragg rise)
MC_COLOR, PRED_COLOR = "#2a78d6", "#1f1f1e"
TPCS = [   # name, selection on x, colour
    ("east TPC (x < 0)", lambda x: x < 0, "#2a78d6"),
    ("west TPC (x > 0)", lambda x: x > 0, "#eb6834"),
]


def _angle_drift(dx):
    return np.degrees(np.arccos(np.clip(np.abs(dx), 0, 1)))


def _angle_yz(dy, dz):
    """Angle in the wire (y-z) plane to the beam axis, the line folded to [-90, 90] deg."""
    return np.degrees(np.arctan(dy / np.where(dz == 0, 1e-9, dz)))


def _segments(seg, m):
    s = seg[m(seg)]
    return pd.DataFrame({"dqdx": s.dqdx / DQDX_MIP, "x": s.x, "y": s.y, "z": s.z, "angle_drift": _angle_drift(s.dir_x),
                         "angle_yz": _angle_yz(s.dir_y, s.dir_z), "rr": s.rr}, index=s.index)


def _trk1_vtx(w):
    s = w[w.trk1_dqdx.notna() & (w.trk1_len > 10)]
    return pd.DataFrame({"dqdx": s.trk1_dqdx / DQDX_MIP, "x": s.nu_x, "y": s.nu_y, "z": s.nu_z}, index=s.index)


def _mu1_angle(w):
    s = w[w.mu1_dqdx.notna() & (w.cosmict_2_filled == 1)]
    th, ph = np.radians(s.cosmict_2_theta), np.radians(s.cosmict_2_phi)
    dx, dy, dz = np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)
    return pd.DataFrame({"dqdx": s.mu1_dqdx / DQDX_MIP, "x": s.nu_x, "y": s.nu_y, "z": s.nu_z,
                         "angle_drift": _angle_drift(dx), "angle_yz": _angle_yz(dy, dz)}, index=s.index)


SAMPLES = {   # name: (needs the trkseg table, function of (candidates, segments) -> points, title, dQ/dx bins)
    "mu_mip": (True, lambda w, s: _segments(s, lambda d: (d.pdg == PID["muon"]) & (d.rr > RR_MIP_MIN)),
               "muon-labelled tracks, residual range > %g cm (4 cm pieces)" % RR_MIP_MIN, np.arange(0, 3.01, 0.05)),
    "proton": (True, lambda w, s: _segments(s, lambda d: d.pdg == PID["proton"]),
               "proton-labelled tracks (4 cm pieces)", np.arange(0, 8.01, 0.2)),
    "all_tracks": (True, lambda w, s: _segments(s, lambda d: d.pdg > 0), "all tracks (4 cm pieces)", np.arange(0, 6.01, 0.1)),
    "trk1_vtx": (False, lambda w, s: _trk1_vtx(w), "longest track (> 10 cm) median dQ/dx, at the reco vertex",
                 np.arange(0, 6.01, 0.1)),
    "mu1_angle": (False, lambda w, s: _mu1_angle(w), "longest μ-labelled track median dQ/dx, cosmic-tagger track direction",
                  np.arange(0, 3.01, 0.05)),
}

COORDS = [   # column, bins, label
    ("x", np.linspace(-200, 200, 21), "x [cm] (drift; cathode at 0, anodes at ±200)"),
    ("y", np.linspace(-200, 200, 21), "y [cm]"),
    ("z", np.linspace(0, 500, 26), "z [cm]"),
    ("angle_drift", np.linspace(0, 90, 19), "angle to the drift (x) direction [deg]"),
    ("angle_yz", np.linspace(-90, 90, 19), "angle in the y-z (wire) plane to the beam axis [deg]"),
    ("rr", np.linspace(0, 150, 31), "residual range [cm]"),
]


def wquantile(v, w, q):
    o = np.argsort(v)
    c = np.cumsum(w[o])
    return v[o][np.minimum(np.searchsorted(c, np.asarray(q) * c[-1]), len(v) - 1)]


def median_err(v, w):
    """(weighted median, its error), NaN for fewer than 5 entries."""
    if len(v) < 5:
        return np.nan, np.nan
    q16, q50, q84 = wquantile(v, w, [0.16, 0.5, 0.84])
    neff = w.sum() ** 2 / (w ** 2).sum()
    return q50, 1.253 * (q84 - q16) / 2 / np.sqrt(neff)


def profile(p, coord, bins):
    """Median dQ/dx and its error per coordinate bin."""
    b = np.digitize(p[coord].to_numpy(float), bins) - 1
    v, w = p.dqdx.to_numpy(float), p.w.to_numpy(float)
    out = np.full((len(bins) - 1, 2), np.nan)
    for i in range(len(bins) - 1):
        m = (b == i) & np.isfinite(v)
        out[i] = median_err(v[m], w[m])
    return out


def plot_dist(name, title, bins, on, mc, off, out):
    hd = np.histogram(on.dqdx, bins)[0]
    hm = np.histogram(mc.dqdx, bins, weights=mc.w)[0]
    ho = np.histogram(off.dqdx, bins, weights=off.w)[0]
    k = hd.sum() / max(hm.sum() + ho.sum(), 1e-12)   # shape comparison
    c = 0.5 * (bins[1:] + bins[:-1])
    fig, (a, r) = plt.subplots(2, 1, figsize=(10, 7), sharex=True, gridspec_kw=dict(height_ratios=[3, 1], hspace=0.05))
    a.hist([c, c], bins=bins, weights=[ho * k, hm * k], stacked=True, color=[OFFBEAM_COLOR, MC_COLOR],
           label=["off-beam data", "MC (ν + cosmics in the spill)"], edgecolor="white", linewidth=0.4)
    a.errorbar(c, hd, np.sqrt(hd), fmt="o", color="k", ms=4, label="on-beam data (%d)" % hd.sum())
    pred = (hm + ho) * k
    ok = pred > 0
    r.errorbar(c[ok], hd[ok] / pred[ok], np.sqrt(hd[ok]) / pred[ok], fmt="o", color="k", ms=3)
    r.axhline(1, color="k", lw=1)
    r.set_ylim(0, 2)
    r.set_ylabel("data / pred.")
    r.set_xlabel("dQ/dx [MIP]")
    a.set_ylabel("entries (pred. scaled to data)")
    md = median_err(on.dqdx.to_numpy(float), on.w.to_numpy(float))[0]
    mp = median_err(np.r_[mc.dqdx, off.dqdx], np.r_[mc.w, off.w])[0]
    a.set_title("%s\nmedian: data %.3f, prediction %.3f, data / pred. %.3f" % (title, md, mp, md / mp), loc="left", fontsize=11)
    a.legend(frameon=False, fontsize=10)
    a.set_ylim(bottom=0)
    fig.savefig(os.path.join(out, "dqdx_%s_dist.png" % name), dpi=130, bbox_inches="tight")
    plt.close(fig)
    return md, mp


def plot_profile(name, title, coord, bins, label, on, mc, off, out):
    pred = pd.concat([mc, off])
    P = {k: profile(d, coord, bins) for k, d in [("on", on), ("pred", pred), ("mc", mc), ("off", off)]}
    c = 0.5 * (bins[1:] + bins[:-1])
    xe = np.diff(bins) / 2
    fig, (a, r) = plt.subplots(2, 1, figsize=(10, 7), sharex=True, gridspec_kw=dict(height_ratios=[3, 1], hspace=0.05))
    a.stairs(P["mc"][:, 0], bins, baseline=None, color=MC_COLOR, lw=1.2, label="MC only")
    a.stairs(P["off"][:, 0], bins, baseline=None, color=OFFBEAM_COLOR, lw=1.2, label="off-beam only")
    lo, hi = P["pred"][:, 0] - P["pred"][:, 1], P["pred"][:, 0] + P["pred"][:, 1]
    a.fill_between(bins, np.r_[lo, lo[-1]], np.r_[hi, hi[-1]], step="post", color="#c9c8c3", alpha=0.7, lw=0)
    a.stairs(P["pred"][:, 0], bins, baseline=None, color=PRED_COLOR, lw=2, label="prediction (MC + off-beam)")
    a.errorbar(c, P["on"][:, 0], P["on"][:, 1], xerr=xe, fmt="o", color="k", ms=4, label="on-beam data")
    ratio = P["on"][:, 0] / P["pred"][:, 0]
    rerr = ratio * np.hypot(P["on"][:, 1] / P["on"][:, 0], P["pred"][:, 1] / P["pred"][:, 0])
    r.errorbar(c, ratio, rerr, xerr=xe, fmt="o", color="k", ms=3)
    r.axhline(1, color="k", lw=1)
    lim = np.nanmax(np.abs(ratio - 1) + np.nan_to_num(rerr)) if np.isfinite(ratio).any() else 0.1
    r.set_ylim(1 - max(0.1, 1.2 * lim), 1 + max(0.1, 1.2 * lim))
    r.set_ylabel("data / pred.")
    r.set_xlabel(label)
    r.set_xlim(bins[0], bins[-1])
    a.set_ylabel("median dQ/dx [MIP]")
    a.set_title(title, loc="left", fontsize=11)
    a.legend(frameon=False, fontsize=9.5)
    fig.savefig(os.path.join(out, "dqdx_%s_%s.png" % (name, coord)), dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_profile_tpc(name, title, coord, bins, label, on, mc, off, out):
    """The median profile and data / prediction ratio for each TPC."""
    pred = pd.concat([mc, off])
    c = 0.5 * (bins[1:] + bins[:-1])
    xe = np.diff(bins) / 2
    fig, (a, r) = plt.subplots(2, 1, figsize=(10, 7), sharex=True, gridspec_kw=dict(height_ratios=[3, 1], hspace=0.05))
    lims = [0.1]
    for i, (tpc, sel, col) in enumerate(TPCS):
        d, p = on[sel(on.x)], pred[sel(pred.x)]
        Pd, Pp = profile(d, coord, bins), profile(p, coord, bins)
        lo, hi = Pp[:, 0] - Pp[:, 1], Pp[:, 0] + Pp[:, 1]
        a.fill_between(bins, np.r_[lo, lo[-1]], np.r_[hi, hi[-1]], step="post", color=col, alpha=0.2, lw=0)
        a.stairs(Pp[:, 0], bins, baseline=None, color=col, lw=2, label="%s: prediction" % tpc)
        a.errorbar(c, Pd[:, 0], Pd[:, 1], xerr=xe, fmt="o", color=col, mec="k", ms=4.5, label="%s: on-beam data" % tpc)
        ratio = Pd[:, 0] / Pp[:, 0]
        rerr = ratio * np.hypot(Pd[:, 1] / Pd[:, 0], Pp[:, 1] / Pp[:, 0])
        sh = (i - 0.5) * 0.3 * xe   # side by side
        r.errorbar(c + sh, ratio, rerr, fmt="o", color=col, mec="k", ms=3.5)
        if np.isfinite(ratio).any():
            lims.append(1.2 * np.nanmax(np.abs(ratio - 1) + np.nan_to_num(rerr)))
    r.axhline(1, color="k", lw=1)
    r.set_ylim(1 - max(lims), 1 + max(lims))
    r.set_ylabel("data / pred.")
    r.set_xlabel(label)
    r.set_xlim(bins[0], bins[-1])
    a.set_ylabel("median dQ/dx [MIP]")
    a.set_title(title, loc="left", fontsize=11)
    a.legend(frameon=False, fontsize=9.5)
    fig.savefig(os.path.join(out, "dqdx_%s_%s_bytpc.png" % (name, coord)), dpi=130, bbox_inches="tight")
    plt.close(fig)


def tpc_medians(on, mc, off):
    """'east data/pred r, west data/pred r' of the median dQ/dx."""
    pred = pd.concat([mc, off])
    out = []
    for tpc, sel, _ in TPCS:
        d, p = on[sel(on.x)], pred[sel(pred.x)]
        md = median_err(d.dqdx.to_numpy(float), d.w.to_numpy(float))[0]
        mp = median_err(p.dqdx.to_numpy(float), p.w.to_numpy(float))[0]
        out.append("%s %.3f" % (tpc.split()[0], md / mp))
    return ", ".join(out)


def load_segments(inputs, keep, nproc=None):
    """The trkseg table for the candidates in keep (an index over CAND_LEVELS);
    None if the inputs have none."""
    files = dfio.expand_inputs(inputs)
    if SEG_KEY not in dfio.file_keys(files[0]):
        return None

    def select(k, d):
        if k != SEG_KEY or not len(d):
            return d
        return d[d.index.droplevel(["kine_index", "seg"]).isin(keep)]
    frames, _ = dfio.load(files, keys=[SEG_KEY], nproc=nproc, select=select)
    return frames[SEG_KEY]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mc", nargs="+", required=True, help="MC .df files, directories or globs")
    ap.add_argument("--onbeam", nargs="+", required=True, help="on-beam data .df files, directories or globs")
    ap.add_argument("--offbeam", nargs="+", required=True, help="off-beam data .df files, directories or globs")
    ap.add_argument("-o", "--outdir", default="dqdx_plots")
    ap.add_argument("--stage", default="muon", choices=list(STAGES), help="candidates used (default: muon)")
    ap.add_argument("--area-norm", action="store_true", help="without POT / gate information: fit the off-beam share "
                    "(as plot_data_mc.py); only the MC / off-beam mix matters for the medians")
    ap.add_argument("--offbeam-frac", type=float, help="off-beam share of the prediction at preselection, with --area-norm")
    ap.add_argument("--samples", nargs="+", default=list(SAMPLES), choices=list(SAMPLES))
    ap.add_argument("-j", "--nproc", type=int, help="worker processes (default: up to 8)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    bdt = args.stage not in ("presel", "muon", "muon_nofv", "cosmic")   # the later stages cut on the uB scores
    mc, _, pot_mc = load_mc(args.mc, UBOONE_BDT_DIR if bdt else None, nproc=args.nproc)
    on, pot_on, _ = load_data(args.onbeam, args.nproc, uboone_bdt=bdt)
    off, _, gates_off = load_data(args.offbeam, args.nproc, uboone_bdt=bdt)
    sc_mc, sc_off, norm_line = normalisation(mc, on, off, pot_mc, pot_on, gates_off, args.area_norm, args.offbeam_frac)
    chain, stage_title = STAGES[args.stage]
    cands = {k: d[cuts.passes(d, chain)] for k, d in [("mc", mc), ("on", on), ("off", off)]}
    weight = {"mc": sc_mc, "on": 1.0, "off": sc_off}
    segs = {}
    if any(SAMPLES[s][0] for s in args.samples):
        inputs = {"mc": args.mc, "on": args.onbeam, "off": args.offbeam}
        segs = {k: load_segments(inputs[k], cands[k].index, args.nproc) for k in cands}
        if any(v is None for v in segs.values()):
            print("no %s table in some inputs (dataframes made before make_wc_trkseg_df): track-piece samples skipped" % SEG_KEY)
            segs = {}
    comp = cands["mc"].cat.map(cuts.CATEGORY_LABELS).value_counts(normalize=True)
    lines = [norm_line, "stage: %s (%d MC candidates; MC truth: %s)" % (
        stage_title, len(cands["mc"]), ", ".join("%s %.0f%%" % (k, 100 * v) for k, v in comp.head(5).items())), ""]
    for name in args.samples:
        need_seg, fn, title, bins = SAMPLES[name]
        if need_seg and not segs:
            continue
        pts = {k: fn(cands[k], segs.get(k)).assign(w=weight[k]) for k in cands}
        pts = {k: p[np.isfinite(p.dqdx)] for k, p in pts.items()}
        full = "%s, %s (%d data / %d MC / %d off-beam entries)" % (title, stage_title, len(pts["on"]), len(pts["mc"]), len(pts["off"]))
        md, mp = plot_dist(name, full, bins, pts["on"], pts["mc"], pts["off"], args.outdir)
        lines.append("%-12s median dQ/dx data %.3f, prediction %.3f, ratio %.3f (by TPC: %s)"
                     % (name, md, mp, md / mp, tpc_medians(pts["on"], pts["mc"], pts["off"])))
        for coord, cbins, label in COORDS:
            if coord in pts["on"]:
                plot_profile(name, full, coord, cbins, label, pts["on"], pts["mc"], pts["off"], args.outdir)
                if coord != "x":
                    plot_profile_tpc(name, full, coord, cbins, label, pts["on"], pts["mc"], pts["off"], args.outdir)
    open(os.path.join(args.outdir, "summary.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("wrote the plots to %s" % args.outdir)


if __name__ == "__main__":
    main()
