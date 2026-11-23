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
  from the pandtrk table (make_pandora_hits.py: Pandora tracks > 30 cm,
  per plane; every event, the --stage cut does not apply): Pandora
  muon-like tracks (track score > 0.5, chi2_mu < chi2_p), per plane p0-p2:
    pandora_dqdx_p<n>         median hit dQ/dx (the CAF's calibrated value)
    pandora_width_p<n>        median hit width [ticks]
    pandora_mult_p<n>         mean hit multiplicity
    pandora_int_over_sum_p<n> median hit integral / ADC sum
  all over hits with residual range > 30 cm; Pandora's hits come from the
  same signal-processing output as Wire-Cell's imaging
  from the candidate table (any dataframes):
    trk1_vtx     longest track median dQ/dx at the reco vertex (no direction)
    mu1_angle    longest muon-labelled track median dQ/dx, with the direction
                 of the cosmic tagger's main track (cosmict_2_theta / _phi:
                 polar angle to z, azimuth from x; the same track as the
                 trkseg direction to |cos| > 0.99 for 3/4 of the candidates),
                 position at the reco vertex

Angles: --angle-bin wide bins (default 5 deg) in the profiles; the
statistics at small angles go as sin(angle). Each angle profile also has
the data and MC entries per bin in a third panel.

Recombination: for the samples with a track angle, the MC is also shown
with its dQ/dx re-derived under the SBND data recombination parameters
(wc_calo.swap_recombination: measured dQ/dx -> dE/dx with the MC set ->
dQ/dx with the data set, at the track's angle to the drift field), as a
second prediction ("MC, data recombination"; off-beam unchanged). What is
left of the data / MC difference after it is not recombination. Off with
--no-recomb-swap.

--fit-recomb: for the samples with a track angle, fits the median dQ/dx vs
angle to drift (--angle-bin bins) with the MC swapped to (alpha, beta_90)
of the data set and a free ellipsoid ratio R, times a free MC charge scale
k (off-beam unchanged): chi2 on a grid of R x k, minimum and the Delta chi2
= 1 range of R (profiled over k; the bin errors ignore correlations, so the
range is indicative). dqdx_<sample>_recomb_fit.png shows the data, the
nominal MC, the MC with the data set, the best fit and chi2(R); the numbers
go to summary.txt. A best R far from the measured 1.57, or a poor best chi2,
says the angle dependence is not (only) recombination.

Plots, per sample: dqdx_<sample>_<coordinate>.png, the median profile with
the MC-only median as a thin line, for --coords (default angle_drift only;
also x y z angle_yz rr). Optional: --dist adds
dqdx_<sample>_dist.png (dQ/dx distribution, prediction scaled to the data),
--by-tpc adds dqdx_<sample>_<coordinate>_bytpc.png (profile and ratio for
each TPC: east x < 0, west x > 0; by the piece's x for trkseg samples, by
the reco vertex x for the candidate-level ones, so a track crossing the
cathode counts in its vertex TPC). summary.txt always has the overall
medians. Errors on a median: 1.253 x (16-84% half width) /
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
from analysis_village.SinglePhotonWC import cuts, dfio, wc_calo  # noqa: E402
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
PAND_KEY = "pandtrk"
SEG_LEVELS = ["file", "__ntuple", "entry", "nu_index", "kine_index", "seg"]   # after dfio.load
PAND_LEVELS = ["file", "__ntuple", "entry", "rec.slc..index", "rec.slc.reco.pfp..index", "plane"]
CAND_LEVELS = ["file", "__ntuple", "entry", "nu_index"]
RR_MIP_MIN = 30.0   # cm, residual range above which a muon is a MIP (no Bragg rise)
MC_COLOR, PRED_COLOR, ALT_COLOR = "#2a78d6", "#1f1f1e", "#1baf7a"
ALT_LABEL = "prediction, MC with data recombination"
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


def _pandora(t, plane, q):
    """Pandora muon-like tracks (track score > 0.5, collection-plane chi2_mu
    < chi2_p) on one plane: quantity q (column of the pandtrk table)."""
    s = t[(t.index.get_level_values("plane") == plane) & (t.track_score > 0.5) & (t.chi2_mu < t.chi2_p)]
    return pd.DataFrame({"dqdx": s[q], "x": s.x, "y": 0.5 * (s.start_y + s.end_y), "z": 0.5 * (s.start_z + s.end_z),
                         "angle_drift": _angle_drift(s.dir_x), "angle_yz": _angle_yz(s.dir_y, s.dir_z)}, index=s.index)


# name: table needed (None: the candidates only), function of (candidates, table) -> points (the plotted
# quantity in column "dqdx"), title, bins of its distribution, y label, Wire-Cell dQ/dx in MIPs (recombination
# swap and fit apply)
SAMPLES = {
    "mu_mip": dict(table=SEG_KEY, fn=lambda w, s: _segments(s, lambda d: (d.pdg == PID["muon"]) & (d.rr > RR_MIP_MIN)),
                   title="muon-labelled tracks, residual range > %g cm (4 cm pieces)" % RR_MIP_MIN,
                   bins=np.arange(0, 3.01, 0.05)),
    "proton": dict(table=SEG_KEY, fn=lambda w, s: _segments(s, lambda d: d.pdg == PID["proton"]),
                   title="proton-labelled tracks (4 cm pieces)", bins=np.arange(0, 8.01, 0.2)),
    "all_tracks": dict(table=SEG_KEY, fn=lambda w, s: _segments(s, lambda d: d.pdg > 0), title="all tracks (4 cm pieces)",
                       bins=np.arange(0, 6.01, 0.1)),
    "trk1_vtx": dict(table=None, fn=lambda w, s: _trk1_vtx(w),
                     title="longest track (> 10 cm) median dQ/dx, at the reco vertex", bins=np.arange(0, 6.01, 0.1)),
    "mu1_angle": dict(table=None, fn=lambda w, s: _mu1_angle(w),
                      title="longest μ-labelled track median dQ/dx, cosmic-tagger track direction",
                      bins=np.arange(0, 3.01, 0.05)),
}
for _pl in range(3):   # Pandora hits (pandtrk table), every event (not --stage), per plane
    for _q, _lab, _bins in [("dqdx_mip", "hit dQ/dx [CAF units]", np.arange(0, 3001, 50)),
                            ("width_mip", "hit width [ticks]", np.arange(0, 15.01, 0.25)),
                            ("mult_mip", "hit multiplicity", np.arange(1, 3.01, 0.05)),
                            ("int_over_sum_mip", "hit integral / ADC sum", np.arange(0.8, 1.201, 0.01))]:
        SAMPLES["pandora_%s_p%d" % (_q.replace("_mip", ""), _pl)] = dict(
            table=PAND_KEY, fn=(lambda pl, q: lambda w, t: _pandora(t, pl, q))(_pl, _q),
            title="Pandora muon-like tracks > 30 cm, plane %d: per-track %s %s (hits with rr > 30 cm)"
                  % (_pl, "mean" if _q == "mult_mip" else "median", _lab),
            bins=_bins, ylabel=("mean " if _q == "mult_mip" else "median ") + _lab, wc=False)
DEFAULT_SAMPLES = ["mu_mip", "proton", "all_tracks", "mu1_angle"] + [   # trk1_vtx has no angle: --samples trk1_vtx --coords x
    "pandora_%s_p%d" % (q, pl) for q in ("dqdx", "width", "mult", "int_over_sum") for pl in range(3)]

COORDS = [   # column, bins, label
    ("x", np.linspace(-200, 200, 21), "x [cm] (drift; cathode at 0, anodes at ±200)"),
    ("y", np.linspace(-200, 200, 21), "y [cm]"),
    ("z", np.linspace(0, 500, 26), "z [cm]"),
    ("angle_drift", None, "angle to the drift (x) direction [deg]"),   # angle bins: --angle-bin
    ("angle_yz", None, "angle in the y-z (wire) plane to the beam axis [deg]"),
    ("rr", np.linspace(0, 150, 31), "residual range [cm]"),
]
ANGLE_RANGE = {"angle_drift": (0, 90), "angle_yz": (-90, 90)}


def angle_bins(coord, width):
    lo, hi = ANGLE_RANGE[coord]
    return np.linspace(lo, hi, int(round((hi - lo) / width)) + 1)


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


def plot_dist(name, title, bins, on, mc, off, out, mc_alt=None):
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
    alt = ""
    if mc_alt is not None:
        ha = np.histogram(mc_alt.dqdx, bins, weights=mc_alt.w)[0]
        pa = (ha + ho) * hd.sum() / max(ha.sum() + ho.sum(), 1e-12)
        a.stairs(pa, bins, color=ALT_COLOR, lw=2, label=ALT_LABEL)
        oka = pa > 0
        r.errorbar(c[oka], hd[oka] / pa[oka], np.sqrt(hd[oka]) / pa[oka], fmt="s", color=ALT_COLOR, ms=3)
        ma = median_err(np.r_[mc_alt.dqdx, off.dqdx], np.r_[mc_alt.w, off.w])[0]
        alt = "; with data recombination %.3f, data / pred. %.3f" % (ma, median_err(on.dqdx.to_numpy(float), on.w.to_numpy(float))[0] / ma)
    r.axhline(1, color="k", lw=1)
    r.set_ylim(0, 2)
    r.set_ylabel("data / pred.")
    r.set_xlabel("dQ/dx [MIP]")
    a.set_ylabel("entries (pred. scaled to data)")
    md = median_err(on.dqdx.to_numpy(float), on.w.to_numpy(float))[0]
    mp = median_err(np.r_[mc.dqdx, off.dqdx], np.r_[mc.w, off.w])[0]
    a.set_title("%s\nmedian: data %.3f, prediction %.3f, data / pred. %.3f%s" % (title, md, mp, md / mp, alt), loc="left", fontsize=11)
    a.legend(frameon=False, fontsize=10)
    a.set_ylim(bottom=0)
    fig.savefig(os.path.join(out, "dqdx_%s_dist.png" % name), dpi=130, bbox_inches="tight")
    plt.close(fig)
    return md, mp


def plot_profile(name, title, coord, bins, label, on, mc, off, out, mc_alt=None, ylabel="median dQ/dx [MIP]"):
    pred = pd.concat([mc, off])
    P = {k: profile(d, coord, bins) for k, d in [("on", on), ("pred", pred), ("mc", mc)]}
    c = 0.5 * (bins[1:] + bins[:-1])
    xe = np.diff(bins) / 2
    counts = coord.startswith("angle")
    fig, axs = plt.subplots(3 if counts else 2, 1, figsize=(10, 8.6 if counts else 7), sharex=True,
                            gridspec_kw=dict(height_ratios=[3, 1, 1] if counts else [3, 1], hspace=0.05))
    a, r = axs[0], axs[1]
    if counts:   # entries per bin: how far a feature in one bin can be trusted
        n = axs[2]
        nd = np.histogram(on[coord], bins)[0]
        nm = np.histogram(mc[coord], bins)[0]
        n.stairs(nd, bins, color="k", lw=1.5, label="data entries")
        n.stairs(nm * nd.sum() / max(nm.sum(), 1), bins, color=MC_COLOR, lw=1.5, label="MC entries (scaled to data)")
        n.set_ylabel("entries")
        n.set_ylim(bottom=0)
        n.legend(frameon=False, fontsize=8.5)
        n.set_xlabel(label)
    a.stairs(P["mc"][:, 0], bins, baseline=None, color=MC_COLOR, lw=1.2, label="MC only")
    lo, hi = P["pred"][:, 0] - P["pred"][:, 1], P["pred"][:, 0] + P["pred"][:, 1]
    a.fill_between(bins, np.r_[lo, lo[-1]], np.r_[hi, hi[-1]], step="post", color="#c9c8c3", alpha=0.7, lw=0)
    a.stairs(P["pred"][:, 0], bins, baseline=None, color=PRED_COLOR, lw=2, label="prediction (MC + off-beam)")
    a.errorbar(c, P["on"][:, 0], P["on"][:, 1], xerr=xe, fmt="o", color="k", ms=4, label="on-beam data")
    ratio = P["on"][:, 0] / P["pred"][:, 0]
    rerr = ratio * np.hypot(P["on"][:, 1] / P["on"][:, 0], P["pred"][:, 1] / P["pred"][:, 0])
    r.errorbar(c, ratio, rerr, xerr=xe, fmt="o", color="k", ms=3)
    lims = [np.nanmax(np.abs(ratio - 1) + np.nan_to_num(rerr)) if np.isfinite(ratio).any() else 0.1]
    if mc_alt is not None:
        Pa = profile(pd.concat([mc_alt, off]), coord, bins)
        a.stairs(Pa[:, 0], bins, baseline=None, color=ALT_COLOR, lw=2, ls="--", label=ALT_LABEL)
        ra = P["on"][:, 0] / Pa[:, 0]
        rae = ra * np.hypot(P["on"][:, 1] / P["on"][:, 0], Pa[:, 1] / Pa[:, 0])
        r.errorbar(c + 0.25 * xe, ra, rae, fmt="s", color=ALT_COLOR, ms=3)
        if np.isfinite(ra).any():
            lims.append(np.nanmax(np.abs(ra - 1) + np.nan_to_num(rae)))
    r.axhline(1, color="k", lw=1)
    lim = max(lims)
    r.set_ylim(1 - max(0.1, 1.2 * lim), 1 + max(0.1, 1.2 * lim))
    r.set_ylabel("data / pred.")
    if not counts:
        r.set_xlabel(label)
    r.set_xlim(bins[0], bins[-1])
    a.set_ylabel(ylabel)
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


def _fast_profile(v, w, b, nb):
    """Weighted median and its error per bin index b (0..nb-1)."""
    out = np.full((nb, 2), np.nan)
    o = np.lexsort((v, b))
    v, w, b = v[o], w[o], b[o]
    edges = np.searchsorted(b, np.arange(nb + 1))
    for i in range(nb):
        lo, hi = edges[i], edges[i + 1]
        if hi - lo < 5:
            continue
        vi, wi = v[lo:hi], w[lo:hi]
        c = np.cumsum(wi)
        q16, q50, q84 = vi[np.minimum(np.searchsorted(c, np.array([0.16, 0.5, 0.84]) * c[-1]), len(vi) - 1)]
        out[i] = q50, 1.253 * (q84 - q16) / 2 / np.sqrt(c[-1] ** 2 / (wi ** 2).sum())
    return out


R_GRID = np.round(np.arange(1.0, 3.001, 0.05), 3)
K_GRID = np.round(np.arange(0.90, 1.1001, 0.005), 4)


def fit_recombination(on, mc, off, bins):
    """chi2 of the median dQ/dx vs angle to drift over (R, k): the MC
    swapped to the data (alpha, beta_90) with ellipsoid ratio R, times k.
    Returns a dict of the grid, chi2 and the best point."""
    nb = len(bins) - 1
    idx = lambda d: np.clip(np.digitize(d.angle_drift.to_numpy(float), bins) - 1, -1, nb)
    bd, bm, bo = idx(on), idx(mc), idx(off)
    Pd = _fast_profile(on.dqdx.to_numpy(float), on.w.to_numpy(float), bd, nb)
    phi = np.radians(mc.angle_drift.to_numpy(float))
    q = mc.dqdx.to_numpy(float) * DQDX_MIP
    a, b90, _ = wc_calo.RECOMB["data"]
    vo, wo = off.dqdx.to_numpy(float), off.w.to_numpy(float)
    wm = mc.w.to_numpy(float)
    keep = lambda b: (b >= 0) & (b < nb)
    chi2 = np.full((len(R_GRID), len(K_GRID)), np.nan)
    for i, R in enumerate(R_GRID):
        vR = wc_calo.dedx_to_dqdx(wc_calo.dqdx_to_dedx(q, phi, recomb=wc_calo.RECOMB["mc"]), phi, recomb=(a, b90, R)) / DQDX_MIP
        for j, k in enumerate(K_GRID):
            v = np.r_[vR * k, vo]
            w = np.r_[wm, wo]
            bb = np.r_[bm, bo]
            m = keep(bb) & np.isfinite(v)
            Pp = _fast_profile(v[m], w[m], bb[m], nb)
            ok = np.isfinite(Pd[:, 0]) & np.isfinite(Pp[:, 0])
            chi2[i, j] = np.sum(((Pd[ok, 0] - Pp[ok, 0]) / np.hypot(Pd[ok, 1], Pp[ok, 1])) ** 2)
    i, j = np.unravel_index(np.nanargmin(chi2), chi2.shape)
    prof = np.nanmin(chi2, axis=1)   # profiled over k
    inside = R_GRID[prof <= prof[i] + 1]
    return dict(chi2=chi2, prof=prof, R=R_GRID[i], k=K_GRID[j], chi2_min=chi2[i, j], ndf=int(ok.sum()) - 2,
                R_range=(inside.min(), inside.max()), Pd=Pd)


def plot_recomb_fit(name, title, bins, on, mc, off, out, fit):
    a, b90, _ = wc_calo.RECOMB["data"]
    c = 0.5 * (bins[1:] + bins[:-1])
    swap = lambda R, k=1.0: mc.assign(dqdx=k * wc_calo.dedx_to_dqdx(
        wc_calo.dqdx_to_dedx(mc.dqdx * DQDX_MIP, np.radians(mc.angle_drift), recomb=wc_calo.RECOMB["mc"]),
        np.radians(mc.angle_drift), recomb=(a, b90, R)) / DQDX_MIP)
    lines = [("nominal MC", mc, PRED_COLOR, "-"),
             ("MC with the data recombination (R = %.2f)" % wc_calo.RECOMB["data"][2], swap(wc_calo.RECOMB["data"][2]), ALT_COLOR, "--"),
             ("best fit: R = %.2f, k = %.3f" % (fit["R"], fit["k"]), swap(fit["R"], fit["k"]), "#eb6834", "-")]
    fig = plt.figure(figsize=(16, 7))
    a1 = fig.add_axes([0.06, 0.36, 0.55, 0.56])
    r1 = fig.add_axes([0.06, 0.1, 0.55, 0.24], sharex=a1)
    a2 = fig.add_axes([0.7, 0.1, 0.28, 0.82])
    Pd = fit["Pd"]
    for lab, m_, col, ls in lines:
        P = profile(pd.concat([m_, off]), "angle_drift", bins)
        a1.stairs(P[:, 0], bins, baseline=None, color=col, ls=ls, lw=2, label=lab)
        r1.plot(c, Pd[:, 0] / P[:, 0], "o", color=col, ms=3.5)
    a1.errorbar(c, Pd[:, 0], Pd[:, 1], fmt="o", color="k", ms=4, label="on-beam data")
    a1.set_ylabel("median dQ/dx [MIP]")
    a1.legend(frameon=False, fontsize=9.5)
    a1.set_title(title, loc="left", fontsize=10)
    plt.setp(a1.get_xticklabels(), visible=False)
    r1.axhline(1, color="k", lw=1)
    r1.set_ylabel("data / pred.")
    r1.set_xlabel("angle to the drift (x) direction [deg]")
    r1.set_xlim(bins[0], bins[-1])
    a2.plot(R_GRID, fit["prof"] - fit["chi2_min"], color="#eb6834", lw=2)
    a2.axhline(1, color="#888", ls=":")
    a2.axvline(wc_calo.RECOMB["data"][2], color=ALT_COLOR, ls="--", label="SBND data R")
    a2.axvline(wc_calo.RECOMB["mc"][2], color=PRED_COLOR, ls="--", label="SBND MC R")
    a2.set_ylim(0, min(50, np.nanmax(fit["prof"] - fit["chi2_min"])))
    a2.set_xlabel("ellipsoid ratio R")
    a2.set_ylabel("$\Delta\chi^2$ (k profiled)")
    a2.legend(frameon=False, fontsize=9.5)
    a2.set_title("best R %.2f (%.2f-%.2f), $\chi^2$/ndf %.1f/%d" % (fit["R"], *fit["R_range"], fit["chi2_min"], fit["ndf"]),
                 fontsize=10, loc="left")
    fig.savefig(os.path.join(out, "dqdx_%s_recomb_fit.png" % name), dpi=130, bbox_inches="tight")
    plt.close(fig)


def load_table(inputs, key, keep=None, nproc=None):
    """The trkseg (rows of the candidates in keep, an index over CAND_LEVELS)
    or pandtrk table; None if the inputs have none."""
    files = dfio.expand_inputs(inputs)
    if key not in dfio.file_keys(files[0]):
        return None

    names = {SEG_KEY: SEG_LEVELS, PAND_KEY: PAND_LEVELS}[key]
    lost = []

    def select(k, d):
        if k != key or not len(d):
            return d
        # a dataframe split whose first CAF had no rows was saved without the level names (or, rarely, with
        # the levels collapsed): name them by position, drop the unrecoverable ones
        if d.index.nlevels == len(names):
            d.index = d.index.set_names(names)
        else:
            lost.append(len(d))
            return pd.DataFrame(columns=d.columns, index=pd.MultiIndex.from_arrays([[]] * len(names), names=names))
        if k != SEG_KEY or keep is None:
            return d
        return d[d.index.droplevel(["kine_index", "seg"]).isin(keep)]
    frames, _ = dfio.load(files, keys=[key], nproc=nproc, select=select)
    if lost:
        print("%s: dropped %d rows in %d splits with a broken index" % (key, sum(lost), len(lost)))
    return frames[key]


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
    ap.add_argument("--samples", nargs="+", default=DEFAULT_SAMPLES, choices=list(SAMPLES))
    ap.add_argument("--angle-bin", type=float, default=5.0, help="angle bin width in the profiles [deg] (default 5)")
    ap.add_argument("--coords", nargs="+", default=["angle_drift"], choices=[c[0] for c in COORDS],
                    help="profile coordinates (default: angle_drift)")
    ap.add_argument("--dist", action="store_true", help="also the dQ/dx distribution plots")
    ap.add_argument("--by-tpc", action="store_true", help="also the profiles split by TPC")
    ap.add_argument("--no-recomb-swap", action="store_true", help="no MC prediction with the data recombination")
    ap.add_argument("--fit-recomb", action="store_true", help="fit the ellipsoid ratio R and an MC scale to the angle profile")
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
    inputs = {"mc": args.mc, "on": args.onbeam, "off": args.offbeam}
    tables = {}
    for key in {SAMPLES[n]["table"] for n in args.samples} - {None}:
        t = {k: load_table(inputs[k], key, cands[k].index if key == SEG_KEY else None, args.nproc) for k in cands}
        if any(v is None for v in t.values()):
            print("no %s table in some inputs (older dataframes): its samples are skipped" % key)
        else:
            tables[key] = t
    comp = cands["mc"].cat.map(cuts.CATEGORY_LABELS).value_counts(normalize=True)
    lines = [norm_line, "stage: %s (%d MC candidates; MC truth: %s)" % (
        stage_title, len(cands["mc"]), ", ".join("%s %.0f%%" % (k, 100 * v) for k, v in comp.head(5).items())), ""]
    for name in args.samples:
        spec = SAMPLES[name]
        key, title, bins = spec["table"], spec["title"], spec["bins"]
        if key is not None and key not in tables:
            continue
        pts = {k: spec["fn"](cands[k], tables[key][k] if key else None).assign(w=weight[k]) for k in cands}
        pts = {k: p[np.isfinite(p.dqdx)] for k, p in pts.items()}
        full = "%s, %s (%d data / %d MC / %d off-beam entries)" % (
            title, stage_title if key != PAND_KEY else "all events", len(pts["on"]), len(pts["mc"]), len(pts["off"]))
        on, mc, off = pts["on"], pts["mc"], pts["off"]
        wcq = spec.get("wc", True)
        alt = None
        if wcq and "angle_drift" in mc and not args.no_recomb_swap:
            alt = mc.assign(dqdx=wc_calo.swap_recombination(mc.dqdx * DQDX_MIP, np.radians(mc.angle_drift)) / DQDX_MIP)
        if args.dist:
            plot_dist(name, full, bins, on, mc, off, args.outdir, alt)
        md = median_err(on.dqdx.to_numpy(float), on.w.to_numpy(float))[0]
        mp = median_err(np.r_[mc.dqdx, off.dqdx], np.r_[mc.w, off.w])[0]
        line = "%-22s median data %.4g, prediction %.4g, ratio %.3f (by TPC: %s)" % (
            name, md, mp, md / mp, tpc_medians(on, mc, off))
        if alt is not None:
            ma = median_err(np.r_[alt.dqdx, off.dqdx], np.r_[alt.w, off.w])[0]
            line += "; MC with data recombination: prediction %.3f, ratio %.3f" % (ma, md / ma)
        lines.append(line)
        for coord, cbins, label in COORDS:
            if coord in on and coord in args.coords:
                cbins = angle_bins(coord, args.angle_bin) if cbins is None else cbins
                plot_profile(name, full, coord, cbins, label, on, mc, off, args.outdir, alt,
                             ylabel=spec.get("ylabel", "median dQ/dx [MIP]"))
                if args.by_tpc and coord != "x":
                    plot_profile_tpc(name, full, coord, cbins, label, on, mc, off, args.outdir)
        if "angle_drift" in on and wcq:
            if args.fit_recomb:
                fb = angle_bins("angle_drift", args.angle_bin)
                fit = fit_recombination(on, mc, off, fb)
                plot_recomb_fit(name, full, fb, on, mc, off, args.outdir, fit)
                lines.append("%-12s recombination fit: R = %.2f (Delta chi2 = 1: %.2f-%.2f), MC scale k = %.3f, chi2/ndf = %.1f/%d"
                             % ("", fit["R"], *fit["R_range"], fit["k"], fit["chi2_min"], fit["ndf"]))
    open(os.path.join(args.outdir, "summary.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("wrote the plots to %s" % args.outdir)


if __name__ == "__main__":
    main()
