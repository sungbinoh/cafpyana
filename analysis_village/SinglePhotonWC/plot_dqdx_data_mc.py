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
    muon         muon-labelled tracks, every residual range (for --coords rr;
                 not in the defaults)
    proton       proton-labelled tracks
    all_tracks   every track
  from the pandtrk table (make_pandora_hits.py: Pandora tracks > 30 cm,
  per plane; every event, the --stage cut does not apply): Pandora
  muon-like tracks (track score > 0.5, chi2_mu < chi2_p), per plane p0-p2:
    pandora_dqdx_p<n>         median hit dQ/dx (the CAF's calibrated value,
                              ADC/cm); with the MC also re-derived with the
                              data recombination (to electrons with the MC
                              plane gain of chi2pid, the track's median hit
                              angle to the field, 0.5 kV/cm)
    pandora_width_p<n>        median hit width [ticks]
    pandora_mult_p<n>         mean hit multiplicity
    pandora_int_over_sum_p<n> median hit integral / ADC sum
    pandora_integral_p<n>     median hit integral (all hits)
    t_wire (Pandora samples): how long the track's signal lasts on one wire,
    pitch x |cos(angle to the field)| / 1.563 mm/us; an AC-coupling (RC)
    response mismatch between data and simulation acts on long signals, so a
    ratio falling with t_wire the same way on every plane points to it
    pandora_pitch_p<n>        median hit pitch (all hits):
                              dQ/dx = integral / pitch (x calibration), so
                              these split a dQ/dx difference into charge and
                              length
  profiled in angle to drift and in angle_wire, the angle of the track's
  y-z projection to the plane's wires (WIRE_ANGLE_Z: induction at +-60 deg
  to the vertical, mirrored between the TPCs; collection vertical), plus
  dqdx_<sample>_drift_wire_map.png: the data / prediction ratio of the
  per-track median in cells of both angles (the MC with the data
  recombination for dQ/dx)
  all over hits with residual range > 30 cm; Pandora's hits come from the
  same signal-processing output as Wire-Cell's imaging
  from the pandwire, pandtrk and hdr tables, through the cafpyana dE/dx
  calibration chain (pandora_calib.py, no CVMFS needed: YZ maps, electron lifetime on data and
  MC, gains, then recombination; the CAF's own dQ/dx leaves the MC without its
  35 ms lifetime correction), the longest Pandora muon-like track of each
  event with a --stage candidate, per plane, median over its wires at
  residual range > 30 cm; the Pandora counterparts of mu1_angle (not in the
  defaults: --samples pandcal_dedx_p0 pandcal_dedx_p1 pandcal_dedx_p2):
    pandcal_dedx_p<n>  dE/dx [MeV/cm], data with the SBND data recombination,
                       MC with ICARUS (the cafpyana standard); the dashed
                       prediction converts the MC with the data recombination
    pandcal_dqdx_p<n>  dQ/dx [ke/cm] (before recombination); the dashed
                       prediction is the MC with the data recombination
  from the candidate table (any dataframes):
    trk1_vtx     longest track median dQ/dx at the reco vertex (no direction)
    mu1_angle    longest muon-labelled track median dQ/dx, with the direction
                 of the cosmic tagger's main track (cosmict_2_theta / _phi:
                 polar angle to z, azimuth from x; the same track as the
                 trkseg direction to |cos| > 0.99 for 3/4 of the candidates),
                 position at the reco vertex

Residual range (--coords rr, trkseg samples): also
dqdx_<sample>_rr_angle_slices.png, the profile in slices of angle to drift
(--angle-slices), ratios to the nominal MC and to the MC with the data
recombination per slice: a ratio changing with residual range points to the
dE/dx dependence of recombination, one changing between slices to its angle
dependence (or another angle effect), a common flat offset to a scale.

Pandora grids: dqdx_pandora_grid_<coordinate>.png (angle_drift, angle_wire,
t_wire): the data / prediction ratio for every Pandora quantity (rows: dQ/dx,
width, multiplicity, integral / ADC sum, integral, pitch) and plane
(columns), rows on a shared y range: a change of the charge before readout
(recombination, lifetime, a common gain) moves the planes together, signal
processing / field response / electronics do not.

profiles.csv: every plotted bin (sample, coordinate, angle slice for the
residual-range and wire-angle cells, entries, data / prediction medians and
errors, and the MC with the data recombination): aggregates only, no
events, small enough to share.

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
from analysis_village.SinglePhotonWC import cuts, dfio, pandora_calib, wc_calo  # noqa: E402
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
TRAIN_KEY = "pandtrain"
TRAIN_LEVELS = PAND_LEVELS + ["wire", "train"]
WIRE_KEY = "pandwire"
WIRE_LEVELS = PAND_LEVELS + ["wire"]
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


def _wc_swap(mc):
    """Wire-Cell dQ/dx (MIPs) of the MC with the data recombination, at the track's angle to drift."""
    return wc_calo.swap_recombination(mc.dqdx * DQDX_MIP, np.radians(mc.angle_drift)) / DQDX_MIP


def _pandora_swap(plane):
    """Pandora MC hit dQ/dx (ADC/cm) with the data recombination: to electrons with the MC gain of the plane,
    at the track's median hit angle to the field (phi), nominal 0.5 kV/cm field."""
    scale = 1.0 / wc_calo.PANDORA_GAIN["mc"][plane]
    return lambda mc: wc_calo.swap_recombination(mc.dqdx, np.radians(mc.phi), scale=scale)


def _pandora(t, plane, q):
    """Pandora muon-like tracks (track score > 0.5, collection-plane chi2_mu
    < chi2_p) on one plane: quantity q (column of the pandtrk table)."""
    s = t[(t.index.get_level_values("plane") == plane) & (t.track_score > 0.5) & (t.chi2_mu < t.chi2_p)]
    # drift length one wire sees from the track: the hit pitch times |cos(angle to the field)|; / drift speed
    # = how long the signal on one wire lasts, what an AC-coupling (RC) mismatch would act on
    dx_wire = (s.pitch * np.abs(np.cos(np.radians(s.phi)))).to_numpy(float)
    return pd.DataFrame({"dqdx": s[q], "phi": s.phi, "x": s.x, "y": 0.5 * (s.start_y + s.end_y), "z": 0.5 * (s.start_z + s.end_z),
                         "angle_drift": _angle_drift(s.dir_x), "angle_yz": _angle_yz(s.dir_y, s.dir_z),
                         "angle_wire": _angle_wire(plane, s.x, s.dir_y, s.dir_z),
                         "t_wire": dx_wire / DRIFT_SPEED_CM_US}, index=s.index)


# name: table needed (None: the candidates only), function of (candidates, table) -> points (the plotted
# quantity in column "dqdx"), title, bins of its distribution, y label, Wire-Cell dQ/dx in MIPs (recombination
# swap and fit apply)
SAMPLES = {
    "mu_mip": dict(table=SEG_KEY, fn=lambda w, s: _segments(s, lambda d: (d.pdg == PID["muon"]) & (d.rr > RR_MIP_MIN)),
                   title="muon-labelled tracks, residual range > %g cm (4 cm pieces)" % RR_MIP_MIN,
                   bins=np.arange(0, 3.01, 0.05)),
    "muon": dict(table=SEG_KEY, fn=lambda w, s: _segments(s, lambda d: d.pdg == PID["muon"]),
                 title="muon-labelled tracks, all residual ranges (4 cm pieces)", bins=np.arange(0, 4.01, 0.05)),
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
                            ("int_over_sum_mip", "hit integral / ADC sum", np.arange(0.8, 1.201, 0.01)),
                            ("integral", "hit integral [ADC tick]", np.arange(0, 2001, 25)),
                            ("pitch", "hit pitch [cm]", np.arange(0, 3.01, 0.05))]:
        SAMPLES["pandora_%s_p%d" % (_q.replace("_mip", ""), _pl)] = dict(
            table=PAND_KEY, fn=(lambda pl, q: lambda w, t: _pandora(t, pl, q))(_pl, _q),
            title="Pandora muon-like tracks > 30 cm, plane %d: per-track %s %s (%s)"
                  % (_pl, "mean" if _q == "mult_mip" else "median", _lab, "hits with rr > 30 cm" if _q.endswith("_mip") else "all hits"),
            bins=_bins, ylabel=("mean " if _q == "mult_mip" else "median ") + _lab, wc=False,
            swap=(lambda pl: _pandora_swap(pl))(_pl) if _q == "dqdx_mip" else None)
DEFAULT_SAMPLES = ["mu_mip", "proton", "all_tracks", "mu1_angle"] + [   # trk1_vtx has no angle: --samples trk1_vtx --coords x
    "pandora_%s_p%d" % (q, pl) for q in ("dqdx", "width", "mult", "int_over_sum", "integral", "pitch") for pl in range(3)]

COORDS = [   # column, bins, label
    ("x", np.linspace(-200, 200, 21), "x [cm] (drift; cathode at 0, anodes at ±200)"),
    ("y", np.linspace(-200, 200, 21), "y [cm]"),
    ("z", np.linspace(0, 500, 26), "z [cm]"),
    ("angle_drift", None, "angle to the drift (x) direction [deg]"),   # angle bins: --angle-bin
    ("angle_yz", None, "angle in the y-z (wire) plane to the beam axis [deg]"),
    ("angle_wire", None, "angle of the track's y-z projection to the plane's wires [deg]"),
    ("rr", np.r_[np.arange(0, 30, 2), np.arange(30, 151, 10)], "residual range [cm] (4 cm pieces: mean of the piece)"),
    ("t_wire", np.r_[0.01, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64],
     "signal duration per wire, pitch |cos φ| / v_drift [µs] (log axis)"),
]
ANGLE_RANGE = {"angle_drift": (0, 90), "angle_yz": (-90, 90), "angle_wire": (0, 90)}
# wire direction in the y-z plane, angle to z [deg], per (plane, TPC; 0: x < 0, 1: x > 0); fitted to the
# Pandora hit wire numbers vs (y, z) on local MC CAFs (0.30 cm pitch): induction wires at +-60 deg to the
# vertical, mirrored between the TPCs, collection wires vertical
WIRE_ANGLE_Z = {(0, 0): 30.0, (0, 1): -30.0, (1, 0): -30.0, (1, 1): 30.0, (2, 0): 90.0, (2, 1): 90.0}
DRIFT_SPEED_CM_US = 0.1563   # Wire-Cell simulation / reconstruction (simparams.jsonnet)
MAP_BINS = np.linspace(0, 90, 10)
LOG_COORDS = {"t_wire"}   # most tracks sit at 1-3 us, the tail reaches 60 us: unreadable on a linear axis   # both angles of the Pandora drift x wire-angle maps
MAP_MIN_ENTRIES = 30


def _angle_wire(plane, x, dy, dz):
    """Angle [deg] between the track's y-z projection and the wires of the plane (in the track's TPC),
    folded to 0-90: 0 along the wires, 90 across them."""
    a = np.radians(np.where(np.asarray(x) < 0, WIRE_ANGLE_Z[(plane, 0)], WIRE_ANGLE_Z[(plane, 1)]))
    dy, dz = np.asarray(dy, float), np.asarray(dz, float)
    c = np.abs(dy * np.sin(a) + dz * np.cos(a)) / np.maximum(np.hypot(dy, dz), 1e-9)
    return np.degrees(np.arccos(np.clip(c, 0, 1)))


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


PROFILE_ROWS = []   # binned medians of every plot, written to profiles.csv (aggregates only, shareable)


def _record(sample, coord, lo, hi, on, pred, alt=None, extra=None):
    """Append one row per bin: entries and median +- error of data, prediction (and the MC with the data
    recombination)."""
    nd = np.histogram(on[coord], np.r_[lo, hi[-1]])[0] if len(on) else np.zeros(len(lo))
    npred = np.histogram(pred[coord], np.r_[lo, hi[-1]])[0] if len(pred) else np.zeros(len(lo))
    Pd, Pp = profile(on, coord, np.r_[lo, hi[-1]]), profile(pred, coord, np.r_[lo, hi[-1]])
    Pa = profile(alt, coord, np.r_[lo, hi[-1]]) if alt is not None else np.full_like(Pd, np.nan)
    for i in range(len(lo)):
        PROFILE_ROWS.append(dict(sample=sample, coord=coord, **(extra or {}), lo=lo[i], hi=hi[i], n_data=nd[i],
                                 data=Pd[i, 0], data_err=Pd[i, 1], n_pred=npred[i], pred=Pp[i, 0], pred_err=Pp[i, 1],
                                 pred_datarecomb=Pa[i, 0], pred_datarecomb_err=Pa[i, 1]))


def plot_profile(name, title, coord, bins, label, on, mc, off, out, mc_alt=None, ylabel="median dQ/dx [MIP]"):
    pred = pd.concat([mc, off])
    _record(name, coord, bins[:-1], bins[1:], on, pred, pd.concat([mc_alt, off]) if mc_alt is not None else None)
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
    if coord in LOG_COORDS:
        r.set_xscale("log")
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


def plot_wire_map(name, title, on, mc, off, out, mc_alt=None, qlabel="median"):
    """Data / prediction of the per-cell median in angle to drift x angle to the plane's wires (cell values
    printed; blank: fewer than MAP_MIN_ENTRIES data or MC entries). Against the MC with the data
    recombination when given, else the nominal MC."""
    def cells(p):
        ix = np.digitize(p.angle_drift.to_numpy(float), MAP_BINS) - 1
        iy = np.digitize(p.angle_wire.to_numpy(float), MAP_BINS) - 1
        n = len(MAP_BINS) - 1
        ok = (ix >= 0) & (ix < n) & (iy >= 0) & (iy < n) & np.isfinite(p.dqdx.to_numpy(float))
        b = (ix * n + iy)[ok]
        P = _fast_profile(p.dqdx.to_numpy(float)[ok], p.w.to_numpy(float)[ok], b, n * n)
        return P[:, 0].reshape(n, n), P[:, 1].reshape(n, n), np.bincount(b, minlength=n * n).reshape(n, n)
    m_ = mc if mc_alt is None else mc_alt
    for lo, hi in zip(MAP_BINS[:-1], MAP_BINS[1:]):   # the cells, as rows: angle_wire bins per angle_drift slice
        sl = lambda d: d[(d.angle_drift >= lo) & (d.angle_drift < hi)]
        _record(name, "angle_wire", MAP_BINS[:-1], MAP_BINS[1:], sl(on), pd.concat([sl(mc), sl(off)]),
                pd.concat([sl(mc_alt), sl(off)]) if mc_alt is not None else None,
                extra=dict(angle_drift_lo=lo, angle_drift_hi=hi))
    Md, Ed, Nd = cells(on)
    Mp, Ep, _ = cells(pd.concat([m_, off]))
    Nm = cells(m_)[2]
    R = np.where((Nd >= MAP_MIN_ENTRIES) & (Nm >= MAP_MIN_ENTRIES), Md / Mp, np.nan)
    Re = R * np.hypot(Ed / Md, Ep / Mp)
    fig, a = plt.subplots(figsize=(10, 7.5))
    lim = max(0.05, min(0.3, np.nanmax(np.abs(R - 1)) if np.isfinite(R).any() else 0.1))
    im = a.pcolormesh(MAP_BINS, MAP_BINS, R.T, cmap="RdBu_r", vmin=1 - lim, vmax=1 + lim, shading="flat")
    c = 0.5 * (MAP_BINS[1:] + MAP_BINS[:-1])
    for i in range(len(c)):
        for j in range(len(c)):
            if np.isfinite(R[i, j]):
                a.text(c[i], c[j], "%.3f\n±%.3f" % (R[i, j], Re[i, j]), ha="center", va="center", fontsize=7.5)
    fig.colorbar(im, ax=a, label="data / prediction, %s" % qlabel)
    a.set_xlabel("angle to the drift (x) direction [deg]")
    a.set_ylabel("angle to the plane's wires (y-z projection) [deg]")
    a.grid(False)
    a.set_title("%s\nprediction: %s + off-beam; blank: < %d data or MC entries"
                % (title, "MC with the data recombination" if mc_alt is not None else "nominal MC", MAP_MIN_ENTRIES),
                loc="left", fontsize=9.5)
    fig.savefig(os.path.join(out, "dqdx_%s_drift_wire_map.png" % name), dpi=130, bbox_inches="tight")
    plt.close(fig)


SLICE_COLORS = ["#2a78d6", "#1baf7a", "#eb6834", "#7f5bd5", "#a05a2c"]


def plot_rr_angle_slices(name, title, bins, label, on, mc, off, out, slices, mc_alt=None):
    """Median dQ/dx vs residual range in slices of angle to drift: data (points) and prediction (lines) per
    slice; ratio panel data / nominal prediction (filled) and / the MC with the data recombination (open)."""
    c = 0.5 * (bins[1:] + bins[:-1])
    xe = np.diff(bins) / 2
    fig, (a, r) = plt.subplots(2, 1, figsize=(10, 7.5), sharex=True, gridspec_kw=dict(height_ratios=[3, 1.4], hspace=0.05))
    lims = [0.05]
    sel = lambda d, lo, hi: d[(d.angle_drift >= lo) & (d.angle_drift < hi)]
    for k, (lo, hi) in enumerate(zip(slices[:-1], slices[1:])):
        col = SLICE_COLORS[k % len(SLICE_COLORS)]
        Pd = profile(sel(on, lo, hi), "rr", bins)
        Pp = profile(pd.concat([sel(mc, lo, hi), sel(off, lo, hi)]), "rr", bins)
        _record(name, "rr", bins[:-1], bins[1:], sel(on, lo, hi), pd.concat([sel(mc, lo, hi), sel(off, lo, hi)]),
                pd.concat([sel(mc_alt, lo, hi), sel(off, lo, hi)]) if mc_alt is not None else None,
                extra=dict(angle_drift_lo=lo, angle_drift_hi=hi))
        a.stairs(Pp[:, 0], bins, baseline=None, color=col, lw=1.8, label="%g-%g°: prediction" % (lo, hi))
        a.errorbar(c, Pd[:, 0], Pd[:, 1], fmt="o", color=col, mec="k", ms=4, label="%g-%g°: data" % (lo, hi))
        sh = (k - (len(slices) - 2) / 2) * 0.25 * xe
        ratio = Pd[:, 0] / Pp[:, 0]
        rerr = ratio * np.hypot(Pd[:, 1] / Pd[:, 0], Pp[:, 1] / Pp[:, 0])
        r.errorbar(c + sh, ratio, rerr, fmt="o", color=col, ms=4)
        if mc_alt is not None:
            Pa = profile(pd.concat([sel(mc_alt, lo, hi), sel(off, lo, hi)]), "rr", bins)
            ra = Pd[:, 0] / Pa[:, 0]
            r.errorbar(c + sh, ra, ra * np.hypot(Pd[:, 1] / Pd[:, 0], Pa[:, 1] / Pa[:, 0]), fmt="o", mfc="white",
                       color=col, ms=4)
            if np.isfinite(ra).any():
                lims.append(1.2 * np.nanmax(np.abs(ra - 1)))
        if np.isfinite(ratio).any():
            lims.append(1.2 * np.nanmax(np.abs(ratio - 1) + np.nan_to_num(rerr)))
    r.axhline(1, color="k", lw=1)
    lim = min(max(lims), 0.3)
    r.set_ylim(1 - lim, 1 + lim)
    r.set_ylabel("data / pred.")
    r.text(0.99, 0.04, "filled: vs nominal MC" + ("; open: vs MC with data recombination" if mc_alt is not None else ""),
           transform=r.transAxes, ha="right", va="bottom", fontsize=8.5)
    r.set_xlabel(label)
    r.set_xlim(bins[0], bins[-1])
    a.set_ylabel("median dQ/dx [MIP]")
    a.set_title(title + "\nslices of angle to the drift direction", loc="left", fontsize=10.5)
    a.legend(frameon=False, fontsize=8.5, ncol=2)
    fig.savefig(os.path.join(out, "dqdx_%s_rr_angle_slices.png" % name), dpi=130, bbox_inches="tight")
    plt.close(fig)


PANDORA_QUANTITIES = [("dqdx", "dQ/dx"), ("width", "hit width"), ("mult", "hit multiplicity"),
                      ("int_over_sum", "integral / ADC sum"), ("integral", "hit integral"), ("pitch", "hit pitch")]


def plot_pandora_grid(coord, label, bins, results, out):
    """One figure: data / prediction vs coord for every Pandora quantity (rows) and plane (columns). For dQ/dx,
    open markers compare with the MC with the data recombination. Rows share a y range (set by the
    well-measured bins), so planes compare at a glance. A change of the charge before readout moves the three
    columns together; signal processing, field response and electronics act per plane."""
    rows = [(q, lab) for q, lab in PANDORA_QUANTITIES if any("pandora_%s_p%d" % (q, pl) in results for pl in range(3))]
    if not rows:
        return
    c = 0.5 * (bins[1:] + bins[:-1])
    fig, axs = plt.subplots(len(rows), 3, figsize=(16, 2.7 * len(rows) + 1), sharex=True, squeeze=False)
    for i, (q, lab) in enumerate(rows):
        lim = 0.02
        for pl in range(3):
            a = axs[i, pl]
            key = "pandora_%s_p%d" % (q, pl)
            if key in results:
                on, mc, off, alt = results[key]
                Pd = profile(on, coord, bins)
                for m_, mfc, name in [(mc, None, "vs nominal MC"), (alt, "white", "vs MC with data recombination")]:
                    if m_ is None:
                        continue
                    Pp = profile(pd.concat([m_, off]), coord, bins)
                    ratio = Pd[:, 0] / Pp[:, 0]
                    rerr = ratio * np.hypot(Pd[:, 1] / Pd[:, 0], Pp[:, 1] / Pp[:, 0])
                    a.errorbar(c, ratio, rerr, fmt="o", color="#1f1f1e", mfc=mfc, ms=4, label=name)
                    good = np.isfinite(ratio) & (rerr < 0.03)
                    if good.any():
                        lim = max(lim, 1.25 * np.max(np.abs(ratio[good] - 1) + rerr[good]))
            a.axhline(1, color="#1f1f1e", lw=1)
            if i == 0:
                a.set_title("plane %d%s" % (pl, " (collection)" if pl == 2 else " (induction)"), fontsize=11)
            if pl == 0:
                a.set_ylabel("%s\ndata / pred." % lab, fontsize=10)
            if i == len(rows) - 1:
                a.set_xlabel(label, fontsize=10)
            if i == 0 and pl == 2 and q == "dqdx":
                a.legend(frameon=False, fontsize=8.5)
        lim = min(lim, 0.4)
        for pl in range(3):
            axs[i, pl].set_ylim(1 - lim, 1 + lim)
            axs[i, pl].set_xlim(bins[0], bins[-1])
            if coord in LOG_COORDS:
                axs[i, pl].set_xscale("log")
    fig.suptitle("Pandora muon-like tracks: data / prediction by quantity (rows) and plane (columns) vs %s" % coord,
                 x=0.01, ha="left", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(os.path.join(out, "dqdx_pandora_grid_%s.png" % coord), dpi=120, bbox_inches="tight")
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
    """The trkseg (rows of the candidates in keep, an index over CAND_LEVELS),
    pandtrk, pandtrain or pandwire table; None if the inputs have none."""
    files = dfio.expand_inputs(inputs)
    if key not in dfio.file_keys(files[0]):
        return None

    names = {SEG_KEY: SEG_LEVELS, PAND_KEY: PAND_LEVELS, TRAIN_KEY: TRAIN_LEVELS, WIRE_KEY: WIRE_LEVELS}[key]
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


CAL_KEY = "pandcal"   # pandwire + pandtrk (+ hdr for data) through pandora_calib, per (track, plane)


def load_pandcal(inputs, cands, isdata, cal, nproc=None):
    """Per (track, plane): the longest Pandora muon-like track (track score > 0.5, chi2_mu < chi2_p) of each
    event with a candidate in cands, medians over its wires at residual range > RR_MIP_MIN of the calibrated
    dQ/dx and dE/dx (pandora_calib.calibrate); None without the pandwire table."""
    wires = load_table(inputs, WIRE_KEY, nproc=nproc)
    trk = load_table(inputs, PAND_KEY, nproc=nproc)
    if wires is None or trk is None:
        return None
    events = cands.index.droplevel([l for l in cands.index.names if l not in ("file", "__ntuple", "entry")]).unique()
    evlev = ["file", "__ntuple", "entry"]
    wires = wires[(wires.track_score > 0.5) & (wires.chi2_mu < wires.chi2_p) & (wires.rr_min > RR_MIP_MIN) & (wires.pitch > 0)
                  & wires.index.droplevel([l for l in wires.index.names if l not in evlev]).isin(events)]
    key = pandora_calib.TRK_LEVELS + ["plane"]
    cols = ["dqdx_e", "dedx", "dedx_mcrecomb", "dedx_datarecomb", "x", "phi", "len", "dir_x", "dir_y", "dir_z"]
    if not len(wires):
        return pd.DataFrame(columns=cols, index=pd.MultiIndex.from_arrays([[]] * len(key), names=key), dtype=float)
    runs = None
    if isdata:
        frames, _ = dfio.load(dfio.expand_inputs(inputs), keys=["hdr"], nproc=nproc)
        runs = frames["hdr"].run
        runs = runs[~runs.index.duplicated()]
    w = pandora_calib.calibrate(wires, trk, cal, isdata, runs)
    g = w.groupby(key)
    t = pd.DataFrame({c: g[c].median() for c in ("dqdx_e", "dedx", "dedx_mcrecomb", "dedx_datarecomb", "x", "phi")})
    t = t.join(g[["len", "dir_x", "dir_y", "dir_z"]].first())
    longest = t.groupby(evlev + ["plane"]).len.transform("max")
    return t[t.len == longest]


def _pandcal(t, plane, q):
    s = t[t.index.get_level_values("plane") == plane]
    v = s.dqdx_e / 1e3 if q == "dqdx" else s.dedx
    return pd.DataFrame({"dqdx": v, "phi": s.phi, "x": s.x, "angle_drift": _angle_drift(s.dir_x),
                         "angle_yz": _angle_yz(s.dir_y, s.dir_z), "_dedx_datarecomb": s.dedx_datarecomb,
                         "_dqdx_e": s.dqdx_e}, index=s.index)


def _pandcal_swap(q):
    """The MC with the data recombination: dE/dx converted with it, or dQ/dx re-derived with it (e/cm)."""
    if q == "dedx":
        return lambda mc: mc._dedx_datarecomb
    return lambda mc: wc_calo.swap_recombination(mc._dqdx_e, np.radians(mc.phi), scale=1.0) / 1e3


for _pl in range(3):
    for _q, _lab, _bins in [("dedx", "dE/dx [MeV/cm]", np.arange(0, 6.01, 0.05)), ("dqdx", "dQ/dx [ke/cm]", np.arange(0, 120.1, 1))]:
        SAMPLES["pandcal_%s_p%d" % (_q, _pl)] = dict(
            table=CAL_KEY, fn=(lambda pl, q: lambda w, t: _pandcal(t, pl, q))(_pl, _q),
            title="longest Pandora muon-like track, plane %d: median calibrated %s (wires at rr > %g cm)" % (_pl, _lab, RR_MIP_MIN),
            bins=_bins, ylabel="median " + _lab, wc=False, swap=_pandcal_swap(_q))


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
    ap.add_argument("--coords", nargs="+", default=["angle_drift", "angle_wire", "t_wire"], choices=[c[0] for c in COORDS],
                    help="profile coordinates (default: angle_drift, and angle_wire and t_wire for the Pandora samples)")
    ap.add_argument("--dist", action="store_true", help="also the dQ/dx distribution plots")
    ap.add_argument("--angle-slices", type=float, nargs="+", default=[0, 40, 60, 90],
                    help="angle-to-drift slice edges of the residual-range profiles [deg] (default 0 40 60 90)")
    ap.add_argument("--by-tpc", action="store_true", help="also the profiles split by TPC")
    ap.add_argument("--no-recomb-swap", action="store_true", help="no MC prediction with the data recombination")
    ap.add_argument("--fit-recomb", action="store_true", help="fit the ellipsoid ratio R and an MC scale to the angle profile")
    ap.add_argument("--sbnd-data", help="sbnd_data directory for the pandcal samples (default: the copies in pandora_calib)")
    ap.add_argument("--lifetime-db", default="v01_41_00", choices=list(pandora_calib.ETAU_TABLES),
                    help="embedded sbnd_data lifetime database version (default v01_41_00, as chi2pid)")
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
    cal = pandora_calib.Calibration(args.sbnd_data, args.lifetime_db) if any(SAMPLES[n]["table"] == CAL_KEY for n in args.samples) else None
    for key in {SAMPLES[n]["table"] for n in args.samples} - {None}:
        if key == CAL_KEY:
            t = {k: load_pandcal(inputs[k], cands[k], k != "mc", cal, args.nproc) for k in cands}
        else:
            t = {k: load_table(inputs[k], key, cands[k].index if key == SEG_KEY else None, args.nproc) for k in cands}
        if any(v is None for v in t.values()):
            print("no %s table in some inputs (older dataframes): its samples are skipped" % key)
        else:
            tables[key] = t
    comp = cands["mc"].cat.map(cuts.CATEGORY_LABELS).value_counts(normalize=True)
    lines = [norm_line, "stage: %s (%d MC candidates; MC truth: %s)" % (
        stage_title, len(cands["mc"]), ", ".join("%s %.0f%%" % (k, 100 * v) for k, v in comp.head(5).items())), ""]
    results = {}
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
        swap = spec.get("swap", _wc_swap if wcq else None)
        alt = None
        if swap is not None and "angle_drift" in mc and not args.no_recomb_swap:
            alt = mc.assign(dqdx=swap(mc))
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
        results[name] = (on, mc, off, alt)
        for coord, cbins, label in COORDS:
            if coord in on and coord in args.coords:
                cbins = angle_bins(coord, args.angle_bin) if cbins is None else cbins
                plot_profile(name, full, coord, cbins, label, on, mc, off, args.outdir, alt,
                             ylabel=spec.get("ylabel", "median dQ/dx [MIP]"))
                if args.by_tpc and coord != "x":
                    plot_profile_tpc(name, full, coord, cbins, label, on, mc, off, args.outdir)
                if coord == "rr" and "angle_drift" in on:
                    plot_rr_angle_slices(name, full, cbins, label, on, mc, off, args.outdir, args.angle_slices, alt)
        if "angle_wire" in on:
            plot_wire_map(name, full, on, mc, off, args.outdir, alt, spec.get("ylabel", "median dQ/dx"))
        if "angle_drift" in on and wcq:
            if args.fit_recomb:
                fb = angle_bins("angle_drift", args.angle_bin)
                fit = fit_recombination(on, mc, off, fb)
                plot_recomb_fit(name, full, fb, on, mc, off, args.outdir, fit)
                lines.append("%-12s recombination fit: R = %.2f (Delta chi2 = 1: %.2f-%.2f), MC scale k = %.3f, chi2/ndf = %.1f/%d"
                             % ("", fit["R"], *fit["R_range"], fit["k"], fit["chi2_min"], fit["ndf"]))
    for coord, cbins, label in COORDS:   # Pandora quantity x plane grids, for every Pandora coordinate requested
        if coord in args.coords and any(k.startswith("pandora_") and coord in v[0] for k, v in results.items()):
            plot_pandora_grid(coord, label, angle_bins(coord, args.angle_bin) if cbins is None else cbins, results,
                              args.outdir)
    pd.DataFrame(PROFILE_ROWS).to_csv(os.path.join(args.outdir, "profiles.csv"), index=False, float_format="%.6g")
    open(os.path.join(args.outdir, "summary.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("wrote the plots to %s" % args.outdir)


if __name__ == "__main__":
    main()
