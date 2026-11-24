#!/usr/bin/env python3
"""Data / MC of the charge one Pandora track leaves on one wire (pandwire
table, make_pandora_hits.make_pandora_wire_df), to tell whether the dQ/dx
discrepancy of tracks near the drift direction is in dQ or in dx.

SBND's calorimetry (Gnocchi, ChargeMethod 3) makes one calo point per ROI
snippet a track crosses on a wire, dQ/dx = the snippet's summed track-hit
integral / the full per-wire pitch. A long signal the signal processing splits
into several ROIs gives several points on one wire, each with the full pitch,
so the per-point dQ/dx falls with the signal duration (in MC too) and measures
the ROI splitting. Summed per wire (muon-like tracks: track score > 0.5,
chi2_mu < chi2_p; every point at residual range > --rr-min):

  wire_dqdx    sum of the points' dQ/dx x pitch / pitch [ADC/cm, calibrated]:
               the wire's track charge over the track length the wire sees
  roi_dqdx     sum of the snippets' ADC sums / pitch [ADC/cm, uncalibrated]:
               all the charge in the ROIs (independent of the hit fits)
  amp_dqdx     the ROIs' ADC per tick / (v_drift x tick x |cos phi|) [ADC/cm]:
               dQ/dx from the amplitude alone, independent of how much of the
               signal's length the ROIs keep, and of the pitch; only for
               t_wire >= 8 us (shorter ROIs are mostly padding)
  coverage     summed ROI length / expected duration (pitch |cos phi| / v_drift)
  n_points     ROI snippets on the wire (mean): the splitting
  point_dqdx   wire_dqdx / n_points: what the CAF's per-point dQ/dx shows

Reading: a wire_dqdx / roi_dqdx deficit that amp_dqdx shares is in the
amplitude (dQ: gain, recombination, the deconvolved amplitude); one that
amp_dqdx does not share, with coverage low, is signal length lost at the ROI
edges; a point_dqdx deficit with wire_dqdx fine is more ROI splitting in data
(bookkeeping, not charge). The pitch itself is geometric, from the track.

Coordinates: t_wire (the expected duration, log axis) and angle_drift (from
the track direction). Every event, no candidate selection. Wires with t_wire
< 4 us are a 1-in-10 sample in the table (weighted back).

Outputs (-o): wirecharge_<quantity>_<coord>.png (planes in columns, ratio
below), wirecharge_dq_vs_dx_<coord>.png (data / prediction of wire_dqdx,
amp_dqdx, coverage and n_points together per plane: the dQ-or-dx plot),
wirecharge_profiles.csv (binned medians only, shareable) and
wirecharge_summary.txt (per plane, a fit of data / pred. = a (1 - b t_wire)
and the RC time constant a flat-pulse loss of b would need: loss kappa T / 2,
kappa = 1/tau_true - 1/tau_assumed). The MC with the data recombination
(dashed) multiplies each MC wire by the MIP-charge ratio of the two
recombination models at its angle to the field.

  python analysis_village/SinglePhotonWC/plot_wirecharge_data_mc.py \\
      --mc <MC> --onbeam <on-beam> --offbeam <off-beam> -o wirecharge_plots/ [--area-norm]
"""

import argparse
import csv
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from analysis_village.SinglePhotonWC import wc_calo  # noqa: E402
from analysis_village.SinglePhotonWC.cutflow import load as load_mc  # noqa: E402
from analysis_village.SinglePhotonWC.plot_data_mc import load_data, normalisation  # noqa: E402
from analysis_village.SinglePhotonWC.plot_dqdx_data_mc import (  # noqa: E402
    DRIFT_SPEED_CM_US, WIRE_KEY, _angle_drift, load_table, median_err)

plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), "dune.mplstyle"))

TICK_US = 0.5
TAU_ASSUMED_US = 1100.0      # RC time constant of the SBND Wire-Cell configs
AMP_MIN_US = 8.0           # amp_dqdx only above this t_wire: shorter ROIs are mostly padding and response tails
MIP_ADC_CM = 1100.0          # Pandora MIP hit dQ/dx [ADC/cm] (short signals), only sets the recombination factor
COORDS = {
    "t_wire": (np.r_[0.1, 0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128],
               "expected signal duration per wire, pitch |cos φ| / v_drift [µs] (log axis)", True),
    "angle_drift": (np.r_[0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 90], "track angle to the drift (x) direction [deg]", False),
}
QUANTITIES = [   # column, label, statistic, ratio range
    ("wire_dqdx", "wire track charge / pitch [ADC/cm]", "median", (0.6, 1.4)),
    ("roi_dqdx", "wire ROI ADC sum / pitch [ADC/cm]", "median", (0.6, 1.4)),
    ("amp_dqdx", "ROI ADC per tick → dQ/dx [ADC/cm]", "median", (0.6, 1.4)),
    ("coverage", "ROI length / expected duration", "median", (0.6, 1.4)),
    ("n_points", "ROI snippets on the wire", "mean", (0.8, 1.2)),
    ("point_dqdx", "per-point dQ/dx (wire / snippets) [ADC/cm]", "median", (0.6, 1.4)),
]
CHARGE = ("wire_dqdx", "roi_dqdx", "amp_dqdx", "point_dqdx")   # get the data-recombination variant
INK, BLUE, GREEN, ORANGE, PURPLE = "#1f1f1e", "#2a78d6", "#1baf7a", "#eb6834", "#7f5bd5"


def recomb_factor(phi_deg):
    """MIP charge under the data recombination / under the MC one, at the angle to the field."""
    phi = np.radians(np.asarray(phi_deg, float))
    x = np.full(phi.shape, MIP_ADC_CM)
    return wc_calo.swap_recombination(x, phi, scale=1.0 / wc_calo.PANDORA_GAIN["mc"][2]) / x


def prepare(t, rr_min):
    s = t[(t.track_score > 0.5) & (t.chi2_mu < t.chi2_p) & (t.rr_min > rr_min) & (t.pitch > 0)].copy()
    cphi = np.abs(np.cos(np.radians(s.phi.to_numpy(float))))
    s["angle_drift"] = _angle_drift(s.dir_x.to_numpy(float))
    s["wire_dqdx"] = s.q / s.pitch
    s["roi_dqdx"] = s.sumadc / s.pitch
    s["point_dqdx"] = s.wire_dqdx / s.n_points
    with np.errstate(divide="ignore", invalid="ignore"):
        s["amp_dqdx"] = np.where((cphi > 0.02) & (s.t_wire >= AMP_MIN_US), s.sumadc / s.roi_ticks / (DRIFT_SPEED_CM_US * TICK_US * cphi), np.nan)
        s["coverage"] = np.where(s.t_wire > 0, s.roi_ticks * TICK_US / s.t_wire, np.nan)
    return s


def stat(v, w, kind):
    if kind == "median":
        return median_err(v, w)
    if len(v) < 5:
        return np.nan, np.nan
    m = np.average(v, weights=w)
    neff = w.sum() ** 2 / (w ** 2).sum()
    return m, np.sqrt(np.average((v - m) ** 2, weights=w) / neff)


def profile(d, w, coord, bins, col, kind, factor=None):
    x = d[coord].to_numpy(float)
    v = d[col].to_numpy(float) * (1.0 if factor is None else factor)
    b = np.digitize(x, bins) - 1
    out = np.full((len(bins) - 1, 3), np.nan)
    for i in range(len(bins) - 1):
        m = (b == i) & np.isfinite(v)
        out[i, :2] = stat(v[m], w[m], kind) if m.any() else (np.nan, np.nan)
        out[i, 2] = m.sum()
    return out


def ratio(a, b):
    r = a[:, 0] / b[:, 0]
    return r, np.abs(r) * np.hypot(a[:, 1] / a[:, 0], b[:, 1] / b[:, 0])


def linear_fit(x, r, e):
    """Weighted fit r = a (1 - b x): (a, b, err a, err b, chi2, ndf)."""
    ok = np.isfinite(x) & np.isfinite(r) & np.isfinite(e) & (e > 0)
    if ok.sum() < 3:
        return (np.nan,) * 6
    A = np.vstack([np.ones(ok.sum()), x[ok]]).T
    W = 1 / e[ok] ** 2
    cov = np.linalg.inv(A.T @ (A * W[:, None]))
    p0, p1 = cov @ (A.T @ (W * r[ok]))
    chi2 = float(np.sum(W * (r[ok] - p0 - p1 * x[ok]) ** 2))
    b = -p1 / p0
    eb = abs(b) * np.hypot(np.sqrt(cov[1, 1]) / p1, np.sqrt(cov[0, 0]) / p0) if p1 != 0 else np.nan
    return p0, b, np.sqrt(cov[0, 0]), eb, chi2, int(ok.sum()) - 2


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mc", nargs="+", required=True, help="MC .df files, directories or globs")
    ap.add_argument("--onbeam", nargs="+", required=True, help="on-beam data .df files, directories or globs")
    ap.add_argument("--offbeam", nargs="+", required=True, help="off-beam data .df files, directories or globs")
    ap.add_argument("-o", "--outdir", default="wirecharge_plots")
    ap.add_argument("--area-norm", action="store_true", help="without POT / gate information (as plot_data_mc.py)")
    ap.add_argument("--offbeam-frac", type=float, help="off-beam share of the prediction, with --area-norm")
    ap.add_argument("--rr-min", type=float, default=30.0, help="residual range of every hit on the wire > this [cm]")
    ap.add_argument("--fit-min", type=float, default=4.0, help="lowest t_wire in the linear fit [us] (default 4)")
    ap.add_argument("-j", "--nproc", type=int)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    mc, _, pot_mc = load_mc(args.mc, None, nproc=args.nproc)
    on, pot_on, _ = load_data(args.onbeam, args.nproc, uboone_bdt=False)
    off, _, gates_off = load_data(args.offbeam, args.nproc, uboone_bdt=False)
    sc_mc, sc_off, norm_line = normalisation(mc, on, off, pot_mc, pot_on, gates_off, args.area_norm, args.offbeam_frac)
    del mc, on, off
    tab = {k: load_table(v, WIRE_KEY, nproc=args.nproc) for k, v in
           {"mc": args.mc, "on": args.onbeam, "off": args.offbeam}.items()}
    if any(v is None for v in tab.values()):
        sys.exit("no %s table in some inputs: make the dataframes with configs/sp_mc.py / sp_data.py" % WIRE_KEY)
    tab = {k: prepare(v, args.rr_min) for k, v in tab.items()}

    rows = []
    lines = [norm_line, "wires: muon-like Pandora tracks (track score > 0.5, chi2_mu < chi2_p), every point at rr > %g cm; fit of data / pred. = a (1 - b t_wire) for t_wire >= %g us" % (args.rr_min, args.fit_min), ""]
    for coord, (bins, xlabel, logx) in COORDS.items():
        c = np.sqrt(bins[1:] * bins[:-1]) if logx else 0.5 * (bins[1:] + bins[:-1])
        prof = {}
        for col, ylab, kind, rlim in QUANTITIES:
            fig, axs = plt.subplots(2, 3, figsize=(17, 7.5), sharex=True, gridspec_kw=dict(height_ratios=[3, 1.4], hspace=0.05))
            for pl in range(3):
                d = {k: v[v.index.get_level_values("plane") == pl] for k, v in tab.items()}
                pred = pd.concat([d["mc"], d["off"]])
                wp = np.r_[sc_mc * d["mc"].prescale.to_numpy(float), sc_off * d["off"].prescale.to_numpy(float)]
                Pd = profile(d["on"], d["on"].prescale.to_numpy(float), coord, bins, col, kind)
                Pp = profile(pred, wp, coord, bins, col, kind)
                Pm = profile(d["mc"], d["mc"].prescale.to_numpy(float), coord, bins, col, kind)
                Pa = None
                if col in CHARGE:
                    f = np.r_[recomb_factor(d["mc"].phi), np.ones(len(d["off"]))]
                    Pa = profile(pred, wp, coord, bins, col, kind, factor=f)
                prof[(col, pl)] = (Pd, Pp, Pa)
                a, r = axs[0, pl], axs[1, pl]
                a.plot(c, Pp[:, 0], "-", color=INK, lw=2, drawstyle="default", label="prediction (MC + off-beam)")
                a.plot(c, Pm[:, 0], "-", color=BLUE, lw=1.2, label="MC only")
                if Pa is not None:
                    a.plot(c, Pa[:, 0], "--", color=GREEN, lw=1.5, label="prediction, MC with data recombination")
                a.errorbar(c, Pd[:, 0], Pd[:, 1], fmt="o", color=INK, ms=4, label="on-beam data")
                rr_, re_ = ratio(Pd, Pp)
                r.errorbar(c, rr_, re_, fmt="o", color=INK, ms=4)
                if Pa is not None:
                    ra, rae = ratio(Pd, Pa)
                    r.errorbar(c * 1.04 if logx else c + 1, ra, rae, fmt="s", color=GREEN, ms=3)
                r.axhline(1, color=INK, lw=1)
                r.set_ylim(*rlim)
                a.set_ylim(bottom=0)
                if logx:
                    a.set_xscale("log")
                a.set_title("plane %d: %d data / %d MC wires" % (pl, len(d["on"]), len(d["mc"])), fontsize=10.5, loc="left")
                r.set_xlabel(xlabel, fontsize=10)
                if pl == 0:
                    a.set_ylabel("%s %s" % (kind, ylab), fontsize=10)
                    r.set_ylabel("data / pred.")
                    a.legend(frameon=False, fontsize=8.5)
                for i in range(len(c)):
                    rows.append(dict(coord=coord, quantity=col, stat=kind, plane=pl, lo=bins[i], hi=bins[i + 1],
                                     n_data=int(Pd[i, 2]), data=Pd[i, 0], data_err=Pd[i, 1],
                                     n_pred=int(Pp[i, 2]), pred=Pp[i, 0], pred_err=Pp[i, 1],
                                     pred_datarecomb=np.nan if Pa is None else Pa[i, 0]))
            fig.suptitle("Charge of one Pandora track on one wire (its ROI snippets summed): %s %s, data vs prediction" % (kind, ylab),
                         x=0.01, ha="left", fontsize=11)
            fig.savefig(os.path.join(args.outdir, "wirecharge_%s_%s.png" % (col, coord)), dpi=120, bbox_inches="tight")
            plt.close(fig)

        # the dQ-or-dx plot: ratios of the three together
        fig, axs = plt.subplots(1, 3, figsize=(17, 5), sharey=True)
        for pl in range(3):
            a = axs[pl]
            for col, lab, colr, off_ in [("wire_dqdx", "wire charge / pitch (amplitude and length)", INK, 1.0),
                                         ("amp_dqdx", "ADC per tick (amplitude only)", ORANGE, 1.03),
                                         ("coverage", "ROI length / expected (length only)", PURPLE, 0.97),
                                         ("n_points", "ROI snippets per wire (splitting)", BLUE, 1.06)]:
                Pd, Pp, _ = prof[(col, pl)]
                rr_, re_ = ratio(Pd, Pp)
                a.errorbar(c * off_ if logx else c + (off_ - 1) * 100, rr_, re_, fmt="o-", color=colr, ms=4, lw=1,
                           label=lab)
            a.axhline(1, color="#8a8984", lw=0.8)
            a.set_ylim(0.6, 1.4)
            if logx:
                a.set_xscale("log")
            a.set_xlabel(xlabel, fontsize=10)
            a.set_title("plane %d" % pl, fontsize=11, loc="left")
            if pl == 0:
                a.set_ylabel("data / prediction (medians)")
                a.legend(frameon=False, fontsize=9)
        fig.suptitle("Is the deficit in dQ or in dx? per-wire charge, data / prediction", x=0.01, ha="left", fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        fig.savefig(os.path.join(args.outdir, "wirecharge_dq_vs_dx_%s.png" % coord), dpi=120, bbox_inches="tight")
        plt.close(fig)

        if coord == "t_wire":
            sel = c >= args.fit_min
            for pl in range(3):
                for col in ("wire_dqdx", "roi_dqdx", "amp_dqdx", "coverage", "n_points"):
                    Pd, Pp, _ = prof[(col, pl)]
                    rr_, re_ = ratio(Pd, Pp)
                    a0, b, ea, eb, chi2, ndf = linear_fit(np.where(sel, c, np.nan), rr_, re_)
                    line = "plane %d %-12s a = %.3f +- %.3f, b = %+.3f +- %.3f %%/us, chi2/ndf %.1f/%d" % (
                        pl, col, a0, ea, 100 * b, 100 * eb, chi2, ndf)
                    k = 2 * b + 1 / TAU_ASSUMED_US
                    if col in CHARGE and np.isfinite(b) and b > 2 * eb:
                        line += "; as a flat-pulse RC loss: kappa = %.2f /ms, tau_true = %s" % (
                            2000 * b, "%.2f ms" % (1 / k / 1000) if k > 0 else "none")
                    lines.append(line)
                lines.append("")

    with open(os.path.join(args.outdir, "wirecharge_profiles.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    open(os.path.join(args.outdir, "wirecharge_summary.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("wrote the plots, wirecharge_profiles.csv and wirecharge_summary.txt to %s" % args.outdir)


if __name__ == "__main__":
    main()
