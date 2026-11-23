#!/usr/bin/env python3
"""Data / MC of Pandora hit trains (pandtrain table, make_pandora_hits.py):
the consecutive hits one track leaves on one wire, i.e. how the charge of a
single wire's long signal evolves with time. A residual AC-coupling (RC
high-pass) effect after the reconstruction's correction makes a flat signal
x(t) come out as x(t) (1 - kappa t) for kappa t << 1, with
kappa = 1/tau_true - 1/tau_assumed: the charge per tick of the later hits of a
train falls (kappa > 0: the correction too weak) or rises (kappa < 0). The
simulation closes by construction, so the data - MC difference of the
fractional slope is -kappa, the same for every train span. Gain, recombination
and the track direction cancel in a fractional slope within one train.

Trains: muon-like Pandora tracks (track score > 0.5, chi2_mu < chi2_p), all
hits of the train at residual range > --rr-min (MIP region), >= --min-hits
hits (default 2; most trains are 2 hits). The track-score cut matters: in
local MC most trains with track score < 0.5 come from showers (true photons
and electrons). Every event (no candidate selection).

Plots (planes in columns):
  hittrain_slope_vs_span.png   median fractional slope (per us) vs train span,
                               data, prediction (MC + off-beam) and MC only;
                               below, data - prediction with the fitted -kappa
  hittrain_lastfirst_vs_span.png
                               median charge per tick of the last / first hit
                               (model-free), data / prediction
  hittrain_shapes.png          train span and hits per train, data vs
                               prediction (unit area): Pandora splits a long
                               signal into hits, and a different split in data
                               would bias the slopes
hittrain_summary.txt: per plane, kappa (the weighted mean of -(data - pred.)
over the span bins) and the tau_true it implies for tau_assumed = 1.1 ms,
and the lifetime slope for comparison. Errors of medians: 1.253 x half the
16-84% width / sqrt(effective entries); trains of one track are correlated.

  python analysis_village/SinglePhotonWC/plot_hittrain_data_mc.py \\
      --mc <MC> --onbeam <on-beam> --offbeam <off-beam> -o hittrain_plots/ [--area-norm]
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2 as chi2_dist

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from analysis_village.SinglePhotonWC.cutflow import load as load_mc  # noqa: E402
from analysis_village.SinglePhotonWC.plot_data_mc import load_data, normalisation  # noqa: E402
from analysis_village.SinglePhotonWC.plot_dqdx_data_mc import TRAIN_KEY, load_table, median_err  # noqa: E402

plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), "dune.mplstyle"))

SPAN_BINS = np.array([5, 10, 15, 20, 30, 40, 60, 80, 120, 200, 400])   # us
TAU_ASSUMED_US = 1100.0       # the RC time constant the SBND Wire-Cell configs assume (rc_resp, chndb rcrc)
LIFETIME_US = 35000.0         # the simulated electron lifetime, for the size of its slope
INK, BLUE, GREY = "#1f1f1e", "#2a78d6", "#a3a29d"


def interpret(kappa, kerr, p_const):
    """What a fitted kappa says under the single-RC-stage model (kappa in 1/us)."""
    notes = []
    if np.isfinite(p_const) and p_const < 0.01:
        notes.append("the data - MC difference changes with span (p < 0.01), which a single RC mismatch does not "
                     "do: another effect, or errors too small (trains of one track are correlated)")
    if abs(kappa) < 2 * kerr:
        notes.append("kappa consistent with 0 (< 2 sigma): no RC mismatch seen at this precision")
    elif kappa + 1 / TAU_ASSUMED_US > 0:
        notes.append("if a single-stage RC mismatch: tau_true = %.2f ms (%.2f-%.2f ms at 1 sigma) for tau_assumed %.1f ms"
                     % tuple([1 / (k + 1 / TAU_ASSUMED_US) / 1000 if k + 1 / TAU_ASSUMED_US > 0 else np.inf
                              for k in (kappa, kappa + kerr, kappa - kerr)] + [TAU_ASSUMED_US / 1000]))
    else:
        notes.append("kappa below -1/tau_assumed = %.2f /ms: later hits gain more than any single-stage tau_true "
                     "(even infinite) allows; not a single RC mismatch" % (-1000 / TAU_ASSUMED_US))
    return "; ".join(notes)


def select(t, rr_min, min_hits):
    s = t[(t.track_score > 0.5) & (t.chi2_mu < t.chi2_p) & (t.rr_min > rr_min) & (t.n_hits >= min_hits)
          & np.isfinite(t.slope_frac_per_us)]
    return s.assign(lastfirst=s.q_last / s.q_first.replace(0, np.nan))


def binned_median(d, col, w):
    b = np.digitize(d.span_us.to_numpy(float), SPAN_BINS) - 1
    v = d[col].to_numpy(float)
    out = np.full((len(SPAN_BINS) - 1, 2), np.nan)
    for i in range(len(SPAN_BINS) - 1):
        m = (b == i) & np.isfinite(v)
        out[i] = median_err(v[m], w[m]) if m.any() else (np.nan, np.nan)
    return out


def per_plane(t, plane):
    return t[t.index.get_level_values("plane") == plane]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mc", nargs="+", required=True)
    ap.add_argument("--onbeam", nargs="+", required=True)
    ap.add_argument("--offbeam", nargs="+", required=True)
    ap.add_argument("-o", "--outdir", default="hittrain_plots")
    ap.add_argument("--area-norm", action="store_true", help="without POT / gate information (as plot_data_mc.py)")
    ap.add_argument("--offbeam-frac", type=float, help="off-beam share of the prediction, with --area-norm")
    ap.add_argument("--rr-min", type=float, default=30.0, help="residual range of every hit in the train > this [cm]")
    ap.add_argument("--min-hits", type=int, default=2, help="hits per train (default 2: most trains of MIP tracks have 2)")
    ap.add_argument("-j", "--nproc", type=int)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    mc, _, pot_mc = load_mc(args.mc, None, nproc=args.nproc)
    on, pot_on, _ = load_data(args.onbeam, args.nproc, uboone_bdt=False)
    off, _, gates_off = load_data(args.offbeam, args.nproc, uboone_bdt=False)
    sc_mc, sc_off, norm_line = normalisation(mc, on, off, pot_mc, pot_on, gates_off, args.area_norm, args.offbeam_frac)
    inputs = {"mc": args.mc, "on": args.onbeam, "off": args.offbeam}
    tab = {k: load_table(v, TRAIN_KEY, nproc=args.nproc) for k, v in inputs.items()}
    if any(v is None for v in tab.values()):
        sys.exit("no %s table in some inputs: remake the dataframes with the current configs" % TRAIN_KEY)
    tab = {k: select(v, args.rr_min, args.min_hits) for k, v in tab.items()}
    weight = {"mc": sc_mc, "on": 1.0, "off": sc_off}

    c = np.sqrt(SPAN_BINS[1:] * SPAN_BINS[:-1])
    lines = [norm_line, "trains: muon-like Pandora tracks, every hit at rr > %g cm, >= %d hits; slopes in %%/us"
             % (args.rr_min, args.min_hits), "lifetime slope for comparison: %.4f %%/us at %.0f ms"
             % (-100 / LIFETIME_US, LIFETIME_US / 1000), ""]
    figs = {k: plt.subplots(2, 3, figsize=(17, 8.5), sharex=True, gridspec_kw=dict(height_ratios=[3, 2]))
            for k in ("slope", "lastfirst")}
    fsh, ash = plt.subplots(2, 3, figsize=(17, 8))
    for pl in range(3):
        d = {k: per_plane(v, pl) for k, v in tab.items()}
        pred = pd.concat([d["mc"], d["off"]])
        wpred = np.r_[np.full(len(d["mc"]), weight["mc"]), np.full(len(d["off"]), weight["off"])]
        for key, col, scale, lab in [("slope", "slope_frac_per_us", 100.0, "median slope [%/µs]"),
                                     ("lastfirst", "lastfirst", 1.0, "median last / first hit")]:
            fig, axs = figs[key]
            a, r = axs[0, pl], axs[1, pl]
            Pd = binned_median(d["on"], col, np.ones(len(d["on"])))
            Pp = binned_median(pred, col, wpred)
            Pm = binned_median(d["mc"], col, np.full(len(d["mc"]), 1.0))
            a.plot(c, scale * Pp[:, 0], "-", color=INK, lw=2, label="prediction (MC + off-beam)")
            a.plot(c, scale * Pm[:, 0], "-", color=BLUE, lw=1.2, label="MC only")
            a.errorbar(c, scale * Pd[:, 0], scale * Pd[:, 1], fmt="o", color=INK, ms=4, label="on-beam data")
            a.set_xscale("log")
            a.set_title("plane %d%s: %d data / %d MC trains" % (pl, " (collection)" if pl == 2 else " (induction)",
                                                                len(d["on"]), len(d["mc"])), fontsize=10.5, loc="left")
            if pl == 0:
                a.set_ylabel(lab)
                a.legend(frameon=False, fontsize=9)
            if key == "slope":
                diff = scale * (Pd[:, 0] - Pp[:, 0])
                derr = scale * np.hypot(Pd[:, 1], Pp[:, 1])
                r.errorbar(c, diff, derr, fmt="o", color=INK, ms=4)
                r.axhline(0, color=INK, lw=1)
                ok = np.isfinite(diff) & np.isfinite(derr) & (derr > 0)
                if ok.any():
                    w = 1 / derr[ok] ** 2
                    m, e = np.sum(w * diff[ok]) / np.sum(w), 1 / np.sqrt(np.sum(w))
                    chi2 = np.sum(w * (diff[ok] - m) ** 2)
                    r.axhline(m, color=BLUE, lw=1.5, ls="--", label="mean = −κ = %.4f ± %.4f %%/µs" % (m, e))
                    r.legend(frameon=False, fontsize=9)
                    kappa, kerr = -m / 100.0, e / 100.0      # per us
                    ndf = int(ok.sum()) - 1
                    p_const = chi2_dist.sf(chi2, ndf) if ndf > 0 else np.nan
                    lines.append("plane %d: kappa = 1/tau_true - 1/tau_assumed = %+.2f +- %.2f /ms (%.1f sigma); "
                                 "chi2/ndf of a span-independent difference %.1f/%d (p = %.2g)"
                                 % (pl, 1000 * kappa, 1000 * kerr, abs(kappa) / kerr, chi2, ndf, p_const))
                    lines.append("         " + interpret(kappa, kerr, p_const))
                r.set_ylabel("data − pred. [%/µs]" if pl == 0 else "")
            else:
                ratio = Pd[:, 0] / Pp[:, 0]
                r.errorbar(c, ratio, ratio * np.hypot(Pd[:, 1] / Pd[:, 0], Pp[:, 1] / Pp[:, 0]), fmt="o", color=INK, ms=4)
                r.axhline(1, color=INK, lw=1)
                r.set_ylabel("data / pred." if pl == 0 else "")
            r.set_xlabel("train span [µs]")
            r.set_xlim(SPAN_BINS[0], SPAN_BINS[-1])
        # shapes: span and hits per train
        for k2, (col, bins, xl) in enumerate([("span_us", SPAN_BINS, "train span [µs]"),
                                               ("n_hits", np.arange(2.5, 21.5), "hits per train")]):
            a = ash[k2, pl]
            hd = np.histogram(d["on"][col], bins)[0].astype(float)
            hp = np.histogram(pred[col], bins, weights=wpred)[0]
            if hd.sum() > 0 and hp.sum() > 0:
                a.stairs(hp / hp.sum(), bins, color=INK, lw=1.8, label="prediction")
                cc = np.sqrt(bins[1:] * bins[:-1]) if col == "span_us" else 0.5 * (bins[1:] + bins[:-1])
                a.errorbar(cc, hd / hd.sum(), np.sqrt(hd) / hd.sum(), fmt="o", color=INK, ms=3, label="on-beam data")
            if col == "span_us":
                a.set_xscale("log")
            a.set_xlabel(xl)
            a.set_ylim(bottom=0)
            if pl == 0:
                a.set_ylabel("fraction of trains")
                a.legend(frameon=False, fontsize=9)
            if k2 == 0:
                a.set_title("plane %d" % pl, fontsize=11, loc="left")
    for key, (fig, axs) in figs.items():
        fig.suptitle("Pandora hit trains (consecutive hits of one track on one wire): %s vs train span (log axis)\n"
                     "an RC mismatch gives a span-independent data − MC slope" % (
                         "fractional slope of the charge per tick" if key == "slope" else "charge per tick, last / first hit"),
                     x=0.01, ha="left", fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        fig.savefig(os.path.join(args.outdir, "hittrain_%s_vs_span.png" % key), dpi=120, bbox_inches="tight")
        plt.close(fig)
    fsh.suptitle("Hit-train shapes, data vs prediction (unit area)", x=0.01, ha="left", fontsize=11)
    fsh.tight_layout(rect=(0, 0, 1, 0.96))
    fsh.savefig(os.path.join(args.outdir, "hittrain_shapes.png"), dpi=120, bbox_inches="tight")
    plt.close(fsh)
    open(os.path.join(args.outdir, "hittrain_summary.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("wrote the plots to %s" % args.outdir)


if __name__ == "__main__":
    main()
