#!/usr/bin/env python3
"""Data / MC overlay of the average in-signal shape of long single-wire
signals, from dump_waveforms.py's shapes_plane<p>.npy (each row: one
signal's deconvolved waveform, from its start, NaN after its end).

The data (BNB + light) and MC (CV) samples hold different track mixes, and
the shape of a wire's signal depends on the track geometry as well as on the
electronics, so the signals are compared like for like:
  * in classes of signal duration (DURATION_CLASSES, from where each row
    turns NaN), so the two samples' different length mixes do not shape the
    average;
  * each signal normalised to its own mean over REF_US (20-40 us: inside
    every signal, all are >= 40 us), so the curves agree there by
    construction and the comparison is of how the signal evolves;
  * against absolute time since the start (an RC sag is a fall with time)
    and against the fraction of the signal's duration (geometry scales with
    it).
Rows of the figures are planes, columns duration classes; the lower panels
are data / MC.

  python analysis_village/SinglePhotonWC/plot_waveform_shapes.py --data <dir> --mc <dir> -o <outdir>
"""

import argparse
import glob
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TICK_US = 0.5
REF_US = (20.0, 40.0)
DURATION_CLASSES = [(40, 60), (60, 100), (100, 200), (200, 400)]   # us; 400 = the stored length
NREL = 20
DATA, MC = "#1f1f1e", "#2a78d6"


def load(d):
    out = {}
    for f in sorted(glob.glob(os.path.join(d, "shapes_plane*.npy"))):
        out[int(os.path.basename(f)[len("shapes_plane"):-4])] = np.load(f)
    if not out:
        raise SystemExit("no shapes_plane*.npy in %s (dump_waveforms.py writes them next to the plots)" % d)
    return out


def prepare(S):
    """(signals renormalised to their REF_US mean, durations in us)."""
    valid = np.isfinite(S)
    n = np.where(valid.all(axis=1), S.shape[1], np.argmin(valid, axis=1))
    dur = n * TICK_US
    a, b = int(REF_US[0] / TICK_US), int(REF_US[1] / TICK_US)
    ref = np.nanmean(S[:, a:b], axis=1)
    ok = (ref > 0) & (dur >= REF_US[1])
    return S[ok] / ref[ok, None], dur[ok]


def average(S):
    nn = np.sum(np.isfinite(S), axis=0)
    m = np.where(nn >= 5, np.nanmean(S, axis=0), np.nan)
    e = np.where(nn >= 5, np.nanstd(S, axis=0) / np.sqrt(nn.clip(1)), np.nan)
    return m, e


def relative(S, dur):
    """Each signal resampled to NREL bins of the fraction of its duration."""
    out = np.full((len(S), NREL), np.nan)
    for i, (s, d) in enumerate(zip(S, dur)):
        n = int(round(d / TICK_US))
        edges = np.linspace(0, n, NREL + 1).astype(int)
        out[i] = [np.nanmean(s[edges[k]:max(edges[k + 1], edges[k] + 1)]) for k in range(NREL)]
    return out


def figure(kind, data, mc, out):
    planes = sorted(set(data) & set(mc))
    nc = len(DURATION_CLASSES)
    fig = plt.figure(figsize=(4.2 * nc, 4.8 * len(planes)))
    gs = fig.add_gridspec(3 * len(planes), nc, height_ratios=[3, 1.2, 0.55] * len(planes), hspace=0.08, wspace=0.25)
    rows = []
    for i, p in enumerate(planes):
        (Sd, dd), (Sm, dm) = data[p], mc[p]
        for j, (lo, hi) in enumerate(DURATION_CLASSES):
            a = fig.add_subplot(gs[3 * i, j])
            r = fig.add_subplot(gs[3 * i + 1, j], sharex=a)
            sd, sm = (dd >= lo) & (dd < hi), (dm >= lo) & (dm < hi)
            if kind == "abs":
                md, ed = average(Sd[sd][:, :int(hi / TICK_US)])
                mm, em = average(Sm[sm][:, :int(hi / TICK_US)])
                x = np.arange(len(md)) * TICK_US
            else:
                md, ed = average(relative(Sd[sd], dd[sd]))
                mm, em = average(relative(Sm[sm], dm[sm]))
                x = (np.arange(NREL) + 0.5) / NREL
            for m, e, col, lab, n in [(mm, em, MC, "MC (CV)", sm.sum()), (md, ed, DATA, "data", sd.sum())]:
                a.fill_between(x, m - e, m + e, color=col, alpha=0.2, lw=0)
                a.plot(x, m, color=col, lw=1.5, label="%s: %d" % (lab, n))
            ratio = md / mm
            rerr = ratio * np.hypot(ed / md, em / mm)
            r.fill_between(x, ratio - rerr, ratio + rerr, color=DATA, alpha=0.25, lw=0)
            r.plot(x, ratio, color=DATA, lw=1.2)
            r.axhline(1, color="#8a8984", lw=0.8)
            r.set_ylim(0.6, 1.4)
            a.axhline(1, color="#8a8984", lw=0.6)
            if kind == "abs":
                a.axvspan(*REF_US, color="#eb6834", alpha=0.12, lw=0)
            a.legend(frameon=False, fontsize=8, loc="upper right")
            plt.setp(a.get_xticklabels(), visible=False)
            if i == 0:
                a.set_title("signals %d-%d µs" % (lo, hi), fontsize=10)
            if j == 0:
                a.set_ylabel("plane %d\nsignal / its %g-%g µs mean" % (p, *REF_US), fontsize=9)
                r.set_ylabel("data / MC", fontsize=9)
            if True:
                r.set_xlabel("time since the start [µs]" if kind == "abs" else "fraction of the signal's duration", fontsize=9)
            for k in range(len(x)):
                rows.append(dict(kind=kind, plane=p, dur_lo=lo, dur_hi=hi, x=x[k], data=md[k], data_err=ed[k],
                                 mc=mm[k], mc_err=em[k]))
    fig.suptitle("Average shape of long single-wire signals (DNN-ROI deconvolved), data vs MC, by plane and duration; "
                 "%s" % ("normalised in the shaded window" if kind == "abs" else "against the fraction of each signal"),
                 x=0.01, ha="left", fontsize=11)
    fig.savefig(os.path.join(out, "waveform_shapes_%s.png" % ("time" if kind == "abs" else "fraction")), dpi=110,
                bbox_inches="tight")
    plt.close(fig)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="dump_waveforms.py output directory of the data")
    ap.add_argument("--mc", required=True, help="the same for MC")
    ap.add_argument("-o", "--outdir", default="waveform_shapes")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    data = {p: prepare(S) for p, S in load(args.data).items()}
    mc = {p: prepare(S) for p, S in load(args.mc).items()}
    rows = figure("abs", data, mc, args.outdir) + figure("rel", data, mc, args.outdir)
    import csv
    with open(os.path.join(args.outdir, "waveform_shapes.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for p in sorted(set(data) & set(mc)):
        print("plane %d: %d data / %d MC signals" % (p, len(data[p][0]), len(mc[p][0])))
    print("wrote %s/waveform_shapes_time.png, _fraction.png and .csv" % args.outdir)


if __name__ == "__main__":
    main()
