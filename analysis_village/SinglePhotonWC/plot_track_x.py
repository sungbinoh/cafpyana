#!/usr/bin/env python3
"""Data / MC x (drift coordinate) of track ends, to check the time -> x
conversion (drift speed, t0, space-charge position correction): the anode
planes (x = +-200 cm) and the cathode (x = 0) are edges in the end-point
distributions, and the drift length (200 cm) is the edge of the x extent of
tracks crossing from the cathode to an anode. A drift-speed mismatch shifts
the anode edges and the 200 cm edge, scaling with the drift distance; a t0
offset shifts the anode edges the same way in both TPCs (towards / away from
the cathode); a cathode-edge shift alone points to the cathode-side position
correction.

Inputs as plot_dqdx_data_mc.py; uses the pandtrk table (Pandora tracks
>= 30 cm: start_x / end_x) and the trkseg table (Wire-Cell tracks: the
lowest and highest x of their 4 cm pieces, so ~2 cm inside the true ends).
Shapes: each distribution is normalised to unit area (data) and the
prediction (MC + off-beam, mixed as plot_data_mc.normalisation) likewise.
--nu-like keeps the Pandora tracks of slices with nu_score > 0.4 that are not
clear cosmics (in time with the beam, so x is t0-correct).

Plots (one figure per table): end-point x over the full drift and zoomed on
the anodes and the cathode, and the x extent of tracks ending within 10 cm
of the cathode, each with a data / prediction ratio; edges.txt lists the
half-height edge positions from data and prediction.

  python analysis_village/SinglePhotonWC/plot_track_x.py \\
      --mc <MC> --onbeam <on-beam> --offbeam <off-beam> -o track_x_plots/ [--area-norm] [--nu-like]
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
from analysis_village.SinglePhotonWC.cutflow import load as load_mc  # noqa: E402
from analysis_village.SinglePhotonWC.plot_data_mc import load_data, normalisation  # noqa: E402
from analysis_village.SinglePhotonWC.plot_dqdx_data_mc import PAND_KEY, SEG_KEY, load_table  # noqa: E402

plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), "dune.mplstyle"))

PANELS = [   # name, bins, title
    ("full", np.arange(-210, 210.1, 2), "all x"),
    ("anode_east", np.arange(-215, -184.9, 0.5), "east anode (x = -200 cm)"),
    ("cathode", np.arange(-15, 15.01, 0.5), "cathode (x = 0)"),
    ("anode_west", np.arange(185, 215.01, 0.5), "west anode (x = +200 cm)"),
]
SPAN_BINS = np.arange(150, 220.1, 2)
CATHODE_END = 10.0   # cm, a track end this close to the cathode starts the drift-length measurement


def pandora_ends(t, nu_like=False):
    """(end-point x, x extent of cathode-ending tracks) of the Pandora tracks, one row per track."""
    t = t[t.index.get_level_values("plane") == 2]   # the track variables repeat on each plane
    if nu_like:
        t = t[(t.nu_score > 0.4) & (t.is_clear_cosmic == 0)]
    ends = np.r_[t.start_x.to_numpy(float), t.end_x.to_numpy(float)]
    return ends, _cathode_spans(t.start_x.to_numpy(float), t.end_x.to_numpy(float))


def wirecell_ends(s):
    """The same for Wire-Cell tracks, from the lowest and highest x of their pieces."""
    g = s.x.groupby(level=["file", "__ntuple", "entry", "nu_index", "kine_index"])
    lo, hi = g.min().to_numpy(float), g.max().to_numpy(float)
    return np.r_[lo, hi], _cathode_spans(lo, hi)


def _cathode_spans(a, b):
    near = (np.abs(a) < CATHODE_END) | (np.abs(b) < CATHODE_END)
    return np.abs(a - b)[near]


def half_height_edge(x, bins, rising):
    """Position where the histogram crosses half its plateau (the mean of the upper third of bins),
    from the first rising / last falling crossing; NaN if none."""
    h, _ = np.histogram(x, bins)
    if not h.any():
        return np.nan
    top = np.sort(h)[-max(1, len(h) // 3):].mean()
    c = 0.5 * (bins[1:] + bins[:-1])
    idx = np.flatnonzero(h >= 0.5 * top)
    if not len(idx):
        return np.nan
    i = idx[0] if rising else idx[-1]
    j = i - 1 if rising else i + 1
    if j < 0 or j >= len(h):
        return c[i]
    return np.interp(0.5 * top, sorted([h[j], h[i]]), [c[j], c[i]] if h[j] < h[i] else [c[i], c[j]])


def _shape(x, bins, w=None):
    h, _ = np.histogram(x, bins, weights=w)
    h2, _ = np.histogram(x, bins, weights=None if w is None else w ** 2)
    s = h.sum()
    return (h / s, np.sqrt(h2) / s) if s > 0 else (h * np.nan, h * np.nan)


def plot_table(name, label, data, mc, off, sc_mc, sc_off, out):
    """data, mc, off: (end-point x, cathode-ending x extents)."""
    fig, axs = plt.subplots(2, len(PANELS) + 1, figsize=(6 * (len(PANELS) + 1), 7.5), sharex="col",
                            gridspec_kw=dict(height_ratios=[3, 1], hspace=0.05))
    lines = []
    for k, (key, bins, title) in enumerate(PANELS + [("span", SPAN_BINS, "x extent, tracks ending at the cathode")]):
        i = 1 if key == "span" else 0
        xd, xm, xo = data[i], mc[i], off[i]
        xp = np.r_[xm, xo]
        wp = np.r_[np.full(len(xm), sc_mc), np.full(len(xo), sc_off)]
        hd, ed = _shape(xd, bins)
        hp, ep = _shape(xp, bins, wp)
        c = 0.5 * (bins[1:] + bins[:-1])
        a, r = axs[0, k], axs[1, k]
        a.stairs(hp, bins, color="#1f1f1e", lw=1.8, label="prediction (MC + off-beam)")
        a.errorbar(c, hd, ed, fmt="o", color="k", ms=2.5, label="on-beam data")
        a.set_title(title, fontsize=10.5, loc="left")
        a.set_ylim(bottom=0)
        ok = hp > 0
        r.errorbar(c[ok], hd[ok] / hp[ok], np.hypot(ed[ok] / hp[ok], hd[ok] * ep[ok] / hp[ok] ** 2), fmt="o", color="k", ms=2)
        r.axhline(1, color="k", lw=1)
        r.set_ylim(0.5, 1.5)
        r.set_xlabel("x extent [cm]" if key == "span" else "track end x [cm]")
        if k == 0:
            a.set_ylabel("fraction / bin")
            r.set_ylabel("data / pred.")
            a.legend(frameon=False, fontsize=9)
        if key != "full":   # edge positions: the anodes fall off outwards, the cathode is a dip
            pairs = {"anode_east": [("east anode", True)], "anode_west": [("west anode", False)],
                     "cathode": [], "span": [("drift length", False)]}[key]
            for what, rising in pairs:
                e_d = half_height_edge(xd, bins, rising)
                e_p = half_height_edge(xp, bins, rising)
                for v, col in [(e_d, "k"), (e_p, "#2a78d6")]:
                    if np.isfinite(v):
                        a.axvline(v, color=col, ls="--", lw=1)
                lines.append("%-14s %-12s data %7.2f cm  prediction %7.2f cm  data - pred. %+5.2f cm  (%d data / %d MC entries)"
                             % (label, what, e_d, e_p, e_d - e_p, ((xd >= bins[0]) & (xd < bins[-1])).sum(),
                                ((xm >= bins[0]) & (xm < bins[-1])).sum()))
    fig.suptitle("%s: track end-point x, shapes (unit area); dashed: half-height edges (black data, blue prediction)"
                 % label, x=0.01, ha="left", fontsize=11)
    fig.savefig(os.path.join(out, "track_x_%s.png" % name), dpi=120, bbox_inches="tight")
    plt.close(fig)
    return lines


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mc", nargs="+", required=True)
    ap.add_argument("--onbeam", nargs="+", required=True)
    ap.add_argument("--offbeam", nargs="+", required=True)
    ap.add_argument("-o", "--outdir", default="track_x_plots")
    ap.add_argument("--area-norm", action="store_true", help="without POT / gate information (as plot_data_mc.py)")
    ap.add_argument("--offbeam-frac", type=float, help="off-beam share of the prediction, with --area-norm")
    ap.add_argument("--nu-like", action="store_true", help="Pandora: only nu-like, not clear-cosmic slices")
    ap.add_argument("-j", "--nproc", type=int)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    mc, _, pot_mc = load_mc(args.mc, None, nproc=args.nproc)
    on, pot_on, _ = load_data(args.onbeam, args.nproc, uboone_bdt=False)
    off, _, gates_off = load_data(args.offbeam, args.nproc, uboone_bdt=False)
    sc_mc, sc_off, norm_line = normalisation(mc, on, off, pot_mc, pot_on, gates_off, args.area_norm, args.offbeam_frac)
    inputs = {"mc": args.mc, "on": args.onbeam, "off": args.offbeam}
    lines = [norm_line, ""]
    for key, name, label, fn in [(PAND_KEY, "pandora", "Pandora", lambda t: pandora_ends(t, args.nu_like)),
                                 (SEG_KEY, "wirecell", "Wire-Cell", wirecell_ends)]:
        t = {k: load_table(v, key, nproc=args.nproc) for k, v in inputs.items()}
        if any(v is None for v in t.values()):
            print("no %s table in some inputs: %s skipped" % (key, label))
            continue
        lines += plot_table(name, label, fn(t["on"]), fn(t["mc"]), fn(t["off"]), sc_mc, sc_off, args.outdir)
    open(os.path.join(args.outdir, "edges.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("wrote the plots to %s" % args.outdir)


if __name__ == "__main__":
    main()
