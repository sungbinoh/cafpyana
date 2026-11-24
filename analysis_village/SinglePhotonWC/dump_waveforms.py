#!/usr/bin/env python3
"""Raw and deconvolved waveforms of chosen channels from SBND reco1 art-ROOT
files, to look at the AC-coupling (RC) tail after long single-wire signals in
data and MC. Needs the SBND software environment (gallery + lardataobj
dictionaries, e.g. after `setup sbndcode ...`): run it there, not in the
cafpyana venv.

Input: reco1 files and the CSV of select_waveform_channels.py (run, subrun,
event, channel, plane, ...). For each listed (run, subrun, event, channel)
found in the files it reads
  raw::RawDigit   --raw-label     (MC: simtpc2d:daq; data: the DAQ decoder's
                                   label, check with `lar -c eventdump.fcl`)
  recob::Wire     --wire-labels   (MC: simtpc2d:gauss simtpc2d:dnnsp;
                                   data: sptpc2d:gauss sptpc2d:dnnsp)
and writes
  waveforms.npz   (np.load(..., allow_pickle=True)) per channel: raw ADC - pedestal, each wire label's
                  deconvolved signal (zeros outside its ROIs), and the CSV row
  wf_<run>_<sub>_<evt>_ch<channel>.png
                  raw (top) and deconvolved (below) over the whole readout,
                  the signal window shaded, plus a zoom of the raw baseline
                  after the signal (where an RC undershoot sits)
  tail_average_plane<p>.png
                  raw ADC - pedestal averaged over the listed channels of a
                  plane, aligned at the end of the signal (the longest ROI of
                  the first wire label) and normalised to the signal's
                  deconvolved charge (induction-plane raw signals are bipolar,
                  their raw sum is ~0): an undershoot that does not return to 0 within ~1 ms
                  is the RC tail; compare data and MC by running on each

The raw ADC includes noise and the electronics' own baseline; the
pedestal is RawDigit::GetPedestal() (0 in some productions: then the median
of the waveform is used). Untested here (no LArSoft in the cafpyana
container): if the gallery calls need adjusting for your release, the reading
is confined to read_event().

  python dump_waveforms.py -i reco1_*.root --channels channels.csv -o waveforms/ \\
      --raw-label simtpc2d:daq --wire-labels simtpc2d:gauss simtpc2d:dnnsp
"""

import argparse
import csv
import os

import numpy as np

TICK_US = 0.5
TAIL_TICKS = 3000      # 1.5 ms after the signal
PRE_TICKS = 200


def load_gallery():
    import ROOT
    ROOT.gROOT.SetBatch(True)
    for h in ["gallery/Event.h", "lardataobj/RawData/RawDigit.h", "lardataobj/RawData/raw.h",
              "lardataobj/RecoBase/Wire.h", "canvas/Utilities/InputTag.h"]:
        ROOT.gInterpreter.Declare('#include "%s"' % h)
    return ROOT


def handle(ROOT, ev, cls, label):
    """ev.getValidHandle<cls>(tag), for both PyROOT template spellings."""
    tag = ROOT.art.InputTag(label)
    vec = ROOT.std.vector(cls)
    try:
        return ev.getValidHandle[vec](tag).product()
    except TypeError:
        return ev.getValidHandle(vec)(tag).product()


def read_event(ROOT, ev, channels, raw_label, wire_labels):
    """{channel: {"raw": adc - pedestal, label: deconvolved}} for the requested channels of the current event."""
    out = {c: {} for c in channels}
    for d in handle(ROOT, ev, ROOT.raw.RawDigit, raw_label):
        c = int(d.Channel())
        if c not in out:
            continue
        adc = ROOT.std.vector("short")(int(d.Samples()))
        ROOT.raw.Uncompress(d.ADCs(), adc, d.Compression())
        a = np.array(adc, dtype=float)
        ped = float(d.GetPedestal())
        out[c]["raw"] = a - (ped if ped != 0 else np.median(a))
    for lab in wire_labels:
        for w in handle(ROOT, ev, ROOT.recob.Wire, lab):
            c = int(w.Channel())
            if c in out:
                out[c][lab] = np.array(w.Signal(), dtype=float)
    return out


def signal_window(sig, hint):
    """(start, end) tick of the longest stretch of nonzero deconvolved signal; the CSV hint if none."""
    nz = np.flatnonzero(sig > 0)
    if not len(nz):
        return hint
    breaks = np.flatnonzero(np.diff(nz) > 1)
    starts, ends = np.r_[nz[0], nz[breaks + 1]], np.r_[nz[breaks], nz[-1]]
    i = np.argmax(ends - starts)
    return int(starts[i]), int(ends[i])


def plot_channel(path, row, wf, labels, win):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n = 2 + len(labels)
    fig, axs = plt.subplots(n, 1, figsize=(12, 2.4 * n), sharex=False)
    t = np.arange(len(wf["raw"])) * TICK_US
    s0, s1 = win
    axs[0].plot(t, wf["raw"], lw=0.7, color="#1f1f1e")
    axs[0].set_ylabel("raw ADC − ped.")
    for k, lab in enumerate(labels):
        if lab in wf:
            axs[k + 1].plot(np.arange(len(wf[lab])) * TICK_US, wf[lab], lw=0.8, color="#2a78d6")
        axs[k + 1].set_ylabel(lab, fontsize=9)
    for a in axs[:-1]:
        a.axvspan(s0 * TICK_US, s1 * TICK_US, color="#eb6834", alpha=0.15, lw=0)
        a.axhline(0, color="#8a8984", lw=0.8)
    z = axs[-1]   # zoom on the baseline after the signal
    lo, hi = max(0, s0 - PRE_TICKS), min(len(wf["raw"]), s1 + TAIL_TICKS)
    z.plot(t[lo:hi], wf["raw"][lo:hi], lw=0.7, color="#1f1f1e")
    k = 50   # running mean, 25 us
    if hi - lo > k:
        z.plot(t[lo:hi][k - 1:], np.convolve(wf["raw"][lo:hi], np.ones(k) / k, "valid"), lw=1.8, color="#eb6834",
               label="25 µs running mean")
    z.axhline(0, color="#8a8984", lw=0.8)
    z.axvspan(s0 * TICK_US, s1 * TICK_US, color="#eb6834", alpha=0.15, lw=0)
    after = wf["raw"][s1 + 20:hi]
    ymax = 5 * np.std(after) + 1 if len(after) else 10
    z.set_ylim(-ymax, ymax)
    z.set_xlabel("time [µs]")
    z.set_ylabel("raw, zoom")
    z.legend(frameon=False, fontsize=8)
    axs[0].set_title("run %s subrun %s event %s, channel %s (plane %s); shaded: the signal (%.0f µs)"
                     % (row["run"], row["subrun"], row["event"], row["channel"], row["plane"], (s1 - s0) * TICK_US),
                     loc="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def plot_tail_average(path, plane, tails):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    T = np.array(tails)
    m, e = np.nanmean(T, axis=0), np.nanstd(T, axis=0) / np.sqrt(np.sum(np.isfinite(T), axis=0).clip(1))
    t = (np.arange(T.shape[1]) - PRE_TICKS) * TICK_US
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.fill_between(t, m - e, m + e, color="#2a78d6", alpha=0.25, lw=0)
    ax.plot(t, m, color="#2a78d6", lw=0.6, alpha=0.6, label="average")
    k = 50   # running mean, 25 us
    mm = np.convolve(np.nan_to_num(m), np.ones(k) / k, "same")
    mm[~np.isfinite(m)] = np.nan
    ax.plot(t, mm, color="#eb6834", lw=2, label="25 µs running mean")
    ax.legend(frameon=False, fontsize=9)
    ax.axhline(0, color="#8a8984", lw=0.8)
    ax.axvline(0, color="#1f1f1e", lw=0.8, ls="--")
    ax.set_xlabel("time after the end of the signal [µs]")
    ax.set_ylabel("raw ADC − ped., per unit\nof the signal's deconvolved charge")
    lim = 6 * np.nanmax(e[PRE_TICKS + 40:]) if T.shape[1] > PRE_TICKS + 40 else 1
    ax.set_ylim(-lim, lim)
    ax.set_title("plane %d: average after %d long signals (aligned at their end)" % (plane, len(T)), loc="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--inputs", nargs="+", required=True, help="reco1 art-ROOT files")
    ap.add_argument("--channels", required=True, help="CSV from select_waveform_channels.py")
    ap.add_argument("-o", "--outdir", default="waveforms")
    ap.add_argument("--raw-label", default="simtpc2d:daq")
    ap.add_argument("--wire-labels", nargs="+", default=["simtpc2d:gauss", "simtpc2d:dnnsp"])
    ap.add_argument("--max-plots", type=int, default=30, help="single-channel plots written (all go in the npz and the average)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    rows = list(csv.DictReader(open(args.channels)))
    want = {}
    for r in rows:
        want.setdefault((int(r["run"]), int(r["subrun"]), int(r["event"])), []).append(r)
    ROOT = load_gallery()
    ev = ROOT.gallery.Event(ROOT.std.vector("std::string")(args.inputs))
    saved, tails, nplot = {}, {}, 0
    while not ev.atEnd():
        aux = ev.eventAuxiliary()
        key = (int(aux.run()), int(aux.subRun()), int(aux.event()))
        if key in want:
            chans = {int(r["channel"]): r for r in want[key]}
            data = read_event(ROOT, ev, list(chans), args.raw_label, args.wire_labels)
            for c, wf in data.items():
                if "raw" not in wf:
                    print("run %d subrun %d event %d: channel %d not in %s" % (*key, c, args.raw_label))
                    continue
                r = chans[c]
                lab0 = next((l for l in args.wire_labels if l in wf), None)
                win = signal_window(wf[lab0], (int(r["tick_start"]), int(r["tick_end"]))) if lab0 else \
                    (int(r["tick_start"]), int(r["tick_end"]))
                name = "wf_%d_%d_%d_ch%d" % (*key, c)
                saved[name] = dict(row=r, window=win, **wf)
                if nplot < args.max_plots:
                    plot_channel(os.path.join(args.outdir, name + ".png"), r, wf, args.wire_labels, win)
                    nplot += 1
                s0, s1 = win
                # normalise by the deconvolved charge: induction-plane raw signals are bipolar, their sum ~0
                q = wf[lab0][s0:s1 + 1].sum() if lab0 else wf["raw"][s0:s1 + 1].sum()
                seg = np.full(PRE_TICKS + TAIL_TICKS, np.nan)
                a, b = max(0, s1 - PRE_TICKS), min(len(wf["raw"]), s1 + TAIL_TICKS)
                seg[a - (s1 - PRE_TICKS):b - (s1 - PRE_TICKS)] = wf["raw"][a:b]
                if q > 0:
                    tails.setdefault(int(r["plane"]), []).append(seg / q)
        ev.next()
    np.savez_compressed(os.path.join(args.outdir, "waveforms.npz"),
                        **{k: np.array(v, dtype=object) for k, v in saved.items()})
    for p, tl in sorted(tails.items()):
        plot_tail_average(os.path.join(args.outdir, "tail_average_plane%d.png" % p), p, tl)
    print("found %d of %d listed channels; wrote %d plots and the averages to %s"
          % (len(saved), len(rows), nplot, args.outdir))


if __name__ == "__main__":
    main()
