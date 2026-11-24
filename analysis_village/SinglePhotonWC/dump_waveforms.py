#!/usr/bin/env python3
"""Raw and deconvolved waveforms of chosen channels from SBND reco1 art-ROOT
files, to look at the AC-coupling (RC) tail after long single-wire signals in
data and MC. Needs the SBND software environment (gallery + lardataobj
dictionaries, e.g. after `setup sbndcode ...`): run it there, not in the
cafpyana venv.

Input: reco1 files, and either the CSV of select_waveform_channels.py (run,
subrun, event, channel, plane, ...) or --auto N: then, in every event, the N
channels per plane with the longest single ROI (>= --min-us) of the first
wire label are taken, so any reco1 file works, data or MC. For each channel
it reads
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
  shape_average_plane<p>.png (+ shapes_plane<p>.npy)
                  the deconvolved signal / its mean, aligned at the start of
                  the signal and averaged per plane: the in-signal sag, the
                  part visible without raw digits (the Gen2 MC reco1 files
                  keep only simtpc2d:dnnsp)
  tail_average_plane<p>.png (only with raw digits)
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

Parallel: -j N runs N worker processes on slices of the file list (each into
<outdir>/part_<k>; their output is discarded, a failing worker's error
output is printed) and merges their averages into <outdir>; Ctrl-C or
killing the main process stops all workers; --merge DIR... merges existing output directories. Most of the time
is reading the files over xrootd, so a few workers per machine help;
--events-per-file K spreads the same number of events over more files.

  python dump_waveforms.py -i $(cat reco1_urls_mc.txt) -o waveforms_mc/ --auto 3 \\
      --raw-label simtpc2d:daq --wire-labels simtpc2d:gauss simtpc2d:dnnsp
  python dump_waveforms.py -i $(cat reco1_urls_data.txt) -o waveforms_data/ --auto 3 \\
      --raw-label daq --wire-labels sptpc2d:gauss sptpc2d:dnnsp

Data: the signal processing reads raw::RawDigit "daq" (wcsp_data_sbnd.fcl);
whether reco1 files keep it depends on the output's drop list. Without raw
digits only the deconvolved signals can be shown (no tails): check one file
with `lar -c eventdump.fcl -s <file> -n 1`.
"""

import argparse
import csv
import os
import time

try:
    from tqdm import tqdm
except ImportError:   # a minimal text progress bar (the LArSoft python has no pip to install tqdm)
    import sys as _sys

    class tqdm:
        def __init__(self, total=None, unit="it", desc="", disable=False, **_):
            self.total, self.unit, self.desc, self.disable = total, unit, desc, disable
            self.n, self.post, self.t0 = 0, "", time.time()
            self._draw()

        def _draw(self):
            if self.disable:
                return
            dt = time.time() - self.t0
            frac = self.n / self.total if self.total else 0
            bar = "#" * int(30 * frac) + "-" * (30 - int(30 * frac))
            eta = " eta %ds" % (dt / self.n * (self.total - self.n)) if self.n and self.total else ""
            _sys.stderr.write("\r%s [%s] %d/%s %s %ds%s %s" % (self.desc, bar, self.n, self.total, self.unit, dt, eta, self.post))
            _sys.stderr.flush()

        def update(self, n=1):
            self.n += n
            self._draw()

        def set_postfix_str(self, s, **_):
            self.post = s
            self._draw()

        def close(self):
            if not self.disable:
                _sys.stderr.write("\n")

import numpy as np

TICK_US = 0.5
SHAPE_TICKS = 800      # 400 us of signal in the shape average
TAIL_TICKS = 3000      # 1.5 ms after the signal
PRE_TICKS = 200


def load_gallery():
    import ROOT
    ROOT.gROOT.SetBatch(True)
    for h in ["gallery/Event.h", "lardataobj/RawData/RawDigit.h", "lardataobj/RawData/raw.h",
              "lardataobj/RecoBase/Wire.h", "canvas/Utilities/InputTag.h"]:
        ROOT.gInterpreter.Declare('#include "%s"' % h)
    return ROOT


_WARNED = set()


class ProductNotFound(Exception):
    pass


def _template(ROOT, tmpl, arg):
    """tmpl<arg> for both PyROOT spellings (tmpl[arg] in cppyy, tmpl(arg) in older PyROOT)."""
    try:
        return tmpl[arg]
    except TypeError:
        return tmpl(arg)


def handle(ROOT, ev, cls, label):
    """The std::vector<cls> product `label` of the current event, raising ProductNotFound when absent.
    gallery versions differ: getHandle (newer), getByLabel into a gallery::Handle (older), and as a last
    resort getValidHandle, whose exception on a missing product is turned into ProductNotFound."""
    tag = ROOT.art.InputTag(label)
    vec = ROOT.std.vector(cls)
    what = "no std::vector<%s> with label '%s'" % (getattr(cls, "__cpp_name__", getattr(cls, "__name__", cls)), label)
    if hasattr(ev, "getHandle"):
        h = _template(ROOT, ev.getHandle, vec)(tag)
        if not h.isValid():
            raise ProductNotFound(what)
        return h.product()
    if hasattr(ev, "getByLabel"):
        h = _template(ROOT, ROOT.gallery.Handle, vec)()
        try:
            ok = ev.getByLabel(tag, h)
        except TypeError:   # explicit template argument needed
            ok = _template(ROOT, ev.getByLabel, vec)(tag, h)
        if not ok:
            raise ProductNotFound(what)
        return h.product()
    try:
        return _template(ROOT, ev.getValidHandle, vec)(tag).product()
    except Exception as e:
        raise ProductNotFound("%s (%s)" % (what, type(e).__name__))


def list_products(ROOT, path, patterns=("recob::Wire", "raw::RawDigit")):
    """The module:instance labels of the art products of `patterns` in a file, from its Events tree branches
    (<class>s_<module>_<instance>_<process>.)."""
    f = ROOT.TFile.Open(path)
    t = f.Get("Events")
    out = []
    for b in t.GetListOfBranches():
        name = b.GetName()
        for p in patterns:
            if name.startswith(p.replace("::", "::") + "s_"):
                parts = name.rstrip(".").split("_")
                out.append("%-15s label %s:%s  (process %s)" % (p, parts[1], parts[2], parts[3] if len(parts) > 3 else "?"))
    f.Close()
    return out


def read_event(ROOT, ev, channels, raw_label, wire_labels):
    """{channel: {"raw": adc - pedestal, label: deconvolved}} for the requested channels of the current event.
    A missing product (e.g. raw digits dropped from the reco1 output) is skipped with a warning."""
    out = {c: {} for c in channels}
    try:
        digits = handle(ROOT, ev, ROOT.raw.RawDigit, raw_label)
    except Exception as e:
        if raw_label not in _WARNED:
            print("no raw::RawDigit '%s' (%s): deconvolved signals only (reported once)" % (raw_label, type(e).__name__))
            _WARNED.add(raw_label)
        digits = []
    for d in digits:
        c = int(d.Channel())
        if c not in out:
            continue
        adc = ROOT.std.vector("short")(int(d.Samples()))
        ROOT.raw.Uncompress(d.ADCs(), adc, d.Compression())
        a = np.array(adc, dtype=float)
        ped = float(d.GetPedestal())
        out[c]["raw"] = a - (ped if ped != 0 else np.median(a))
    for lab in wire_labels:
        try:
            wires = handle(ROOT, ev, ROOT.recob.Wire, lab)
        except Exception as e:
            if lab not in _WARNED:
                print("no recob::Wire '%s' (%s; reported once)" % (lab, type(e).__name__))
                _WARNED.add(lab)
            continue
        for w in wires:
            c = int(w.Channel())
            if c in out:
                out[c][lab] = np.array(w.Signal(), dtype=float)
    return out


# (first channel, plane) per block of the SBND channel map (make_wc_events.CHANNEL_BLOCKS; TPC 0 then TPC 1)
CHANNEL_BLOCKS = [(0, 0), (1984, 1), (3968, 2), (5638, 0), (7622, 1), (9606, 2), (11276, -1)]


def channel_plane(c):
    return [p for start, p in CHANNEL_BLOCKS if start <= c][-1]


def longest_rois(ROOT, ev, label, min_ticks, per_plane):
    """The channels with the longest single ROI of recob::Wire `label` in this event, per plane: {channel: row}."""
    best = {}
    try:
        wires = handle(ROOT, ev, ROOT.recob.Wire, label)
    except ProductNotFound as e:
        print("%s in this event: check the labels with --list-products" % e)
        return {}
    for w in wires:
        c = int(w.Channel())
        sig = np.array(w.Signal(), dtype=float)
        nz = np.flatnonzero(sig != 0)
        if len(nz) < min_ticks:
            continue
        brk = np.flatnonzero(np.diff(nz) > 1)
        st, en = np.r_[nz[0], nz[brk + 1]], np.r_[nz[brk], nz[-1]]
        i = int(np.argmax(en - st))
        if en[i] - st[i] + 1 >= min_ticks:
            best.setdefault(channel_plane(c), []).append((en[i] - st[i], c, int(st[i]), int(en[i])))
    out = {}
    for plane, lst in best.items():
        for n, c, a, b in sorted(lst, reverse=True)[:per_plane]:
            out[c] = dict(channel=c, plane=plane, tick_start=a, tick_end=b, duration_us=(b - a) * TICK_US)
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


def plot_channel_wires(path, row, wf, labels, win):
    """The deconvolved signals only (files without raw digits): the whole readout and a zoom on the signal."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    labs = [l for l in labels if l in wf]
    fig, axs = plt.subplots(len(labs) + 1, 1, figsize=(12, 2.6 * (len(labs) + 1)))
    s0, s1 = win
    for a, lab in zip(axs, labs):
        a.plot(np.arange(len(wf[lab])) * TICK_US, wf[lab], lw=0.8, color="#2a78d6")
        a.axvspan(s0 * TICK_US, s1 * TICK_US, color="#eb6834", alpha=0.15, lw=0)
        a.set_ylabel(lab, fontsize=9)
    z = axs[-1]
    lo, hi = max(0, s0 - 100), min(len(wf[labs[0]]), s1 + 100)
    z.plot(np.arange(lo, hi) * TICK_US, wf[labs[0]][lo:hi], lw=1.0, color="#2a78d6")
    z.set_xlabel("time [µs]")
    z.set_ylabel("zoom")
    axs[0].set_title("run %s subrun %s event %s, channel %s (plane %s); deconvolved only (no raw digits); signal %.0f µs"
                     % (row["run"], row["subrun"], row["event"], row["channel"], row["plane"], (s1 - s0) * TICK_US),
                     loc="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def plot_shape_average(path, plane, shapes):
    """Deconvolved signal / its mean over the signal, aligned at the start, averaged over channels: the
    in-signal sag (a fall with time), comparable between data and MC."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    S = np.array(shapes)
    nn = np.sum(np.isfinite(S), axis=0)
    m = np.where(nn >= 3, np.nanmean(S, axis=0), np.nan)
    e = np.where(nn >= 3, np.nanstd(S, axis=0) / np.sqrt(nn.clip(1)), np.nan)
    t = np.arange(S.shape[1]) * TICK_US
    fig, (a, b) = plt.subplots(2, 1, figsize=(10, 6), sharex=True, gridspec_kw=dict(height_ratios=[3, 1]))
    a.fill_between(t, m - e, m + e, color="#2a78d6", alpha=0.25, lw=0)
    a.plot(t, m, color="#2a78d6", lw=1.2)
    a.axhline(1, color="#8a8984", lw=0.8)
    a.set_ylabel("deconvolved / its mean\nover the signal")
    a.set_title("plane %d: average shape of %d long signals, aligned at their start" % (plane, len(S)), loc="left", fontsize=10)
    b.plot(t, nn, color="#1f1f1e", lw=1.2)
    b.set_ylabel("signals")
    b.set_xlabel("time since the start of the signal [µs]")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
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
    ap.add_argument("-i", "--inputs", nargs="+", required=True,
                    help="reco1 art-ROOT files or xrootd URLs, and / or text files listing them one per line "
                         "(a .txt / .list argument is read as a list: no shell argument-length limit)")
    ap.add_argument("--chunk", type=int, default=1, help="files opened per gallery pass (default 1: per-file progress)")
    ap.add_argument("--progress-file", help=argparse.SUPPRESS)   # set by -j: the worker reports files done there
    ap.add_argument("--channels", help="CSV from select_waveform_channels.py; without it, --auto")
    ap.add_argument("--auto", type=int, default=3, metavar="N",
                    help="without --channels: the N channels per plane with the longest ROI of the first wire label, "
                         "in every event (default 3)")
    ap.add_argument("--min-us", type=float, default=40.0, help="--auto: shortest ROI kept [us] (default 40)")
    ap.add_argument("--max-events", type=int, default=0, help="stop after this many events (0: all)")
    ap.add_argument("--events-per-file", type=int, default=0,
                    help="at most this many events from each file (0: all; opens the files one at a time): "
                         "more files, fewer events each, for the same time")
    ap.add_argument("-j", "--jobs", type=int, default=1,
                    help="split the files over this many worker processes (outputs in <outdir>/part_<k>), then "
                         "merge their averages into <outdir>")
    ap.add_argument("--merge", nargs="+", metavar="DIR",
                    help="only merge the averages (shapes_/tails_plane*.npy) of these output directories into -o")
    ap.add_argument("--list-products", action="store_true",
                    help="print the recob::Wire and raw::RawDigit labels in the first input file and stop")
    ap.add_argument("-o", "--outdir", default="waveforms")
    ap.add_argument("--raw-label", default="simtpc2d:daq")
    ap.add_argument("--wire-labels", nargs="+", default=["simtpc2d:dnnsp"],
                    help="recob::Wire labels (the Gen2 MC reco1 keeps simtpc2d:dnnsp only; see --list-products)")
    ap.add_argument("--max-plots", type=int, default=30, help="single-channel plots written (all go in the npz and the average)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    rows, want = [], {}
    if args.channels:
        rows = list(csv.DictReader(open(args.channels)))
        for r in rows:
            want.setdefault((int(r["run"]), int(r["subrun"]), int(r["event"])), []).append(r)
    if args.merge:
        merge(args.merge, args.outdir)
        return
    files = expand_inputs(args.inputs)
    print("%d input files" % len(files))
    if args.jobs > 1 and len(files) > 1:
        run_parallel(args, files)
        return
    ROOT = load_gallery()
    if args.list_products:
        print("\n".join(list_products(ROOT, files[0])) or "no recob::Wire / raw::RawDigit branches found")
        return
    saved, tails, shapes, nplot, nev = {}, {}, {}, 0, 0
    if args.events_per_file:
        args.chunk = 1
    bar = tqdm(total=len(files), unit="file", desc="reading", disable=bool(args.progress_file))
    for k in range(0, len(files), args.chunk):
        if args.max_events and nev >= args.max_events:
            break
        if k:
            bar.update(args.chunk)
            bar.set_postfix_str("%d events, %d signals" % (nev, sum(len(v) for v in shapes.values())))
            if args.progress_file:
                open(args.progress_file, "w").write("%d\n" % k)
        chunk = files[k:k + args.chunk]
        try:
            ev = ROOT.gallery.Event(ROOT.std.vector("std::string")(chunk))
        except Exception as e:   # one unreadable file: try the chunk's files one by one
            print("chunk %d: %s; opening its files one at a time" % (k // args.chunk, e))
            ev = None
        for evt in ([ev] if ev is not None else []) or [_open_one(ROOT, f) for f in chunk]:
            if evt is None:
                continue
            nev, nplot = loop_events(ROOT, evt, args, want, rows, saved, tails, shapes, nev, nplot)
    bar.update(len(files) - bar.n)
    bar.close()
    if args.progress_file:
        open(args.progress_file, "w").write("%d\n" % len(files))
    _finish(args, rows, saved, tails, shapes, nplot)


def _open_one(ROOT, f):
    try:
        return ROOT.gallery.Event(ROOT.std.vector("std::string")([f]))
    except Exception as e:
        print("skip %s (%s)" % (f, e))
        return None


def expand_inputs(items):
    """Files / URLs, with .txt / .list arguments read as one entry per line (# comments, blank lines skipped)."""
    out = []
    for it in items:
        if it.endswith((".txt", ".list")) and os.path.isfile(it):
            out += [l.strip() for l in open(it) if l.strip() and not l.startswith("#")]
        else:
            out.append(it)
    return list(dict.fromkeys(out))


def loop_events(ROOT, ev, args, want, rows, saved, tails, shapes, nev, nplot):
    here = 0
    while not ev.atEnd() and (not args.max_events or nev < args.max_events) and \
            (not args.events_per_file or here < args.events_per_file):
        nev += 1
        here += 1
        aux = ev.eventAuxiliary()
        key = (int(aux.run()), int(aux.subRun()), int(aux.event()))
        chans = None
        if args.channels and key in want:
            chans = {int(r["channel"]): r for r in want[key]}
        elif not args.channels:
            found = longest_rois(ROOT, ev, args.wire_labels[0], int(args.min_us / TICK_US), args.auto)
            chans = {c: dict(r, run=key[0], subrun=key[1], event=key[2]) for c, r in found.items()}
            rows += list(chans.values())
        if chans:
            data = read_event(ROOT, ev, list(chans), args.raw_label, args.wire_labels)
            for c, wf in data.items():
                r = chans[c]
                lab0 = next((l for l in args.wire_labels if l in wf), None)
                if lab0 is None and "raw" not in wf:
                    continue
                if lab0 is not None:   # the in-signal shape, from the deconvolved signal
                    a0, a1 = signal_window(wf[lab0], (int(r["tick_start"]), int(r["tick_end"])))
                    seg = np.full(SHAPE_TICKS, np.nan)
                    n = min(SHAPE_TICKS, a1 - a0 + 1)
                    mean = wf[lab0][a0:a1 + 1].mean()
                    if mean > 0:
                        seg[:n] = wf[lab0][a0:a0 + n] / mean
                        shapes.setdefault(int(r["plane"]), []).append(seg)
                if "raw" not in wf:   # deconvolved only (raw digits not in the file): no tail
                    name = "wf_%d_%d_%d_ch%d" % (*key, c)
                    win = signal_window(wf[lab0], (int(r["tick_start"]), int(r["tick_end"])))
                    saved[name] = dict(row=r, window=win, **wf)
                    if nplot < args.max_plots:
                        plot_channel_wires(os.path.join(args.outdir, name + ".png"), r, wf, args.wire_labels, win)
                        nplot += 1
                    continue
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
    return nev, nplot


def _finish(args, rows, saved, tails, shapes, nplot):
    np.savez_compressed(os.path.join(args.outdir, "waveforms.npz"),
                        **{k: np.array(v, dtype=object) for k, v in saved.items()})
    for p, tl in sorted(tails.items()):
        plot_tail_average(os.path.join(args.outdir, "tail_average_plane%d.png" % p), p, tl)
        np.save(os.path.join(args.outdir, "tails_plane%d.npy" % p), np.array(tl))
    for p, sh in sorted(shapes.items()):
        plot_shape_average(os.path.join(args.outdir, "shape_average_plane%d.png" % p), p, sh)
        np.save(os.path.join(args.outdir, "shapes_plane%d.npy" % p), np.array(sh))
    print("found %d of %d listed channels; wrote %d plots and the averages to %s"
          % (len(saved), len(rows), nplot, args.outdir))


def _die_with_parent():
    """In a worker (Linux): get SIGTERM when the parent process dies, however it dies (even kill -9)."""
    try:
        import ctypes
        import signal
        ctypes.CDLL("libc.so.6", use_errno=True).prctl(1, signal.SIGTERM)   # PR_SET_PDEATHSIG
    except Exception:
        pass


def run_parallel(args, files):
    """Run this script on args.jobs slices of the files (round robin) in parallel, each into
    <outdir>/part_<k>, then merge the parts' averages into <outdir>. The workers' output is discarded
    (a failing worker's error output is shown); Ctrl-C, kill or the parent dying stops every worker."""
    import signal
    import subprocess
    import sys
    import tempfile
    passthrough, skip = [], False
    drop = {"-i", "--inputs", "-o", "--outdir", "-j", "--jobs"}
    for a in sys.argv[1:]:   # every option except the inputs, outdir and jobs
        if skip:
            if a.startswith("-"):
                skip = False
            else:
                continue
        if a in drop:
            skip = True
            continue
        if any(a.startswith(d + "=") for d in drop):
            continue
        passthrough.append(a)
    procs = []   # (process, part directory, temporary file holding its error output)
    for k in range(args.jobs):
        mine = files[k::args.jobs]
        if not mine:
            continue
        part = os.path.join(args.outdir, "part_%d" % k)
        os.makedirs(part, exist_ok=True)
        lst = os.path.join(part, "inputs.txt")
        open(lst, "w").write("\n".join(mine) + "\n")
        err = tempfile.TemporaryFile(mode="w+")
        cmd = [sys.executable, os.path.abspath(__file__), "-i", lst, "-o", part,
               "--progress-file", os.path.join(part, ".progress")] + passthrough
        procs.append((subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=err, preexec_fn=_die_with_parent), part, err))

    def stop(*_):
        for p, _, _ in procs:
            if p.poll() is None:
                p.terminate()
        deadline = time.time() + 10
        for p, _, _ in procs:
            try:
                p.wait(max(0.1, deadline - time.time()))
            except subprocess.TimeoutExpired:
                p.kill()
        raise SystemExit("stopped: %d workers terminated" % len(procs))
    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGHUP, stop)

    def done():
        n = 0
        for _, part, _ in procs:
            try:
                n += int(open(os.path.join(part, ".progress")).read().strip() or 0)
            except (OSError, ValueError):
                pass
        return n
    bar = tqdm(total=len(files), unit="file", desc="reading (%d workers)" % len(procs))
    try:
        while any(p.poll() is None for p, _, _ in procs):
            bar.update(max(0, done() - bar.n))
            time.sleep(2)
    except KeyboardInterrupt:
        bar.close()
        stop()
    bar.update(max(0, done() - bar.n))
    bar.close()
    for p, part, err in procs:
        if p.returncode != 0:
            err.seek(0)
            tail = err.read().strip().splitlines()[-15:]
            print("worker %s exited with %d:\n  %s" % (part, p.returncode, "\n  ".join(tail)))
        err.close()
        for f in (".progress", "inputs.txt"):
            try:
                os.remove(os.path.join(part, f))
            except OSError:
                pass
    merge([part for _, part, _ in procs], args.outdir)


def merge(dirs, outdir):
    """Combine shapes_plane*.npy and tails_plane*.npy of several output directories and replot the averages."""
    import glob
    os.makedirs(outdir, exist_ok=True)
    for kind, plot in [("shapes", plot_shape_average), ("tails", plot_tail_average)]:
        per = {}
        for d in dirs:
            for f in glob.glob(os.path.join(d, "%s_plane*.npy" % kind)):
                p = int(os.path.basename(f).split("plane")[1].split(".")[0])
                per.setdefault(p, []).append(np.load(f))
        for p, arrs in sorted(per.items()):
            a = np.concatenate(arrs)
            np.save(os.path.join(outdir, "%s_plane%d.npy" % (kind, p)), a)
            plot(os.path.join(outdir, "%s_average_plane%d.png" % ("shape" if kind == "shapes" else "tail", p)), p, list(a))
            print("%s plane %d: %d signals from %d directories" % (kind, p, len(a), len(arrs)))


if __name__ == "__main__":
    main()
