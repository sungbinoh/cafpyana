#!/usr/bin/env python3
"""Pick (run, subrun, event, channel) with long single-wire signals for a
waveform look (dump_waveforms.py): tracks near the drift direction keep one
wire busy for tens to hundreds of us, where an AC-coupling (RC) mismatch
between the hardware and the reconstruction's model shows up.

Runs here (uproot, no LArSoft). Inputs: Wire-Cell flat CAFs (*.flat.caf.root;
the runs come from wc.T_proj_data via make_wc_events.make_wc_chanrun_df), or
dataframes with the wcchanrun and wcevt tables. Keeps the runs of at least
--min-us on one channel, the --per-plane longest per plane, and writes a CSV:
run, subrun, event, channel, plane, tpc, tick_start / tick_end (time slice x
TICKS_PER_SLICE: approximate, the dump finds the signal itself), duration_us,
q (summed charge), on_candidate (the main cluster of a Wire-Cell candidate),
source (the flat CAF or dataframe it came from), reco1 (the reco1 file name
the flat CAF was made from, by the production's pairing rule
reco2-<reco1 stem>[.wc].flat.caf.root; empty for dataframe inputs).
reco1_files.txt (next to the CSV) lists those names once each, e.g. for

  for f in $(cat reco1_files.txt); do samweb -e sbnd locate-file $f; done

  python analysis_village/SinglePhotonWC/select_waveform_channels.py -i <flat CAFs or .df> -o channels.csv \\
      [--min-us 40] [--per-plane 20]
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import uproot

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from analysis_village.SinglePhotonWC import dfio  # noqa: E402
from analysis_village.SinglePhotonWC.make_wc_events import make_wc_chanrun_df  # noqa: E402
from makedf import wirecell as wc  # noqa: E402

TICKS_PER_SLICE = 4     # Wire-Cell imaging time slices (tick_per_slice in wcsimsp_sbnd.fcl)
SLICE_US = 2.0


def from_caf(path):
    f = uproot.open(path)
    runs = make_wc_chanrun_df(f)
    if not len(runs):
        return pd.DataFrame()
    ev = wc.make_wc_evtdf(f)[["run", "subrun", "evt"]]
    runs = runs.join(ev, on="entry")
    return runs.assign(source=os.path.basename(path))


def from_df(path):
    frames, _ = dfio.load([path], keys=["wcchanrun", "wcevt"], nproc=1)
    runs, ev = frames["wcchanrun"], frames["wcevt"]
    if not len(runs):
        return pd.DataFrame()
    key = [l for l in runs.index.names if l != "run"]
    runs = runs.join(ev[["run", "subrun", "evt"]].rename(columns={"run": "run_number"}), on=key)
    return runs.rename(columns={"run_number": "run"}).assign(source=os.path.basename(path))


def reco1_name(caf):
    """reco2-<reco1 stem>[.wc].flat.caf.root -> <reco1 stem>.root (the production's pairing rule)."""
    if not caf.startswith("reco2-") or ".flat.caf" not in caf:
        return ""
    stem = caf[len("reco2-"):].split(".flat.caf")[0]
    return (stem[:-3] if stem.endswith(".wc") else stem) + ".root"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--inputs", nargs="+", required=True, help="flat CAFs (.root) and / or dataframes (.df)")
    ap.add_argument("-o", "--output", default="waveform_channels.csv")
    ap.add_argument("--min-us", type=float, default=40.0, help="shortest single-channel signal kept [us] (default 40)")
    ap.add_argument("--per-plane", type=int, default=20, help="longest runs kept per plane (default 20)")
    ap.add_argument("--candidates-only", action="store_true", help="only the main clusters of Wire-Cell candidates")
    args = ap.parse_args()

    parts = []
    for p in args.inputs:
        try:
            parts.append(from_df(p) if p.endswith(".df") else from_caf(p))
        except Exception as e:
            print("skip %s (%s: %s)" % (p, type(e).__name__, e))
    parts = [d for d in parts if len(d)]
    if not parts:
        sys.exit("no long single-channel signals found")
    d = pd.concat(parts, ignore_index=True)
    d["duration_us"] = d.n_slices * SLICE_US
    d = d[d.duration_us >= args.min_us]
    if args.candidates_only:
        d = d[d.nu_index >= 0]
    d = d.sort_values("duration_us", ascending=False).groupby("plane").head(args.per_plane)
    out = pd.DataFrame({
        "run": d.run.astype(int), "subrun": d.subrun.astype(int), "event": d.evt.astype(int),
        "channel": d.channel.astype(int), "plane": d.plane.astype(int), "tpc": d.tpc.astype(int),
        "tick_start": (d.ts_start * TICKS_PER_SLICE).astype(int),
        "tick_end": ((d.ts_start + d.n_slices) * TICKS_PER_SLICE).astype(int),
        "duration_us": d.duration_us, "q": d.q, "on_candidate": (d.nu_index >= 0).astype(int), "source": d.source,
        "reco1": d.source.map(reco1_name)})
    out = out.sort_values(["plane", "duration_us"], ascending=[True, False])
    out.to_csv(args.output, index=False, float_format="%.6g")
    names = sorted(n for n in out.reco1.unique() if n)
    if names:
        open(os.path.join(os.path.dirname(os.path.abspath(args.output)), "reco1_files.txt"), "w").write("\n".join(names) + "\n")
        print("%d reco1 files hold these channels: reco1_files.txt" % len(names))
    print(out.groupby("plane").agg(n=("channel", "size"), longest_us=("duration_us", "max"),
                                   median_us=("duration_us", "median")).to_string())
    print("wrote %d channels to %s" % (len(out), args.output))


if __name__ == "__main__":
    main()
