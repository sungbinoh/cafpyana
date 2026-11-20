"""Random file-level subsample of an EXTERNAL_BEAM on-beam .df, written as an
ordinary on-beam .df carrying its own spill information.

The input (e.g. ICARUS_SpringRun4BNB_ReCAF2026.df, the full ICARUS Run 4
dataset) has no spill information of its own (empty bnb table, NaN evt spill_*
columns, hdr POT 0; see dataquality.EXTERNAL_BEAM). A random fraction of its
CAF files (__ntuple) is kept -- the unit the existing *_unblind samples are cut
in -- and every table is restricted to the kept (__ntuple, entry) records.
The spills come from the beam-quality ROOT file of the input:

  * bnb_<i>: the bnbinfo spills attributed (bnbinfo run, subrun, event) to a
    kept event, hung off the first header record of the event's file, in the
    native bnb layout (spill_time_sec/nsec, event, TOR860, TOR875, LM875A/B/C,
    THCURR, FOM, good_spill). The FOM is the production FOM; dataquality
    replaces it by the corrected FOM_best on load (correct_fom).
  * evt_<i> spill_*: the beam readings of each event's own spill (spill tree).
  * hdr_<i>.pot: on each file's first header record, the sum of the valid
    TOR875 of the file's spills.
  * flashevt_<i> (flash-only productions, one row per event) and trig_<i>:
    restricted to the kept events.

So the POT, the spill-level beam-quality fraction and the gate counts all scale
with the kept files, and dataquality / loaddf treat the output like any other
on-beam file (it must NOT be listed in dataquality.EXTERNAL_BEAM).

Usage:
    python make_unblind_subsample.py --input DF --beamfile ROOT --output DF [--fraction 0.1] [--seed N]
"""
import argparse
import os

import numpy as np
import pandas as pd
import uproot

BNB_COLS = ["spill_time_sec", "spill_time_nsec", "event", "TOR860", "TOR875",
            "LM875A", "LM875B", "LM875C", "THCURR", "FOM", "good_spill"]


def keys(store, prefix):
    return sorted([k[1:] for k in store.keys() if k[1:].startswith(prefix + "_")],
                  key=lambda k: int(k.split("_")[-1]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--beamfile", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--fraction", type=float, default=0.1)
    ap.add_argument("--seed", type=int, default=20261020)
    args = ap.parse_args()
    assert not os.path.exists(args.output), "refusing to overwrite %s" % args.output

    rng = np.random.default_rng(args.seed)
    beam = uproot.open(args.beamfile)
    bi = pd.DataFrame(beam["bnbinfo"].arrays(["run", "subrun", "event", "spill_time_sec", "spill_time_nsec",
                                              "TOR860", "TOR875", "THCURR", "FOM"], library="np"))
    # duplicate spill records (same time stamp) kept once, as dataquality._external_spills
    bi = bi.drop_duplicates(["spill_time_sec", "spill_time_nsec"])
    for c in ("run", "subrun", "event"):
        bi[c] = bi[c].astype(np.int64)
    sp = pd.DataFrame(beam["spill"].arrays(["run", "subrun", "evt", "TOR860", "TOR875", "THCURR", "FOM"],
                                           library="np"))
    for c in ("run", "subrun", "evt"):
        sp[c] = sp[c].astype(np.int64)
    sp = sp.drop_duplicates(["run", "subrun", "evt"]).set_index(["run", "subrun", "evt"])

    with pd.HDFStore(args.input, "r") as inp, pd.HDFStore(args.output, "w") as out:
        tot = dict(files=0, kept=0, ev=0, ev_kept=0, pot_all=0., pot_kept=0., spill_kept=0)
        for hk in keys(inp, "hdr"):
            i = int(hk.split("_")[-1])
            hdr = inp[hk]
            files = np.unique(hdr.index.get_level_values("__ntuple"))
            keep_files = np.sort(rng.choice(files, size=int(round(args.fraction*len(files))), replace=False))
            hkeep = hdr[hdr.index.get_level_values("__ntuple").isin(keep_files)].copy()
            kept_idx = hkeep.index

            # spills of the kept events, hung off each file's first header record
            ekey = hkeep[["run", "subrun", "evt"]].astype(np.int64)
            first = hkeep.reset_index().groupby("__ntuple").entry.min()
            ev2file = pd.Series(hkeep.index.get_level_values("__ntuple"),
                                index=pd.MultiIndex.from_frame(ekey)).groupby(level=[0, 1, 2]).first()
            bkey = pd.MultiIndex.from_frame(bi[["run", "subrun", "event"]])
            m = bkey.isin(ev2file.index)
            b = bi[m].copy()
            b["__ntuple"] = ev2file.reindex(pd.MultiIndex.from_frame(b[["run", "subrun", "event"]])).to_numpy()
            b["entry"] = first.reindex(b["__ntuple"]).to_numpy()
            b = b.sort_values(["__ntuple", "spill_time_sec", "spill_time_nsec"])
            b["rec.hdr.bnbinfo..index"] = b.groupby("__ntuple").cumcount()
            for c in ("LM875A", "LM875B", "LM875C"):
                b[c] = np.nan
            b["good_spill"] = False
            bnb = b.set_index(["__ntuple", "entry", "rec.hdr.bnbinfo..index"])[BNB_COLS]
            bnb.index = bnb.index.set_levels([bnb.index.levels[0].astype(np.int64),
                                              bnb.index.levels[1].astype(np.int64),
                                              bnb.index.levels[2]])

            # hdr POT on each file's first record
            valid = np.isfinite(b.TOR875) & (b.TOR875 > 0)
            fpot = b.TOR875.where(valid, 0.).groupby(b["__ntuple"]).sum()
            hkeep["pot"] = 0.
            fi = pd.MultiIndex.from_arrays([first.index, first.to_numpy()], names=["__ntuple", "entry"])
            hkeep.loc[fi, "pot"] = fpot.reindex(first.index).fillna(0.).to_numpy()

            out.put(hk, hkeep, format="fixed")
            out.put("bnb_%i" % i, bnb, format="fixed")
            for pref in ("evt", "flashevt", "trig", "histgenevtdf", "histpotdf"):
                k = "%s_%i" % (pref, i)
                if "/" + k not in inp.keys():
                    continue
                d = inp[k]
                if pref == "evt":
                    d = d[d.index.droplevel(-1).isin(kept_idx)].copy()
                    eh = hkeep[["run", "subrun", "evt"]].reindex(d.index.droplevel(-1)).astype(np.int64)
                    r = sp.reindex(pd.MultiIndex.from_frame(eh))
                    for c in ("TOR860", "TOR875", "THCURR", "FOM"):
                        d["spill_" + c] = r[c].to_numpy(dtype=np.float64)
                elif pref in ("trig", "flashevt"):
                    d = d[d.index.isin(kept_idx)]
                elif len(d) and "__ntuple" in (d.index.names or []):
                    d = d[d.index.get_level_values("__ntuple").isin(keep_files)]
                out.put(k, d, format="fixed")

            tot["files"] += len(files); tot["kept"] += len(keep_files)
            tot["ev"] += len(hdr); tot["ev_kept"] += len(hkeep)
            allv = np.isfinite(bi.TOR875) & (bi.TOR875 > 0)
            tot["pot_all"] += float(bi.TOR875[allv & bkey.isin(
                pd.MultiIndex.from_frame(hdr[["run", "subrun", "evt"]].astype(np.int64)))].sum())
            tot["pot_kept"] += float(fpot.sum()); tot["spill_kept"] += len(bnb)
        out.put("split", inp["split"], format="fixed")

    print("input  %s" % args.input)
    print("output %s  (fraction %.3f, seed %d)" % (args.output, args.fraction, args.seed))
    print("files  %d -> %d (%.2f%%)" % (tot["files"], tot["kept"], 100*tot["kept"]/tot["files"]))
    print("events %d -> %d (%.2f%%)" % (tot["ev"], tot["ev_kept"], 100*tot["ev_kept"]/tot["ev"]))
    print("beamfile spills: %d (%.4e valid TOR875); attributed to an input event: %.4e POT (%.2f%%)"
          % (len(bi), bi.TOR875[np.isfinite(bi.TOR875) & (bi.TOR875 > 0)].sum(), tot["pot_all"],
             100*tot["pot_all"]/bi.TOR875[np.isfinite(bi.TOR875) & (bi.TOR875 > 0)].sum()))
    print("kept spills %d, POT %.4e (%.2f%% of the attributed POT)"
          % (tot["spill_kept"], tot["pot_kept"], 100*tot["pot_kept"]/tot["pot_all"]))


if __name__ == "__main__":
    main()
