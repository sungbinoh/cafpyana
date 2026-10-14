"""Apply a preselection to a dataframe (.df) file.

Keeps only candidates passing the configured cut in the candidate table, and
only events (identified by the first two index levels: __ntuple, entry) with
at least one passing candidate in the other per-event tables. The index
structure matching tables within a file is preserved. Tables carrying
file-global information (POT, gate counts, normalization histograms) are
copied through unfiltered, as is the "split" metadata key. Tables listed in
DROP are omitted from the output entirely.

Usage:
  python preselct.py -c <config.py> -o <output.df> <input.df>

The config file is python code (exec'd, as in run_df_maker.py) that must define:
  CANDIDATE: name of the candidate table the cut is applied to (e.g. "evt")
  CUT: callable mapping the candidate dataframe to a boolean Series
  KEEP_ALL (optional): list of tables to copy through unfiltered
  DROP (optional): list of tables to drop from the output entirely.
"""
import argparse
import sys

import h5py
import pandas as pd

# Tables that hold file-global information needed for normalization and so
# must never be filtered by an event selection:
#  hdr           -- per-event POT (hdr.pot) and SBND offbeam gates (hdr.noffbeambnb)
#  trig          -- ICARUS offbeam gate counting (trig.gate_delta)
#  bnb           -- per-spill beam info
#  histpotdf     -- TotalPOT histogram
#  histgenevtdf  -- TotalGenEvents histogram
DEFAULT_KEEP_ALL = ["hdr", "trig", "bnb", "histpotdf", "histgenevtdf"]

SPLIT_KEY = "split"


def load_config(path):
    cfg = {"__file__": path}
    exec(open(path).read(), cfg)
    if "CANDIDATE" not in cfg or "CUT" not in cfg:
        raise ValueError("Config %s must define CANDIDATE and CUT" % path)
    return cfg["CANDIDATE"], cfg["CUT"], cfg.get("KEEP_ALL", DEFAULT_KEEP_ALL), cfg.get("DROP", [])


def event_index(df):
    """Reduce a dataframe index to its (__ntuple, entry) event levels."""
    idx = df.index
    if idx.nlevels > 2:
        idx = idx.droplevel(list(range(2, idx.nlevels)))
    return idx


def main(fname, fout, candidate, cut, keep_all, drop):
    with h5py.File(fname, "r") as f:
        keys = [k for k in f.keys() if k != SPLIT_KEY]
        has_split_key = SPLIT_KEY in f.keys()

    tbls = [(k.rsplit("_", 1)[0], int(k.rsplit("_", 1)[1])) for k in keys]
    dropped = sorted(set(t for t, i in tbls if t in drop))
    if dropped:
        print("dropping tables: %s" % ", ".join(dropped))
    tbls = [(t, i) for t, i in tbls if t not in drop]
    isplits = sorted(set(i for t, i in tbls if t == candidate))
    if not isplits:
        raise ValueError("No '%s' tables found in %s" % (candidate, fname))

    orphans = [(t, i) for t, i in tbls if i not in isplits]
    if orphans:
        print("WARNING: keys %s have no matching %s split, copying unfiltered" %
              (["%s_%i" % o for o in orphans], candidate))

    first = True

    def write(df, t, i):
        nonlocal first
        df.to_hdf(fout, key="%s_%i" % (t, i), mode="w" if first else "a")
        first = False

    nevt_in = nevt_out = 0
    for i in isplits:
        cand = pd.read_hdf(fname, "%s_%i" % (candidate, i))
        passed = cut(cand).values.astype(bool)
        cand_out = cand[passed]
        select = event_index(cand_out).unique()

        n_in = event_index(cand).nunique()
        nevt_in += n_in
        nevt_out += len(select)
        print("split %i: selected %i of %i candidate-events (%i of %i candidates)" %
              (i, len(select), n_in, len(cand_out), len(cand)))

        for t, j in tbls:
            if j != i:
                continue
            if t == candidate:
                df = cand_out
            else:
                df = pd.read_hdf(fname, "%s_%i" % (t, i))
                if t not in keep_all:
                    df = df[event_index(df).isin(select)]
            write(df, t, i)

    for t, i in orphans:
        write(pd.read_hdf(fname, "%s_%i" % (t, i)), t, i)

    if has_split_key:
        pd.read_hdf(fname, SPLIT_KEY).to_hdf(fout, key=SPLIT_KEY, mode="a")

    print("total: selected %i of %i candidate-events" % (nevt_out, nevt_in))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-c", dest="config", required=True,
        help="Path to preselection configuration file, i.e.) -c ./analysis_village/gump/configs/preselect_fv.py")
    parser.add_argument("-o", dest="output", required=True, help="Output dataframe file.")
    parser.add_argument("input", help="Input dataframe file.")
    args = parser.parse_args()

    candidate, cut, keep_all, drop = load_config(args.config)
    main(args.input, args.output, candidate, cut, keep_all, drop)
