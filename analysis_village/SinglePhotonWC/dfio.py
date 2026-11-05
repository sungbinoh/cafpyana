"""Read and merge run_df_maker .df files, e.g. from the grid.

A .df file is an HDF5 store holding each dataframe in splits, under keys
"<name>_0", "<name>_1", ..., with the number of splits in key "split". Every
file also has histpotdf / histgenevtdf (TotalPOT, TotalGenEvents per CAF).

Each grid job numbers its own CAFs from 0 in the __ntuple level, so frames
from different jobs collide. load() adds an outer "file" level (position in
the input list) to keep them apart.

Merge a grid output directory into one .df:
  python -m analysis_village.SinglePhotonWC.dfio -i $CAFPYANA_GRID_OUT_DIR/dfs/<dir> -o merged.df
"""

import argparse
import glob
import os
import re
import sys

import pandas as pd

POT_KEY = "histpotdf"
GENEVT_KEY = "histgenevtdf"


def get_n_split(path):
    try:
        return int(pd.read_hdf(path, key="split").n_split.iloc[0])
    except KeyError:
        return 1


def _natural_key(path):
    """Sort x_2.df before x_10.df."""
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", os.path.basename(path))]


def expand_inputs(patterns):
    """Directories, .df files, glob patterns, or .list/.txt files -> list of .df paths."""
    files = []
    for p in patterns:
        if os.path.isdir(p):  # checked first: grid output directories end in .df too
            found = [f for f in glob.glob(os.path.join(p, "*.df")) if os.path.isfile(f)]
            if not found:
                print("Warning: no .df files in directory %s" % p)
            files += sorted(found, key=_natural_key)
        elif p.endswith(".list") or p.endswith(".txt"):
            with open(p) as fl:
                files += [l.strip() for l in fl if l.strip() and not l.startswith("#")]
        else:
            matched = sorted((f for f in glob.glob(p) if os.path.isfile(f)), key=_natural_key)
            files += matched if matched else [p]
    files = list(dict.fromkeys(files))  # drop duplicates, keep order
    missing = [f for f in files if not os.path.exists(f)]
    if missing:
        sys.exit("Input files not found: %s" % missing[:5])
    if not files:
        sys.exit("No input files.")
    return files


def file_keys(path):
    """Dataframe names stored in a file, e.g. ["histgenevtdf", "histpotdf", "wc"]."""
    with pd.HDFStore(path, mode="r") as st:
        names = {re.sub(r"_\d+$", "", k.lstrip("/")) for k in st.keys()}
    return sorted(names - {"split"})


def _read_file(path, keys):
    """{key: [split dfs]} for one file. Raises if unreadable."""
    out = {k: [] for k in keys}
    with pd.HDFStore(path, mode="r") as st:
        stored = set(st.keys())
        for i in range(get_n_split(path)):
            for k in keys:
                key = "/%s_%d" % (k, i)
                if key in stored:  # absent when no rows passed in this split
                    out[k].append(st[key])
    return out


def _tag_file(df, offset):
    """Add the "file" level, or shift it when the input is already merged."""
    if "file" in df.index.names:
        return df.set_index(df.index.set_levels(df.index.levels[df.index.names.index("file")] + offset, level="file"))
    return pd.concat({offset: df}, names=["file"])


def _n_files(got):
    """File ids one input takes up: 1, or all of them for a merged input."""
    ids = [df.index.get_level_values("file").max() for dfs in got.values() for df in dfs
           if "file" in df.index.names and len(df)]
    return int(max(ids)) + 1 if ids else 1


def load(files, keys=("wc",)):
    """Concatenate dataframes over all readable files.

    Returns ({key: DataFrame}, totals). Each frame gets an outer "file" index
    level (position in `files`; merged inputs keep their own file ids,
    shifted to stay unique). totals is a DataFrame indexed by "file" with
    the POT and generated events of each file (use .sum() for the sample).
    A file counts only if every key in it reads cleanly, so the frames and
    totals come from the same files.
    """
    read_keys = list(dict.fromkeys(list(keys) + [POT_KEY, GENEVT_KEY]))
    parts = {k: [] for k in read_keys}
    bad = []
    offset = 0
    for ifile, path in enumerate(files):
        try:
            got = _read_file(path, read_keys)
        except Exception as e:  # truncated/corrupt HDF5 from a crashed job, etc.
            lines = [l.strip() for l in str(e).splitlines() if l.strip()]
            bad.append((path, "%s: %s" % (type(e).__name__, lines[-1] if lines else "")))
            continue
        for k in read_keys:
            parts[k] += [_tag_file(d, offset) for d in got[k]]
        offset += _n_files(got)
        if (ifile + 1) % 50 == 0:
            print("  read %d / %d files" % (ifile + 1, len(files)))

    if bad:
        print("Warning: skipped %d unreadable file(s):" % len(bad))
        for path, err in bad[:10]:
            print("  %s (%s)" % (path, err))
        if len(bad) > 10:
            print("  ... and %d more" % (len(bad) - 10))
    print("Used %d of %d file(s)." % (len(files) - len(bad), len(files)))

    frames = {k: pd.concat(parts[k]) if parts[k] else pd.DataFrame() for k in read_keys}
    for k in keys:
        if not len(frames[k]):
            print("Warning: no '%s' rows in the inputs." % k)

    totals = pd.DataFrame({
        "pot": frames[POT_KEY].TotalPOT.groupby(level="file").sum() if len(frames[POT_KEY]) else pd.Series(dtype=float),
        "ngen": frames[GENEVT_KEY].TotalGenEvents.groupby(level="file").sum() if len(frames[GENEVT_KEY]) else pd.Series(dtype=float),
    }).rename_axis("file")
    return frames, totals


def write(path, frames):
    """Write {key: DataFrame} as a single-split .df (keys <name>_0)."""
    with pd.HDFStore(path, mode="w") as st:
        for k, df in frames.items():
            if len(df):
                st.put("%s_0" % k, df, format="fixed")
        st.put("split", pd.DataFrame({"n_split": [1]}), format="fixed")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--inputs", nargs="+", required=True,
                    help=".df files, directories, glob patterns or .list files")
    ap.add_argument("-o", "--output", required=True, help="merged .df to write")
    ap.add_argument("-k", "--keys", nargs="+", default=None,
                    help="dataframes to merge (default: every one in the first input)")
    args = ap.parse_args()

    files = expand_inputs(args.inputs)
    keys = args.keys or [k for k in file_keys(files[0]) if k not in (POT_KEY, GENEVT_KEY)]
    print("Merging %s from %d file(s)" % (keys, len(files)))
    frames, totals = load(files, keys)
    write(args.output, frames)
    for k, df in frames.items():
        print("  %-14s %d rows" % (k, len(df)))
    print("POT %.4g, generated events %d -> %s" % (totals.pot.sum(), totals.ngen.sum(), args.output))


if __name__ == "__main__":
    main()
