"""Read NCDelta .df files (cafpyana run_df_maker output). numpy/pandas only.

A .df file is an HDF5 store holding each dataframe in splits, under keys
"<name>_0", "<name>_1", ..., with the number of splits in key "split".
"""

import glob
import os
import re
import sys

import pandas as pd


def get_n_split(path):
    """Number of splits in a file (same as pyanalib.split_df_helpers.get_n_split)."""
    try:
        return int(pd.read_hdf(path, key="split").n_split.iloc[0])
    except KeyError:
        return 1


def _natural_key(path):
    """Sort ncdelta_2.df before ncdelta_10.df."""
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


def _read_file(path, keys):
    """{key: [split dfs]} and total hdr POT for one file. Raises if unreadable."""
    out = {k: [] for k in keys}
    pot = 0.
    for i in range(get_n_split(path)):
        for k in keys:
            try:
                df = pd.read_hdf(path, key="%s_%d" % (k, i))
            except KeyError:
                df = None  # nothing stored for this split (e.g. no slices passed)
            if df is not None and len(df):
                out[k].append(df)
        try:
            pot += pd.read_hdf(path, key="hdr_%d" % i).pot.sum()
        except KeyError:
            print("Warning: %s has no hdr_%d; its POT is not counted." % (path, i))
    return out, pot


def load(files, keys=("slc",), require_truth=True):
    """Concatenate dataframes over all readable files.

    Returns ({key: DataFrame}, POT per file). The POT is a Series indexed by
    the "file" level (use .sum() for the total). A file counts only if every
    key in it reads cleanly, so all keys and the POT come from the same files.
    Each frame gets an outer "file" index level (position in `files`), which
    keeps indices unique and can be used to split samples, e.g. for training
    and validation.
    """
    parts = {k: [] for k in keys}
    pot, bad = {}, []
    for ifile, path in enumerate(files):
        try:
            got, file_pot = _read_file(path, keys)
        except Exception as e:  # truncated/corrupt HDF5 from a crashed job, etc.
            lines = [l.strip() for l in str(e).splitlines() if l.strip()]
            bad.append((path, "%s: %s" % (type(e).__name__, lines[-1] if lines else "")))
            continue
        for k in keys:
            parts[k] += [pd.concat({ifile: d}, names=["file"]) for d in got[k]]
        pot[ifile] = file_pot
        if (ifile + 1) % 50 == 0:
            print("  read %d / %d files" % (ifile + 1, len(files)))

    if bad:
        print("Warning: skipped %d unreadable file(s):" % len(bad))
        for path, err in bad[:10]:
            print("  %s (%s)" % (path, err))
        if len(bad) > 10:
            print("  ... and %d more" % (len(bad) - 10))
    print("Used %d of %d file(s)." % (len(files) - len(bad), len(files)))

    frames = {}
    for k in keys:
        if not parts[k]:
            sys.exit("No '%s' rows found in the inputs." % k)
        frames[k] = pd.concat(parts[k])
    if require_truth and "slc" in frames and "true_pdg" not in frames["slc"].columns:
        sys.exit("No truth columns: these look like data dataframes. Use the MC ones.")
    return frames, pd.Series(pot, name="pot", dtype=float).rename_axis("file")
