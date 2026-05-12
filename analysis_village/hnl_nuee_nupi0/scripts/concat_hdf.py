#!/usr/bin/env python3
"""Concatenate split HDF5 files in a directory into a single output file."""
import argparse
import glob
import os
import re
import sys
import warnings

import math

import pandas as pd
import tables
from tqdm import tqdm

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=tables.exceptions.NaturalNameWarning)


def get_n_split(file):
    return pd.read_hdf(file, key="split").n_split.iloc[0]


def load_dfs(file, keys2load, n_max_concat=None, start_split=0):
    n_splits = get_n_split(file) - start_split
    n_concat = n_splits if n_max_concat is None else min(n_max_concat, n_splits)
    out = {}
    for key in keys2load:
        dfs = [pd.read_hdf(file, key=f"{key}_{i}") for i in range(start_split, start_split + n_concat)]
        out[key] = pd.concat(dfs, ignore_index=False)
    return out


def get_keys(file):
    """Return unique base key names from an HDF5 file, excluding 'split'."""
    with pd.HDFStore(file, mode='r') as store:
        raw_keys = store.keys()
    base_keys = {re.sub(r'_\d+$', '', k.lstrip('/')) for k in raw_keys}
    base_keys.discard('split')
    return sorted(base_keys)


def concat_hdf_files(directory, keys2load, output_dir=None, pattern="*.df", n_max_concat=None, split_size=1.0):
    """Concatenate DataFrames across multiple HDF5 files and write to a new file.

    Output preserves the split-file structure (keys stored as {key}_N, plus a
    split key with n_split) so it is readable by load_dfs. Each row gets a
    file_idx index level indicating which source file it came from.
    """
    all_files = sorted(glob.glob(os.path.join(directory, pattern)))
    if not all_files:
        raise FileNotFoundError(f"No files matching '{pattern}' found in {directory}")

    first_basename = os.path.basename(all_files[0])
    out_basename = re.sub(r'_\d+(\.[^.]+)$', r'\1', first_basename)
    out_dir = output_dir if output_dir is not None else directory
    out_path = os.path.join(out_dir, out_basename)

    # exclude the output file in case it already exists in the same directory
    files = [f for f in all_files if os.path.abspath(f) != os.path.abspath(out_path)]
    if not files:
        raise FileNotFoundError(f"No source files found after excluding output file {out_path}")

    merged = {key: [] for key in keys2load}

    for file_idx, file in enumerate(tqdm(files, desc="Concatenating files")):
        n_splits = get_n_split(file)
        n_load = n_splits if n_max_concat is None else min(n_max_concat, n_splits)
        file_dfs = load_dfs(file, keys2load, n_max_concat=n_load)
        for key in keys2load:
            df = file_dfs[key].copy()
            df["file_idx"] = file_idx
            df = df.set_index("file_idx", append=True)
            merged[key].append(df)

    # concatenate all data and determine number of splits from total memory
    full = {key: pd.concat(dfs, ignore_index=False) for key, dfs in merged.items()}
    total_gb = sum(df.memory_usage(deep=True).sum() for df in full.values()) / (1024 ** 3)
    n_splits = max(1, math.ceil(total_gb / split_size))

    with pd.HDFStore(out_path, mode='w') as store:
        for key in keys2load:
            n_rows = len(full[key])
            chunk_size = math.ceil(n_rows / n_splits)
            for i in range(n_splits):
                chunk = full[key].iloc[i * chunk_size:(i + 1) * chunk_size]
                store.put(f"{key}_{i}", chunk, format="fixed")
        store.put("split", pd.DataFrame({"n_split": [n_splits]}), format="fixed")

    print(f"Wrote {out_path} ({n_splits} split(s))")
    return out_path


def main():
    parser = argparse.ArgumentParser(
        description="Concatenate split HDF5 files and write a merged output file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "directory",
        help="Directory containing the HDF5 files to concatenate.",
    )
    parser.add_argument(
        "keys",
        nargs="*",
        help="HDF5 keys to load and concatenate (e.g. nu slc). Defaults to all keys.",
    )
    parser.add_argument(
        "--output-dir", "-o",
        default=None,
        help="Directory to write the output file. Defaults to the input directory.",
    )
    parser.add_argument(
        "--pattern", "-p",
        default="*.df",
        help="Glob pattern for matching input files.",
    )
    parser.add_argument(
        "--n-max-concat", "-n",
        type=int,
        default=None,
        help="Maximum number of splits to load per file. Loads all splits by default.",
    )
    parser.add_argument(
        "--split", "-s",
        type=float,
        default=1.0,
        help="Target split size in GB for the output file.",
    )

    args = parser.parse_args()

    keys2load = args.keys
    if not keys2load:
        files = sorted(glob.glob(os.path.join(args.directory, args.pattern)))
        if not files:
            print(f"Error: no files matching '{args.pattern}' in {args.directory}", file=sys.stderr)
            sys.exit(1)
        keys2load = get_keys(files[0])
        print(f"Auto-detected keys: {keys2load}")

    concat_hdf_files(
        directory=args.directory,
        keys2load=keys2load,
        output_dir=args.output_dir,
        pattern=args.pattern,
        n_max_concat=args.n_max_concat,
        split_size=args.split,
    )


if __name__ == "__main__":
    main()
