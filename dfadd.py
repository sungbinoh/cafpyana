#!/usr/bin/env python3
import os, sys, re
import shutil
import argparse
import warnings
import tables
import pandas as pd
from tqdm import tqdm

from pyanalib.split_df_helpers import get_n_split

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=tables.exceptions.NaturalNameWarning)

GB = 1024**3

## Arguments
parser = argparse.ArgumentParser(
    description="hadd for cafpyana dataframe (.df) files: merge the outputs of grid run_df_maker.py jobs into consolidated files.",
    epilog="""\
Examples:

  $ python dfadd.py -o merged.df input_0.df input_1.df input_2.df

  -- Stage output locally before copying to dCache (required for /pnfs destinations,
     since HDF5 cannot write directly to dCache):
  $ python dfadd.py -o /pnfs/sbn/persistent/users/$USER/out.df -tmpdir /exp/sbnd/data/users/$USER/tmp input_*.df

Outputs follow the run_df_maker.py conventions: each dataframe is stored in splits
<name>_<k> kept below the -split size, all dataframes flush to the same split index
together so they stay consistent on the same events, and a "split" key records the
number of splits. If the output exceeds the -filesplit size it is broken into
multiple files named <output>_<i>.df.
""",
    formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument('inputs', nargs='*', help="input .df files to merge")
parser.add_argument('-o', dest='output', default="", help="output data frame file, i.e.) -o merged.df")
parser.add_argument('-l', dest='inputfilelist', default="", help="a file of list for input df files")
parser.add_argument('-split', dest='SplitSize', default=1.0, type=float, help="Max in-memory size in GB of each dataframe split. Default = 1.0 GB.")
parser.add_argument('-filesplit', dest='FileSplitSize', default=30.0, type=float, help="Max output file size in GB before splitting into multiple files. Default = 30.0 GB.")
parser.add_argument('-tmpdir', dest='tmpdir', default="", help="Local staging directory. Each output file is written here, then copied to the -o destination and the staged copy deleted. Required when -o is on dCache (/pnfs).")
parser.add_argument('--no-remap', dest='remap', action='store_false', help="Do not offset the __ntuple index level per input file. By default __ntuple is offset so it stays unique across the merged inputs.")

args = parser.parse_args()

def natural_key(path):
    # sort NAME_2.df before NAME_10.df
    b = os.path.basename(path)
    m = re.match(r'^(.*)_(\d+)\.df$', b)
    if m:
        return (m.group(1), int(m.group(2)))
    return (b, -1)

class OutputWriter(object):
    """Manages the output HDF5 file(s): writes batches of dataframe splits,
    rolls over to a new file when the size limit is hit, and (optionally)
    stages each file locally before copying it to the destination."""
    def __init__(self, output, tmpdir, filesplit_gb):
        output = str(output)
        if not output.endswith(".df"):
            output += ".df"
        self.outdir = os.path.dirname(output) or "."
        self.stem = os.path.basename(output)[:-len(".df")]
        self.tmpdir = tmpdir if tmpdir else None
        self.filesplit_bytes = filesplit_gb * GB
        self.ifile = 0
        self.k_idx = 0
        self.store = None
        self.part_path = None
        self.written_files = []

    def _open(self):
        stagedir = self.tmpdir if self.tmpdir else self.outdir
        self.part_path = os.path.join(stagedir, "%s_part%d.df" % (self.stem, self.ifile))
        if os.path.exists(self.part_path):
            os.remove(self.part_path)
        self.store = pd.HDFStore(self.part_path)

    def write_batch(self, batch, batch_bytes):
        """Write one batch (dict of key name -> concatenated dataframe) as the next
        split index, all keys together. Rolls to a new file first if needed."""
        if self.store is None:
            self._open()
        elif self.k_idx > 0 and os.path.getsize(self.part_path) + batch_bytes > self.filesplit_bytes:
            self._finalize(last=False)
            self.ifile += 1
            self.k_idx = 0
            self._open()

        for k in sorted(batch):
            this_key = k + "_" + str(self.k_idx)
            try:
                self.store.put(key=this_key, value=batch[k], format="fixed")
                print(f"Saved {this_key}: {batch[k].memory_usage(deep=True).sum() / GB:.4f} GB")
            except Exception as e:
                print(f"Table {this_key} failed to save, skipping. Exception: {str(e)}")
        self.store.flush()
        self.k_idx += 1

    def _finalize(self, last):
        # Save the split count metadata
        split_df = pd.DataFrame({"n_split": [self.k_idx]})
        self.store.put(key="split", value=split_df, format="fixed")
        self.store.close()
        self.store = None

        # Only file ever produced -> no _<i> suffix
        if last and self.ifile == 0:
            final_name = "%s.df" % self.stem
        else:
            final_name = "%s_%d.df" % (self.stem, self.ifile)
        final_path = os.path.join(self.outdir, final_name)

        if self.tmpdir:
            print(f"Copying {self.part_path} -> {final_path}")
            shutil.copyfile(self.part_path, final_path)
            if os.path.getsize(final_path) != os.path.getsize(self.part_path):
                raise IOError("Size mismatch copying %s to %s!" % (self.part_path, final_path))
            os.remove(self.part_path)
        else:
            os.replace(self.part_path, final_path)

        print(f"Finalized {final_path}: {self.k_idx} splits, {os.path.getsize(final_path) / GB:.2f} GB")
        self.written_files.append(final_path)

    def close(self):
        if self.store is None:
            self._open()
        self._finalize(last=True)

def main():
    inputs = list(args.inputs)
    if args.inputfilelist != "":
        with open(args.inputfilelist) as f:
            inputs += [line.strip() for line in f if line.strip() and not line.startswith("#")]

    if not inputs or args.output == "":
        parser.print_help()
        print(parser.epilog)
        sys.exit(1)

    inputs.sort(key=natural_key)

    split_bytes = args.SplitSize * GB
    writer = OutputWriter(args.output, args.tmpdir, args.FileSplitSize)

    buffers = {}  # key name -> list of dataframes
    sizes = {}    # key name -> accumulated in-memory bytes
    ntuple_offset = 0
    skipped = []

    def flush():
        if not any(buffers.values()):
            return
        batch = {}
        batch_bytes = 0
        for k, buffer in buffers.items():
            if buffer:
                batch[k] = pd.concat(buffer, ignore_index=False)
                batch_bytes += sizes[k]
        writer.write_batch(batch, batch_bytes)
        del batch
        for k in buffers:
            buffers[k] = []
            sizes[k] = 0

    for fname in tqdm(inputs, unit="file"):
        try:
            n_split = int(get_n_split(fname))
            with pd.HDFStore(fname, mode='r') as store:
                keys = [k.lstrip('/') for k in store.keys()]

            # map split index -> key names present at that split
            key_map = {}
            for k in keys:
                if k == "split":
                    continue
                m = re.match(r'^(.*)_(\d+)$', k)
                if m is None:
                    print(f"Key {k} in {fname} does not follow the <name>_<split> convention, skipping.")
                    continue
                key_map.setdefault(int(m.group(2)), []).append(m.group(1))

            file_max_ntuple = -1
            for i in range(n_split):
                for base in sorted(key_map.get(i, [])):
                    df = pd.read_hdf(fname, key=f"{base}_{i}")

                    if df.index.names[0] == "__ntuple":
                        if args.remap and ntuple_offset > 0:
                            if isinstance(df.index, pd.MultiIndex):
                                df.index = df.index.set_levels(df.index.levels[0] + ntuple_offset, level=0)
                            else:
                                df.index = df.index + ntuple_offset
                        if len(df):
                            file_max_ntuple = max(file_max_ntuple, int(df.index.get_level_values(0).max()))
                    else:
                        print(f"Key {base}_{i} in {fname} does not have __ntuple as the first index level.")

                    buffers.setdefault(base, []).append(df)
                    sizes[base] = sizes.get(base, 0) + df.memory_usage(deep=True).sum()

                if any(v > split_bytes for v in sizes.values()):
                    flush()

            ntuple_offset = max(ntuple_offset, file_max_ntuple + 1)

        except (Exception,) as e:
            print(f"Could not process file ({fname}). Skipping. Exception: {str(e)}")
            skipped.append(fname)

    flush()
    writer.close()

    print(f"\nMerged {len(inputs) - len(skipped)}/{len(inputs)} input files into:")
    for f in writer.written_files:
        print(f"  {f}")
    if skipped:
        print(f"WARNING: skipped {len(skipped)} input files:")
        for f in skipped:
            print(f"  {f}")

if __name__ == "__main__":
    main()
