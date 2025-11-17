#import glob
import XRootD.client.glob_funcs as glob
import numpy as np
import uproot
import pandas as pd
from tqdm.auto import tqdm
import subprocess
from multiprocessing import Pool
import multiprocessing
import os
import dill
import sys
from functools import partial
import time
import uuid
import tempfile

from makedf.makedf import make_histpotdf

CPU_COUNT = multiprocessing.cpu_count()

if CPU_COUNT == 0:
    CPU_COUNT = os.cpu_count()

class NTupleProc(object):
    def __init__(self, f=None, name="None"):
        self.name = name
        self.f = f

    def __call__(self, df):
        return self.f(df)

    def __bool__(self):
        return self.f is not None

    # make work with pickle for multiprocessing
    def __getstate__(self):
        return dill.dumps({"f": self.f, "name": self.name})

    def __setstate__(self, state):
        data = dill.loads(state)
        self.f = data["f"]
        self.name = data["name"]

def _open_with_retries(path, attempts=5, sleep=2.0):
    last_exc = None
    for k in range(attempts):
        try:
            return uproot.open(path, timeout=120)
        except (OSError, ValueError) as e:
            last_exc = e
            if k + 1 < attempts:
                time.sleep(sleep * (k + 1))
    raise last_exc

def _loaddf(applyfs, preprocess, g):
    # fname, index, applyfs = inp
    index, fname = g
    # Convert pnfs to xroot URL's
    if fname.startswith("/pnfs"):
        fname = fname.replace("/pnfs", "root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr")
    # fix xroot URL's
    elif fname.startswith("xroot"):
        fname = fname[1:]

    madef = False

    # run any preprocess-ing commands
    tempfiles = []
    if preprocess is not None:
        for i, p in enumerate(preprocess):
            temp_directory = tempfile.gettempdir()
            temp_file_name = os.path.join(temp_directory, "temp%i_%s.flat.caf.root" % (i, str(uuid.uuid4()))) 
            p.run(fname, temp_file_name)
            tempfiles.append(temp_file_name)
            fname = temp_file_name

    try:
        # Open AND close strictly within the context manager
        with _open_with_retries(fname) as f:

            dfs = []
            if "recTree" not in f:
                print("File (%s) missing recTree. Try only histpotdf and skipping other dfs..." % fname)
            else:
                for applyf in applyfs:
                    df = applyf(f)  # must fully read from 'f' here
                    if df is None:
                        dfs.append(None)
                        continue

                    # --- CRITICAL: detach from file-backed/lazy data ---
                    # If it's a pandas obj, deep-copy; if not, try to materialize.
                    if isinstance(df, pd.DataFrame):
                        df = df.copy(deep=True)
                    elif hasattr(df, "to_numpy"):  # Series / array-like
                        df = pd.DataFrame(df.to_numpy()).copy(deep=True)
                    # ---------------------------------------------------

                    # Tag with __ntuple and move it to front of MultiIndex
                    df["__ntuple"] = index
                    df.set_index("__ntuple", append=True, inplace=True)
                    new_order = [df.index.nlevels - 1] + list(range(df.index.nlevels - 1))
                    df = df.reorder_levels(new_order)

                    dfs.append(df)

            if "TotalPOT" in f:
                df = make_histpotdf(f)
                df["__ntuple"] = index
                df.set_index("__ntuple", append=True, inplace=True)
                new_order = [df.index.nlevels - 1] + list(range(df.index.nlevels - 1))
                df = df.reorder_levels(new_order)

                dfs.append(df)
            else:
                print("File (%s) missing TotalPOT histgoram. Cannot make histpotdf..." % fname)

            if not dfs:
                return None

            return dfs
    except (OSError, ValueError) as e:
        print(f"Could not open file ({fname}). Skipping...")
        print(e)
        dfs = None


    if madef:
        os.remove(fname)

    for f in tempfiles:
        os.remove(f)

    return dfs

class NTupleGlob(object):
    def __init__(self, g, branches):
        if isinstance(g, list) and len(g) == 1:
            g = g[0]
        if isinstance(g, list):
            self.glob = g
        elif g.endswith(".list"):
            with open(g) as f:
                self.glob = [line.rstrip("\n") for line in f]
        else:
            self.glob = glob.glob(g, raise_error=True)
        self.branches = branches

    def dataframes(self, fs, maxfile=None, nproc=1, savemeta=False, preprocess=None):
        if not isinstance(fs, list):
            fs = [fs]

        thisglob = self.glob 
        if maxfile:
            thisglob = thisglob[:maxfile]

        if nproc == "auto":
            CPU_COUNT_use = int(CPU_COUNT * 0.8)
            nproc = min(CPU_COUNT_use, len(thisglob))
            print("CPU_COUNT : " + str(CPU_COUNT) + ", len(thisglob): " + str(len(thisglob)) + ", nproc: " + str(nproc))

        ret = []

        try:
            with Pool(processes=nproc) as pool:
                for i, dfs in enumerate(tqdm(pool.imap_unordered(partial(_loaddf, fs, preprocess), enumerate(thisglob)), total=len(thisglob), unit="file", delay=5, smoothing=0.2)):
                    if dfs is not None:
                        ret.append(dfs)
        # Ctrl-C handling
        except KeyboardInterrupt:
            print('Received Ctrl-C. Returning dataframes collected so far.')

        return ret


