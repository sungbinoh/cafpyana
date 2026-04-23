# Standard library imports
import multiprocessing
import os
import random
import subprocess
import sys
import tempfile
import time
import traceback
import uuid
from functools import partial
from multiprocessing import Pool

# Third-party imports
import dill
import numpy as np
import pandas as pd
import uproot
import XRootD.client.glob_funcs as glob
from tqdm.auto import tqdm

# Local application/custom imports
from makedf.makedf import make_histgenevtdf, make_histpotdf

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

def _open_with_retries(path, attempts=3, sleep=5.0):
    last_exc = None

    # --- Try streaming with a "Pre-flight" check ---
    for k in range(attempts):
        try:
            # 1. Open the file handle
            f = uproot.open(path, timeout=300)

            # 2. PRE-FLIGHT: Force a small read (getting keys)
            # If the door is timing out, this will raise the 'Operation expired' error HERE.
            _ = f.keys()

            return f
        except Exception as e:
            last_exc = e
            print(f"[Attempt {k+1}] Network lag for {os.path.basename(path)}. Error: {e}", flush=True)
            time.sleep(sleep * (k + 1) + random.uniform(0, 3))

    # If streaming fails after all attempts, raise the exception
    # so _loaddf can trigger the local copy failover.
    raise last_exc

def _execute_load(f, applyfs, index, fname):
    """Internal helper to process an open file handle."""
    results = []

    # 1. Read TotalEvents immediately to check file validity
    totevt = f['TotalEvents'].values()[0]

    # --- PHASE A: Main Trees ---
    if "recTree" not in f:
        print(f"File ({fname}) missing recTree. Skipping main DFS.", flush=True)
    elif totevt < 1e-6:
        print(f"File ({fname}) has 0 in TotalEvents. Skipping main DFS.", flush=True)
    else:
        for i, applyf in enumerate(applyfs):
            df = applyf(f)
            if df is not None:
                # CRITICAL: Deep copy ensures data is in RAM, not 'linked' to the file
                df = df.copy(deep=True)

                # Metadata tagging
                df["__ntuple"] = index
                df.set_index("__ntuple", append=True, inplace=True)
                new_order = [df.index.nlevels - 1] + list(range(df.index.nlevels - 1))
                df = df.reorder_levels(new_order)
                results.append(df)

    # --- PHASE B: Metadata Histograms ---
    # POT Histogram
    try:
        df_histpot = make_histpotdf(f).copy(deep=True)
        df_histpot["__ntuple"] = index
        df_histpot.set_index("__ntuple", append=True, inplace=True)
        new_order = [df_histpot.index.nlevels - 1] + list(range(df_histpot.index.nlevels - 1))
        results.append(df_histpot.reorder_levels(new_order))
    except Exception as e:
        print(f"Warning: Could not read TotalPOT from {os.path.basename(fname)}: {e}", flush=True)

    # GenEvents Histogram
    try:
        df_histgenevt = make_histgenevtdf(f).copy(deep=True)
        df_histgenevt["__ntuple"] = index
        df_histgenevt.set_index("__ntuple", append=True, inplace=True)
        new_order = [df_histgenevt.index.nlevels - 1] + list(range(df_histgenevt.index.nlevels - 1))
        results.append(df_histgenevt.reorder_levels(new_order))
    except Exception as e:
        print(f"Warning: Could not read TotalGenEvents from {os.path.basename(fname)}: {e}", flush=True)

    return results

def _loaddf(applyfs, preprocess, g):
    index, fname = g
    original_fname = fname

    # Convert pnfs to xroot URL's
    if fname.startswith("/pnfs"):
        fname = fname.replace("/pnfs", "root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr")
    elif fname.startswith("xroot"):
        fname = fname[1:]

    dfs = []
    tempfiles = []
    local_copy_path = None

    # Run any preprocess-ing commands
    if preprocess is not None:
        for i, p in enumerate(preprocess):
            temp_directory = tempfile.gettempdir()
            temp_file_name = os.path.join(temp_directory, f"temp{i}_{uuid.uuid4()}.flat.caf.root")
            p.run(fname, temp_file_name)
            tempfiles.append(temp_file_name)
            fname = temp_file_name

    try:
        # ATTEMPT 1: Normal Streaming
        try:
            with _open_with_retries(fname) as f:
                dfs = _execute_load(f, applyfs, index, fname)
        except Exception as e:
            # If we hit an error (like XRootD Expired or exhausted retries), trigger failover
            if "Operation expired" in str(e) or "vector_read" in str(e) or "I/O operation on closed file" in str(e):
                print(f"Streaming timeout for {os.path.basename(fname)}. Triggering local copy failover...", flush=True)

                # Setup local scratch path using your environment variable!
                scratch_dir = os.getenv("CAFPYANA_TMP_SCRATCH", "/scratch/7DayLifetime/sungbino/tmp_failover")
                os.makedirs(scratch_dir, exist_ok=True)
                local_copy_path = os.path.join(scratch_dir, f"{uuid.uuid4()}_{os.path.basename(original_fname)}")

                # XRDCP the file (Robust)
                subprocess.run(["xrdcp", "-s", "-f", fname, local_copy_path], check=True)

                # ATTEMPT 2: Process the local file
                with uproot.open(local_copy_path) as f_local:
                    dfs = _execute_load(f_local, applyfs, index, local_copy_path)
                    print(f"Successfully recovered {os.path.basename(fname)} via local copy.", flush=True)
            else:
                # If it's a completely different error, don't bother copying; just raise it.
                raise e

    except Exception as e:
        print(f"PERMANENT FAILURE: Could not process file ({fname}). Skipping...", flush=True)
        print(f"Error Detail: {e}", flush=True)
        traceback.print_exc()
        dfs = None

    # --- CLEANUP ---
    # 1. Preprocessed files
    for tmp in tempfiles:
        if os.path.exists(tmp):
            os.remove(tmp)

    # 2. Local failover copy
    if local_copy_path and os.path.exists(local_copy_path):
        try:
            os.remove(local_copy_path)
        except Exception as e:
            print(f"Cleanup error for {local_copy_path}: {e}")

    return dfs if (dfs and len(dfs) > 0) else None

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
