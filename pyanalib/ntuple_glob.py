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
import traceback
from makedf.makedf import make_histpotdf, make_histgenevtdf

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

def _open_with_retries(path, attempts=10, sleep=2.0, timeout=300):
    last_exc = None
    for k in range(attempts):
        try:
            return uproot.open(path, timeout=timeout)
        except Exception as e:
            last_exc = e
            if k + 1 < attempts:
                print(f"Attempt {k+1}/{attempts} failed for {path}: {type(e).__name__}: {e}", file=sys.stderr, flush=True)
                time.sleep(sleep * (k + 1))
            else:
                print(f"All {attempts} attempts failed for {path}: {type(e).__name__}: {e}", file=sys.stderr, flush=True)
    raise last_exc

def _loaddf(applyfs, preprocess, g):
    # fname, index, applyfs = inp
    index, fname = g
    og_fname = fname
    # Convert pnfs to xroot URL's
    if fname.startswith("/pnfs"):
        fname = fname.replace("/pnfs", "root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr")
    # fix xroot URL's
    elif fname.startswith("xroot"):
        fname = fname[1:]
    # Revert to original filename if it was an XRootD URL
    elif fname.startswith("root://"):
        og_fname = fname.replace("root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr","/pnfs")
    madef = False

    # run any preprocess-ing commands
    tempfiles = []
    if preprocess is not None:
        for i, p in enumerate(preprocess):
            temp_directory = '/exp/sbnd/data/users/brindenc/analyze_sbnd/numu/v10_06_00_validation/pandora/tmp'
            temp_file_name = os.path.join(temp_directory, "temp%i_%s.flat.caf.root" % (i, str(uuid.uuid4()))) 
            p.run(og_fname, temp_file_name)
            print(f"Preprocess command: {p.script} {og_fname} {temp_file_name} completed")
            #After running, make sure it's visible to EAF
            temp_file_name = temp_file_name.replace("/pnfs", "root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr")
            temp_directory = temp_directory.replace("/pnfs", "root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr")
            fname = temp_file_name
            tempfiles.append(temp_file_name)
            print(f"Temp file: {temp_file_name} added to tempfiles list")
            print(f'Files in {temp_directory} after preprocess: {os.listdir(temp_directory)}')
    # Retry the entire file operation (open, read, close)
    attempts = 1
    sleep = 1.0
    timeout = 10
    dfs = None
    
    for k in range(attempts):
        try:
            # Open AND close strictly within the context manager
            with _open_with_retries(fname, attempts=1, sleep=0, timeout=timeout) as f:

                if "recTree" not in f:
                    print("File (%s) missing recTree. Skipping..." % fname)
                    print(f'Keys in file: {f.keys()}')
                    return None

                dfs = []
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
                    dfs.append(None)

                if "TotalGenEvents" in f:
                    df = make_histgenevtdf(f)
                    df["__ntuple"] = index
                    df.set_index("__ntuple", append=True, inplace=True)
                    new_order = [df.index.nlevels - 1] + list(range(df.index.nlevels - 1))
                    df = df.reorder_levels(new_order)

                    dfs.append(df)
                else:
                    print("File (%s) missing TotalGenEvents histgoram. Cannot make histgenevtdf..." % fname)
                    dfs.append(None)
            
            # Success - break out of retry loop
            break
            
        except Exception as e:
            if k + 1 < attempts:
                print(f"Attempt {k+1}/{attempts} failed for {fname}: {type(e).__name__}: {e}", file=sys.stderr, flush=True)
                print(f"Retrying in {sleep * (k + 1):.1f} seconds...", file=sys.stderr, flush=True)
                time.sleep(sleep * (k + 1))
                print(f"Starting attempt {k+2}/{attempts}...", file=sys.stderr, flush=True)
                continue  # Explicitly continue to next iteration
            else:
                print(f"All {attempts} attempts failed for {fname}: {type(e).__name__}: {e}", file=sys.stderr, flush=True)
                print(f"Could not open file ({fname}). Skipping...", file=sys.stderr, flush=True)
                traceback.print_exc(file=sys.stderr)
                dfs = None

    if tempfiles:
        os.system("source /exp/sbnd/app/users/brindenc/develop/cafpyana/pyanalib/remove_temp.sh")
        
    if not dfs:
        return None
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
            # Convert /pnfs to xroot URL before globbing
            if g.startswith("/pnfs"):
                g_xroot = g.replace("/pnfs", "root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr")
                self.glob = glob.glob(g_xroot)
            else:
                self.glob = glob.glob(g)
            if len(self.glob) == 0:
                print(f"Warning: No files matched pattern: {g}")
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
            nproc = 1 if nproc < 1 else nproc
            print("CPU_COUNT : " + str(CPU_COUNT) + ", len(thisglob): " + str(len(thisglob)) + ", nproc: " + str(nproc))

        ret = []

        try:
            with Pool(processes=nproc) as pool:
                for i, dfs in enumerate(tqdm(pool.imap(partial(_loaddf, fs, preprocess), enumerate(thisglob)), total=len(thisglob), unit="file", delay=5, smoothing=0.2)):
                    if dfs is not None:
                        ret.append(dfs)
        # Ctrl-C handling
        except KeyboardInterrupt:
            print('Received Ctrl-C. Returning dataframes collected so far.')

        return ret


