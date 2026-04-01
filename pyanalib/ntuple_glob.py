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
import random

from makedf.makedf import make_histpotdf
from makedf.makedf import make_histgenevtdf

CPU_COUNT = multiprocessing.cpu_count()

if CPU_COUNT == 0:
    CPU_COUNT = os.cpu_count()

NTUPLE_OPEN_TIMEOUT = int(os.getenv("CAF_NTUPLE_OPEN_TIMEOUT", "300"))
NTUPLE_READ_RETRIES = int(os.getenv("CAF_NTUPLE_READ_RETRIES", "3"))
NTUPLE_RETRY_SLEEP = float(os.getenv("CAF_NTUPLE_RETRY_SLEEP", "4.0"))
NTUPLE_MAX_REMOTE_NPROC = int(os.getenv("CAF_NTUPLE_MAX_REMOTE_NPROC", "8"))
NTUPLE_REDIRECTORS = [r.strip() for r in os.getenv("CAF_XROOTD_REDIRECTORS", "").split(",") if r.strip()]


def _is_remote_path(path):
    return path.startswith("/pnfs") or path.startswith("root://") or path.startswith("xroot")


def _normalize_remote_path(path):
    # Convert pnfs to xroot URL's
    if path.startswith("/pnfs"):
        return path.replace("/pnfs", "root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr")
    # fix xroot URL's
    if path.startswith("xroot"):
        return path[1:]
    return path


def _candidate_paths(path):
    path = _normalize_remote_path(path)
    if not path.startswith("root://"):
        return [path]

    if not NTUPLE_REDIRECTORS:
        return [path]

    url_no_proto = path[len("root://"):]
    host_sep = url_no_proto.find("/")
    if host_sep < 0:
        return [path]

    suffix = url_no_proto[host_sep:]
    candidates = [path]
    for host in NTUPLE_REDIRECTORS:
        if not host:
            continue
        candidates.append("root://%s%s" % (host, suffix))

    # Preserve order and remove duplicates.
    seen = set()
    ordered = []
    for c in candidates:
        if c not in seen:
            seen.add(c)
            ordered.append(c)
    return ordered


def _is_transient_io_error(exc):
    msg = str(exc).lower()
    transient_markers = (
        "operation expired",
        "vector_read",
        "timed out",
        "timeout",
        "pool unavailable",
        "[1010]",
        "[3012]",
        "failed to open file",
        "temporary failure",
        "connection reset",
        "socket",
        "xrootd",
    )
    return any(marker in msg for marker in transient_markers)

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

def _open_with_retries(path, attempts=5, sleep=2.0, timeout=NTUPLE_OPEN_TIMEOUT):
    last_exc = None
    candidates = _candidate_paths(path)
    for k in range(attempts):
        for candidate in candidates:
            try:
                return uproot.open(candidate, timeout=timeout)
            except (OSError, ValueError, RuntimeError, TimeoutError) as e:
                last_exc = e
                continue

        if k + 1 < attempts:
            # Add jitter so many workers do not stampede back at once.
            wait_s = sleep * (k + 1) * (1.0 + 0.25 * random.random())
            time.sleep(wait_s)

    raise last_exc

def _loaddf(applyfs, preprocess, g):
    # fname, index, applyfs = inp
    index, fname = g
    fname = _normalize_remote_path(fname)

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

    dfs = None
    try:
        for attempt in range(max(1, NTUPLE_READ_RETRIES)):
            try:
                # Open AND close strictly within the context manager
                with _open_with_retries(fname, sleep=NTUPLE_RETRY_SLEEP) as f:
                    dfs = []
                    # totevt = f["TotalEvents"].values()[0] if "TotalEvents" in f else 0.0
                    totevt=1.0

                    if "recTree" not in f:
                        print("File (%s) missing recTree. Try only histpotdf & histgenevtdf and skipping other dfs..." % fname)
                    elif totevt < 1e-6:
                        print("File (%s) has 0 in TotalEvents. Try only histpotdf & histgenevtdf and skipping other dfs..." % fname)
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

                    df_histpot = make_histpotdf(f)
                    if "TotalPOT" not in f:
                        print(f"File ({fname}) missing TotalPOT histogram. Using empty DataFrame.")
                    df_histpot["__ntuple"] = index
                    df_histpot.set_index("__ntuple", append=True, inplace=True)
                    new_order = [df_histpot.index.nlevels - 1] + list(range(df_histpot.index.nlevels - 1))
                    df_histpot = df_histpot.reorder_levels(new_order)
                    dfs.append(df_histpot)

                    df_histgenevt = make_histgenevtdf(f)
                    if "TotalGenEvents" not in f:
                        print(f"File ({fname}) missing TotalGenEvents histogram. Using empty DataFrame.")
                    df_histgenevt["__ntuple"] = index
                    df_histgenevt.set_index("__ntuple", append=True, inplace=True)
                    new_order = [df_histgenevt.index.nlevels - 1] + list(range(df_histgenevt.index.nlevels - 1))
                    df_histgenevt = df_histgenevt.reorder_levels(new_order)
                    dfs.append(df_histgenevt)

                break
            except (OSError, ValueError, RuntimeError, TimeoutError, KeyError) as e:
                dfs = None
                if attempt + 1 < NTUPLE_READ_RETRIES and _is_transient_io_error(e):
                    wait_s = NTUPLE_RETRY_SLEEP * (attempt + 1)
                    print(f"Transient read error for ({fname}) [attempt {attempt + 1}/{NTUPLE_READ_RETRIES}]. Retrying in {wait_s:.1f}s...")
                    print(e)
                    time.sleep(wait_s)
                    continue

                print(f"Could not process file ({fname}). Skipping...")
                print(e)
                break
    finally:
        for f in tempfiles:
            if os.path.exists(f):
                os.remove(f)
            
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
            remote_inputs = sum(1 for path in thisglob if _is_remote_path(str(path)))
            if len(thisglob) > 0 and remote_inputs == len(thisglob):
                CPU_COUNT_use = min(CPU_COUNT_use, NTUPLE_MAX_REMOTE_NPROC)
            nproc = max(1, min(CPU_COUNT_use, len(thisglob)))
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


