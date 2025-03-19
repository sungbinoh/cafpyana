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

def _loaddf(hdfout, applyfs, g):
    # fname, index, applyfs = inp
    index, fname = g
    #print("index: %d, fname: %s" %(index, fname))
    # Convert pnfs to xroot URL's
    if fname.startswith("/pnfs"):
        fname = fname.replace("/pnfs", "root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr")
    # fix xroot URL's
    elif fname.startswith("xroot"):
        fname = fname[1:]

    madef = False

    # TODO: make available?
    # 
    # Flatten non-flat cafs
    # if "flat" not in fname.split("/")[-1].split("."):
    #     flatcaf = "/tmp/" + fname.split("/")[-1].split(".")[0] + ".flat.root"
    #     subprocess.run(["flatten_caf", fname, flatcaf],  stdout=subprocess.DEVNULL)
    #     fname = flatcaf
    #     madef = True
   
    try:
        f = uproot.open(fname, timeout=120)
    except (OSError, ValueError) as e:
        print("Could not open file (%s) due to exception. Skipping..." % fname) 
        print(e)
        return None
    with f:
        try:
            dfs = [applyf(f) for applyf in applyfs]
        except Exception as e:
            if True:
                raise
            print("Error processing file (%s). Skipping..." % fname)
            print(e)
            return None

        # Set an index on the NTuple number to make sure we keep track of what is where
        for i in range(len(dfs)):
            if dfs[i] is not None:
                dfs[i]["__ntuple"] = index
                dfs[i].set_index("__ntuple", append=True, inplace=True)
                dfs[i] = dfs[i].reorder_levels([dfs[i].index.nlevels-1] + list(range(0, dfs[i].index.nlevels-1)))
            else:
                dfs[i] = []

    # Save DataFrame to HDF5 file
    if dfs:
        #with pd.HDFStore(hdfout, mode="a", complevel=5, complib="blosc") as store:
        with pd.HDFStore(hdfout) as store:
            for i, df in enumerate(dfs):
                if not df.empty:
                    store.put(f"tempdf_{index}_{i}", df, format="fixed")
    if madef:
        os.remove(fname)

    return [f"tempdf_{index}_{i}" for i in range(len(dfs))]  # Return dataset keys

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
            self.glob = glob.glob(g)
        self.branches = branches

    def dataframes(self, fs, maxfile=None, nproc=1, savemeta=False):
        if not isinstance(fs, list):
            fs = [fs]

        thisglob = self.glob 
        if maxfile:
            thisglob = thisglob[:maxfile]

        if nproc == "auto":
            print("CPU_COUNT : " + str(CPU_COUNT) + ", len(thisglob): " + str(len(thisglob)))
            nproc = min(CPU_COUNT, len(thisglob))

        saved_keys = []  # Store dataset keys
        hdf5_file = os.path.join("output_data", "data.h5")
        os.makedirs(os.path.dirname(hdf5_file), exist_ok=True)

        try:
            with Pool(processes=nproc) as pool:
                for keys in tqdm(
                    pool.imap_unordered(partial(_loaddf, hdf5_file, fs), enumerate(thisglob)),
                    total=len(thisglob),
                    unit="file",
                    delay=5,
                    smoothing=0.2
                ):
                    if keys is not None:
                        saved_keys.extend(keys)

        # Ctrl-C handling
        except KeyboardInterrupt:
            print('Received Ctrl-C. Returning dataframes collected so far.')
            return hdf5_file, saved_keys

        return hdf5_file, saved_keys
