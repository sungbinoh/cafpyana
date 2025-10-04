import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from plot_tools import *
# Add the head direcoty to sys.path
workspace_root = os.getcwd()
sys.path.insert(0, workspace_root + "/../../")

# import this repo's classes
import pyanalib.pandas_helpers as ph
from makedf.util import *

import kinematics

def load_data(file, nfiles=1):
    """Load event, header, and mcnu data from HDF file."""

    for s in range(nfiles):
        print("df index:"+str(s))
        df_evt = pd.read_hdf(file, "evt_"+str(s))
        df_hdr = pd.read_hdf(file, "hdr_"+str(s))
        df_mcnu = pd.read_hdf(file, "mcnu_"+str(s))
        df_stub = pd.read_hdf(file, "stub_"+str(s))

        matchdf = df_evt.copy()
        matchdf.columns = pd.MultiIndex.from_tuples([(col, '') for col in matchdf.columns])
        df_evt = ph.multicol_merge(matchdf.reset_index(), df_mcnu.reset_index(),
                                    left_on=[("__ntuple", ""), ("entry", ""), ("tmatch_idx", "")],
                                    right_on=[("__ntuple", ""), ("entry", ""), ("rec.mc.nu..index", "")],
                                    how="left") ## -- save all sllices
        cols_to_drop = []
        for c in df_evt.columns:
            if 'GENIE' in c[0] or 'Flux' in c[0]:
                cols_to_drop.append(c)

        df_evt.drop(columns=cols_to_drop, inplace=True)
        del df_mcnu

        if s == 0:
            res_df_evt = df_evt
            res_df_hdr = df_hdr
            res_df_stub = df_stub
        else:
            res_df_evt = pd.concat([res_df_evt, df_evt])
            res_df_hdr = pd.concat([res_df_hdr, df_hdr])
            res_df_stub = pd.concat([res_df_stub, df_stub])

        del df_evt
        del df_hdr
        del df_stub

    return res_df_evt, res_df_hdr, res_df_stub

def scale_pot(df, df_hdr, desired_pot):
    """Scale DataFrame by desired POT."""
    pot = sum(df_hdr.pot.tolist())
    print(f"POT: {pot}\nScaling to: {desired_pot}")
    scale = desired_pot / pot
    df['glob_scale'] = scale
    return pot, scale
