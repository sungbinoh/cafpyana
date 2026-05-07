import os
import sys
import pandas as pd
import uproot
import pyanalib.pandas_helpers as ph
import awkward as ak
import numpy as np


# Add the head directly to sys.path
sys.path.insert(0, "analysis_village/gump")

from gump_cuts import *
from loaddf import *
from rwt_map import *

def apply_truth(recodf, mcnudf):
    matchdf = recodf.copy()
    matchdf.columns = pd.MultiIndex.from_tuples([(col, '') for col in matchdf.columns])
    matchdf = tmatch(matchdf, mcnudf)

    recodf = recodf.reset_index()
    recodf = pd.concat([recodf, matchdf['is_nc'].reset_index(), matchdf['genie_mode'].reset_index()], axis = 1)
    recodf['is_nc'] = recodf['is_nc'].fillna(0)
    recodf['genie_mode'] = recodf['genie_mode'].fillna(0)
    return recodf

def apply_systs(recodf, mcnuwgtdf, DETECTOR, det_run):
    # 1. Setup and Matching
    matchdf = recodf.copy()
    matchdf.columns = pd.MultiIndex.from_tuples([(col, '') for col in matchdf.columns])
    matchdf = tmatch(matchdf, mcnuwgtdf)

    recodf_wgt_out = pd.DataFrame({},index=matchdf.index)

    # 2. Process existing weight columns
    wgt_cols = [c for c in mcnuwgtdf.columns.get_level_values(0).unique()
                if any(k in c for k in ["GENIE", "Flux", "SBNNuSyst", "InterpWeighting"])]

    for col in wgt_cols:
        newcol = f"multisigma{col}" if "MECq0q3InterpWeighting" in col else col
        recodf_wgt_out[newcol] = np.array([matchdf[col][u].values for u in matchdf[col].columns]).T.tolist()

    recodf_wgt_out = apply_det_syst(recodf, recodf_wgt_out, DETECTOR)
    recodf_wgt_out = apply_pot_syst(recodf, recodf_wgt_out, det_run)

    return recodf_wgt_out

def apply_pot_syst(recodf, recodf_wgt_out, det_run):
    # POT Systematics
    pot_scales = {'ps1': 0.5, 'ms1': -0.525805, 'ps2': 1.0, 'ms2': -1.12726, 'ps3': 1.5, 'ms3': -1.7286, 'cv': 0.0}

    for det in ["SBNDRun1", "ICARUSRun2", "ICARUSRun4"]:
        if int(det[-1]) == det_run:
            recodf_wgt_out["POT_multisigma_"+det] = [[1.0 + (v/100.0) for v in pot_scales.values()] for _ in range(len(recodf))]
        else:
            recodf_wgt_out["POT_multisigma_"+det] = [[1.0 for _ in pot_scales.values()] for _ in range(len(recodf))]

    return recodf_wgt_out

def apply_det_syst(recodf, recodf_wgt_out, DETECTOR):
    # 1. Load Mappings
    # These maps are always loaded, but used differently based on DETECTOR
    sce_df = apply_double_map(recodf, 'analysis_village/gump/rwt_outputs/min_SCE.txt', 'analysis_village/gump/rwt_outputs/pls_SCE.txt', 'CAFPYANA_SBN_v1_multisigma_SCE')
    wmxthetaxw_df = apply_map(recodf, 'analysis_village/gump/rwt_outputs/XThetaXW.txt', 'WireMod_SBN_v1_multisigma_XThetaXW')
    wmyz_df = apply_map(recodf, 'analysis_village/gump/rwt_outputs/YZ.txt', 'WireMod_SBN_v1_multisigma_YZ')

    # 2. Apply Detector-Specific Logic
    is_sbnd = (DETECTOR == "SBND")
    not_DETECTOR = "ICARUS" if is_sbnd else "SBND"

    # Common Logic: Use values if SBND, else Use ones
    sys_map = {
        'CAFPYANA_SBN_v1_multisigma_SCE': sce_df,
        'WireMod_SBN_v1_multisigma_XThetaXW': wmxthetaxw_df,
        'WireMod_SBN_v1_multisigma_YZ': wmyz_df
    }

    for key, df in sys_map.items():
        if is_sbnd:
            recodf_wgt_out[key] = np.array([df[key][u].values for u in df[key].columns]).T.tolist()
        else:
            # Create a list of 'ones' with the same shape as the data
            recodf_wgt_out[key] = [np.ones(df[key].shape[1]) for _ in range(len(df))]

    # 3. Handle Smear and Gain (Active for current DETECTOR, 1.0 for the 'other')
    for feat in ['smear', 'gain']:
        col_name = f'CAFPYANA_SBN_v1_multisigma_{DETECTOR}_{feat}'
        other_col = f'CAFPYANA_SBN_v1_multisigma_{not_DETECTOR}_{feat}'

        # Load the map dynamically
        feat_df = apply_double_map(recodf, f'analysis_village/gump/rwt_outputs/min_{DETECTOR}_{feat}.txt', f'analysis_village/gump/rwt_outputs/pls_{DETECTOR}_{feat}.txt', col_name)

        recodf_wgt_out[col_name] = feat_df[col_name].values.tolist()
        recodf_wgt_out[other_col] = [([1.0] * feat_df[col_name].shape[1]) for _ in range(len(feat_df))]

    return recodf_wgt_out

def apply_flash(df, detector, det_run, fname, idf, ismc):
    if "flash_maxpe" in df.columns:
      del df["flash_maxpe"]

    flashes = pd.read_hdf(fname, "flash_%i" % idf)

    time_name = "firsttime" if detector == "SBND" else "time"
    if ismc: # Scale PE for MC-only
        if detector == "SBND": pe_scale = 0.66
        elif detector == "ICARUS" and det_run == 2: pe_scale = 0.6
        elif detector == "ICARUS" and det_run == 4: pe_scale = 0.4
    else:
        pe_scale = 1.0

    intime = (flashes[time_name] > -5) & (flashes[time_name] < 5)
    maxpe = (flashes.totalpe*intime).groupby(level=[0, 1]).max().rename("flash_maxpe")*pe_scale
    df = df.join(maxpe)
    print(df.flash_maxpe)
    return df

def make_gump_ttree_mc(dfname, split):

    # This should be replaced with reading from df later
    if 'ICARUSRun4' in dfname:
        det_run = 4
    if 'ICARUSRun2' in dfname:
        det_run = 2
    elif 'SBND' in dfname:
        det_run = 1

    with pd.HDFStore(dfname, mode='r') as store:
        keys = store.keys()

    if len(keys) > 0:
        if '/evt_0' in keys:
            suffix = '_' + str(split)
        else:
            suffix = ''

        recodf_key = 'evt' + suffix
        hdrdf_key = 'hdr' + suffix
        if '/wgt' + suffix in keys:
            mcnuwgtdf_key = 'wgt' + suffix
        else:
            mcnuwgtdf_key = 'slimwgt' + suffix
        mcnudf_key = 'mcnu' + suffix

        recodf = pd.read_hdf(dfname, key=recodf_key)
        hdrdf = pd.read_hdf(dfname, key=hdrdf_key)
        mcnuwgtdf = pd.read_hdf(dfname, key=mcnuwgtdf_key)
    else:
        recodf = pd.read_hdf(dfname)
        mcnuwgtdf  = recodf.copy()

    mcnudf = pd.read_hdf(dfname, key=mcnudf_key)

    ## Figure out which detector this is
    DETECTOR = recodf.detector.iloc[0]

    recodf = apply_flash(recodf, DETECTOR, det_run, dfname, split, ismc=True)

    ## apply selection
    recodf = all_cuts(recodf, DETECTOR, det_run)

    ## add in systematics
    recodf_wgt_out = apply_systs(recodf, mcnuwgtdf, DETECTOR, det_run)

    ## grab some truth level information
    recodf = apply_truth(recodf, mcnudf)

    ## combine
    recodf_wgt_out = recodf_wgt_out.reset_index()
    recodf = pd.concat([recodf, recodf_wgt_out], axis = 1)

    return recodf

def make_gump_ttree_data(dfname, split):

    # This should be replaced with reading from df later
    if 'ICARUSRun4' in dfname:
        det_run = 4
    if 'ICARUSRun2' in dfname:
        det_run = 2
    elif 'SBND' in dfname:
        det_run = 1
    elif 'ICARUS' in dfname:
        det_run = 2

    with pd.HDFStore(dfname, mode='r') as store:
        keys = store.keys()

    if '/evt_0' in keys:
        suffix = '_' + str(split)
    else:
        suffix = ''

    recodf_key = 'evt' + suffix
    hdrdf_key = 'hdr' + suffix

    recodf = pd.read_hdf(dfname, key=recodf_key)
    hdrdf = pd.read_hdf(dfname, key=hdrdf_key)

    ## Figure out which detector this is
    DETECTOR = recodf.detector.iloc[0]

    recodf = apply_flash(recodf, DETECTOR, det_run, dfname, split, ismc=True)
    ## apply selection
    recodf = all_cuts(recodf, DETECTOR, det_run)

    return recodf
