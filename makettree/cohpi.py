import os
import sys
import pandas as pd
import uproot
import pyanalib.pandas_helpers as ph
import awkward as ak
import numpy as np

def make_cohpi_ttree_mc(dfname, split):
    recodf_key = 'cohpi_' + str(split)
    hdrdf_key = 'hdr_' + str(split)
    mcnuwgtdf_key = 'mcnu_' + str(split)

    recodf = pd.read_hdf(dfname, key=recodf_key)
    hdrdf = pd.read_hdf(dfname, key=hdrdf_key)
    mcnuwgtdf = pd.read_hdf(dfname, key=mcnuwgtdf_key)

    ## Collect POT and scale factor to the target POT
    this_pot = sum(hdrdf.pot)
    target_POT = 4.58e18
    POT_scale = target_POT / this_pot

    ## Work for the reco df
    matchdf = recodf.copy()
    matchdf.columns = pd.MultiIndex.from_tuples([(col, '') for col in matchdf.columns])
    matchdf = ph.multicol_merge(matchdf.reset_index(), mcnuwgtdf.reset_index(),
                               left_on=[("__ntuple", ""), ("entry", ""), ("tmatch_idx", "")],
                               right_on=[("__ntuple", ""), ("entry", ""), ("rec.mc.nu..index", "")],
                               how="left") ## -- save all sllices

    wgt_columns = [c for c in list(set(mcnuwgtdf.columns.get_level_values(0)))if (c.startswith("GENIE") or c.startswith("Flux"))]
    #print(wgt_columns)
    recodf_wgt_out = pd.DataFrame({}, index=matchdf.index)
    for col in wgt_columns:
        recodf_wgt_out[col] = np.array([matchdf[col][u].values for u in matchdf[col].columns]).T.tolist()

    recodf = recodf.reset_index()
    recodf = pd.concat([recodf, recodf_wgt_out], axis = 1)
    
    ## Work for the true df
    mcnuwgtdf = mcnuwgtdf[mcnuwgtdf.nuint_categ == 1]
    mcnuwgtdf = mcnuwgtdf.reset_index()

    truedf_wgt_out = pd.DataFrame({}, index=mcnuwgtdf.index)
    for col in wgt_columns:
        truedf_wgt_out[col] = np.array([mcnuwgtdf[col][u].values for u in mcnuwgtdf[col].columns]).T.tolist()

    non_syst_columns = [col for col in mcnuwgtdf.columns if not (col[1].startswith("univ") or col[1].startswith("ms") or col[1].startswith("ps") or col[1].startswith("cv") or col[1].startswith("morph"))]
    #print(list(non_syst_columns))
    truedf_out = mcnuwgtdf[non_syst_columns]
    truedf_out.columns = truedf_out.columns.get_level_values(0)
    truedf_out = pd.concat([truedf_out, truedf_wgt_out], axis = 1)
    
    return recodf, truedf_out

def make_cohpi_ttree_data(dfname, split):
    recodf_key = 'cohpi_' + str(split)
    hdrdf_key = 'hdr_' + str(split)

    recodf = pd.read_hdf(dfname, key=recodf_key)
    hdrdf = pd.read_hdf(dfname, key=hdrdf_key)

    ## Collect POT and scale factor to the target POT
    POT_scale = sum(hdrdf.pot)

    return recodf
