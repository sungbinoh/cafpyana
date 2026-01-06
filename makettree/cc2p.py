import os
import sys
import pandas as pd
import uproot
import pyanalib.pandas_helpers as ph
import awkward as ak
import numpy as np

def make_cc2p_ttree_mc(dfname):
    
    recodf = pd.read_hdf(dfname, key='cc2p')
    hdrdf = pd.read_hdf(dfname, key='hdr')
    mcnuwgtdf = pd.read_hdf(dfname, key='mcnu')

    ## Collect POT and scale factor to the target POT
    target_POT = sum(hdrdf.pot)

    ## Work for the reco df
    matchdf = recodf.copy()
    matchdf.columns = pd.MultiIndex.from_tuples([(col, '') for col in matchdf.columns])
    matchdf = ph.multicol_merge(matchdf.reset_index(), mcnuwgtdf.reset_index(),
                               left_on=[("__ntuple", ""), ("entry", ""), ("tmatch_idx", "")],
                               right_on=[("__ntuple", ""), ("entry", ""), ("rec.mc.nu..index", "")],
                               how="left") ## -- save all slices

    wgt_columns = [c for c in list(set(mcnuwgtdf.columns.get_level_values(0))) if c.startswith("GENIEReWeight")]
    
    recodf_wgt_out = pd.DataFrame({}, index=matchdf.index)
    for col in wgt_columns:
        recodf_wgt_out[col] = np.array([matchdf[col][u].values for u in matchdf[col].columns]).T.tolist()

    recodf = recodf.reset_index()
    recodf = pd.concat([recodf, recodf_wgt_out], axis = 1)
    
    ## Work for the true df
    print(mcnuwgtdf.keys())
    print("mcnuwgtdf.nuint_categ")
    print(mcnuwgtdf.nuint_categ.value_counts())
    mcnuwgtdf = mcnuwgtdf[mcnuwgtdf.nuint_categ == 1]
    mcnuwgtdf = mcnuwgtdf.reset_index()

    print("mcnuwgtdf")
    print(mcnuwgtdf)
    truedf_wgt_out = pd.DataFrame({}, index=mcnuwgtdf.index)
    for col in wgt_columns:
        print(col)
        print(np.array([mcnuwgtdf[col][u].values for u in mcnuwgtdf[col].columns]).T.tolist())
        truedf_wgt_out[col] = np.array([mcnuwgtdf[col][u].values for u in mcnuwgtdf[col].columns]).T.tolist()

    print("truedf_wgt_out")
    print(truedf_wgt_out)

    non_syst_columns = [col for col in mcnuwgtdf.columns if not col[0].startswith("GENIEReWeight")]
    print(non_syst_columns)

    truedf_out = mcnuwgtdf[non_syst_columns]
    truedf_out.columns = truedf_out.columns.get_level_values(0)
    truedf_out = pd.concat([truedf_out, truedf_wgt_out], axis = 1)
    
    return recodf, truedf_out


def make_cc2p_ttree_data(dfname):
    
    recodf = pd.read_hdf(dfname, key='cc2p')
    hdrdf = pd.read_hdf(dfname, key='hdr')

    ## Collect POT and scale factor to the target POT
    target_POT = sum(hdrdf.pot)

    recodf = recodf.reset_index()
    return recodf