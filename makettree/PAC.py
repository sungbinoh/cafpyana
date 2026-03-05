import os
import sys
import pandas as pd
import uproot
import pyanalib.pandas_helpers as ph
import awkward as ak
import numpy as np

def make_ttree_mc(dfname, split):

    with pd.HDFStore(dfname, mode='r') as store:
        keys = store.keys()

    print(keys)

    if len(keys) > 0:
        if '/mcnu_0' in keys:
            suffix = '_' + str(split)
        else:
            suffix = ''

        mcnudf_key = 'mcnu' + suffix
        mcnudf = pd.read_hdf(dfname, key=mcnudf_key)

    else:
        mcnudf = pd.read_hdf(dfname)

    #wgt_columns = [c for c in list(set(mcnudf.columns.get_level_values(0)))if (c.startswith("GENIEReWeight_SBN_v1_multisigma_EtaNCEL"))]
    wgt_columns = [c for c in list(set(mcnudf.columns.get_level_values(0)))if (c.startswith("GENIE") or "Flux" in c or "SBNNuSyst" in c or "InterpWeighting" in c)]

    mcnudf_ret = mcnudf.copy()

    to_drop = []
    for col in wgt_columns:
        print(col)
        mcnudf_ret.drop(columns=[col], inplace=True)
        #if "horncurrent_Flux" in col:
        #    print(mcnudf[col])
        #    sys.exit()
        if "MECq0q3InterpWeighting" in col:
            newcol =  "multisigma" + col
        else:
            newcol = col
            if len(mcnudf[col].columns) == 1:
                print(f"Will drop {col}!")
                to_drop.append(col)

        if len(mcnudf[col].columns) > 100:
            keep_univs = [f'univ_{i}' for i in range(100)]
            mcnudf_ret[newcol] = np.array([mcnudf[col][u].values for u in keep_univs]).T.tolist()
        else:
            mcnudf_ret[newcol] = np.array([ mcnudf[col][u].values for u in mcnudf[col].columns]).T.tolist()


    ## just get NC from here
    #mcnudf = pd.read_hdf(dfname, key=mcnudf_key)
    #mcnudf_ret.drop(columns=["GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi", "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi", "GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad"])
    mcnudf_ret.drop(columns=to_drop)
    mcnudf_ret.columns = mcnudf_ret.columns.get_level_values(0)
    return mcnudf_ret.reset_index()
