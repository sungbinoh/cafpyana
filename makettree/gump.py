import os
import sys
import pandas as pd
import uproot
import pyanalib.pandas_helpers as ph
import awkward as ak
import numpy as np

from analysis_village.gump.gump_cuts import *
from analysis_village.gump.SCE_map import *

def make_gump_ttree_mc(dfname, split):
    recodf_key = 'evt_' + str(split)
    hdrdf_key = 'hdr_' + str(split)
    mcnuwgtdf_key = 'wgt_' + str(split)
    mcnudf_key = 'mcnu_' + str(split)

    recodf = pd.read_hdf(dfname, key=recodf_key)
    hdrdf = pd.read_hdf(dfname, key=hdrdf_key)
    mcnuwgtdf = pd.read_hdf(dfname, key=mcnuwgtdf_key)

    ## Figure out which detector this is
    DETECTOR = recodf.detector.iloc[0]

    ## Apply cuts
    recodf = recodf[slcfv_cut(recodf, DETECTOR)]

    ### NuScore cut
    recodf = recodf[cosmic_cut(recodf)]

    ### Two prong cut
    recodf = recodf[twoprong_cut(recodf)]

    ### containment cut
    recodf = recodf[mufv_cut(recodf, DETECTOR)]
    recodf = recodf[pfv_cut(recodf, DETECTOR)]

    ### PID cut
    recodf = recodf[pid_cut(recodf.mu_chi2_of_mu_cand, recodf.mu_chi2_of_prot_cand, 
                            recodf.prot_chi2_of_mu_cand, recodf.prot_chi2_of_prot_cand, 
                            recodf.mu_len)]

    ### crthitveto cut
    if DETECTOR == "ICARUS":
        recodf = recodf[crthitveto_cut(recodf)]

    ## Work for the reco df
    matchdf = recodf.copy()
    matchdf.columns = pd.MultiIndex.from_tuples([(col, '') for col in matchdf.columns])
    matchdf = ph.multicol_merge(matchdf.reset_index(), mcnuwgtdf.reset_index(),
                               left_on=[("__ntuple", ""), ("entry", ""), ("tmatch_idx", "")],
                               right_on=[("__ntuple", ""), ("entry", ""), ("rec.mc.nu..index", "")],
                               how="left") ## -- save all sllices

    wgt_columns = [c for c in list(set(mcnuwgtdf.columns.get_level_values(0)))if (c.startswith("GENIE") or "Flux" in c)]
    recodf_wgt_out = pd.DataFrame({}, index=matchdf.index)

    for col in wgt_columns:
            recodf_wgt_out[col] = np.array([matchdf[col][u].values for u in matchdf[col].columns]).T.tolist()

    sce_df = apply_sce_map(recodf, 'analysis_village/gump/min_SCE.txt', 'analysis_village/gump/pls_SCE.txt')
    recodf_wgt_out['CAFPYANA_SBN_v1_multisigma_SCE'] = np.array([sce_df['CAFPYANA_SBN_v1_multisigma_SCE'][u].values for u in sce_df['CAFPYANA_SBN_v1_multisigma_SCE'].columns]).T.tolist()

    ## just get NC from here
    mcnudf = pd.read_hdf(dfname, key=mcnudf_key)
    matchdf = recodf.copy()
    matchdf.columns = pd.MultiIndex.from_tuples([(col, '') for col in matchdf.columns])
    matchdf = ph.multicol_merge(matchdf.reset_index(), mcnudf.reset_index(),
                               left_on=[("__ntuple", ""), ("entry", ""), ("tmatch_idx", "")],
                               right_on=[("__ntuple", ""), ("entry", ""), ("rec.mc.nu..index", "")],
                               how="left") ## -- save all sllices

    recodf = recodf.reset_index()
    recodf = pd.concat([recodf, matchdf['is_nc'].reset_index()], axis = 1)
    recodf['is_nc'] = recodf['is_nc'].fillna(0)
    recodf_wgt_out = recodf_wgt_out.reset_index()
    recodf = pd.concat([recodf, recodf_wgt_out], axis = 1)

    recodf = recodf.reset_index()
    recodf = pd.concat([recodf, matchdf['genie_mode'].reset_index()], axis = 1)
    recodf['genie_mode'] = recodf['genie_mode'].fillna(0)
    recodf_wgt_out = recodf_wgt_out.reset_index()
    recodf = pd.concat([recodf, recodf_wgt_out], axis = 1)

    return recodf

def make_gump_ttree_data(dfname, split):
    recodf_key = 'evt_' + str(split)
    hdrdf_key = 'hdr_' + str(split)

    recodf = pd.read_hdf(dfname, key=recodf_key)
    hdrdf = pd.read_hdf(dfname, key=hdrdf_key)

    ## Figure out which detector this is
    DETECTOR = recodf.detector.iloc[0]

    ## Apply cuts
    recodf = recodf[slcfv_cut(recodf, DETECTOR)]

    ### NuScore cut
    recodf = recodf[cosmic_cut(recodf)]

    ### Two prong cut
    recodf = recodf[twoprong_cut(recodf)]

    ### containment cut
    recodf = recodf[mufv_cut(recodf, DETECTOR)]
    recodf = recodf[pfv_cut(recodf, DETECTOR)]

    ### PID cut
    recodf = recodf[pid_cut(recodf.mu_chi2_of_mu_cand, recodf.mu_chi2_of_prot_cand, 
                            recodf.prot_chi2_of_mu_cand, recodf.prot_chi2_of_prot_cand, 
                            recodf.mu_len)]
    ### crthitveto cut
    if DETECTOR == "ICARUS":
        recodf = recodf[crthitveto_cut(recodf)]

    return recodf
