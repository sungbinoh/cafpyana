import os
import sys
import pandas as pd
import uproot
import pyanalib.pandas_helpers as ph
import awkward as ak
import numpy as np

from analysis_village.gump.gump_cuts import *

def make_gump_ttree_mc(dfname, split):
    recodf_key = 'evt_' + str(split)
    hdrdf_key = 'hdr_' + str(split)
    mcnuwgtdf_key = 'wgt_' + str(split)
    # mcnuwgtdf_key = 'mcnu_' + str(split)

    recodf = pd.read_hdf(dfname, key=recodf_key)
    hdrdf = pd.read_hdf(dfname, key=hdrdf_key)
    mcnuwgtdf = pd.read_hdf(dfname, key=mcnuwgtdf_key)

    ## Figure out which detector this is
    DETECTOR = recodf.detector.iloc[0]

    ## Work for the reco df
    matchdf = recodf.copy()
    matchdf.columns = pd.MultiIndex.from_tuples([(col, '') for col in matchdf.columns])
    matchdf = ph.multicol_merge(matchdf.reset_index(), mcnuwgtdf.reset_index(),
                               left_on=[("__ntuple", ""), ("entry", ""), ("tmatch_idx", "")],
                               right_on=[("__ntuple", ""), ("entry", ""), ("rec.mc.nu..index", "")],
                               how="left") ## -- save all sllices
    wgt_columns = [c for c in list(set(mcnuwgtdf.columns.get_level_values(0)))if (c.startswith("GENIE") or c.startswith("Flux"))]
    recodf_wgt_out = pd.DataFrame({}, index=matchdf.index)
    for col in wgt_columns:
        recodf_wgt_out[col] = np.array([matchdf[col][u].values for u in matchdf[col].columns]).T.tolist()

    recodf = recodf.reset_index()
    recodf = pd.concat([recodf, recodf_wgt_out], axis = 1)

    ## Apply cuts
    slc_vtx = pd.DataFrame({'x':recodf.slc_vtx_x, 
                            'y':recodf.slc_vtx_y,
                            'z':recodf.slc_vtx_z})

    recodf = recodf[fv_cut(slc_vtx, DETECTOR)]

    ### NuScore cut
    recodf = recodf[cosmic_cut(recodf)]

    ### Two prong cut
    recodf = recodf[twoprong_cut(recodf)]

    ### PID cut
    recodf = recodf[pid_cut(recodf.mu_chi2_of_mu_cand, recodf.mu_chi2_of_prot_cand, 
                            recodf.prot_chi2_of_mu_cand, recodf.prot_chi2_of_prot_cand, 
                            recodf.mu_len)]

    truedf_wgt_out = pd.DataFrame({}, index=mcnuwgtdf.index)
    for col in wgt_columns:
        truedf_wgt_out[col] = np.array([mcnuwgtdf[col][u].values for u in mcnuwgtdf[col].columns]).T.tolist()

    non_syst_columns = [col for col in mcnuwgtdf.columns if not (col[1].startswith("univ") or col[1].startswith("ms") or col[1].startswith("ps") or col[1].startswith("cv") or col[1].startswith("morph"))]
    truedf_out = mcnuwgtdf[non_syst_columns]
    truedf_out.columns = truedf_out.columns.get_level_values(0)
    truedf_out = pd.concat([truedf_out, truedf_wgt_out], axis = 1)
    
    return recodf, truedf_out

def make_gump_ttree_data(dfname, split):
    recodf_key = 'evt_' + str(split)
    hdrdf_key = 'hdr_' + str(split)

    recodf = pd.read_hdf(dfname, key=recodf_key)
    hdrdf = pd.read_hdf(dfname, key=hdrdf_key)

    ## Figure out which detector this is
    DETECTOR = recodf.detector.iloc[0]

    ## Apply cuts
    slc_vtx = pd.DataFrame({'x':recodf.slc_vtx_x, 
                            'y':recodf.slc_vtx_y,
                            'z':recodf.slc_vtx_z})

    recodf = recodf[fv_cut(slc_vtx, DETECTOR)]

    ### NuScore cut
    recodf = recodf[cosmic_cut(recodf)]

    ### Two prong cut
    recodf = recodf[twoprong_cut(recodf)]

    ### PID cut
    recodf = recodf[pid_cut(recodf.mu_chi2_of_mu_cand, recodf.mu_chi2_of_prot_cand, 
                            recodf.prot_chi2_of_mu_cand, recodf.prot_chi2_of_prot_cand, 
                            recodf.mu_len)]

    return recodf
