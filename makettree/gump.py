import os
import sys
import pandas as pd
import uproot
import pyanalib.pandas_helpers as ph
import awkward as ak
import numpy as np

from analysis_village.gump.gump_cuts import *
from analysis_village.gump.rwt_map import *

def make_gump_ttree_mc(dfname, split):

    with pd.HDFStore(dfname, mode='r') as store:
        keys = store.keys()
    print(keys)
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

    wgt_columns = [c for c in list(set(mcnuwgtdf.columns.get_level_values(0)))if (c.startswith("GENIE") or "Flux" in c or "SBNNuSyst" in c or "InterpWeighting" in c)]
    recodf_wgt_out = pd.DataFrame({}, index=matchdf.index)

    for col in wgt_columns:
        if "MECq0q3InterpWeighting" in col:
            newcol =  "multisigma" + col
        else:
            newcol = col
        recodf_wgt_out[newcol] = np.array([matchdf[col][u].values for u in matchdf[col].columns]).T.tolist()

    sce_df = apply_double_map(recodf, 'analysis_village/gump/rwt_outputs/min_SCE.txt', 'analysis_village/gump/rwt_outputs/pls_SCE.txt', 'CAFPYANA_SBN_v1_multisigma_SCE')
    wmxthetaxw_df = apply_map(recodf, 'analysis_village/gump/rwt_outputs/XThetaXW.txt', 'WireMod_SBN_v1_multisigma_XThetaXW')
    wmyz_df = apply_map(recodf, 'analysis_village/gump/rwt_outputs/YZ.txt', 'WireMod_SBN_v1_multisigma_YZ')

    smear_df = apply_double_map(recodf, f'analysis_village/gump/rwt_outputs/min_{DETECTOR}_smear.txt', f'analysis_village/gump/rwt_outputs/pls_{DETECTOR}_smear.txt', f'CAFPYANA_SBN_v1_multisigma_{DETECTOR}_smear')
    recodf_wgt_out[f'CAFPYANA_SBN_v1_multisigma_{DETECTOR}_smear'] = np.array([smear_df[f'CAFPYANA_SBN_v1_multisigma_{DETECTOR}_smear'][u].values for u in smear_df[f'CAFPYANA_SBN_v1_multisigma_{DETECTOR}_smear'].columns]).T.tolist()

    gain_df = apply_double_map(recodf, f'analysis_village/gump/rwt_outputs/min_{DETECTOR}_gain.txt', f'analysis_village/gump/rwt_outputs/pls_{DETECTOR}_gain.txt', f'CAFPYANA_SBN_v1_multisigma_{DETECTOR}_gain')
    recodf_wgt_out[f'CAFPYANA_SBN_v1_multisigma_{DETECTOR}_gain'] = np.array([gain_df[f'CAFPYANA_SBN_v1_multisigma_{DETECTOR}_gain'][u].values for u in gain_df[f'CAFPYANA_SBN_v1_multisigma_{DETECTOR}_gain'].columns]).T.tolist()

    if DETECTOR == "ICARUS":
        not_DETECTOR = "SBND"
        recodf_wgt_out['CAFPYANA_SBN_v1_multisigma_SCE'] = np.array([np.ones_like(sce_df['CAFPYANA_SBN_v1_multisigma_SCE'][u].values) for u in sce_df['CAFPYANA_SBN_v1_multisigma_SCE'].columns]).T.tolist()
        recodf_wgt_out['WireMod_SBN_v1_multisigma_XThetaXW'] = np.array([np.ones_like(wmxthetaxw_df['WireMod_SBN_v1_multisigma_XThetaXW'][u].values) for u in wmxthetaxw_df['WireMod_SBN_v1_multisigma_XThetaXW'].columns]).T.tolist()
        recodf_wgt_out['WireMod_SBN_v1_multisigma_YZ'] = np.array([np.ones_like(wmyz_df['WireMod_SBN_v1_multisigma_YZ'][u].values) for u in wmyz_df['WireMod_SBN_v1_multisigma_YZ'].columns]).T.tolist()

    elif DETECTOR == "SBND":
        not_DETECTOR = "ICARUS"
        recodf_wgt_out['CAFPYANA_SBN_v1_multisigma_SCE'] = np.array([sce_df['CAFPYANA_SBN_v1_multisigma_SCE'][u].values for u in sce_df['CAFPYANA_SBN_v1_multisigma_SCE'].columns]).T.tolist()
        recodf_wgt_out['WireMod_SBN_v1_multisigma_XThetaXW'] = np.array([wmxthetaxw_df['WireMod_SBN_v1_multisigma_XThetaXW'][u].values for u in wmxthetaxw_df['WireMod_SBN_v1_multisigma_XThetaXW'].columns]).T.tolist()
        recodf_wgt_out['WireMod_SBN_v1_multisigma_YZ'] = np.array([wmyz_df['WireMod_SBN_v1_multisigma_YZ'][u].values for u in wmyz_df['WireMod_SBN_v1_multisigma_YZ'].columns]).T.tolist()


    recodf_wgt_out[f'CAFPYANA_SBN_v1_multisigma_{not_DETECTOR}_smear'] = np.array([np.ones_like(smear_df[f'CAFPYANA_SBN_v1_multisigma_{DETECTOR}_smear'][u].values) for u in smear_df[f'CAFPYANA_SBN_v1_multisigma_{DETECTOR}_smear'].columns]).T.tolist()
    recodf_wgt_out[f'CAFPYANA_SBN_v1_multisigma_{not_DETECTOR}_gain'] = np.array([np.ones_like(gain_df[f'CAFPYANA_SBN_v1_multisigma_{DETECTOR}_gain'][u].values) for u in gain_df[f'CAFPYANA_SBN_v1_multisigma_{DETECTOR}_gain'].columns]).T.tolist()

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
