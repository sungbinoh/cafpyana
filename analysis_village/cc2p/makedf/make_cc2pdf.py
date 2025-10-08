from makedf.makedf import *
from pyanalib.pandas_helpers import *
from makedf.util import *

import sys, os
import pyanalib.cc2p_reco_var as cc2preco

## -- data fram maker for cc2p analysis
def make_cc2pdf(f):
    
    pandora_df = make_pandora_df(f)
    
    #### (1) FV cut
    pandora_df = pandora_df[InFV(df = pandora_df.slc.vertex, inzback = 0, det = "SBND")]

    
    #### (2) Not clear cosmic cut
    pandora_df = pandora_df[pandora_df.slc.is_clear_cosmic == 0]
    #pandora_df = pandora_df[pandora_df[('slc', 'is_clear_cosmic', '', '', '', '')] == 0]
    
    #### The three tracks should satisfy chi2 pid cut 
    #### muon candidate: chi2_muon < 25., and chi2_proton > 100.
    #### proton candidates: if not muon candidates and chi2_proton < 90
        
    #### first classify as muon or proton candidate
    pid_result_series = pandora_df.apply(cc2preco.get_pid_result, axis=1)
    pandora_df[('pfp', 'trk', 'reco_pid', '', '', '')] = pid_result_series
    #### then add the muon and proton counters
    pandora_df = cc2preco.get_n_recopid_per_slc(pandora_df)    
    
    #### (3) Track multiplicity cut
    ######## (3) - a: keep only pfp objects with length > 4 cm and dist_to_vertex < 6 cm
    pandora_df = pandora_df[(pandora_df.pfp.trk.len > 15.) & (pandora_df.pfp.dist_to_vertex < 2.)]
    ######## (3) - b: keep only slices with exactly three pfp object passing the requirement
    pandora_df = cc2preco.pass_slc_with_n_pfps(pandora_df)

    #### (4) The three tracks should have track score > 0.5
    #### since they are track-like objects
    pandora_df = pandora_df[pandora_df.pfp.trackScore > 0.5]
    pandora_df = cc2preco.pass_slc_with_n_pfps(pandora_df)

    #### (5) Nuscore > 0.65
    #pandora_df = pandora_df[pandora_df.slc.nu_score > 0.65]
    
    #### (6) demand track containement
    cc2preco.add_contained_col(pandora_df)
    pandora_df = pandora_df[pandora_df.pfp.contained]
    pandora_df = cc2preco.pass_slc_with_n_pfps(pandora_df)
    
    #### (7) demand reco memomenta above threshold
    pandora_df = pandora_df[(pandora_df.pfp.trk.reco_pid == 13) & (pandora_df.pfp.trk.rangeP.p_muon > 0.1) | 
                            (pandora_df.pfp.trk.reco_pid == 2212) & (pandora_df.pfp.trk.rangeP.p_proton > 0.3)]    
    pandora_df = cc2preco.pass_slc_with_n_pfps(pandora_df)
    
    #### dfs for muons and protons
    muons = pandora_df[pandora_df.pfp.trk.reco_pid == 13]
    protons = pandora_df[pandora_df.pfp.trk.reco_pid == 2212]   
    
    # leading and recoil protons   
    lead_protons = protons.sort_values(by=("pfp", "trk", "len", "", "", ""), ascending=False).groupby(level=[0,1]).nth(0)
    rec_protons = protons.sort_values(by=("pfp", "trk", "len", "", "", ""), ascending=False).groupby(level=[0,1]).nth(1)     
    
    #### (9) collect only slc variables of interest
    range_p_mu_series = muons.pfp.trk.rangeP.p_muon
    range_p_mu_series = range_p_mu_series.reset_index(level='rec.slc.reco.pfp..index', drop=True)

    range_p_lead_proton_series = lead_protons.pfp.trk.rangeP.p_proton
    range_p_lead_proton_series = range_p_lead_proton_series.reset_index(level='rec.slc.reco.pfp..index', drop=True)
    
    range_p_rec_proton_series = rec_protons.pfp.trk.rangeP.p_proton
    range_p_rec_proton_series = range_p_rec_proton_series.reset_index(level='rec.slc.reco.pfp..index', drop=True)    

    cos_theta_mu_series = muons.pfp.trk.dir.z
    cos_theta_mu_series = cos_theta_mu_series.reset_index(level='rec.slc.reco.pfp..index', drop=True)

    cos_theta_lead_proton_series = lead_protons.pfp.trk.dir.z
    cos_theta_lead_proton_series = cos_theta_lead_proton_series.reset_index(level='rec.slc.reco.pfp..index', drop=True)
    
    cos_theta_rec_proton_series = rec_protons.pfp.trk.dir.z
    cos_theta_rec_proton_series = cos_theta_rec_proton_series.reset_index(level='rec.slc.reco.pfp..index', drop=True)    
    
    #### add kinematic imbalance vars
    imbalance_df = pandora_df.groupby(['entry', 'rec.slc..index']).apply(cc2preco.measure_reco_imbalance)
    #imbalance_df = pandora_df.groupby(['entry', 'rec.slc..index']).apply(cc2preco.measure_reco_imbalance)
    print(imbalance_df.keys())

    # reco_deltapt_series = imbalance_df['deltapt']
    # reco_deltaalphat_series = imbalance_df['deltaalphat']
    # reco_deltaphit_series = imbalance_df['deltaphit']
    # reco_costheta_lr_series = imbalance_df['costheta_lr']
    # reco_costheta_mu_sum_series = imbalance_df['costheta_mu_sum']
    # reco_e_cal_series = imbalance_df['e_cal']
    # reco_pn_series = imbalance_df['pn']
    # reco_phi_3d_series = imbalance_df['phi_3d']
    # reco_alpha_3d_series = imbalance_df['alpha_3d']  
    
    # ######## (8) - output data frame
    if pandora_df.empty:
        empty_index = pd.MultiIndex(
        levels=[[], []],
        codes=[[], []],
        names=['entry', 'rec.slc..index']
        )
        # reco_deltapt_series = pd.Series(dtype='float', name='reco_deltapt', index=empty_index)
        # reco_deltaalphat_series = pd.Series(dtype='float', name='reco_deltaalphat', index=empty_index)
        # reco_deltaphit_series = pd.Series(dtype='float', name='reco_deltaphit', index=empty_index)
        # reco_costheta_lr_series = pd.Series(dtype='float', name='reco_costheta_lr', index=empty_index)
        # reco_costheta_mu_sum_series = pd.Series(dtype='float', name='reco_costheta_mu_sum', index=empty_index)
        # reco_e_cal_series = pd.Series(dtype='float', name='reco_e_cal', index=empty_index)
        # reco_pn_series = pd.Series(dtype='float', name='reco_pn', index=empty_index)
        # reco_phi_3d_series = pd.Series(dtype='float', name='reco_phi_3d', index=empty_index)
        # reco_alpha_3d_series = pd.Series(dtype='float', name='reco_alpha_3d', index=empty_index)                             

    ######## (8) - c: slc.tmatch.idx for truth matching
    tmatch_idx_series = pandora_df.slc.tmatch.idx
    tmatch_idx_series = tmatch_idx_series.reset_index(level='rec.slc.reco.pfp..index', drop=True)
    
    ## (9) creat a slice-based reco df
    slcdf = pd.DataFrame({
        'range_p_mu': range_p_mu_series,
        'range_p_lead_proton': range_p_lead_proton_series,
        'range_p_rec_proton': range_p_rec_proton_series,        
        'cos_theta_mu': cos_theta_mu_series,
        'cos_theta_lead_proton': cos_theta_lead_proton_series,    
        'cos_theta_rec_proton': cos_theta_rec_proton_series,            
        'tmatch_idx': tmatch_idx_series#,
        # 'reco_deltapt': reco_deltapt_series,
        # 'reco_deltaalphat': reco_deltaalphat_series,
        # 'reco_deltaphit': reco_deltaphit_series,
        # 'reco_costheta_lr': reco_costheta_lr_series,
        # 'reco_costheta_mu_sum': reco_costheta_mu_sum_series,
        # 'reco_e_cal': reco_e_cal_series,
        # 'reco_pn': reco_pn_series,
        # 'reco_phi_3d': reco_phi_3d_series,
        # 'reco_alpha_3d': reco_alpha_3d_series        
    })

    return slcdf
    
    
def make_cc2p_nudf(f):
    
    nudf = make_mcdf(f)
    nudf["ind"] = nudf.index.get_level_values(1)
    wgtdf = geniesyst.geniesyst_sbnd(f, nudf.ind)
    #print("done with genie weights")
    
    is_fv = InFV(df = nudf.position, inzback = 0, det = "SBND")
    is_signal = cc2preco.Signal(nudf)
    is_cc = nudf.iscc
    genie_mode = nudf.genie_mode
    w = nudf.w
    
    #print(type(nudf), nudf.index, len(nudf), nudf.head())

    try :
        nuint_categ = pd.Series(8, index=nudf.index)
        #print(f"done init nuint_categ")
    except Exception as e:
        print(f"Error init nuint_categ")
        return

    nuint_categ[~is_fv] = -1  # Out of FV
    nuint_categ[is_fv & ~is_cc] = 0  # NC
    nuint_categ[is_fv & is_cc & is_signal] = 1  # Signal
    nuint_categ[is_fv & is_cc & ~is_signal & (genie_mode == 3)] = 2  # Non-signal CCCOH
    nuint_categ[is_fv & is_cc & ~is_signal & (genie_mode == 0)] = 3  # Non-signal CCQE
    nuint_categ[is_fv & is_cc & ~is_signal & (genie_mode == 10)] = 4  # Non-signal 2p2h
    nuint_categ[is_fv & is_cc & ~is_signal & (genie_mode == 1)] = 5  # Non-signal RES
    nuint_categ[is_fv & is_cc & ~is_signal & (genie_mode == 2)] = 6  # Non-signal DIS

    #nudf = pd.DataFrame(index=nudf.index)
    #print(nuint_categ)
    nudf['nuint_categ'] = nuint_categ
    
    #true_deltapt_series = nudf.groupby(['entry', 'rec.mc.nu..index']).apply(get_true_deltapt)

    mu_p_series = magdf(nudf.mu.genp)
    p_p_series = magdf(nudf.p.genp)

    mu_cos_theta_series = nudf.mu.genp.z / mu_p_series
    p_cos_theta_series = nudf.p.genp.z / p_p_series

    this_nudf = pd.DataFrame({
        'p_mu': mu_p_series,
        'p_proton': p_p_series,
        'cos_theta_mu': mu_cos_theta_series,
        'cos_theta_p': p_cos_theta_series,
        'nuint_categ': nuint_categ
    })
    this_nudf.columns = pd.MultiIndex.from_tuples([(col, '') for col in this_nudf.columns])
    
    #print(this_nudf.columns)
    #print(wgtdf.columns)
    
    this_nudf = multicol_concat(this_nudf, wgtdf)

    return this_nudf
