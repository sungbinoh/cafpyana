import pandas as pd
import numpy as np
from data.dictionary.chi2pid.chi2pid_template import *
## Functions for updating dE/dx and chi2_pid

## LAr properties
temp = 88.4 ## Kelvin, https://github.com/SBNSoftware/sbndcode/blob/v10_06_00_05/sbndcode/Utilities/detectorproperties_sbnd.fcl
rho = 1.928 - 0.00615 * temp ## g/cm^3, https://github.com/LArSoft/lardataalg/blob/v10_00_06/lardataalg/DetectorInfo/DetectorPropertiesStandard.cxx#L132
W_ion = 23.6e-6 ## MeV
etau_cv = 100. # ms

## EMB central parameters
alpha_emb_cv = 0.904
beta_90_cv = 0.204
R_emb_cv = 1.25

## C_cal central parameters
c_cal_data_cv = [0.0223037, 0.0219534, 0.0215156]
c_cal_mc_cv = [0.0203521, 0.0202351, 0.0200727]

def etau_corr(x, etau = etau_cv):
    ## == corr to update etau correction to etau_new from etau_cv
    x_max = 200.
    v_drift = 156.267 ## cm/ms
    x_clipped = x.clip(lower=-1. * x_max, upper=x_max)    
    t_drift = (x_max - np.abs(x_clipped)) / v_drift
    corr = np.exp(t_drift / etau)
    return corr

def new_dedx(hit_df, c_cal_frac = 1.0, plane = 2, alpha_emb = alpha_emb_cv, beta_90 = beta_90_cv, R_emb = R_emb_cv, etau = 100., is_mc = True):
    cos_phi = np.cos(hit_df.phi)
    beta_emb = beta_90 / (hit_df.efield * rho * np.sqrt(1 - cos_phi * cos_phi + cos_phi*cos_phi/(R_emb * R_emb)))

    this_c_cal = 0.02
    if is_mc:
        this_c_cal = c_cal_mc_cv[plane]
    else:
        this_c_cal = c_cal_data_cv[plane]
    this_c_cal *= c_cal_frac
    
    dqdx_e = hit_df.dqdx / this_c_cal
    this_etau_corr = 1.
    this_etau_corr = etau_corr(hit_df.x, etau)

    dqdx_e = dqdx_e * this_etau_corr
    dedx = np.exp(beta_emb * W_ion *dqdx_e) / beta_emb - alpha_emb/beta_emb
    return dedx

def particle_chi2(dEdx, ResRange, particle, dedx_a0, dedx_a1):
    #print(f'len(dEdx) = {len(dEdx)}, len(ResRange) = {len(ResRange)}')
    if particle != "kaon" and particle != "proton" and particle != "muon" and particle != "pion" :
        #print("Not a valid particle input")
        #return 99999., -1
        return 0., 0
    if len(dEdx) < 1 or len(ResRange) < 1:
        #print("No dEdx or ResRange data")
        #return 88888.0, 0
        return 0., 0
    
    N_max_hits = 1000
    this_N_calo = len(dEdx)
    this_N_hits = min(N_max_hits, this_N_calo)
    N_skip = 1
    dEdx_truncate_upper = 1000.0
    dEdx_truncate_bellow = 0.0
    max_res_range = 25.
    npt = 0

    dedx_exp = pd.cut(ResRange, chi2pid_temp[particle]["chi2_rr_arrs"], labels=chi2pid_temp[particle]["chi2_dEdx_arrs"]).astype(float)
    dedx_err = pd.cut(ResRange, chi2pid_temp[particle]["chi2_rr_arrs"], labels=chi2pid_temp[particle]["chi2_yerr_arrs"]).astype(float)
    dedx_res = (dedx_a0 + dedx_a1*dEdx**2)*dEdx

    v_chi2 = (dEdx - dedx_exp)**2 / (dedx_err**2 + dedx_res**2)

    # print(f'ResRange (min/max) = {ResRange.min()}, {ResRange.max()}')
    # print(f'dEdx (min/max) = {dEdx.min()}, {dEdx.max()}')
    when_chi2 = (ResRange < max_res_range) & (ResRange > 0.) & (dEdx < dEdx_truncate_upper) & (dEdx > dEdx_truncate_bellow)

    v_chi2= v_chi2.iloc[N_skip:-1 * N_skip]
    chi2_series = v_chi2[when_chi2]
    len_all = len(chi2_series)

    if(len_all) < 1:
        #return 77777.0, 0
        return 0., 0

    return chi2_series.sum() / len_all, len_all

def calculate_chi2_for_entry(group, particle, dedx_a0 = 0.04231, dedx_a1 = 0.0001783):
    chi2, ndof = particle_chi2(group[('dedx', '')], group[('rr', '')], particle, dedx_a0, dedx_a1)
    return chi2, ndof