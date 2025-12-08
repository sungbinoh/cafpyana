import numpy as np

PDG = {
    "muon": [13, "muon", 0.105,],
    "proton": [2212, "proton", 0.938272,],
    "neutron": [2112, "neutron", 0.9395654,],
    "pizero": [111, "pizero", 0.1349768],
    "pipm": [211, "piplus", 0.13957039],
    "argon": [1000180400, "argon", (18*0.938272 + 22*0.9395654)],
    "gamma": [22, "gamma", 0 ],
    "lambda": [3122, "lambda", 1.115683],
    "kaon_p": [321, "kaon_p",  0.493677],
    "sigma_p": [3222, "sigma_p", 1.18936],
    "kaon_0": [311, "kaon_0", 0.497648],
    "sigma_0": [3212, "sigma_0", 1.19246],
    "lambda_p_c": [4122, "lambda_p_c", 2.28646],
    "sigma_pp_c": [4222, "sigma_pp_c", 2.45397],
    "electron": [11, "electron", 0.510998950],
    "sigma_p_c": [4212, "sigma_p_c", 2.4529],
}

def InFV_trk(data): # cm
    xmin = -190.
    ymin = -190.
    zmin = 10.
    xmax = 190.
    ymax =  190.
    zmax =  490.
    return (data.x > xmin) & (data.x < xmax) & (data.y > ymin) & (data.y < ymax) & (data.z > zmin) & (data.z < zmax)


def InFV(data): # cm
    xmin = -190.
    ymin = -190.
    zmin = 10.
    xmax = 190.
    ymax =  190.
    zmax =  450.
    return (np.abs(data.x) > 10) & (np.abs(data.x) < 190) & (data.y > ymin) & (data.y < ymax) & (data.z > zmin) & (data.z < zmax)

## -- truth level flags
def Signal(df): # definition
    
    is_fv = InFV(df.position)
    is_numu = (df.pdg == 14)
    is_cc = (df.iscc == 1)
    is_2p0pi = (df.nmu_27MeV == 1) & (df.npi == 0) & (df.np_50MeV == 2) & (df.npi0 == 0)
    return is_fv & is_numu & is_cc & is_2p0pi

def cc2pNpi(df): # definition
    
    is_fv = InFV(df.position)
    is_numu = (df.pdg == 14)
    is_cc = (df.iscc == 1)
    is_2pNpi = (df.nmu_27MeV == 1) & ( df.npi + df.npi0 > 0 ) & (df.np_50MeV == 2) 
    return is_fv & is_numu & is_cc & is_2pNpi

def cc1p0pi(df): # definition
    
    is_fv = InFV(df.position)
    is_numu = (df.pdg == 14)
    is_cc = (df.iscc == 1)
    is_1p0pi = (df.nmu_27MeV == 1) & (df.npi == 0) & (df.np_50MeV == 1) & (df.npi0 == 0)
    return is_fv & is_numu & is_cc & is_1p0pi

def cc0p0pi(df): # definition
    
    is_fv = InFV(df.position)
    is_numu = (df.pdg == 14)
    is_cc = (df.iscc == 1)
    is_0p0pi = (df.nmu_27MeV == 1) & (df.npi == 0) & (df.np_50MeV == 0) & (df.npi0 == 0)
    return is_fv & is_numu & is_cc & is_0p0pi

def cc0pNpi(df): # definition
    
    is_fv = InFV(df.position)
    is_numu = (df.pdg == 14)
    is_cc = (df.iscc == 1)
    is_0pNpi = (df.nmu_27MeV == 1) & (df.npi + df.npi0 > 0) & (df.np_50MeV == 0)
    return is_fv & is_numu & is_cc & is_0pNpi

def cc1pMpi(df): # definition
    
    is_fv = InFV(df.position)
    is_numu = (df.pdg == 14)
    is_cc = (df.iscc == 1)
    is_1pMpi = (df.nmu_27MeV == 1) & (df.npi + df.npi0 > 1) & (df.np_50MeV == 1)
    return is_fv & is_numu & is_cc & is_1pMpi

def cc1p1pi(df): # definition
    
    is_fv = InFV(df.position)
    is_numu = (df.pdg == 14)
    is_cc = (df.iscc == 1)
    is_1p1pi = (df.nmu_27MeV == 1) & (df.npi + df.npi0 == 1) & (df.np_50MeV == 1)
    return is_fv & is_numu & is_cc & is_1p1pi

def out_range(df): # definition
    
    is_fv = InFV(df.position)
    is_numu = (df.pdg == 14)
    is_cc = (df.iscc == 1)
    is_out_range = (df.nmu_27MeV == 0) & (df.nmu == 1) & (df.np_50MeV == 0) & (df.np == 1)
    return is_fv & is_numu & is_cc & is_out_range

def bkg_other(df): # definition
    
    is_fv = InFV(df.position)
    is_numu = (df.pdg == 14)
    is_cc = (df.iscc == 1)
    is_1p0pi = (df.nmu_27MeV == 1) & (df.npi == 0) & (df.np_50MeV == 1) & (df.npi0 == 0)
    is_2p0pi = (df.nmu_27MeV == 1) & (df.npi == 0) & (df.np_50MeV == 2) & (df.npi0 == 0)    
    is_1p1pi = (df.nmu_27MeV == 1) & (df.npi + df.npi0 == 1) & (df.np_50MeV == 1)    
    is_2pNpi = (df.nmu_27MeV == 1) & (df.npi + df.npi0 > 0) & (df.np_50MeV == 2)   
    is_1pMpi = (df.nmu_27MeV == 1) & (df.npi + df.npi0 > 1) & (df.np_50MeV == 1)
    #is_out_range = (df.nmu_27MeV == 0) & (df.nmu == 1) & (df.np_50MeV == 0) & (df.np == 1)    

    mask = is_fv & is_numu & is_cc & ~is_1p0pi & ~is_2p0pi & ~is_2pNpi & ~is_1p1pi & ~is_1pMpi
    # if mask.any():
    #     print(df.loc[mask, ['nmu_27MeV', 'np_50MeV', 'npi', 'npi0']])        
        
    return mask

## -- reco level flags
def pass_slc_with_n_pfps(df, n = 3):
    group_levels = ['entry', 'rec.slc..index']
    
    # Count how many pfps per slice
    pfp_counts = df.groupby(level=group_levels).size()

    # Get only slices with exactly 3 pfps
    valid_slices = pfp_counts[pfp_counts == n].index

    # Apply the mask to original DataFrame
    filtered_df = df.loc[df.index.droplevel('rec.slc.reco.pfp..index').isin(valid_slices)]
    
    filtered_df = get_n_recopid_per_slc(filtered_df)
    
    filtered_df = filtered_df[ (filtered_df["muon_counter"] == 1) & (filtered_df["proton_counter"] == 2) & (filtered_df["pion_counter"] == 0)]    

    return filtered_df

def Avg(df, pid, drop_0=True):  # exclude value if 0
    if drop_0:
        df = df.replace(0, np.nan)
    # let's just use only the collectin planes
    average = df[("chi2pid", "I2", "chi2_"+pid)]
    return average

def get_pid_result(row):

    chi2_muon = row[('pfp', 'trk', 'chi2pid', 'I2', 'chi2_muon', '')]    
    chi2_proton = row[('pfp', 'trk', 'chi2pid', 'I2', 'chi2_proton', '')]
    len = row[('pfp', 'trk', 'len', '', '', '')]

    if chi2_muon < 25. and chi2_proton > 90.:
        return 13  # muon
    elif chi2_proton < 90.:
        return 2212  # proton
    else:
        return 211  # charged pion   
        
def add_n_slice_col(reco_df):
    df_reset = reco_df.reset_index()
    slc_counts = (
        df_reset[['entry', 'rec.slc..index']]
        .drop_duplicates()
        .groupby(['entry'])
        .size()
        .reset_index(name='n_slc_per_entry')
    )

    slc_counts.columns = pd.MultiIndex.from_tuples([
        ('entry', '', '', '', '', ''),
        ('slc', 'n_slc_per_entry', '', '', '', '')
    ])
    df_reset = df_reset.merge(slc_counts, on=[('entry', '', '', '', '', '')])
    df_reset = df_reset.set_index(["entry", "rec.slc..index", "rec.slc.reco.pfp..index"], verify_integrity=True)
    return df_reset     

def get_n_recopid_per_slc(df):
    
    pid_series = df.pfp.trk.reco_pid
    this_df = pid_series.reset_index()

    muons = this_df[this_df["reco_pid"] == 13]
    protons = this_df[this_df["reco_pid"] == 2212]
    pions = this_df[this_df["reco_pid"] == 211]    

    muon_counts = muons.groupby(["entry", "rec.slc..index"]).size().rename("n_mu")
    proton_counts = protons.groupby(["entry", "rec.slc..index"]).size().rename("n_proton")
    pion_counts = pions.groupby(["entry", "rec.slc..index"]).size().rename("n_pion")    

    this_df = this_df.merge(muon_counts, on=["entry", "rec.slc..index"], how="left")
    this_df = this_df.merge(proton_counts, on=["entry", "rec.slc..index"], how="left")
    this_df = this_df.merge(pion_counts, on=["entry", "rec.slc..index"], how="left")    

    this_df["n_mu"] = this_df["n_mu"].fillna(0).astype(int)
    this_df["n_proton"] = this_df["n_proton"].fillna(0).astype(int)
    this_df["n_pion"] = this_df["n_pion"].fillna(0).astype(int)    

    this_df.set_index(["entry", "rec.slc..index", "rec.slc.reco.pfp..index"], inplace=True)

    df.loc[:, ('muon_counter', '', '', '', '', '')] = this_df['n_mu']
    df.loc[:, ('proton_counter', '', '', '', '', '')] = this_df['n_proton']
    df.loc[:, ('pion_counter', '', '', '', '', '')] = this_df['n_pion']

    #df[('muon_counter', '', '', '', '', '')] = this_df.n_mu
    #df[('proton_counter', '', '', '', '', '')] = this_df.n_proton
    #df[('pion_counter', '', '', '', '', '')] = this_df.n_pion    

    return df

def add_contained_col(df):
    contained = InFV(df.pfp.trk.start) & InFV(df.pfp.trk.end)
    df[('pfp', 'contained', '', '', '', '')] = contained

def reco_imbalance(muon_dir_x, muon_dir_y, muon_dir_z, range_P_muon, 
                 lead_proton_dir_x, lead_proton_dir_y, lead_proton_dir_z, lead_range_P_proton,
                 rec_proton_dir_x, rec_proton_dir_y, rec_proton_dir_z, rec_range_P_proton):

    # deltapt
    deltapt_x = muon_dir_x.iloc[0] * range_P_muon.iloc[0] + lead_proton_dir_x.iloc[0] * lead_range_P_proton.iloc[0] + rec_proton_dir_x.iloc[0] * rec_range_P_proton.iloc[0]
    px_sq = np.power(deltapt_x, 2.)
    deltapt_y = muon_dir_y.iloc[0] * range_P_muon.iloc[0] + lead_proton_dir_y.iloc[0] * lead_range_P_proton.iloc[0] + rec_proton_dir_y.iloc[0] * rec_range_P_proton.iloc[0]
    py_sq = np.power(deltapt_y, 2.)
    deltapt = np.sqrt(px_sq + py_sq)

    # deltaalphat
    muon_px = muon_dir_x.iloc[0] * range_P_muon.iloc[0]
    muon_py = muon_dir_y.iloc[0] * range_P_muon.iloc[0]
    muon_pz = muon_dir_z.iloc[0] * range_P_muon.iloc[0]    
    muon_pt = np.sqrt( np.power(muon_px,2.) + np.power(muon_py,2.) ) 
    muon_p = range_P_muon.iloc[0]
    deltaalphat_denom = muon_pt * deltapt
    deltaalphat_num = - ( muon_px * deltapt_x + muon_py * deltapt_y ) 
    deltaalphat = np.arccos( deltaalphat_num / deltaalphat_denom) * 180./np.pi
    
    #deltaphit
    lead_proton_x = lead_proton_dir_x.iloc[0] * lead_range_P_proton.iloc[0] 
    rec_proton_x = rec_proton_dir_x.iloc[0] * rec_range_P_proton.iloc[0]     
    proton_px = lead_proton_x + rec_proton_x
    
    lead_proton_y = lead_proton_dir_y.iloc[0] * lead_range_P_proton.iloc[0] 
    rec_proton_y = rec_proton_dir_y.iloc[0] * rec_range_P_proton.iloc[0]     
    proton_py = lead_proton_y + rec_proton_y
    
    lead_proton_z = lead_proton_dir_z.iloc[0] * lead_range_P_proton.iloc[0] 
    rec_proton_z = rec_proton_dir_z.iloc[0] * rec_range_P_proton.iloc[0]    
    proton_pz = lead_proton_z + rec_proton_z    

    proton_pt = np.sqrt( np.power(proton_px,2.) + np.power(proton_py,2.) )   
    proton_p = np.sqrt( np.power(proton_px,2.) + np.power(proton_py,2.) + np.power(proton_pz,2.) )      
    lead_proton_p = np.sqrt( np.power(lead_proton_x,2.) + np.power(lead_proton_y,2.) + np.power(lead_proton_z,2.) )
    rec_proton_p = np.sqrt( np.power(rec_proton_x,2.) + np.power(rec_proton_y,2.) + np.power(rec_proton_z,2.) )
      
    deltaphit_denom = muon_pt * proton_pt
    deltaphit_num = - ( muon_px * proton_px + muon_py * proton_py )
    deltaphit = np.arccos(deltaphit_num / deltaphit_denom) * 180./np.pi
    
    # cos(theta_LR)
    costheta_lr_num = lead_proton_x * rec_proton_x + lead_proton_y * rec_proton_y + lead_proton_z * rec_proton_z
    costheta_lr_denom = lead_range_P_proton.iloc[0] * rec_range_P_proton.iloc[0]
    costheta_lr = costheta_lr_num / costheta_lr_denom
    
    # costheta_mu_sum
    costheta_mu_sum_num = muon_px * proton_px + muon_py * proton_py + muon_pz * proton_pz
    costheta_mu_sum_denom = muon_p * proton_p
    costheta_mu_sum = costheta_mu_sum_num / costheta_mu_sum_denom
        
    #e_cal
    e_mu = np.sqrt( np.power(muon_p,2.) + np.power(0.105,2) )
    e_lead_p = np.sqrt( np.power(lead_proton_p,2.) + np.power(0.938272,2) )
    ke_lead_p = e_lead_p - 0.938272
    e_rec_p = np.sqrt( np.power(rec_proton_p,2.) + np.power(0.938272,2) )
    ke_rec_p = e_rec_p - 0.938272    
    e_cal = e_mu + ke_lead_p + ke_rec_p + 0.0309 # https://link.springer.com/article/10.1140/epjc/s10052-019-6750-3
    
    #p_l
    p_l = muon_pz + proton_pz - e_cal
    
    #pn
    pn = np.sqrt( np.power(p_l,2.) + np.power(deltapt,2.) )
    pn_x = deltapt_x
    pn_y = deltapt_y
    pn_z = p_l
    
    # q (energy transfer)
    q_x = - muon_px
    q_y = - muon_py
    q_z = - muon_pz + e_cal
    q = np.sqrt( np.power(q_x,2.) + np.power(q_y,2.) + np.power(q_z,2.) )
    
    # phi_3d
    phi_3d_num = q_x * proton_px + q_y * proton_py + q_z * proton_pz
    phi_3d_denom = q * proton_p
    phi_3d = np.arccos(phi_3d_num / phi_3d_denom) * 180./np.pi
    
    #alpha_3d
    alpha_3d_num = q_x * pn_x + q_y * pn_y + q_z * pn_z
    alpha_3d_denom = q * pn
    alpha_3d = np.arccos(alpha_3d_num / alpha_3d_denom) * 180./np.pi
        
    # invariant mass
    invm = np.sqrt(np.power((e_lead_p + e_rec_p), 2) - np.power(proton_p, 2));        
        
    return pd.Series({
                    'deltapt': deltapt,
                    'deltaalphat': deltaalphat,
                    'deltaphit': deltaphit,
                    'costheta_lr':costheta_lr,
                    'costheta_mu_sum':costheta_mu_sum,
                    'e_cal': e_cal,
                    'pn':pn,
                    'phi_3d': phi_3d,
                    'alpha_3d': alpha_3d,
                    'invm': invm       
                    })

def measure_reco_imbalance(group):
    
    muons = group[group.pfp.trk.reco_pid == 13]
    protons = group[group.pfp.trk.reco_pid == 2212]    
    
    lead_protons = protons.sort_values(by=("pfp", "trk", "len", "", "", ""), ascending=False).groupby(level=[0,1]).nth(0)
    rec_protons = protons.sort_values(by=("pfp", "trk", "len", "", "", ""), ascending=False).groupby(level=[0,1]).nth(1)    
    
    muon_dir_x = muons[('pfp', 'trk', 'dir', 'x', '', '')]
    muon_dir_y = muons[('pfp', 'trk', 'dir', 'y', '', '')]
    muon_dir_z = muons[('pfp', 'trk', 'dir', 'z', '', '')]
    muon_range = muons[('pfp', 'trk', 'rangeP', 'p_muon', '', '')]
    
    lead_proton_dir_x = lead_protons[('pfp', 'trk', 'dir', 'x', '', '')]
    lead_proton_dir_y = lead_protons[('pfp', 'trk', 'dir', 'y', '', '')]
    lead_proton_dir_z = lead_protons[('pfp', 'trk', 'dir', 'z', '', '')]
    lead_proton_range = lead_protons[('pfp', 'trk', 'rangeP', 'p_proton', '', '')]
 
    rec_proton_dir_x = rec_protons[('pfp', 'trk', 'dir', 'x', '', '')]
    rec_proton_dir_y = rec_protons[('pfp', 'trk', 'dir', 'y', '', '')]
    rec_proton_dir_z = rec_protons[('pfp', 'trk', 'dir', 'z', '', '')]
    rec_proton_range = rec_protons[('pfp', 'trk', 'rangeP', 'p_proton', '', '')]        
    
    return reco_imbalance(muon_dir_x, muon_dir_y, muon_dir_z, muon_range, 
                        lead_proton_dir_x, lead_proton_dir_y, lead_proton_dir_z, lead_proton_range,
                        rec_proton_dir_x, rec_proton_dir_y, rec_proton_dir_z, rec_proton_range)  