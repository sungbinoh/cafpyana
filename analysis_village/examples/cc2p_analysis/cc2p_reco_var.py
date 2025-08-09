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
    is_2p0pi = (df.nmu_40MeV == 1) & (df.npi_30MeV == 0) & (df.np_50MeV == 2) & (df.npi0 == 0)
    return is_fv & is_numu & is_cc & is_2p0pi

## -- reco level flags
def pass_slc_with_n_pfps(df, n = 3):
    group_levels = ['entry', 'rec.slc..index']
    
    # Count how many pfps per slice
    pfp_counts = df.groupby(level=group_levels).size()

    # Get only slices with exactly 3 pfps
    valid_slices = pfp_counts[pfp_counts == n].index

    # Apply the mask to original DataFrame
    filtered_df = df.loc[df.index.droplevel('rec.slc.reco.pfp..index').isin(valid_slices)]

    return filtered_df

def Avg(df, pid, drop_0=True):  # exclude value if 0
    if drop_0:
        df = df.replace(0, np.nan)
    # let's just use only the collectin planes
    average = df[("chi2pid", "I2", "chi2_"+pid)]
    return average

def add_contained_col(df):
    containd = InFV(df.pfp.trk.start) & InFV(df.pfp.trk.end)
    df[('pfp', 'contained', '', '', '', '')] = contained
    
def get_pid_result(row):
    chi2_muon = row[('pfp', 'trk', 'chi2pid', 'I2', 'chi2_muon', '')]
    chi2_proton = row[('pfp', 'trk', 'chi2pid', 'I2', 'chi2_proton', '')]

    if chi2_muon < 25. and chi2_proton > 100.:
        return 13  # muon
    else:
        return 2212  # proton   
    
def add_n_slice_col(reco_df):
    df_reset = reco_df.reset_index()
    slc_counts = (
        df_reset[['__ntuple', 'entry', 'rec.slc..index']]
        .drop_duplicates()
        .groupby(['__ntuple', 'entry'])
        .size()
        .reset_index(name='n_slc_per_entry')
    )

    slc_counts.columns = pd.MultiIndex.from_tuples([
        ('__ntuple', '', '', '', '', ''),
        ('entry', '', '', '', '', ''),
        ('slc', 'n_slc_per_entry', '', '', '', '')
    ])
    df_reset = df_reset.merge(slc_counts, on=[('__ntuple', '', '', '', '', ''), ('entry', '', '', '', '', '')])
    df_reset = df_reset.set_index(["__ntuple", "entry", "rec.slc..index", "rec.slc.reco.pfp..index"], verify_integrity=True)
    return df_reset     

def get_n_recopid_per_slc(df):
    pid_series = df.pfp.trk.chi2pid.I2.reco_pid
    this_df = pid_series.reset_index()

    muons = this_df[this_df["reco_pid"] == 13]
    protons = this_df[this_df["reco_pid"] == 2212]

    muon_counts = muons.groupby(["__ntuple", "entry", "rec.slc..index"]).size().rename("n_mu")
    proton_counts = protons.groupby(["__ntuple", "entry", "rec.slc..index"]).size().rename("n_proton")

    this_df = this_df.merge(muon_counts, on=["__ntuple", "entry", "rec.slc..index"], how="left")
    this_df = this_df.merge(proton_counts, on=["__ntuple", "entry", "rec.slc..index"], how="left")

    this_df["n_mu"] = this_df["n_mu"].fillna(0).astype(int)
    this_df["n_proton"] = this_df["n_proton"].fillna(0).astype(int)

    this_df.set_index(["__ntuple", "entry", "rec.slc..index", "rec.slc.reco.pfp..index"], inplace=True)

    df[('muon_counter', '', '', '', '', '')] = this_df.n_mu
    df[('proton_counter', '', '', '', '', '')] = this_df.n_proton

    return df

def add_contained_col(df):
    contained = InFV(df.pfp.trk.start) & InFV(df.pfp.trk.end)
    df[('pfp', 'contained', '', '', '', '')] = contained
    
def reco_deltapt(dir_x, dir_y, dir_z, range_P_muon, range_P_proton):
    # -- assume first particle is muon and the others are the protons
    p_mu = range_P_muon.iloc[0]
    # leading proton
    p_p_l = range_P_proton.iloc[0]
    # recoil proton
    p_p_r = range_P_proton.iloc[1]
    # -- if second proton is more energetic, swap the proton momentum assumption
    if(range_P_proton.iloc[0] < range_P_proton.iloc[1]):
        p_l = range_P_proton.iloc[1]
        p_r = range_P_proton.iloc[0]

    # -- each term
    px_sq = np.power(range_P_muon.iloc[0] * dir_x.iloc[0] + range_P_proton.iloc[1] * dir_x.iloc[1] + range_P_proton.iloc[2] * dir_x.iloc[2], 2.)
    py_sq = np.power(range_P_muon.iloc[0] * dir_y.iloc[0] + range_P_proton.iloc[1] * dir_y.iloc[1] + range_P_proton.iloc[2] * dir_y.iloc[2], 2.)
    deltapt = np.sqrt(px_sq + py_sq)
    
    #print(deltapt)
    return deltapt

def measure_reco_deltapt(group):
    dir_x = group[('pfp', 'trk', 'dir', 'x', '', '')]
    dir_y = group[('pfp', 'trk', 'dir', 'y', '', '')]
    dir_z = group[('pfp', 'trk', 'dir', 'z', '', '')]
    range_P_muon = group[('pfp', 'trk', 'rangeP', 'p_muon', '', '')]
    range_P_proton = group[('pfp', 'trk', 'rangeP', 'p_proton', '', '')]

    # Call reco_deltapt function
    return reco_deltapt(dir_x, dir_y, dir_z, range_P_muon, range_P_proton)    