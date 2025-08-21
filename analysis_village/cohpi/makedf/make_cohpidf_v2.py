from makedf.makedf import *
from pyanalib.pandas_helpers import *
from makedf.util import *

def InFV_nohiyz(data):
    xmin = 10.
    xmax = 190.
    zmin = 10.
    zmax = 450.
    ymax_highz = 100.
    pass_xz = (np.abs(data.x) > xmin) & (np.abs(data.x) < xmax) & (data.z > zmin) & (data.z < zmax)
    pass_y = ((data.z < 250) & (np.abs(data.y) < 190.)) | ((data.z > 250) & (data.y > -190.) & (data.y < ymax_highz))
    return pass_xz & pass_y

def InFV_nohiyz_trk(data):
    xmax = 190.
    zmin = 10.
    zmax = 450.
    ymax_highz = 100.
    pass_xz = (np.abs(data.x) < xmax) & (data.z > zmin) & (data.z < zmax)
    pass_y = ((data.z < 250) & (np.abs(data.y) < 190.)) | ((data.z > 250) & (data.y > -190.) & (data.y < ymax_highz))
    return pass_xz & pass_y

## -- truth level flags
def Signal(df): # definition
    is_fv = InFV_nohiyz(df.position)
    is_1pi0p0n = (df.nmu_27MeV == 1) & (df.npi_30MeV == 1) & (df.np_20MeV == 0) & (df.npi0 == 0) & (df.nn_0MeV == 0)
    return is_fv & is_1pi0p0n

## -- reco level flags
def pass_slc_with_n_pfps(df, n = 2):
    group_levels = ['entry', 'rec.slc..index']
    
    # Count how many pfps per slice
    pfp_counts = df.groupby(level=group_levels).size()

    # Get only slices with at least 2 pfps
    valid_slices = pfp_counts[pfp_counts == n].index

    # Apply the mask to original DataFrame
    filtered_df = df.loc[df.index.droplevel('rec.slc.reco.pfp..index').isin(valid_slices)]

    return filtered_df

def apply_dir_z_cut(df):
    if '_trk_rank' in df.columns:
        df = df.drop(columns=[('_trk_rank', '', '', '', '', '')])

    df = df.sort_index()
    # Grouping by MultiIndex index levels
    group_keys = ['entry', 'rec.slc..index']  # passed as strings for MultiIndex levels

    # MultiIndex column names
    trk_dirz_col = ('pfp', 'trk', 'dir', 'z', '', '')
    trk_len_col  = ('pfp', 'trk', 'len', '', '', '')

    # Step 1: Rank track length descending per SLC group
    df['_trk_rank'] = df.groupby(level=group_keys)[[trk_len_col]] \
                        .rank(method='first', ascending=False)

    # Step 2: Apply cuts
    sel_longest  = (df['_trk_rank'] == 1) & (df[trk_dirz_col] > 0.7)
    sel_2nd_long = (df['_trk_rank'] == 2) & (df[trk_dirz_col] > 0.5)

    # Step 3: Keep only desired rows
    df = df[sel_longest | sel_2nd_long].drop(columns=['_trk_rank']).copy()

    if '_trk_rank' in df.columns:
        df = df.drop(columns=[('_trk_rank', '', '', '', '', '')])

    return df

def Avg(df, pid, drop_0=True):  # average score of 3 planes, exclude value if 0
    if drop_0:
        df = df.replace(0, np.nan)
    #average = df[[("chi2pid", "I0", "chi2_"+pid), ("chi2pid", "I1", "chi2_"+pid), ("chi2pid", "I2", "chi2_"+pid)]].mean(skipna=drop_0, axis=1)
    # let's just use only the collectin planes
    average = df[("chi2pid", "I2", "chi2_"+pid)]
    return average

def add_contained_col(df):
    contained = InFV_nohiyz_trk(df.pfp.trk.start) & InFV_nohiyz_trk(df.pfp.trk.end)
    df['contained'] = contained

def reco_t(dir_x, dir_y, dir_z, range_P_muon, range_P_pion):
    # -- assume first particle is muon and the other is pion
    mass_0 = PDG["muon"][2]
    mass_1 = PDG["pipm"][2]
    p_0 = range_P_muon.iloc[0]
    p_1 = range_P_pion.iloc[1]
    # -- if second track is longer, swap the mass assumption
    if(range_P_muon.iloc[0] > range_P_muon.iloc[1]):
        mass_0 = PDG["pipm"][2]
        mass_1 = PDG["muon"][2]
        p_0 = range_P_pion.iloc[0]
        p_1 = range_P_muon.iloc[1]
    E_0 = np.sqrt(mass_0**2 + p_0**2)
    E_1 = np.sqrt(mass_1**2 + p_1**2)

    # -- each term
    px_sq = np.power(p_0 * dir_x.iloc[0] + p_1 * dir_x.iloc[1], 2.)
    py_sq = np.power(p_0 * dir_y.iloc[0] + p_1 * dir_y.iloc[1], 2.)
    pz_sq = np.power(E_0 + E_1 - p_0 * dir_z.iloc[0] - p_1 * dir_z.iloc[1], 2.)
    abs_t = px_sq + py_sq + pz_sq
    
    #print(abs_t)
    return abs_t

def measure_reco_t(group):
    dir_x = group[('pfp', 'trk', 'dir', 'x', '', '')]
    dir_y = group[('pfp', 'trk', 'dir', 'y', '', '')]
    dir_z = group[('pfp', 'trk', 'dir', 'z', '', '')]
    range_P_muon = group[('pfp', 'trk', 'rangeP', 'p_muon', '', '')]
    range_P_pion = group[('pfp', 'trk', 'rangeP', 'p_pion', '', '')]

    # Call reco_t function
    return reco_t(dir_x, dir_y, dir_z, range_P_muon, range_P_pion)

def opening_angle(dir_x, dir_y, dir_z):
    this_cos_theta = dir_x.iloc[0] * dir_x.iloc[1] + dir_y.iloc[0] * dir_y.iloc[1] + dir_z.iloc[0] * dir_z.iloc[1]
    return this_cos_theta

def measure_opening_angle(group):
    dir_x = group[('pfp', 'trk', 'dir', 'x', '', '')]
    dir_y = group[('pfp', 'trk', 'dir', 'y', '', '')]
    dir_z = group[('pfp', 'trk', 'dir', 'z', '', '')]

    # Call reco_t function
    return opening_angle(dir_x, dir_y, dir_z)

def beam_totp_angle(n_trk_mupid, dir_x, dir_y, dir_z, range_P_muon, range_P_pion, mu_pid_pass):
    if n_trk_mupid != 2:
        return -999.  
    dir_x = dir_x[mu_pid_pass]
    dir_y = dir_y[mu_pid_pass]
    dir_z = dir_z[mu_pid_pass]
    range_P_muon = range_P_muon[mu_pid_pass]
    range_P_pion = range_P_pion[mu_pid_pass]
    if(range_P_muon.size != 2):
        print("error, dir_x.len != 2")
        return -888.
    
    # -- assume first particle is muon and the other is pion
    p_0 = range_P_muon.iloc[0]
    p_1 = range_P_pion.iloc[1]
    # -- if second track is longer, swap the mass assumption
    if(range_P_muon.iloc[0] > range_P_muon.iloc[1]):
        p_0 = range_P_pion.iloc[0]
        p_1 = range_P_muon.iloc[1]

    totpx = p_0 * dir_x.iloc[0] + p_1 * dir_x.iloc[1]
    totpy = p_0 * dir_y.iloc[0] + p_1 * dir_y.iloc[1]
    totpz = p_0 * dir_z.iloc[0] + p_1 * dir_z.iloc[1]

    totp_cos = totpz / np.power(np.power(totpx, 2.) + np.power(totpy, 2.) + np.power(totpz, 2.) , 0.5)
    return totp_cos
    
def measure_beam_totp_angle(group):
    n_trk_mupid = group[('n_trk_mupid', '', '')].iloc[0]
    dir_x = group[('trk', 'dir', 'x')]
    dir_y = group[('trk', 'dir', 'y')]
    dir_z = group[('trk', 'dir', 'z')]
    range_P_muon = group[('trk', 'rangeP', 'p_muon')]
    range_P_pion = group[('trk', 'rangeP', 'p_pion')]
    mu_pid_pass = group[('trk', 'mu_pid_pass', '')]

    # Call reco_t function
    return beam_totp_angle(n_trk_mupid, dir_x, dir_y, dir_z, range_P_muon, range_P_pion, mu_pid_pass)

def get_true_t(df):
    t = (df.E - df.mu.genE - df.cpi.genE)**2 - (df.momentum.x - df.mu.genp.x - df.cpi.genp.x)**2 - (df.momentum.y - df.mu.genp.y - df.cpi.genp.y)**2 - (df.momentum.z - df.mu.genp.z - df.cpi.genp.z)**2
    return np.abs(t)

## -- data fram maker for cohpi analysis
def make_cohpidf_v2(f):
    
    pandora_df = make_pandora_df(f)
    
    #### (1) FV cut
    pandora_df = pandora_df[InFV_nohiyz(pandora_df.slc.vertex)]

    #### (2) Not clear cosmic cut
    pandora_df = pandora_df[pandora_df.slc.is_clear_cosmic == 0]
    #pandora_df = pandora_df[pandora_df[('slc', 'is_clear_cosmic', '', '', '', '')] == 0]
    
    #### (3) Track multiplicity cut
    ######## (3) - a: keep only pfp objects with length > 4 cm and dist_to_vertex < 6 cm
    pandora_df = pandora_df[(pandora_df.pfp.trk.len > 4.) & (pandora_df.pfp.dist_to_vertex < 6.)]
    ######## (3) - b: keep onlhy slices with exactly two pfp object passing the requirement
    pandora_df = pass_slc_with_n_pfps(pandora_df)

    #### (4) The two trakcs should have track score > 0.5
    pandora_df = pandora_df[pandora_df.pfp.trackScore > 0.5]
    pandora_df = pass_slc_with_n_pfps(pandora_df)

    #### (5) Tho two tracks should satisfy chi2 pid cut 
    pandora_df = pandora_df[(pandora_df.pfp.trk.chi2pid.I2.chi2_muon < 25.) & (pandora_df.pfp.trk.chi2pid.I2.chi2_proton > 100.)]
    pandora_df = pass_slc_with_n_pfps(pandora_df)

    #### (6) dir Z cut
    pandora_df = apply_dir_z_cut(pandora_df)
    pandora_df = pass_slc_with_n_pfps(pandora_df)

    #### (6) Nuscore > 0.65 -> not applied
    #pandora_df = pandora_df[pandora_df.slc.nu_score > 0.65]

    #### (7) Cosine of the opening angle between the two track should be greater than 0.5
    opening_angle_series = pandora_df.groupby(['entry', 'rec.slc..index']).apply(measure_opening_angle)
    if pandora_df.empty:
        pandora_df[('slc', 'opening_angle', '', '', '', '')] = pd.Series(dtype='float')
    else:
        pandora_df[('slc', 'opening_angle', '', '', '', '')] = opening_angle_series
    pandora_df = pandora_df[pandora_df.slc.opening_angle > 0.5]

    #### (8) Tracks are contained
    add_contained_col(pandora_df)
    pandora_df = pandora_df[pandora_df.contained]
    pandora_df = pass_slc_with_n_pfps(pandora_df)

    #### (9) collect only slc variables of interest
    ######## (9) - a: reco momentum/angle of the two tracks
    long_trk_df = pandora_df.sort_values(by=("pfp", "trk", "len", "", "", ""), ascending=False).groupby(level=[0,1]).nth(0)
    short_trk_df = pandora_df.sort_values(by=("pfp", "trk", "len", "", "", ""), ascending=False).groupby(level=[0,1]).nth(1)

    range_p_mu_series = long_trk_df.pfp.trk.rangeP.p_muon
    range_p_mu_series = range_p_mu_series.reset_index(level='rec.slc.reco.pfp..index', drop=True)

    range_p_pi_series = short_trk_df.pfp.trk.rangeP.p_pion
    range_p_pi_series = range_p_pi_series.reset_index(level='rec.slc.reco.pfp..index', drop=True)

    cos_theta_mu_series = long_trk_df.pfp.trk.dir.z
    cos_theta_mu_series = cos_theta_mu_series.reset_index(level='rec.slc.reco.pfp..index', drop=True)

    cos_theta_pi_series = short_trk_df.pfp.trk.dir.z
    cos_theta_pi_series = cos_theta_pi_series.reset_index(level='rec.slc.reco.pfp..index', drop=True)

    ######## (9) - b: reco |t|
    if pandora_df.empty:
        empty_index = pd.MultiIndex(
            levels=[[], []],
            codes=[[], []],
            names=['entry', 'rec.slc..index']
        )
        reco_t_series = pd.Series(dtype='float', name='reco_t', index=empty_index)
        #reco_t_series = pd.Series(dtype='float', name='reco_t')
    else:
        reco_t_series = pandora_df.groupby(['entry', 'rec.slc..index']).apply(measure_reco_t)
    #reco_t_series = reco_t_series.reset_index(level='rec.slc.reco.pfp..index', drop=True)

    ######## (9) - c: slc.tmatch.idx for truth matching
    tmatch_idx_series = pandora_df.slc.tmatch.idx
    tmatch_idx_series = tmatch_idx_series.reset_index(level='rec.slc.reco.pfp..index', drop=True)
    
    ## (10) creat a slice-based reco df
    slcdf = pd.DataFrame({
        'range_p_mu': range_p_mu_series,
        'range_p_pi': range_p_pi_series,
        'reco_cos_theta_mu': cos_theta_mu_series,
        'reco_cos_theta_pi': cos_theta_pi_series,
        'reco_t': reco_t_series,
        'tmatch_idx': tmatch_idx_series
    })

    #print(slcdf.columns)
    return slcdf
    
def make_cohpi_nudf(f):
    
    nudf = make_mcdf(f)
    nudf["ind"] = nudf.index.get_level_values(1)
    wgtdf = pd.concat([bnbsyst.bnbsyst(f, nudf.ind), geniesyst.geniesyst_sbnd(f, nudf.ind)], axis=1)

    true_t = get_true_t(nudf).fillna(999999)
    nudf['true_t'] = true_t
    
    is_fv = InFV_nohiyz_trk(nudf.position)
    is_signal = Signal(nudf)
    is_cc = nudf.iscc
    genie_mode = nudf.genie_mode
    w = nudf.w

    try :
        nuint_categ = pd.Series(8, index=nudf.index)
    except Exception as e:
        print(f"Error init nuint_categ")
        return

    nuint_categ[~is_fv] = -1  # Out of FV
    nuint_categ[is_fv & is_signal] = 1 # Signal
    nuint_categ[is_fv & ~is_cc & ~is_signal] = 0  # NC
    nuint_categ[is_fv & is_cc & ~is_signal & (genie_mode == 3)] = 2  # Non-signal CCCOH
    nuint_categ[is_fv & is_cc & ~is_signal & (genie_mode == 0)] = 3  # CCQE
    nuint_categ[is_fv & is_cc & ~is_signal & (genie_mode == 10)] = 4  # 2p2h
    nuint_categ[is_fv & is_cc & ~is_signal & (genie_mode != 0) & (genie_mode != 3) & (genie_mode != 10) & ((w < 1.4) | (genie_mode == 1))] = 5  # RES
    nuint_categ[is_fv & is_cc & ~is_signal & (genie_mode != 0) & (genie_mode != 3) & (genie_mode != 10) & ((w > 2.0) | (genie_mode == 2))] = 6  # DIS

    nudf['nuint_categ'] = nuint_categ
    
    true_t_series = nudf.true_t

    mu_p_series = magdf(nudf.mu.genp)
    cpi_p_series = magdf(nudf.cpi.genp)

    mu_cos_theta_series = nudf.mu.genp.z / mu_p_series
    cpi_cos_theta_series = nudf.cpi.genp.z / cpi_p_series

    this_nudf = pd.DataFrame({
        'true_p_mu': mu_p_series,
        'true_p_pi': cpi_p_series,
        'true_cos_theta_mu': mu_cos_theta_series,
        'true_cos_theta_pi': cpi_cos_theta_series,
        'true_t': true_t_series,
        'nuint_categ': nuint_categ
    })
    this_nudf.columns = pd.MultiIndex.from_tuples([(col, '') for col in this_nudf.columns])
    
    this_nudf = multicol_concat(this_nudf, wgtdf)

    return this_nudf
