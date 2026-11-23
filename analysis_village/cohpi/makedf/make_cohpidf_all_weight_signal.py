from makedf.makedf import *
from pyanalib.pandas_helpers import *
from makedf.util import *
from analysis_village.cohpi.makedf.cohpi_selection import *
from makedf import getsyst
from makedf import chi2pid

# Define an empty global variable to hold all systematic weights
all_systs_weights_df = None

def Avg(df, pid, drop_0=True):  # average score of 3 planes, exclude value if 0
    if drop_0:
        df = df.replace(0, np.nan)
    #average = df[[("chi2pid", "I0", "chi2_"+pid), ("chi2pid", "I1", "chi2_"+pid), ("chi2pid", "I2", "chi2_"+pid)]].mean(skipna=drop_0, axis=1)
    # let's just use only the collectin planes
    average = df[("chi2pid", "I2", "chi2_"+pid)]
    return average

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

def make_slc_var(df):
    out = df.groupby(level=['entry', 'rec.slc..index'], sort=False).first()
    return out

def calo_variation_chi2(f, ismc, plane=2):
    """chi2_muon / chi2_proton on `plane` for every chi2pid.CALO_VARIATIONS entry.
    Same recipe as make_pandora_df's CV recalo (no hit FV cut), so 'cv' == *_new.
    Returns {var: {'muon': Series, 'proton': Series}} indexed by (entry, slc, pfp)."""
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det
    if not (det == 1).all():
        raise NotImplementedError("CALO_VARIATIONS are SBND parameter sets")
    trkhitdf = make_trkhitdf(f, plane)
    out = {}
    for var in CALO_VARS:
        trkhitdf["dedx_redo"] = chi2pid.dedx(trkhitdf, gain="SBND", calibrate="SBND", plane=plane,
                                            isMC=ismc, new_calo_params=chi2pid.CALO_VARIATIONS[var])
        out[var] = {par: chi2pid.chi2par(trkhitdf, dedxname="dedx_redo", par=par)[0] for par in ("muon", "proton")}
    return out

## collect syst weights first

## -- data fram maker for cohpi analysis
def make_cohpidf_slc(f):

    hdrdf = make_mchdrdf(f)
    ismc = hdrdf.ismc.iloc[0]

    is_recalo = True
    if is_recalo:
        pandora_df = make_pandora_df_calo_update(f)
    else:
        pandora_df = make_pandora_df(f)

    barycenterFM_df = loadbranches(f["recTree"], barycenterFMbranches).rec
    pandora_df = multicol_merge(barycenterFM_df, pandora_df, left_index=True, right_index=True, how="right", validate="one_to_many")

    #### (0) BCFM and nuscore
    bcfm_score = pandora_df.slc.barycenterFM.score
    bcfm_score = make_slc_var(bcfm_score)
    nu_score = pandora_df.slc.nu_score
    nu_score = make_slc_var(nu_score)

    #### (1) FV cut
    is_fv = InFV_nohiyz(pandora_df.slc.vertex)
    two_trk_df_condition = is_fv
    is_fv = make_slc_var(is_fv)

    #### (2) not clear cosmic
    is_not_clear_cosmic = pandora_df.slc.is_clear_cosmic == 0
    two_trk_df_condition = two_trk_df_condition & is_not_clear_cosmic
    is_not_clear_cosmic = make_slc_var(is_not_clear_cosmic)

    #### (3) assign trk and proton/shower candidates
    pandora_df[('pfp', 'trk', 'chi2pid', 'I2', 'mu_over_p', '')] = pandora_df.pfp.trk.chi2pid.I2.chi2_muon / pandora_df.pfp.trk.chi2pid.I2.chi2_proton
    pandora_df[('pfp', 'trk', 'chi2pid', 'I2', 'arctan_mu_over_p', '')] = ((2 / np.pi ) * np.arctan(pandora_df.pfp.trk.chi2pid.I2.mu_over_p / 0.8)).fillna(-0.2)

    pandora_df[('pfp', 'is_prong_candidate', '', '', '', '')] = is_prong_candidate(pandora_df)
    pandora_df[('pfp', 'is_trk_candidate', '', '', '', '')] = is_trk_candidate(pandora_df)
    pandora_df[('pfp', 'is_proton_candidate', '', '', '', '')] = is_proton_candidate(pandora_df, is_recalo)
    pandora_df[('pfp', 'is_shower_candidate', '', '', '', '')] = is_shower_candidate(pandora_df)

    n_prong = (pandora_df.pfp.is_prong_candidate).groupby(level=[0,1]).sum()
    n_trk = (pandora_df.pfp.is_trk_candidate).groupby(level=[0,1]).sum()
    n_proton = (pandora_df.pfp.is_proton_candidate).groupby(level=[0,1]).sum()
    n_shower = (pandora_df.pfp.is_shower_candidate).groupby(level=[0,1]).sum()

    #### (3) - 1: plane-2 chi2pid per calo variation (MC only), see calo_variation_chi2();
    ####          chi2-based cuts (n_proton, chi2_proton) are NOT applied here, apply offline with calo_var_masks()
    n_proton_calovar = {}
    if ismc:
        calo_chi2 = calo_variation_chi2(f, ismc)
        for var in CALO_VARS:
            pandora_df[('pfp', 'trk', 'calovar', var, 'chi2_mu', '')] = calo_chi2[var]['muon'].reindex(pandora_df.index).fillna(0.)
            pandora_df[('pfp', 'trk', 'calovar', var, 'chi2_p', '')] = calo_chi2[var]['proton'].reindex(pandora_df.index).fillna(0.)
            is_p_v = is_proton_candidate_chi2mu(pandora_df, pandora_df.pfp.trk.calovar[var].chi2_mu)
            n_proton_calovar[var] = make_slc_var(is_p_v.groupby(level=[0,1]).sum())

    pass_n_trk_condition = (n_trk == 2) & (n_proton == 0) & (n_shower == 0)

    ## n_proton is chi2-dependent -> not required here, so two-track kinematics exist for every calo variation
    two_trk_df_condition = (two_trk_df_condition) & (n_prong == 2) & (n_trk == 2) & (n_shower == 0)

    n_prong = make_slc_var(n_prong)
    n_trk = make_slc_var(n_trk)
    n_proton = make_slc_var(n_proton)
    n_shower = make_slc_var(n_shower)

    tmatch_idx = pandora_df.slc.tmatch.idx
    tmatch_eff = pandora_df.slc.tmatch.eff
    tmatch_purity = pandora_df.slc.tmatch.pur
    tmatch_idx = make_slc_var(tmatch_idx)
    tmatch_eff = make_slc_var(tmatch_eff)
    tmatch_purity = make_slc_var(tmatch_purity)

    slcdf = pd.DataFrame({
        'is_fv': is_fv,
        'bcfm_score': bcfm_score,
        'nu_score': nu_score,
        'is_not_clear_cosmic': is_not_clear_cosmic,
        'n_prong': n_prong,
        'n_trk': n_trk,
        'n_proton': n_proton,
        'n_shower': n_shower,
        'tmatch_idx': tmatch_idx,
        'tmatch_eff': tmatch_eff,
        'tmatch_purity': tmatch_purity,
    })

    #### (4) dir Z cut
    longdf = pandora_df.sort_values(by=("pfp", "trk", "len", "", "", ""), ascending=False).groupby(level=[0,1]).nth(0)
    shortdf = pandora_df.sort_values(by=("pfp", "trk", "len", "", "", ""), ascending=False).groupby(level=[0,1]).nth(1)
    is_long_dirz_pass, is_short_dirz_pass = dirz_pass(longdf, shortdf)

    long_start_x = longdf.pfp.trk.start.x
    long_start_y = longdf.pfp.trk.start.y
    long_start_z = longdf.pfp.trk.start.z
    short_start_x = shortdf.pfp.trk.start.x
    short_start_y = shortdf.pfp.trk.start.y
    short_start_z = shortdf.pfp.trk.start.z

    long_end_x = longdf.pfp.trk.end.x
    long_end_y = longdf.pfp.trk.end.y
    long_end_z = longdf.pfp.trk.end.z
    short_end_x = shortdf.pfp.trk.end.x
    short_end_y = shortdf.pfp.trk.end.y
    short_end_z = shortdf.pfp.trk.end.z

    long_start_x = make_slc_var(long_start_x)
    long_start_y = make_slc_var(long_start_y)
    long_start_z = make_slc_var(long_start_z)
    short_start_x = make_slc_var(short_start_x)
    short_start_y = make_slc_var(short_start_y)
    short_start_z = make_slc_var(short_start_z)

    long_end_x = make_slc_var(long_end_x)
    long_end_y = make_slc_var(long_end_y)
    long_end_z = make_slc_var(long_end_z)
    short_end_x = make_slc_var(short_end_x)
    short_end_y = make_slc_var(short_end_y)
    short_end_z = make_slc_var(short_end_z)

    slcdf['long_start_x'] = long_start_x
    slcdf['long_start_y'] = long_start_y
    slcdf['long_start_z'] = long_start_z
    slcdf['short_start_x'] = short_start_x
    slcdf['short_start_y'] = short_start_y
    slcdf['short_start_z'] = short_start_z

    long_dirx = longdf.pfp.trk.dir.x
    long_diry = longdf.pfp.trk.dir.y
    long_dirz = longdf.pfp.trk.dir.z
    short_dirx = shortdf.pfp.trk.dir.x
    short_diry = shortdf.pfp.trk.dir.y
    short_dirz = shortdf.pfp.trk.dir.z

    long_dirx = make_slc_var(long_dirx)
    long_diry = make_slc_var(long_diry)
    long_dirz = make_slc_var(long_dirz)
    short_dirx = make_slc_var(short_dirx)
    short_diry = make_slc_var(short_diry)
    short_dirz = make_slc_var(short_dirz)
    slcdf['long_dirx'] = long_dirx
    slcdf['long_diry'] = long_diry
    slcdf['long_dirz'] = long_dirz
    slcdf['short_dirx'] = short_dirx
    slcdf['short_diry'] = short_diry
    slcdf['short_dirz'] = short_dirz

    is_long_dirz_pass = make_slc_var(is_long_dirz_pass)
    is_short_dirz_pass = make_slc_var(is_short_dirz_pass)

    slcdf['long_dirz_pass'] = is_long_dirz_pass
    slcdf['short_dirz_pass'] = is_short_dirz_pass

    #### (5) add more slc and trk variables for plotting

    ####### (5) - 1: reco vtx
    vtx_x = pandora_df.slc.vertex.x
    vtx_y = pandora_df.slc.vertex.y
    vtx_z = pandora_df.slc.vertex.z
    vtx_x = make_slc_var(vtx_x)
    vtx_y = make_slc_var(vtx_y)
    vtx_z = make_slc_var(vtx_z)
    slcdf['vtx_x'] = vtx_x
    slcdf['vtx_y'] = vtx_y
    slcdf['vtx_z'] = vtx_z

    ###### (5) - 2: trk vars
    long_trk_len = longdf.pfp.trk.len
    short_trk_len = shortdf.pfp.trk.len
    long_trk_len = make_slc_var(long_trk_len)
    short_trk_len = make_slc_var(short_trk_len)
    slcdf['long_trk_len'] = long_trk_len
    slcdf['short_trk_len'] = short_trk_len

    long_trk_pid = longdf.pfp.trk.chi2pid.I2.arctan_mu_over_p
    short_trk_pid = shortdf.pfp.trk.chi2pid.I2.arctan_mu_over_p
    long_trk_pid = make_slc_var(long_trk_pid)
    short_trk_pid = make_slc_var(short_trk_pid)
    slcdf['long_trk_pid'] = long_trk_pid
    slcdf['short_trk_pid'] = short_trk_pid

    long_trk_chi2_mu = longdf.pfp.trk.chi2pid.I2.chi2_muon
    long_trk_chi2_pro = longdf.pfp.trk.chi2pid.I2.chi2_proton
    long_trk_chi2_mu = make_slc_var(long_trk_chi2_mu)
    long_trk_chi2_pro = make_slc_var(long_trk_chi2_pro)
    slcdf['long_trk_chi2_mu'] = long_trk_chi2_mu
    slcdf['long_trk_chi2_pro'] = long_trk_chi2_pro
    if is_recalo:
        long_trk_chi2_mu_new = longdf.pfp.trk.chi2pid.I2.chi2_muon_new
        long_trk_chi2_pro_new = longdf.pfp.trk.chi2pid.I2.chi2_proton_new
        long_trk_chi2_mu_new = make_slc_var(long_trk_chi2_mu_new)
        long_trk_chi2_pro_new = make_slc_var(long_trk_chi2_pro_new)
        slcdf['long_trk_chi2_mu_new'] = long_trk_chi2_mu_new
        slcdf['long_trk_chi2_pro_new'] = long_trk_chi2_pro_new

    short_trk_chi2_mu = shortdf.pfp.trk.chi2pid.I2.chi2_muon
    short_trk_chi2_pro = shortdf.pfp.trk.chi2pid.I2.chi2_proton
    short_trk_chi2_mu = make_slc_var(short_trk_chi2_mu)
    short_trk_chi2_pro = make_slc_var(short_trk_chi2_pro)
    slcdf['short_trk_chi2_mu'] = short_trk_chi2_mu
    slcdf['short_trk_chi2_pro'] = short_trk_chi2_pro
    if is_recalo:
        short_trk_chi2_mu_new = shortdf.pfp.trk.chi2pid.I2.chi2_muon_new
        short_trk_chi2_pro_new = shortdf.pfp.trk.chi2pid.I2.chi2_proton_new
        short_trk_chi2_mu_new = make_slc_var(short_trk_chi2_mu_new)
        short_trk_chi2_pro_new = make_slc_var(short_trk_chi2_pro_new)
        slcdf['short_trk_chi2_mu_new'] = short_trk_chi2_mu_new
        slcdf['short_trk_chi2_pro_new'] = short_trk_chi2_pro_new

    ###### (5) - 3: chi2pid per calo variation
    for var, n_p_v in n_proton_calovar.items():
        slcdf['n_proton_' + var] = n_p_v
        slcdf['long_chi2_mu_' + var] = make_slc_var(longdf.pfp.trk.calovar[var].chi2_mu)
        slcdf['long_chi2_pro_' + var] = make_slc_var(longdf.pfp.trk.calovar[var].chi2_p)
        slcdf['short_chi2_mu_' + var] = make_slc_var(shortdf.pfp.trk.calovar[var].chi2_mu)
        slcdf['short_chi2_pro_' + var] = make_slc_var(shortdf.pfp.trk.calovar[var].chi2_p)

    #### (6) measure two-track variables
    two_trk_df = pandora_df[two_trk_df_condition]
    two_trk_df = two_trk_df[two_trk_df.pfp.is_trk_candidate]

    ###### (6) - 1: reco t
    if two_trk_df.empty:
        empty_index = pd.MultiIndex(
            levels=[[], []],
            codes=[[], []],
            names=['entry', 'rec.slc..index']
        )
        reco_t_series = pd.Series(dtype='float', name='reco_t', index=empty_index)
        #reco_t_series = pd.Series(dtype='float', name='reco_t')
    else:
        reco_t_series = two_trk_df.groupby(['entry', 'rec.slc..index']).apply(measure_reco_t)

    slcdf['reco_t'] = reco_t_series
    #print(reco_t_series.value_counts())

    ###### (6) - 2: two track opening angle
    if two_trk_df.empty:
        empty_index = pd.MultiIndex(
            levels=[[], []],
            codes=[[], []],
            names=['entry', 'rec.slc..index']
        )
        opening_angle_series = pd.Series(dtype='float', name='opening_angle', index=empty_index)
    else:
        opening_angle_series = two_trk_df.groupby(['entry', 'rec.slc..index']).apply(measure_opening_angle)
    slcdf['opening_angle'] = opening_angle_series

    ###### (6) - 3: two tracks are contained
    muondf = two_trk_df.sort_values(by=("pfp", "trk", "len", "", "", ""), ascending=False).groupby(level=[0,1]).nth(0)
    cpidf = two_trk_df.sort_values(by=("pfp", "trk", "len", "", "", ""), ascending=False).groupby(level=[0,1]).nth(1)
    is_muondf_contained = is_contained(muondf)
    is_cpidf_contained = is_contained(cpidf)
    is_muondf_contained = make_slc_var(is_muondf_contained)
    is_cpidf_contained = make_slc_var(is_cpidf_contained)
    slcdf['is_muon_contained'] = is_muondf_contained
    slcdf['is_cpi_contained'] = is_cpidf_contained

    #print("is_muondf_contained")
    #print(is_muondf_contained.value_counts())
    #print("is_cpidf_contained")
    #print(is_cpidf_contained.value_counts())

    ###### (6) - 4: variables for xsec
    range_p_mu = muondf.pfp.trk.rangeP.p_muon
    range_p_cpi = cpidf.pfp.trk.rangeP.p_muon
    range_p_mu = make_slc_var(range_p_mu)
    range_p_cpi = make_slc_var(range_p_cpi)
    slcdf['range_p_mu'] = range_p_mu
    slcdf['range_p_cpi'] = range_p_cpi

    #### (7) PE of the barycenter-matched flash for each slice
    total_pe = pandora_df.slc.barycenterFM.flashPEs
    total_pe = make_slc_var(total_pe)
    slcdf['total_pe'] = total_pe

    #### (8) charge in spheres with centers at vtx
    hitdf = make_trkhitdf_plane2(f)
    new_columns = pd.MultiIndex.from_tuples(
        [('pfp', 'trk', 'hit', col, '') for col in hitdf.columns]
    )
    hitdf.columns = new_columns
    pfp_vtxdist_4cm_df = pandora_df[(pandora_df.pfp.dist_to_vertex < 50.) & (pandora_df.slc.is_clear_cosmic == 0)]
    hittrk_matched_df = multicol_merge(hitdf.reset_index(), pfp_vtxdist_4cm_df.reset_index(),
                                       left_on=[('entry', '', '', '', '', ''), ('rec.slc..index', '', '', '', '', ''), ('rec.slc.reco.pfp..index', '', '', '', '', '')],
                                       right_on=[('entry', '', '', '', '', ''), ('rec.slc..index', '', '', '', '', ''), ('rec.slc.reco.pfp..index', '', '', '', '', '')],
                                       how="right")
    hittrk_matched_df = hittrk_matched_df.set_index(["entry", "rec.slc..index", "rec.slc.reco.pfp..index", "rec.slc.reco.pfp.trk.calo.2.points..index"], verify_integrity=True)
    hittrk_matched_df = multicol_add(hittrk_matched_df, dmagdf(hittrk_matched_df.slc.vertex, hittrk_matched_df.pfp.trk.hit).rename(("pfp", "trk", "hit", "dist_to_vertex", "", "")))
    hittrk_matched_df_4cm = hittrk_matched_df[hittrk_matched_df.pfp.trk.hit.dist_to_vertex < 4.]
    hittrk_matched_df_3cm = hittrk_matched_df[hittrk_matched_df.pfp.trk.hit.dist_to_vertex < 3.]
    hittrk_matched_df_2cm = hittrk_matched_df[hittrk_matched_df.pfp.trk.hit.dist_to_vertex < 2.]
    hittrk_matched_df_1cm = hittrk_matched_df[hittrk_matched_df.pfp.trk.hit.dist_to_vertex < 1.]

    sum_integ_4cm = (hittrk_matched_df_4cm.pfp.trk.hit.integral).groupby(level=[0,1]).sum().fillna(0.)
    sum_integ_3cm = (hittrk_matched_df_3cm.pfp.trk.hit.integral).groupby(level=[0,1]).sum().fillna(0.)
    sum_integ_2cm = (hittrk_matched_df_2cm.pfp.trk.hit.integral).groupby(level=[0,1]).sum().fillna(0.)
    sum_integ_1cm = (hittrk_matched_df_1cm.pfp.trk.hit.integral).groupby(level=[0,1]).sum().fillna(0.)
    sphere_3_over_4m3 = sum_integ_3cm / (sum_integ_4cm - sum_integ_3cm)

    #print(sum_integ_4cm)
    slcdf['sum_integ_4cm'] = sum_integ_4cm
    slcdf['sum_integ_3cm'] = sum_integ_3cm
    slcdf['sum_integ_2cm'] = sum_integ_2cm
    slcdf['sum_integ_1cm'] = sum_integ_1cm
    slcdf['sphere_3_over_4m3'] = sphere_3_over_4m3

    #### (9) last moment evt selection -- cuts defined in cohpi_selection.py
    slcdf = slcdf[preselection_mask(slcdf)]
    ## MC: chi2pid cuts deferred (apply per calo variation with calo_var_masks); data: full selection
    slcdf = apply_final_selection(slcdf, is_recalo, apply_chi2=not ismc)
    if ismc:
        ## pass flag per calo variation ('cv' == nominal selection, use evt[evt.pass_calovar_cv]);
        ## keep only slices selected by at least one variation, so syst weights are not attached to slices no variation selects
        pass_calovar = pd.DataFrame({'pass_calovar_' + var: calo_var_masks(slcdf, var) for var in CALO_VARS}, index=slcdf.index)
        slcdf = pd.concat([slcdf, pass_calovar], axis=1)
        slcdf = slcdf[pass_calovar.any(axis=1)]

    slcdf.columns = pd.MultiIndex.from_tuples(
        [(*col, '') if isinstance(col, tuple) else (col, '')
         for col in slcdf.columns]
    )

    if not ismc:
        return slcdf

    this_nudf = make_cohpi_nudf(f)

    valid_matches = slcdf['tmatch_idx'].dropna()
    match_pairs = pd.MultiIndex.from_arrays([
        valid_matches.index.get_level_values('entry'),
        valid_matches.astype(int)
    ])

    # 3. Filter the second DataFrame where its (entry, rec.mc.nu..index) exists in our pairs
    this_nudf_filtered = this_nudf[this_nudf.index.isin(match_pairs)]
    global all_systs_weights_df
    nudf_filtered_systs = getsyst.filter_systs_nuind(f, all_systs_weights_df, this_nudf_filtered.ind)
    this_nudf_filtered = multicol_concat(this_nudf_filtered, nudf_filtered_systs)

    slcdf = multicol_merge(slcdf.reset_index(), this_nudf_filtered.reset_index(),
                           left_on=[('entry', ''), ('tmatch_idx','')],
                           right_on=[('entry', ''), ('rec.mc.nu..index', '')],
                           how="left") ## -- save all sllices
    #print(slcdf)
    return slcdf

def make_cohpi_nudf(f):

    nudf = make_mcdf(f)
    if nudf.empty:
        cols = [
            ('true_p_mu', ''), ('true_p_pi', ''),
            ('true_cos_theta_mu', ''), ('true_cos_theta_pi', ''),
            ('true_t', ''), ('nuint_categ', '')
        ]
        empty_df = pd.DataFrame(columns=pd.MultiIndex.from_tuples(cols))
        return empty_df

    nudf["ind"] = nudf.index.get_level_values(1)
    true_t = get_true_t(nudf).fillna(999999)
    nudf['true_t'] = true_t

    is_fv = InFV_nohiyz(nudf.position)
    is_signal = Signal(nudf)
    is_cc = nudf.iscc
    genie_mode = nudf.genie_mode
    w = nudf.w
    E = nudf.E
    ind = nudf.ind

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
        'true_E_nu': E,
        'true_p_mu': mu_p_series,
        'true_p_pi': cpi_p_series,
        'true_cos_theta_mu': mu_cos_theta_series,
        'true_cos_theta_pi': cpi_cos_theta_series,
        'true_t': true_t_series,
        'genie_mode': genie_mode,
        'iscc': is_cc,
        'nuint_categ': nuint_categ,
        'ind': ind
    })
    this_nudf.columns = pd.MultiIndex.from_tuples([(col, '') for col in this_nudf.columns])

    return this_nudf

def make_cohpi_signal_nudf(f):
    this_nudf = make_cohpi_nudf(f)

    global all_systs_weights_df
    systs = getsyst.get_all_syst_df(f)
    all_systs_weights_df = systs

    signal_nudf = this_nudf[this_nudf.nuint_categ == 1]
    true_signal_systs = getsyst.filter_systs_nuind(f, all_systs_weights_df, signal_nudf.ind)

    signal_nudf = multicol_concat(signal_nudf, true_signal_systs)
    #print(signal_nudf)

    return signal_nudf
