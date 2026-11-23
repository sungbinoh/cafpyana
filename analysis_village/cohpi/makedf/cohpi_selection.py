"""Single source of truth for the CCCOH (cohpi) event selection.

Read by make_cohpidf_slc_ana.py and make_cohpidf_all_weight_signal.py.
Change cut values / definitions here only.
"""
import numpy as np
from makedf import chi2pid

# ---------------------------------------------------------------------------
# Cut values
# ---------------------------------------------------------------------------
## -- fiducial volume (no high-y region at z > z_split)
FV = dict(
    xmin=10., xmax=190.,
    zmin=10., zmax=450.,
    z_split=250.,
    yabs=190.,          # |y| < yabs for z < z_split
    ymin=-190.,         # ymin < y < ymax_highz for z > z_split
    ymax_highz=100.,
)

## -- truth signal: CC, GENIE coherent, exactly 1 mu (>27 MeV) and 1 charged pi (>30 MeV), in FV
SIGNAL = dict(
    nmu_27MeV=1,
    npi_30MeV=1,
    genie_mode=3,       # CCCOH
)

## -- per-pfp candidate definitions
CAND = dict(
    dist_to_vertex=6.,      # cm
    trk_len=4.,             # cm
    track_score=0.5,        # trk candidate: trackScore > ; shower candidate: trackScore <
    proton_chi2_mu=25.,     # proton candidate: chi2_muon >
    mupi_arctan_mu_over_p=0.4,
)

## -- final slice-level selection (applied in order a -> e)
FINAL_CUTS = dict(
    # (a) light / BCFM
    total_pe=1500.,        # on barycenterFM matched-flash PE (was event opflash sum in +-5 us) -- TODO retune
    bcfm_score=0.05,
    # (b) pfp multiplicities
    n_prong=2, n_trk=2, n_proton=0, n_shower=0,
    # (c) direction and opening angle
    long_dirz=0.7,
    short_dirz=0.5,
    opening_angle=0.5,
    # (d) chi2 proton
    long_chi2_p=100.,
    short_chi2_p=75.,
)

## -- (e) charge in spheres around vertex (not applied by default)
SPHERE_CUTS = dict(
    sum_integ_1cm=4000.,
    sum_integ_2cm=7000.,
    sum_integ_3cm=10000.,
    sum_integ_4cm=13000.,
    sphere_3_over_4m3=6.,
)
APPLY_SPHERE_CUT = False

## -- calorimetry variations (makedf/chi2pid.py); chi2 recomputed per variation, cut offline
CALO_VARS = list(chi2pid.CALO_VARIATIONS.keys())

# ---------------------------------------------------------------------------
# Fiducial volume / truth
# ---------------------------------------------------------------------------
def _pass_y(data):
    low_z = (data.z < FV['z_split']) & (np.abs(data.y) < FV['yabs'])
    high_z = (data.z > FV['z_split']) & (data.y > FV['ymin']) & (data.y < FV['ymax_highz'])
    return low_z | high_z

def InFV_nohiyz(data):
    pass_xz = (np.abs(data.x) > FV['xmin']) & (np.abs(data.x) < FV['xmax']) & (data.z > FV['zmin']) & (data.z < FV['zmax'])
    return pass_xz & _pass_y(data)

def InFV_nohiyz_trk(data):
    # same as InFV_nohiyz but without the cathode (|x| > xmin) requirement
    pass_xz = (np.abs(data.x) < FV['xmax']) & (data.z > FV['zmin']) & (data.z < FV['zmax'])
    return pass_xz & _pass_y(data)

def Signal(df):
    is_fv = InFV_nohiyz(df.position)
    ## signal def used until Nu26:
    ## is_1pi0p0n = (df.nmu_27MeV == 1) & (df.npi_30MeV == 1) & (df.np_20MeV == 0) & (df.npi0 == 0) & (df.nn_0MeV == 0)
    is_1mu1pi_th = (df.nmu_27MeV == SIGNAL['nmu_27MeV']) & (df.npi_30MeV == SIGNAL['npi_30MeV'])
    is_cccoh = (df.iscc) & (df.genie_mode == SIGNAL['genie_mode'])
    return is_fv & is_1mu1pi_th & is_cccoh

# ---------------------------------------------------------------------------
# Reco: per-pfp
# ---------------------------------------------------------------------------
def is_contained(trk_df):
    start_in = InFV_nohiyz_trk(trk_df.pfp.trk.start)
    end_in = InFV_nohiyz_trk(trk_df.pfp.trk.end)
    return start_in & end_in

def is_prong_candidate(df):
    return (df.pfp.dist_to_vertex < CAND['dist_to_vertex']) & (df.pfp.trk.len > CAND['trk_len'])

def is_trk_candidate(df):
    return is_prong_candidate(df) & (df.pfp.trackScore > CAND['track_score'])

def is_proton_candidate_chi2mu(df, chi2_mu):
    return is_trk_candidate(df) & (chi2_mu > CAND['proton_chi2_mu'])

def is_proton_candidate(df, is_recalo):
    chi2_mu = df.pfp.trk.chi2pid.I2.chi2_muon_new if is_recalo else df.pfp.trk.chi2pid.I2.chi2_muon
    return is_proton_candidate_chi2mu(df, chi2_mu)

def is_shower_candidate(df):
    return df.pfp.trackScore < CAND['track_score']

def is_muon_pion_candidate(evt_df):
    return is_trk_candidate(evt_df) & (evt_df.pfp.trk.chi2pid.I2.arctan_mu_over_p < CAND['mupi_arctan_mu_over_p'])

def is_contained_muon_pion_candidate(evt_df):
    return is_muon_pion_candidate(evt_df) & (evt_df.pfp.is_contained)

# ---------------------------------------------------------------------------
# Reco: per-slice
# ---------------------------------------------------------------------------
def dirz_pass(longdf, shortdf, cuts=FINAL_CUTS):
    """dir-z flags for the longest and second-longest pfp."""
    return longdf.pfp.trk.dir.z > cuts['long_dirz'], shortdf.pfp.trk.dir.z > cuts['short_dirz']

def preselection_mask(slcdf):
    return (slcdf.is_not_clear_cosmic) & (slcdf.is_fv)

def final_selection_masks(slcdf, is_recalo, cuts=FINAL_CUTS, apply_sphere=APPLY_SPHERE_CUT, sphere_cuts=SPHERE_CUTS, apply_chi2=True):
    """Ordered dict of per-step masks (cumulative AND gives the final selection).
    Useful for cut-flow tables in notebooks.
    apply_chi2=False drops every chi2pid-dependent cut (n_proton == 0 and the chi2_proton step),
    e.g. to apply them per calo variation offline with calo_var_masks()."""
    sfx = '_new' if is_recalo else ''
    multiplicity = (slcdf.n_prong == cuts['n_prong']) & (slcdf.n_trk == cuts['n_trk']) & (slcdf.n_shower == cuts['n_shower'])
    if apply_chi2:
        multiplicity = multiplicity & (slcdf.n_proton == cuts['n_proton'])
    masks = {
        'pe_bcfm': (slcdf.total_pe > cuts['total_pe']) & (slcdf.bcfm_score > cuts['bcfm_score']),
        'multiplicity': multiplicity,
        'dirz_open_ang': (slcdf.long_dirz > cuts['long_dirz']) & (slcdf.short_dirz > cuts['short_dirz'])
                         & (slcdf.opening_angle > cuts['opening_angle']),
    }
    if apply_chi2:
        masks['chi2_proton'] = ((slcdf['long_trk_chi2_pro' + sfx] > cuts['long_chi2_p'])
                                & (slcdf['short_trk_chi2_pro' + sfx] > cuts['short_chi2_p']))
    if apply_sphere:
        m = None
        for col, val in sphere_cuts.items():
            this = slcdf[col] < val
            m = this if m is None else (m & this)
        masks['sphere'] = m
    return masks

def final_selection_mask(slcdf, is_recalo, **kwargs):
    mask = None
    for m in final_selection_masks(slcdf, is_recalo, **kwargs).values():
        mask = m if mask is None else (mask & m)
    return mask

def apply_final_selection(slcdf, is_recalo, **kwargs):
    return slcdf[final_selection_mask(slcdf, is_recalo, **kwargs)]

def calo_var_masks(slcdf, var, cuts=FINAL_CUTS):
    """chi2pid-dependent part of the selection for one calo variation (see CALO_VARS).
    AND it with final_selection_mask(..., apply_chi2=False)."""
    return ((slcdf['n_proton_' + var] == cuts['n_proton'])
            & (slcdf['long_chi2_pro_' + var] > cuts['long_chi2_p'])
            & (slcdf['short_chi2_pro_' + var] > cuts['short_chi2_p']))
