"""
Make a dataframe that contains the necessary branches for a SPINE and Pandora analysis
"""

from pyanalib.pandas_helpers import *
from pyanalib.variable_calculator import *
from makedf.util import *
import pandas as pd
import numpy as np
from makedf.makedf import *
from makedf.constants import *

#My imports 
import sys
sys.path.append('analysis_village/numuincl') #relative path to my stuff
sys.path.append('/exp/sbnd/app/users/brindenc/develop/cafpyana/analysis_village/numuincl') #relative path to my stuff
from sbnd.cafclasses.slice import CAFSlice
from sbnd.cafclasses.nu import NU
from sbnd.cafclasses.binning import Binning2D
from sbnd.numu.numu_constants import *

DEFAULT_INCLUDE_WEIGHTS = False
DEFAULT_SLIM = False
DEFAULT_UPDATE_RECOMB = False
DEFAULT_FULL_CUTS = False
DEFAULT_PRELIM_CUTS = False
DEFAULT_TRK_CUTS = True
DEFAULT_VERBOSE = False
DEFAULT_ADD_STAT_UNC = False

def _resolve_flag(value, fallback):
    return fallback if value is None else value

def _set_update_recomb(value):
    global UPDATE_RECOMB, MU_KEY_SUFFIXES, CALO_KEEP_COLS
    UPDATE_RECOMB = bool(value)
    if UPDATE_RECOMB:
        MU_KEY_SUFFIXES = ["_c_cal_fracm1","_c_cal_fracp1","_alpha_embm1", "_beta_90m1", "_R_embm1","_alpha_embp1", "_beta_90p1", "_R_embp1","_alpha_emb00"]
        CALO_KEEP_COLS = CALO_KEEP_COLS = [('pfp', 'trackScore', '', '', '', ''),
            ('pfp', 'trk', 'chi2pid', 'I2', 'chi2_muon', ''),
            ('pfp', 'trk', 'chi2pid', 'I2', 'chi2_proton', ''),
            ('pfp', 'trk', 'is_muon', '', '', '')]
    else:
        MU_KEY_SUFFIXES = [""]
        CALO_KEEP_COLS = None

def set_update_recomb(value):
    _set_update_recomb(_resolve_flag(value, DEFAULT_UPDATE_RECOMB))

def apply_setting_dependencies():
    global PRELIM_CUTS
    if FULL_CUTS and not PRELIM_CUTS:
        PRELIM_CUTS = True

INCLUDE_WEIGHTS = DEFAULT_INCLUDE_WEIGHTS
SLIM = DEFAULT_SLIM
FULL_CUTS = DEFAULT_FULL_CUTS
PRELIM_CUTS = DEFAULT_PRELIM_CUTS
TRK_CUTS = DEFAULT_TRK_CUTS
VERBOSE = DEFAULT_VERBOSE
ADD_STAT_UNC = DEFAULT_ADD_STAT_UNC
_set_update_recomb(DEFAULT_UPDATE_RECOMB)
apply_setting_dependencies()

def make_spine_evtdf_wgt(f,include_weights=None, wgt_types=["bnb","genie"],prelim_cuts=None,
                        slim=None):
    include_weights = _resolve_flag(include_weights, INCLUDE_WEIGHTS)
    prelim_cuts = _resolve_flag(prelim_cuts, PRELIM_CUTS)
    slim = _resolve_flag(slim, SLIM)
    multisim_nuniv = 100 if slim else 1000
    # ----- sbnd or icarus? -----
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det
    if (1 == det.unique()):
        DETECTOR = "SBND"
    else:
        DETECTOR = "ICARUS"

    assert DETECTOR == "SBND"

    # load header for run, subrun, and event
    hdrdf = make_hdrdf(f)
    hdrdf = hdrdf.loc[:,["run","subrun","evt"]]
    
    # load slices and particles
    interdf = make_spineinterdf(f,include_weights=include_weights,  wgt_types=wgt_types, slim=slim, multisim_nuniv=multisim_nuniv)
    partdf = make_spinepartdf(f)

    # load the muon candidate
    mudf = partdf[partdf.is_primary & (partdf.pid == 2)].sort_values(partdf.index.names[:2] + [("ke", "", "")]).groupby(level=[0,1]).last()
    mudf.columns = pd.MultiIndex.from_tuples([tuple(["mu"] + list(c)) for c in mudf.columns])

    # Add costheta
    mudf[("mu", "costheta", "", "", "", "")] = mudf.mu.start_dir.z
    mudf[("mu", "tpart", "costheta", "", "", "")] = mudf.mu.tpart.start_dir.I2 #TODO: Fix upstream replacement for I2->z

    #Merge the muon to the interaction
    interdf = multicol_merge(interdf, mudf, left_index=True, right_index=True, how="left", validate="one_to_one")
    
    # Merge the header to the interaction
    interdf = multicol_merge(interdf, hdrdf, left_index=True, right_index=True, how="left", validate="one_to_one")

    if prelim_cuts:
        # Apply pre-selection: Require fiducial vertex, at least one muon, at least one proton
        # require both muon to be present
        interdf = interdf[~np.isnan(interdf.mu.pid)]
        # require fiducial verex
        interdf = interdf[InFV(interdf.vertex, 50,det="SBND")]
    
    # Fix position names (I0, I1, I2) -> (x, y, z) at any depth
    def mappos(s):
        if s == "I0": return "x"
        if s == "I1": return "y"
        if s == "I2": return "z"
        return s
    
    def fixpos_recursive(tuple_or_string):
        """Recursively fix position names at any depth in the MultiIndex"""
        if isinstance(tuple_or_string, tuple):
            return tuple(fixpos_recursive(item) for item in tuple_or_string)
        else:
            return mappos(tuple_or_string)
    
    def fixpos(c):
        """Only apply position fixes to columns starting with specific prefixes"""
        prekeys = ["end_point", "start_point", "start_dir", "vertex", "momentum", "end_dir"]
        find_prekey = False
        for sc in c:
            if sc in prekeys:
                find_prekey = True
                break
        if find_prekey:
            return fixpos_recursive(c)
        else:
            return c

    interdf.columns = pd.MultiIndex.from_tuples([fixpos(c) for c in interdf.columns])
    #Add containment
    #interdf[("mu", "is_contained", "", "")] = (InFV(interdf.mu.start_point, 0, det=DETECTOR+" AV")) & (InFV(interdf.mu.end_point, 0, det=DETECTOR+" AV"))
    #interdf[("mu", "truth", "is_contained", "")] = (InFV(interdf.mu.truth.start_point, 0, det=DETECTOR+" AV")) & (InFV(interdf.mu.truth.end_point, 0, det=DETECTOR+" AV"))

    # mcdf.columns = pd.MultiIndex.from_tuples([tuple(list(c) +["", "", "", "", ""]) for c in mcdf.columns])     # match # of column levels

    # df = multicol_merge(interdf.reset_index(), 
    #               mcdf.reset_index(),
    #               left_on=[("entry", "", "",), 
    #                        ("slc", "tmatch", "idx")], 
    #               right_on=[("entry", "", ""), 
    #                         ("rec.mc.nu..index", "", "")], 
    #               how="left"
    #               )
    #interdf = interdf.set_index(interdf.index.names, verify_integrity=True)
    
    return interdf

def make_custom_trkdf(f, trkScoreCut=None, updaterecomb=None, **trkArgs):
    trkScoreCut = _resolve_flag(trkScoreCut, TRK_CUTS)
    updaterecomb = _resolve_flag(updaterecomb, UPDATE_RECOMB)
    trkdf = make_trkdf(f, scoreCut=trkScoreCut, updaterecomb=updaterecomb, **trkArgs)
    return trkdf

def make_pandora_evtdf(f, include_weights=None, wgt_types=["bnb","genie","g4"], slim=None,
                       trkScoreCut=None, trkDistCut=-1., cutClearCosmic=True, updaterecomb=None,
                       return_mcdf=False, **trkArgs):
    include_weights = _resolve_flag(include_weights, INCLUDE_WEIGHTS)
    slim = _resolve_flag(slim, SLIM)
    multisim_nuniv = 100 if slim else 1000
    updaterecomb = _resolve_flag(updaterecomb, UPDATE_RECOMB)
    trkScoreCut = _resolve_flag(trkScoreCut, TRK_CUTS)
    # ----- sbnd or icarus? -----
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det
    if (1 == det.unique()):
        DETECTOR = "SBND"
    else:
        DETECTOR = "ICARUS"

    assert DETECTOR == "SBND"
    
    mcdf = make_mcnudf(f, include_weights=include_weights,  wgt_types=wgt_types, slim=slim, multisim_nuniv=multisim_nuniv)
    if return_mcdf: ret_mcdf = mcdf.copy() #Make a copy to return
    mcdf.columns = pd.MultiIndex.from_tuples([tuple(["truth"] + list(c)) for c in mcdf.columns])
    trkdf = make_custom_trkdf(f, trkScoreCut=trkScoreCut, updaterecomb=updaterecomb, **trkArgs)
    slcdf = make_slcdf(f)
    hdr = make_hdrdf(f)
    ismc = hdr.ismc.astype(bool).unique()
    if len(ismc) > 1:
        raise ValueError(f'Multiple ismc values: {ismc}')
    else:
        ismc = ismc[0]
    hdr = hdr.loc[:,["run","subrun","evt"]]


    # stubdf = make_stubs(f, det=DETECTOR)
    # load stubs
    # slcdf = multicol_merge(slcdf, stubdf, left_index=True, right_index=True)
    
    # ----- merge dfs -----
    # load pfps
    # slcdf = multicol_merge(slcdf, trkdf, left_index=True, right_index=True, how="right", validate="one_to_many")

    trkdf = multicol_add(trkdf, dmagdf(slcdf.slc.vertex, trkdf.pfp.trk.start).rename(("pfp", "dist_to_vertex")))
    # if trkDistCut > 0:
    #     trkdf = trkdf[trkdf.pfp.dist_to_vertex < trkDistCut]

    # ---- calculate additional info ----
    
    # track containment in active volume
    trkdf[("pfp", "trk", "is_contained", "", "", "")] = (InFV(trkdf.pfp.trk.start, 0, det=DETECTOR+" AV")) & (InFV(trkdf.pfp.trk.end, 0, det=DETECTOR+" AV"))

    # reco momentum -- range for contained, MCS for exiting
    trkdf[("pfp", "trk", "P", "p_muon", "", "")] = np.nan
    trkdf.loc[trkdf.pfp.trk.is_contained, ("pfp", "trk", "P", "p_muon", "", "")]  = trkdf.loc[(trkdf.pfp.trk.is_contained), ("pfp", "trk", "rangeP", "p_muon", "", "")]
    trkdf.loc[np.invert(trkdf.pfp.trk.is_contained), ("pfp", "trk", "P", "p_muon","", "")] = trkdf.loc[np.invert(trkdf.pfp.trk.is_contained), ("pfp", "trk", "mcsP", "fwdP_muon", "", "")]

    trkdf[("pfp", "trk", "P", "p_pion", "", "")] = np.nan
    trkdf.loc[trkdf.pfp.trk.is_contained, ("pfp", "trk", "P", "p_pion", "", "")]  = trkdf.loc[(trkdf.pfp.trk.is_contained), ("pfp", "trk", "rangeP", "p_pion", "", "")]
    trkdf.loc[np.invert(trkdf.pfp.trk.is_contained), ("pfp", "trk", "P", "p_pion", "", "")] = trkdf.loc[np.invert(trkdf.pfp.trk.is_contained), ("pfp", "trk", "mcsP", "fwdP_pion", "", "")]

    trkdf[("pfp", "trk", "P", "p_proton", "", "")] = np.nan
    trkdf.loc[trkdf.pfp.trk.is_contained, ("pfp", "trk", "P", "p_proton", "", "")]  = trkdf.loc[(trkdf.pfp.trk.is_contained), ("pfp", "trk", "rangeP", "p_proton", "", "")]
    trkdf.loc[np.invert(trkdf.pfp.trk.is_contained), ("pfp", "trk", "P", "p_proton", "", "")] = trkdf.loc[np.invert(trkdf.pfp.trk.is_contained), ("pfp", "trk", "mcsP", "fwdP_proton", "", "")]

    # opening angles
    trkdf[("pfp", "trk", "dir", "x", "", "")] = np.nan
    trkdf[("pfp", "trk", "dir", "x", "", "")] = (trkdf.pfp.trk.end.x-trkdf.pfp.trk.start.x)/trkdf.pfp.trk.len
    trkdf[("pfp", "trk", "dir", "y", "", "")] = np.nan
    trkdf[("pfp", "trk", "dir", "y", "", "")] = (trkdf.pfp.trk.end.y-trkdf.pfp.trk.start.y)/trkdf.pfp.trk.len
    trkdf[("pfp", "trk", "dir", "z", "", "")] = np.nan
    trkdf[("pfp", "trk", "dir", "z", "", "")] = (trkdf.pfp.trk.end.z-trkdf.pfp.trk.start.z)/trkdf.pfp.trk.len

    # costheta
    trkdf[("pfp", "trk", "costheta", "", "", "")] = trkdf.pfp.trk.dir.z

    # truth totp and dirz
    if ismc:
        trkdf.loc[:, ("pfp","trk","truth","p","totp","")] = np.sqrt(trkdf.pfp.trk.truth.p.genp.x**2 + trkdf.pfp.trk.truth.p.genp.y**2 + trkdf.pfp.trk.truth.p.genp.z**2)
        trkdf.loc[:, ("pfp","trk","truth","p","dir","x")] = trkdf.pfp.trk.truth.p.genp.x/trkdf.pfp.trk.truth.p.totp
        trkdf.loc[:, ("pfp","trk","truth","p","dir","y")] = trkdf.pfp.trk.truth.p.genp.y/trkdf.pfp.trk.truth.p.totp
        trkdf.loc[:, ("pfp","trk","truth","p","dir","z")] = trkdf.pfp.trk.truth.p.genp.z/trkdf.pfp.trk.truth.p.totp
    # Save the keys in the trkdf, with the nan min and max
    # with open('/exp/sbnd/app/users/brindenc/develop/cafpyana/analysis_village/numuincl/trkdf_keys_inmaker.txt','w') as f:
    #     for k in trkdf.columns:
    #         f.write(f'{k}: (')
    #         f.write(f'min = {trkdf[k].min()}')
    #         f.write(f', max = {trkdf[k].max()})\n')

    mudf_list = []
    for key_suffix in MU_KEY_SUFFIXES:
        trkdf[("pfp", "trk", "is_muon", "", "", "")] = (trkdf[("pfp", "trk", "chi2pid", "I2", f"chi2_muon{key_suffix}","")] < 20) & (trkdf[("pfp", "trk", "chi2pid", "I2", f"chi2_proton{key_suffix}","")] > 90) & (trkdf.pfp.trk.len > 32)
        trkdf.loc[:,("pfp", "trk", "is_muon", "", "", "")] = trkdf.loc[:,("pfp", "trk", "is_muon", "", "", "")].fillna(False) #fillna with False

        # ----- loose PID for candidates ----
        trkdf[("pfp", "trk", "chi2pid", "I2", f"mu_over_p{key_suffix}", "")] = np.nan
        trkdf[("pfp", "trk", "chi2pid", "I2", f"mu_over_p{key_suffix}", "")] = trkdf[("pfp", "trk", "chi2pid", "I2", f"chi2_muon{key_suffix}","")]/trkdf[("pfp", "trk", "chi2pid", "I2", f"chi2_proton{key_suffix}","")]

        # mu candidate is track pfp with smallest chi2_mu/chi2_p
        mudf = trkdf[(trkdf.pfp.trackScore > 0.6) & (trkdf[("pfp", "trk", "is_muon", "", "", "")] ) & ~(trkdf[("pfp", "trk", "is_muon", "", "", "")] .isna())].sort_values(trkdf.pfp.index.names[:-1] + [("pfp", "trk", "chi2pid", "I2", f"mu_over_p{key_suffix}", "")]).groupby(level=[0,1]).head(1)
        
        # Filter columns to only include those matching the current key_suffix or base columns (no variation suffix)
        # This prevents mixing variation suffixes (e.g., mu_R_embm1 containing chi2_proton_alpha_embp1)
        if key_suffix:  # Only filter if we have a variation suffix
            # Get all variation suffixes
            all_suffixes = [s for s in MU_KEY_SUFFIXES if s]
            # Filter columns: keep base columns (no variation suffix) or columns with current key_suffix
            def should_keep_col(col_tuple, contains_keep=True):
                col_str = str(col_tuple)
                # Check if column contains any variation suffix
                contains_any_suffix = any(suffix in col_str for suffix in all_suffixes)
                # Check if the column is in the CALO_KEEP_COLS list
                if (CALO_KEEP_COLS is not None) and contains_keep:
                    contains_keep_key = False
                    # for kk in CALO_KEEP_COLS:
                    #     assert len(kk) == len(col_tuple), f'Column {col_tuple} has length {len(col_tuple[1:])} but CALO_KEEP_COLS has length {len(kk)}'
                    if col_tuple in CALO_KEEP_COLS:
                        contains_keep_key = True
                else:
                    contains_keep_key = True
                
                if not contains_any_suffix and contains_keep_key:
                    # Base column, keep it
                    return True
                # Column has a variation suffix, only keep if it matches current key_suffix
                return (key_suffix in col_str) and contains_keep_key
            # Only apply contains keep for non-null variations
            contains_keep = key_suffix != MU_KEY_SUFFIXES[-1]
            cols_to_keep = [c for c in mudf.columns if should_keep_col(c, contains_keep=contains_keep)]
            mudf = mudf[cols_to_keep]
            if contains_keep:
                assert len(mudf.columns) == len(CALO_KEEP_COLS), f'Number of columns in mudf {len(mudf.columns)} does not match number of columns in CALO_KEEP_COLS {len(CALO_KEEP_COLS)}'
        
        mudf.columns = pd.MultiIndex.from_tuples([tuple(["mu" + key_suffix] + list(c)) for c in mudf.columns])
        mudf_list.append(mudf)

    slcdf = multicol_merge(slcdf, pd.concat(mudf_list, axis=1).droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")
    slcdf = multicol_merge(slcdf, hdr, left_index=True, right_index=True, how="left", validate="one_to_one")

    # ----- apply cuts for lightweight df -----
#     if prelim_cuts:
#         # vertex in FV
#         slcdf = slcdf[InFV(slcdf.slc.vertex, 50, det=DETECTOR)]

#         # neutrino cuts
#         slcdf = slcdf[slcdf.slc.nu_score > 0.5]
#         slcdf = slcdf[slcdf.slc.is_clear_cosmic==0]

#         if full_cuts:

#             # require the muon
#             mask = slcdf[(f'mu{key_suffix}', 'pfp', 'trk', 'is_muon', '', '', '')].astype(bool)
#             slcdf = slcdf[mask]

#             # low z cut
#             mask = (slcdf[(f'mu{key_suffix}', 'pfp', 'trk', 'end', 'z', '', '')].astype(float) > 6)
#             slcdf = slcdf[mask]

#             # require opt0 score
#             mask = (slcdf.slc.opt0.score > 320)
#             slcdf = slcdf[mask]

#             # Find slice with highest opt0 score
#             slcdf = (
#                 slcdf.sort_values([("slc", "opt0", "score")], ascending=False)
#                     .groupby(level=slcdf.index.names[:-1])
#                     .first()
# )

    # ---- truth match ----
    bad_tmatch = np.invert(slcdf.slc.tmatch.eff > 0.5) & (slcdf.slc.tmatch.idx >= 0)
    slcdf.loc[bad_tmatch, ("slc","tmatch","idx", "", "")] = np.nan

    mcdf.columns = pd.MultiIndex.from_tuples([tuple(list(c) +["", "", "", "", ""]) for c in mcdf.columns])     # match # of column levels

    df = multicol_merge(slcdf.reset_index(), 
                  mcdf.reset_index(),
                  left_on=[("entry", "", "",), 
                           ("slc", "tmatch", "idx")], 
                  right_on=[("entry", "", ""), 
                            ("rec.mc.nu..index", "", "")], 
                  how="left"
                  )

    df = df.set_index(slcdf.index.names, verify_integrity=True)
    
    if return_mcdf:
        return df, ret_mcdf
    else:
        return df

def _process_mcnu(mcdf, slc, ismc):
    """
    Helper function to process mcnu dataframe with the same steps used in make_pandora_evtdf_processed.
    Returns the processed mcnu NU object, or None if not MC or empty.
    """
    if ismc and len(mcdf.index.values) != 0:
        mcnu = NU(mcdf)
        mcnu.scale_to_pot(nom_pot=1., sample_pot=1.)

        #Process mcnu
        mcnu.add_fv()
        mcnu.add_av()
        
        mcnu = slc.set_mcnu_containment(mcnu)
        mcnu.cut_muon(cut=False, min_ke=0.1)
        mcnu.cut_fv(cut=False)
        mcnu.cut_cosmic(cut=False)
        mcnu.cut_cont(cut=False)
        mcnu.add_event_type('pandora',)
        
        return mcnu
    return None

def make_pandora_evtdf_processed(f, include_weights=None, wgt_types=["bnb","genie","g4"], slim=None, 
                       trkScoreCut=None, updaterecomb=None, add_stat_unc=None, **trkArgs):
    """
    Utilize my CAF class to add the necessary columns to the dataframe
    """
    include_weights = _resolve_flag(include_weights, INCLUDE_WEIGHTS)
    slim = _resolve_flag(slim, SLIM)
    updaterecomb = _resolve_flag(updaterecomb, UPDATE_RECOMB)
    trkScoreCut = _resolve_flag(trkScoreCut, TRK_CUTS)
    add_stat_unc = _resolve_flag(add_stat_unc, ADD_STAT_UNC)
    df, mcdf = make_pandora_evtdf(f, include_weights=include_weights,  wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, updaterecomb=updaterecomb, return_mcdf=True, **trkArgs)
    hdr = make_hdrdf(f)
    ismc = hdr.ismc.astype(bool).unique()
    if len(ismc) > 1:
        raise ValueError(f'Multiple ismc values: {ismc}')
    else:
        ismc = ismc[0]
    slc = CAFSlice(df) 
    slc.remove_column_suffix(MU_KEY_SUFFIXES[-1]) #Fix null variation
    # with open('/exp/sbnd/app/users/brindenc/develop/cafpyana/analysis_village/numuincl/slc_keys_inmaker.txt','w') as f:
    #     for k in slc.data.keys():
    #         f.write(f'{k}\n')
    #Scale to dummy pot
    slc.scale_to_pot(nom_pot=1.,sample_pot=1.)
    if VERBOSE:
        #print(slc.data.keys())
        #print(slc.data.index)
        qqq = 1
    
    # Process mcnu using the helper function
    mcnu = _process_mcnu(mcdf, slc, ismc)

    slc.clean(dummy_vals=[-9999,-999,999,9999,-5])

    # Need this for calo variations
    for i,key_suffix in enumerate(MU_KEY_SUFFIXES):
        if key_suffix != MU_KEY_SUFFIXES[-1]:
            slc.add_has_muon(suffix=key_suffix)
        else:
            slc.add_has_muon()
    slc.add_in_av()
    slc.add_in_fv()
    slc.add_event_type()
    slc.add_track_flipping()

    if ismc:
        # Fix the costheta and momentum for slices that don't have a true muon
        mask = (slc.data.truth.mu.dir.z == -1) | (np.isnan(slc.data.truth.mu.dir.z))
        dir_col = slc.get_key([f'mu.pfp.trk.truth.p.dir.z'])[0]
        truth_dir_col = slc.get_key([f'truth.mu.dir.z'])[0]
        if mask.sum() > 0:
            slc.data.loc[mask,truth_dir_col] = slc.data.loc[mask,dir_col]

        mask = (slc.data.truth.mu.totp == -1) | (np.isnan(slc.data.truth.mu.totp))
        totp_col = slc.get_key([f'mu.pfp.trk.truth.p.totp'])[0]
        truth_totp_col = slc.get_key([f'truth.mu.totp'])[0]
        if mask.sum() > 0:
            slc.data.loc[mask,truth_totp_col] = slc.data.loc[mask,totp_col]
    slc.add_2d_binning(include_truth=ismc, include_reco=True)


    #Opt0 cuts
    #slc.cut_flashmatch(cut=False)
    #slc.cut_cosmic(cut=False,fmatch_score=320,nu_score=0.5,use_opt0=True,use_isclearcosmic=False)
    #Barycenter FM cuts
    #Prescale flash PE
    if ismc:
        pe_col = slc.get_key([f'slc.barycenterFM.flashPEs'])[0]
        slc.data.loc[:,pe_col] = slc.data.loc[:,pe_col]*0.66 #prescale to match the data
    #Add cuts
    slc.cut_flashpe(cut=False,min_flashpe=2000,prescale=1.)
    slc.cut_cosmic(cut=False,fmatch_score=0.06,use_opt0='barycenterFM',nu_score=None,use_isclearcosmic=False)
    slc.cut_flashmatch(cut=False,method='barycenterFM',use_isclearcosmic=False)
    slc.cut_fv(cut=False)
    slc.cut_muon(cut=False,min_ke=0.1)
    slc.cut_lowz(cut=False,z_max=6,include_start=True)
    slc.cut_is_cont(cut=False) #Don't apply containment cut
    slc.cut_all(cut=False) #Add the all cut column

    if add_stat_unc:
        slc.add_stat_unc()
    

    if VERBOSE:
        from naming import PAND_CUTS
        if mcnu is not None:
            pur,eff,f1,_,_,_ = slc.get_pur_eff_f1(mcnu,PAND_CUTS,categories=[0,1])
            print('Pandora cuts:')
            print(PAND_CUTS)
            print('Pandora pur, eff, f1:')
            print(pur,eff,f1)
    with open('/exp/sbnd/app/users/brindenc/develop/cafpyana/analysis_village/numuincl/slc_keys_inmaker.txt','w') as f:
        for k in slc.data.keys():
            f.write(f'{k}\n')
    return slc.data
    
def make_pandora_evtdf_processed_signal_cut(f, include_weights=None, wgt_types=["bnb","genie","g4"], slim=None, 
    trkScoreCut=None, updaterecomb=None, **trkArgs):
    """
    Utilize my CAF class to add the necessary columns to the dataframe
    """
    include_weights = _resolve_flag(include_weights, INCLUDE_WEIGHTS)
    slim = _resolve_flag(slim, SLIM)
    multisim_nuniv = 100 if slim else 1000
    updaterecomb = _resolve_flag(updaterecomb, UPDATE_RECOMB)
    df = make_pandora_evtdf_processed(f, include_weights=include_weights,  wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, updaterecomb=updaterecomb, **trkArgs)

    if VERBOSE:
        print(f'include weights: {include_weights}')
        print(f'slim: {slim}')
        print(f'updaterecomb: {updaterecomb}')
    slc = CAFSlice(df)
    slc.remove_column_suffix(MU_KEY_SUFFIXES[-1]) #Fix null variation
    is_signal = np.isin(slc.data.truth.event_type,[0,1])
    slc.data = slc.data[is_signal]

    #Apply the cut all for both the reco and truth dataframes (only cut if requested)
    slc.cut_all(cut=True,mode='truth',cont=False)

    with open('/exp/sbnd/app/users/brindenc/develop/cafpyana/analysis_village/numuincl/slc_signal_keys_inmaker.txt','w') as f:
        for k in slc.data.keys():
            f.write(f'{k}\n')

    return slc.data #return just the df

def make_pandora_evtdf_processed_selected_cut(f, include_weights=None, wgt_types=["bnb","genie","g4"], slim=None, 
    trkScoreCut=None, updaterecomb=None, **trkArgs):
    """
    Utilize my CAF class to add the necessary columns to the dataframe
    """
    include_weights = _resolve_flag(include_weights, INCLUDE_WEIGHTS)
    slim = _resolve_flag(slim, SLIM)
    multisim_nuniv = 100 if slim else 1000
    updaterecomb = _resolve_flag(updaterecomb, UPDATE_RECOMB)
    df = make_pandora_evtdf_processed(f, include_weights=include_weights,  wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, updaterecomb=updaterecomb, **trkArgs)

    if VERBOSE:
        print(f'include weights: {include_weights}')
        print(f'slim: {slim}')
        print(f'updaterecomb: {updaterecomb}')
    slc = CAFSlice(df)
    slc.remove_column_suffix(MU_KEY_SUFFIXES[-1]) #Fix null variation
    slc.cut_all(cut=True,mode='reco',cont=False)
    with open('/exp/sbnd/app/users/brindenc/develop/cafpyana/analysis_village/numuincl/slc_selected_keys_inmaker.txt','w') as f:
        for k in slc.data.keys():
            f.write(f'{k}\n')
    return slc.data #return just the df
