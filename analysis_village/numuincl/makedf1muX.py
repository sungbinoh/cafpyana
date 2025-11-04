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

INCLUDE_WEIGHTS = True

def make_spine_evtdf_wgt(f,include_weights=INCLUDE_WEIGHTS, multisim_nuniv=1000, wgt_types=["bnb","genie"],prelim_cuts=False,
                        slim=False):
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
    interdf = make_spineinterdf(f,include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim)
    partdf = make_spinepartdf(f)
    #mcdf = make_mcnudf(f, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim)

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


def make_pandora_evtdf_wgt(f, include_weights=INCLUDE_WEIGHTS, multisim_nuniv=1000, wgt_types=["bnb","genie"], slim=False, 
                       trkScoreCut=False, trkDistCut=-1., cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df


def make_pandora_evtdf(f, include_weights=INCLUDE_WEIGHTS, multisim_nuniv=1000, wgt_types=["bnb","genie"], slim=False, 
                       trkScoreCut=False, trkDistCut=-1., cutClearCosmic=True, prelim_cuts=False, **trkArgs):
    # ----- sbnd or icarus? -----
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det
    if (1 == det.unique()):
        DETECTOR = "SBND"
    else:
        DETECTOR = "ICARUS"

    assert DETECTOR == "SBND"
    
    mcdf = make_mcnudf(f, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim)
    mcdf.columns = pd.MultiIndex.from_tuples([tuple(["truth"] + list(c)) for c in mcdf.columns])
    trkdf = make_trkdf(f, trkScoreCut, trackScoreCut=0.6, **trkArgs)
    slcdf = make_slcdf(f)
    hdr = make_hdrdf(f)
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
    # if cutClearCosmic:
    #     slcdf = slcdf[slcdf.slc.is_clear_cosmic==0]

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

    # Identify if track is a candidate muon
    trkdf[("pfp", "trk", "is_muon", "", "", "")] = (trkdf.pfp.trk.chi2pid.I2.chi2_muon < 18) & (trkdf.pfp.trk.chi2pid.I2.chi2_proton > 87) & (trkdf.pfp.trk.len > 32)
    trkdf.loc[:,("pfp", "trk", "is_muon", "", "", "")] = trkdf.loc[:,("pfp", "trk", "is_muon", "", "", "")].fillna(False) #fillna with False

    # ----- loose PID for candidates ----
    trkdf[("pfp", "trk", "chi2pid", "I2", "mu_over_p", "")] = np.nan
    trkdf[("pfp", "trk", "chi2pid", "I2", "mu_over_p", "")] = trkdf.pfp.trk.chi2pid.I2.chi2_muon/trkdf.pfp.trk.chi2pid.I2.chi2_proton

    # mu candidate is track pfp with smallest chi2_mu/chi2_p
    mudf = trkdf[(trkdf.pfp.trackScore > 0.6) & (trkdf.pfp.trk.is_muon) & ~(trkdf.pfp.trk.is_muon.isna())].sort_values(trkdf.pfp.index.names[:-1] + [("pfp", "trk", "chi2pid", "I2", "mu_over_p", "")]).groupby(level=[0,1]).head(1)
    mudf.columns = pd.MultiIndex.from_tuples([tuple(["mu"] + list(c)) for c in mudf.columns])
    slcdf = multicol_merge(slcdf, mudf.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")
    slcdf = multicol_merge(slcdf, hdr, left_index=True, right_index=True, how="left", validate="one_to_one")
    idx_mu = mudf.index

    # truth
    trkdf.loc[:, ("pfp","trk","truth","p","totp","")] = np.sqrt(trkdf.pfp.trk.truth.p.genp.x**2 + trkdf.pfp.trk.truth.p.genp.y**2 + trkdf.pfp.trk.truth.p.genp.z**2)
    trkdf.loc[:, ("pfp","trk","truth","p","dir","x")] = trkdf.pfp.trk.truth.p.genp.x/trkdf.pfp.trk.truth.p.totp
    trkdf.loc[:, ("pfp","trk","truth","p","dir","y")] = trkdf.pfp.trk.truth.p.genp.y/trkdf.pfp.trk.truth.p.totp
    trkdf.loc[:, ("pfp","trk","truth","p","dir","z")] = trkdf.pfp.trk.truth.p.genp.z/trkdf.pfp.trk.truth.p.totp

    # ----- apply cuts for lightweight df -----
    # if prelim_cuts:
    #     # vertex in FV
    #     slcdf = slcdf[InFV(slcdf.slc.vertex, 50, det=DETECTOR)]

    #     # neutrino cuts
    #     slcdf = slcdf[slcdf.slc.nu_score > 0.5]

    #     # require the muon
    #     mask = (~np.isnan(slcdf.mu.pfp.trk.P.p_muon))
    #     slcdf = slcdf[mask]

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
    
    return df


