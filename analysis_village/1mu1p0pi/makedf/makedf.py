from pyanalib.pandas_helpers import *
from pyanalib.variable_calculator import *
from makedf.util import *
import pandas as pd
import numpy as np
from makedf.makedf import *
from makedf.constants import *

def make_spine_evtdf(f):
    # load slices and particles
    partdf = make_epartdf(f)

    df = make_eslcdf(f)

    # load the proton and muon candidates
    primary = partdf.is_primary
    mudf = partdf[primary & (partdf.pid == 2)].sort_values(partdf.index.names[:2] + [("length", "", "")]).groupby(level=[0,1]).last()
    mudf.columns = pd.MultiIndex.from_tuples([tuple(["mu"] + list(c)) for c in mudf.columns])

    pdf = partdf[primary & (partdf.pid == 4)].sort_values(partdf.index.names[:2] + [("length", "", "")]).groupby(level=[0,1]).last()
    pdf.columns = pd.MultiIndex.from_tuples([tuple(["p"] + list(c)) for c in pdf.columns])

    df = multicol_merge(df, mudf, left_index=True, right_index=True, how="left", validate="one_to_one")
    df = multicol_merge(df, pdf, left_index=True, right_index=True, how="left", validate="one_to_one")

    # in case we want to cut out other objects -- save the highest energy of each other particle
    lead_gamma_energy = partdf.ke[primary & (partdf.pid == 0)].groupby(level=[0,1]).max().rename("lead_gamma_energy")
    df = multicol_add(df, lead_gamma_energy)

    lead_elec_energy = partdf.ke[primary & (partdf.pid == 1)].groupby(level=[0,1]).max().rename("lead_elec_energy")
    df = multicol_add(df, lead_elec_energy)

    lead_pion_length = partdf.length[primary & (partdf.pid == 3)].groupby(level=[0,1]).max().rename("lead_pion_length")
    df = multicol_add(df, lead_pion_length)

    subl_muon_length = partdf[primary & (partdf.pid == 2)].sort_values(partdf.index.names[:2] + [("length", "", "")]).length.groupby(level=[0,1]).nth(-2).rename("subl_muon_length")
    df = multicol_add(df, subl_muon_length)

    subl_proton_length = partdf[primary & (partdf.pid == 4)].sort_values(partdf.index.names[:2] + [("length", "", "")]).length.groupby(level=[0,1]).nth(-2).rename("subl_proton_length")
    df = multicol_add(df, subl_proton_length)

    # Apply pre-selection: Require fiducial vertex, at least one muon, at least one proton

    # require both muon and proton to be present
    df = df[~np.isnan(df.mu.pid) & ~np.isnan(df.p.pid)]

    # require fiducial verex
    df = df[InFV(df.vertex, 50)]

    return df


def make_pandora_evtdf_wgt(f, include_weights=True, multisim_nuniv=1000, wgt_types=["bnb","genie"], slim=True, 
                       trkScoreCut=False, trkDistCut=10., cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df


def make_pandora_evtdf(f, include_weights=True, multisim_nuniv=1000, wgt_types=["bnb","genie"], slim=True, 
                       trkScoreCut=False, trkDistCut=10., cutClearCosmic=True, **trkArgs):
    # ----- sbnd or icarus? -----
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det
    if (1 == det.unique()):
        DETECTOR = "SBND"
    else:
        DETECTOR = "ICARUS"

    assert DETECTOR == "SBND"
    
    mcdf = make_mcnudf(f, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim)
    trkdf = make_trkdf(f, trkScoreCut, **trkArgs)
    slcdf = make_slcdf(f)

    # stubdf = make_stubs(f, det=DETECTOR)
    # load stubs
    # slcdf = multicol_merge(slcdf, stubdf, left_index=True, right_index=True)
    
    # ----- merge dfs -----
    # load pfps
    # slcdf = multicol_merge(slcdf, trkdf, left_index=True, right_index=True, how="right", validate="one_to_many")

    trkdf = multicol_add(trkdf, dmagdf(slcdf.slc.vertex, trkdf.pfp.trk.start).rename(("pfp", "dist_to_vertex")))
    if trkDistCut > 0:
        trkdf = trkdf[trkdf.pfp.dist_to_vertex < trkDistCut]
    if cutClearCosmic:
        slcdf = slcdf[slcdf.slc.is_clear_cosmic==0]

    # ---- calculate additional info ----
    
    # track containment
    trkdf[("pfp", "trk", "is_contained", "", "", "")] = (InFV(trkdf.pfp.trk.start, 0, det=DETECTOR)) & (InFV(trkdf.pfp.trk.end, 0, det=DETECTOR))

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

    # truth
    trkdf.loc[:, ("pfp","trk","truth","p","totp","")] = np.sqrt(trkdf.pfp.trk.truth.p.genp.x**2 + trkdf.pfp.trk.truth.p.genp.y**2 + trkdf.pfp.trk.truth.p.genp.z**2)
    trkdf.loc[:, ("pfp","trk","truth","p","dir","x")] = trkdf.pfp.trk.truth.p.genp.x/trkdf.pfp.trk.truth.p.totp
    trkdf.loc[:, ("pfp","trk","truth","p","dir","y")] = trkdf.pfp.trk.truth.p.genp.y/trkdf.pfp.trk.truth.p.totp
    trkdf.loc[:, ("pfp","trk","truth","p","dir","z")] = trkdf.pfp.trk.truth.p.genp.z/trkdf.pfp.trk.truth.p.totp

    # ----- loose PID for candidates ----
    trkdf[("pfp", "trk", "chi2pid", "I2", "mu_over_p", "")] = np.nan
    trkdf[("pfp", "trk", "chi2pid", "I2", "mu_over_p", "")] = trkdf.pfp.trk.chi2pid.I2.chi2_muon/trkdf.pfp.trk.chi2pid.I2.chi2_proton
    
    # mu candidate is track pfp with smallest chi2_mu/chi2_p
    mudf = trkdf[(trkdf.pfp.trackScore > 0.5)].sort_values(trkdf.pfp.index.names[:-1] + [("pfp", "trk", "chi2pid", "I2", "mu_over_p", "")]).groupby(level=[0,1]).head(1)
    mudf.columns = pd.MultiIndex.from_tuples([tuple(["mu"] + list(c)) for c in mudf.columns])
    slcdf = multicol_merge(slcdf, mudf.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")
    idx_mu = mudf.index
        
    # p candidate is track pfp with largest chi2_mu/chi2_p of remaining pfps
    idx_pfps = trkdf.index
    idx_not_mu = idx_pfps.difference(idx_mu)
    notmudf = trkdf.loc[idx_not_mu]
    pdf = notmudf[(notmudf.pfp.trackScore > 0.5)].sort_values(notmudf.pfp.index.names[:-1] + [("pfp", "trk", "chi2pid", "I2", "mu_over_p", "")]).groupby(level=[0,1]).tail(1)
    pdf.columns = pd.MultiIndex.from_tuples([tuple(["p"] + list(c)) for c in pdf.columns])
    slcdf = multicol_merge(slcdf, pdf.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")
    idx_p = pdf.index
    
    # note if there are any other track/showers
    idx_not_mu_p = idx_not_mu.difference(idx_p)
    otherdf = trkdf.loc[idx_not_mu_p]
    # longest other shower
    othershwdf = otherdf[otherdf.pfp.trackScore < 0.5]
    other_shw_length = othershwdf.pfp.trk.len.groupby(level=[0,1]).max().rename("other_shw_length")
    slcdf = multicol_add(slcdf, other_shw_length)
    # longest other track
    othertrkdf = otherdf[otherdf.pfp.trackScore > 0.5]
    other_trk_length = othertrkdf.pfp.trk.len.groupby(level=[0,1]).max().rename("other_trk_length")
    slcdf = multicol_add(slcdf, other_trk_length)

    # calculate and save transverse kinematic variables for reco slices
    slc_mudf = slcdf.mu.pfp.trk
    slc_pdf = slcdf.p.pfp.trk
    slc_P_mu_col = pad_column_name(("P", "p_muon"), slc_mudf)
    slc_P_p_col = pad_column_name(("P", "p_proton"), slc_pdf)
    tki_reco = get_cc1p0pi_tki(slc_mudf, slc_pdf, slc_P_mu_col, slc_P_p_col)

    slcdf = multicol_add(slcdf, tki_reco["del_alpha"].rename("del_alpha"))
    slcdf = multicol_add(slcdf, tki_reco["del_phi"].rename("del_phi"))
    slcdf = multicol_add(slcdf, tki_reco["del_Tp"].rename("del_Tp"))
    slcdf = multicol_add(slcdf, tki_reco["del_p"].rename("del_p"))

    # calculate and save transverse kinematic variables for MC
    mc_mudf = mcdf.mu
    mc_pdf = mcdf.p
    mc_P_mu_col = pad_column_name(("totp",), mc_mudf)
    mc_P_p_col = pad_column_name(("totp",), mc_pdf)
    tki_mc = get_cc1p0pi_tki(mc_mudf, mc_pdf, mc_P_mu_col, mc_P_p_col)

    mcdf = multicol_add(mcdf, tki_mc["del_alpha"].rename("mc_del_alpha"))
    mcdf = multicol_add(mcdf, tki_mc["del_phi"].rename("mc_del_phi"))
    mcdf = multicol_add(mcdf, tki_mc["del_Tp"].rename("mc_del_Tp"))
    mcdf = multicol_add(mcdf, tki_mc["del_p"].rename("mc_del_p"))

    # ----- apply cuts for lightweight df -----
    # vertex in FV
    slcdf = slcdf[InFV(slcdf.slc.vertex, 0, det=DETECTOR)]

    # neutrino cuts
    slcdf = slcdf[slcdf.slc.nu_score > 0.5]

    # require both muon and proton to be present
    mask = (~np.isnan(slcdf.mu.pfp.trk.P.p_muon)) & (~np.isnan(slcdf.p.pfp.trk.P.p_proton))
    # mask = mask.reindex(slcdf.index, fill_value=False)
    slcdf = slcdf[mask]

    # ---- truth match ----
    bad_tmatch = np.invert(slcdf.slc.tmatch.eff > 0.5) & (slcdf.slc.tmatch.idx >= 0)
    slcdf.loc[bad_tmatch, ("slc","tmatch","idx", "", "", "", "")] = np.nan

    mcdf.columns = pd.MultiIndex.from_tuples([tuple(list(c) +["", "", "", "", ""]) for c in mcdf.columns])     # match # of column levels

    # neutrino cuts
    slcdf = slcdf[InFV(slcdf.slc.vertex, 0, det=DETECTOR)]
    slcdf = slcdf[slcdf.slc.nu_score > 0.5]

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

