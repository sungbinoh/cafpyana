from pyanalib.pandas_helpers import *
from pyanalib.variable_calculator import *
from makedf.util import *
import pandas as pd
import numpy as np
from makedf.makedf import *
from makedf.constants import *
from makedf.mcstat import mcstat

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


def make_pandora_evtdf(f, include_weights=False, multisim_nuniv=1000, wgt_types=["bnb","genie"], slim=True, 
                       trkScoreCut=False, trkDistCut=10., cutClearCosmic=True, **trkArgs):
    # ----- sbnd or icarus? -----
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det
    if (1 == det.unique()):
        DETECTOR = "SBND"
    else:
        DETECTOR = "ICARUS"
    
    hdrdf = loadbranches(f["recTree"], hdrbranches).rec.hdr
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

    # ----- apply some cuts for lightweight df -----
    # ----- use the make_numucc1p0pi_evtdf function to get the events after the full selection -----
    # # vertex in FV
    # slcdf = slcdf[InFV(slcdf.slc.vertex, 0, det=DETECTOR)]

    # # neutrino cuts
    # slcdf = slcdf[slcdf.slc.nu_score > 0.5]

    # # require both muon and proton to be present
    # mask = (~np.isnan(slcdf.mu.pfp.trk.P.p_muon)) & (~np.isnan(slcdf.p.pfp.trk.P.p_proton))
    # slcdf = slcdf[mask]


    # ---- truth match ----
    bad_tmatch = np.invert(slcdf.slc.tmatch.eff > 0.5) & (slcdf.slc.tmatch.idx >= 0)
    slcdf.loc[bad_tmatch, ("slc","tmatch","idx", "", "", "", "")] = np.nan

    mcdf.columns = pd.MultiIndex.from_tuples([tuple(list(c) +["", "", "", "", ""]) for c in mcdf.columns])     # match # of column levels

    df = multicol_merge(slcdf.reset_index(), 
                  mcdf.reset_index(),
                  left_on=[("entry", "", "",), 
                           ("slc", "tmatch", "idx")], 
                  right_on=[("entry", "", ""), 
                            ("rec.mc.nu..index", "", "")], 
                  how="left"
                  ) 

    if include_weights:
        df, _ = mcstat(df, hdrdf, n_universes=multisim_nuniv)

    df = df.set_index(slcdf.index.names, verify_integrity=True)
    
    return df


def make_numucc1p0pi_evtdf(f, include_weights=False, multisim_nuniv=1000, wgt_types=["bnb","genie"], slim=True, 
                       trkScoreCut=False, **trkArgs):

    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det
    if (1 == det.unique()):
        DETECTOR = "SBND"
    else:
        DETECTOR = "ICARUS"
    
    hdrdf = loadbranches(f["recTree"], hdrbranches).rec.hdr
    mcdf = make_mcnudf(f, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim)
    trkdf = make_trkdf(f, trkScoreCut, **trkArgs)
    slcdf = make_slcdf(f)

    # TODO: stubs need to be validated in data
    # stubdf = make_stubs(f, det=DETECTOR)
    # load stubs
    # slcdf = multicol_merge(slcdf, stubdf, left_index=True, right_index=True)

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
    
    # ===== NuInt Selection ===== 

    # not clear cosmic
    slcdf = slcdf[slcdf.slc.is_clear_cosmic == 0]

    # vertex in FV
    slcdf = slcdf[InFV(slcdf.slc.vertex, 0, det=DETECTOR)]

    # neutrino cuts
    slcdf = slcdf[slcdf.slc.nu_score > 0.5]

    # ---- 2-prong cuts ----
    trkdf_index_names = trkdf.index.names
    trkdf = trkdf.reset_index(level=[2]).loc[slcdf.index].reset_index().set_index(trkdf_index_names)
    trkdf = trkdf[trkdf.pfp.trk.producer != 4294967295]
    npfps = trkdf.pfp.id.groupby(level=[0,1]).count()
    assert len(npfps) == len(slcdf)
    # sort tracks by length
    trkdf = trkdf.sort_values(by=('pfp','trk','len'), ascending=False)
    track1 = trkdf.groupby(level=[0,1]).nth(0).reset_index(level=[2])
    track2 = trkdf.groupby(level=[0,1]).nth(1).reset_index(level=[2])
    mask = (npfps == 2)

    # both prongs contained
    mask = mask & InFV(track1.pfp.trk.end, 0, det=DETECTOR)
    mask = mask & InFV(track2.pfp.trk.end, 0, det=DETECTOR)

    # both prongs have track score > 0.5
    mask = mask & (track1.pfp.trackScore > 0.5)
    mask = mask & (track2.pfp.trackScore > 0.5)

    # both (start position - vertex) < 1.2 cm
    def mag3d(df1, df2):
        return np.sqrt((df1.x - df2.x)**2 + (df1.y - df2.y)**2 + (df1.z - df2.z)**2)
    mask = mask & (mag3d(track1.pfp.trk.start, slcdf.slc.vertex) < 1.2)
    mask = mask & (mag3d(track2.pfp.trk.start, slcdf.slc.vertex) < 1.2)

    slcdf = slcdf[mask]
    track1 = track1.loc[slcdf.index].reset_index().set_index(trkdf_index_names)
    track2 = track2.loc[slcdf.index].reset_index().set_index(trkdf_index_names)

    # --- track PID ----
    def avg_chi2(df, var_name):
        planes = ['I0', 'I1', 'I2']
        chi2_vals = []
        for plane in planes:
            chi2 = df['pfp']['trk']['chi2pid'][plane][var_name]
            chi2_vals.append(chi2)
        chi2_df = pd.concat(chi2_vals, axis=1)
        # fill 0 with nan
        chi2_df = chi2_df.replace(0, np.nan)
        avg = chi2_df.mean(axis=1, skipna=True)
        return avg

    trk_2prong = pd.concat([track1.reset_index().set_index(trkdf_index_names), 
                            track2.reset_index().set_index(trkdf_index_names)])

    # chi2 score cut
    chi2mu_avg = avg_chi2(trk_2prong, 'chi2_muon')
    chi2p_avg = avg_chi2(trk_2prong, 'chi2_proton')
    mu_mask = (chi2mu_avg > 0) & (chi2mu_avg < 25)
    mu_mask = mu_mask & (chi2p_avg > 90)

    # length cut
    mu_mask = mu_mask & (trk_2prong.pfp.trk.len > 50)

    # quality cut
    mu_mask = mu_mask & (np.abs(trk_2prong.pfp.trk.rangeP.p_muon - trk_2prong.pfp.trk.mcsP.fwdP_muon)/trk_2prong.pfp.trk.rangeP.p_muon < 0.5)

    mu_candidates = trk_2prong[mu_mask]
    #choose longest
    mu_candidates = mu_candidates.sort_values(by=('pfp','trk','len'), ascending=False)
    mu_candidate = mu_candidates.groupby(level=[0,1]).nth(0)
    
    # choose proton from remaining pfps
    notmu_trk_idx = trk_2prong.index.difference(mu_candidate.index)
    notmu_trk = trk_2prong.loc[notmu_trk_idx]

    # chi2 score cut
    chi2p_avg = avg_chi2(notmu_trk, 'chi2_proton')
    p_mask = (chi2p_avg > 0) & (chi2p_avg < 90)
    p_candidates = notmu_trk[p_mask]

    # --- kinematic phase space cuts ----
    mu_mask = (mu_candidate.pfp.trk.rangeP.p_muon > 0.22)
    mu_mask = mu_mask & (mu_candidate.pfp.trk.rangeP.p_muon < 1)
    mu_candidate = mu_candidate[mu_mask]

    p_mask = (p_candidates.pfp.trk.rangeP.p_proton > 0.3) 
    p_mask = p_mask & (p_candidates.pfp.trk.rangeP.p_proton < 1)
    p_candidates = p_candidates[p_mask]

    # --- select mu-X slices ----
    mu_idx = mu_candidate.reset_index(level=[2]).index.unique()
    slcdf = slcdf.loc[mu_idx]

    # --- select mu-p slices ----
    p_idx = p_candidates.reset_index(level=[2]).index.unique()
    mu_p_idx = mu_idx.intersection(p_idx)

    slcdf = slcdf.loc[mu_p_idx]
    mudf = track1.reset_index(level=[2]).loc[mu_p_idx].reset_index().set_index(trkdf_index_names)
    pdf = track2.reset_index(level=[2]).loc[mu_p_idx].reset_index().set_index(trkdf_index_names)

    # --- kinematic phase space cuts ----
    mudf.columns = pd.MultiIndex.from_tuples([tuple(["mu"] + list(c)) for c in mudf.columns])
    slcdf = multicol_merge(slcdf, mudf.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")
    pdf.columns = pd.MultiIndex.from_tuples([tuple(["p"] + list(c)) for c in pdf.columns])
    slcdf = multicol_merge(slcdf, pdf.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")

    # --- calculate kinematic variables ----
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

    # ---- truth match ----
    bad_tmatch = np.invert(slcdf.slc.tmatch.eff > 0.5) & (slcdf.slc.tmatch.idx >= 0)
    slcdf.loc[bad_tmatch, ("slc","tmatch","idx", "", "", "", "")] = np.nan

    mcdf.columns = pd.MultiIndex.from_tuples([tuple(list(c) +["", "", "", "", ""]) for c in mcdf.columns])     # match # of column levels

    df = multicol_merge(slcdf.reset_index(), 
                  mcdf.reset_index(),
                  left_on=[("entry", "", "",), 
                           ("slc", "tmatch", "idx")], 
                  right_on=[("entry", "", ""), 
                            ("rec.mc.nu..index", "", "")], 
                  how="left"
                  ) 

    # add MCstat weights if requested
    if include_weights:
        df, _ = mcstat(df, hdrdf, n_universes=multisim_nuniv)

    df = df.set_index(slcdf.index.names, verify_integrity=True)
    
    return df

