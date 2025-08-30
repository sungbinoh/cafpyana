from pyanalib.pandas_helpers import *
from pyanalib.variable_calculator import *
from makedf.util import *
import pandas as pd
import numpy as np
from makedf.makedf import *
from makedf.constants import *
from makedf.mcstat import mcstat
# from analysis_village.numucc1p0pi.selection_definitions import *
from analysis_village.numucc1p0pi.makedf.makecov import *
from makedf.branches import *

# def make_geniewgts(df, wgt_types=["bnb","genie"]):
#     df.loc[:,'nuint_categ'] = get_int_category(df)
#     df.loc[:,'genie_categ'] = get_genie_category(df)

#     dict_univ_counts = {}
#     for var_config in var_configs:
#         for topo_mode in mode_list:
#             this_df = df[df.nuint_categ == topo_mode]

#             ## CV
#             n_cv, bins = np.histogram(this_df[var_config.var_nu_col], bins=var_config.bins)
#             cv_col_name = ("nuint_categ_"+str(topo_mode), 
#                         "_".join([str(s) for s in var_config.var_nu_col if s]),
#                         "cv",
#                         "")
#             dict_univ_counts[cv_col_name] = n_cv

#             ## GENIE
#             if "genie" in wgt_types:

#                 # multisim
#                 syst_name = "GENIE"
#                 n_syst_univ = 100

#                 n_univ = np.zeros((n_syst_univ, len(var_config.bin_centers)))
#                 for i in range(n_syst_univ):
#                     weights = this_df[syst_name]["univ_{}".format(i)]
#                     weights[np.isnan(weights)] = 1
#                     n_univ[i, :] = np.histogram(this_df[var_config.var_nu_col], bins=var_config.bins, weights=weights)[0]
#                     col_name = ("nuint_categ_"+str(topo_mode), 
#                                 "_".join([str(s) for s in var_config.var_nu_col if s]),
#                                 syst_name+"_multisim", 
#                                 "univ_{}".format(i))
#                     dict_univ_counts[col_name] = n_univ[i, :]

#                 # multisigma
#                 n_syst_univ = 6
#                 for syst_name in regen_systematics_sbnd_multisigma:
#                     n_univ = np.zeros((n_syst_univ, len(var_config.bin_centers)))
#                     for i in range(n_syst_univ):
#                         if (i < 3):
#                             pm = "ms"
#                         else:
#                             pm = "ps"
#                         sig = i%3 + 1
#                         weights = this_df[syst_name][pm+str(sig)]
#                         weights[np.isnan(weights)] = 1
#                         n_univ[i, :] = np.histogram(this_df[var_config.var_nu_col], bins=var_config.bins, weights=weights)[0]
#                         col_name = ("nuint_categ_"+str(topo_mode), 
#                                     "_".join([str(s) for s in var_config.var_nu_col if s]),
#                                     syst_name,
#                                     "univ_{}".format(i))
#                         dict_univ_counts[col_name] = n_univ[i, :]

#                 # morph
#                 n_syst_univ = 1
#                 for syst_name in regen_systematics_sbnd_morph:
#                     n_univ = np.zeros((n_syst_univ, len(var_config.bin_centers)))
#                     for i in range(n_syst_univ):
#                         weights = this_df[syst_name]["morph"]
#                         weights[np.isnan(weights)] = 1
#                         n_univ[i, :] = np.histogram(this_df[var_config.var_nu_col], bins=var_config.bins, weights=weights)[0]
#                         col_name = ("nuint_categ_"+str(topo_mode), 
#                                     "_".join([str(s) for s in var_config.var_nu_col if s]),
#                                     syst_name,
#                                     "univ_{}".format(i))
#                         dict_univ_counts[col_name] = n_univ[i, :]

#             if "bnb" in wgt_types:

#                 ## Flux
#                 syst_name = "Flux"
#                 n_syst_univ = 500
#                 n_univ = np.zeros((n_syst_univ, len(var_config.bin_centers)))
#                 for i in range(n_syst_univ):
#                     weights = this_df[syst_name]["univ_{}".format(i)]
#                     weights[np.isnan(weights)] = 1
#                     n_univ[i, :] = np.histogram(this_df[var_config.var_nu_col], bins=var_config.bins, weights=weights)[0]
#                     col_name = ("nuint_categ_"+str(topo_mode), 
#                                 "_".join([str(s) for s in var_config.var_nu_col if s]),
#                                 syst_name, 
#                                 "univ_{}".format(i))
#                     dict_univ_counts[col_name] = n_univ[i, :]


#             # if "mcstat" in wgt_types:
#             #     ## MCStat
#             #     n_univ_mcstat = 500
#             #     mc_evt_df, MCstat_univ_events = mcstat(mc_evt_df, mc_hdr_df, n_universes=n_univ_mcstat)

#     df_univ_counts = pd.DataFrame(dict_univ_counts)
#     if not isinstance(df_univ_counts.columns, pd.MultiIndex):
#         df_univ_counts.columns = pd.MultiIndex.from_tuples(df_univ_counts.columns)

#     return df_univ_counts


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
                       trkScoreCut=False, trkDistCut=10., cutClearCosmic=False, **trkArgs):
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
    slcdf = multicol_add(slcdf, tki_reco["del_Tp_x"].rename("del_Tp_x"))
    slcdf = multicol_add(slcdf, tki_reco["del_Tp_y"].rename("del_Tp_y"))
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
    mcdf = multicol_add(mcdf, tki_mc["del_Tp_x"].rename("mc_del_Tp_x"))
    mcdf = multicol_add(mcdf, tki_mc["del_Tp_y"].rename("mc_del_Tp_y"))
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

def make_numucc1p0pi_evtdf_final_wgts(f, mode="df", include_weights=True, multisim_nuniv=500, wgt_types=["bnb","genie"], slim=True, 
                           trkScoreCut=False, **trkArgs):
    if mode == "df":
        return make_numucc1p0pi_evtdf(f, sel_level="final", include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                                    trkScoreCut=trkScoreCut, **trkArgs)
    elif mode == "cov":
        # all weights on mcnudf
        print("all mc nus")
        # mcdf = make_mcnudf(f, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim)

        # print("final selection")
        mcdf = make_numucc1p0pi_evtdf(f, sel_level="final", include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                                    trkScoreCut=trkScoreCut, **trkArgs)

        mcdf_wgts = make_geniewgts(mcdf)

        hdr = loadbranches(f["recTree"], mchdrbranches).rec.hdr
        pot = hdr['pot'].sum()

        return mcdf_wgts, pot


def make_numucc1p0pi_binned_systs(filename):
    
    events = uproot.open(filename)
    try:
        ret = make_numucc1p0pi_evtdf_final_wgts(events, mode="cov")

    except Exception as e:
        print(f"Error in make_numucc1p0pi_binned_systs for file {filename}: {e}")
        return None

    return ret


def make_numucc1p0pi_evtdf_2prong(f, include_weights=False, multisim_nuniv=1000, wgt_types=["bnb","genie"], slim=True, 
                           trkScoreCut=False, **trkArgs):
    return make_numucc1p0pi_evtdf(f, sel_level="2prong", include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                                  trkScoreCut=trkScoreCut, **trkArgs)


def make_numucc1p0pi_evtdf(f, sel_level="all", include_weights=False, multisim_nuniv=1000, wgt_types=["bnb","genie"], slim=True, 
                           trkScoreCut=False, **trkArgs):
    """
    sel_level:
        "all": all slices, no cuts
        "clearcosmic": cosmic rejection
        "fv": vertex in FV
        "nu": n-score cut
        "2prong": 2-prong slices
        "2prong_qual": 2-prong slices with quality cuts on the 2 prongs
        "mux": muon-X cuts
        "final": final selection
    """

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

    # Use slcdf from different places according to sel_level

    # Save slcdf at each selection stage for flexibility
    slcdf_all = slcdf.copy()
    slcdf_clearcosmic = slcdf_all[slcdf_all.slc.is_clear_cosmic == 0]
    slcdf_fv = slcdf_clearcosmic[InFV(slcdf_clearcosmic.slc.vertex, 0, det=DETECTOR)]
    slcdf_nu = slcdf_fv[slcdf_fv.slc.nu_score > 0.5]

    # ---- 2-prong cuts ----
    trkdf_index_names = trkdf.index.names
    trkdf_2prong = trkdf.reset_index(level=[2]).loc[slcdf_nu.index].reset_index().set_index(trkdf_index_names)
    trkdf_2prong = trkdf_2prong[trkdf_2prong.pfp.trk.producer != 4294967295]
    npfps = trkdf_2prong.pfp.id.groupby(level=[0,1]).count()
    assert len(npfps) == len(slcdf_nu)
    # sort tracks by length
    trkdf_2prong = trkdf_2prong.sort_values(by=('pfp','trk','len'), ascending=False)
    track1_2prong = trkdf_2prong.groupby(level=[0,1]).nth(0).reset_index(level=[2])
    track2_2prong = trkdf_2prong.groupby(level=[0,1]).nth(1).reset_index(level=[2])
    mask_2prong = (npfps == 2)
    slcdf_2prong = slcdf_nu[mask_2prong]

    # ---- topo requirements on the 2 prongs ----
    # both prongs contained
    mask_2prong_qual = mask_2prong & InFV(track1_2prong.pfp.trk.end, 0, det=DETECTOR)
    mask_2prong_qual = mask_2prong_qual & InFV(track2_2prong.pfp.trk.end, 0, det=DETECTOR)

    # both prongs have track score > 0.5
    mask_2prong_qual = mask_2prong_qual & (track1_2prong.pfp.trackScore > 0.5)
    mask_2prong_qual = mask_2prong_qual & (track2_2prong.pfp.trackScore > 0.5)

    # both (start position - vertex) < 1.2 cm
    mask_2prong_qual = mask_2prong_qual & (dmagdf(track1_2prong.pfp.trk.start, slcdf_nu.slc.vertex) < 1.2)
    mask_2prong_qual = mask_2prong_qual & (dmagdf(track2_2prong.pfp.trk.start, slcdf_nu.slc.vertex) < 1.2)

    slcdf_2prong_qual = slcdf_nu[mask_2prong_qual]
    track1_2prong_qual = track1_2prong.loc[slcdf_2prong_qual.index].reset_index().set_index(trkdf_index_names)
    track2_2prong_qual = track2_2prong.loc[slcdf_2prong_qual.index].reset_index().set_index(trkdf_index_names)

    # --- track PID ----
    trk_2prong_qual = pd.concat([track1_2prong_qual.reset_index().set_index(trkdf_index_names), 
                            track2_2prong_qual.reset_index().set_index(trkdf_index_names)])

    # chi2 score cut
    chi2mu_avg = avg_chi2(trk_2prong_qual, 'chi2_muon')
    chi2p_avg = avg_chi2(trk_2prong_qual, 'chi2_proton')
    mu_mask = (chi2mu_avg > 0) & (chi2mu_avg < 25)
    mu_mask = mu_mask & (chi2p_avg > 90)

    # length cut
    mu_mask = mu_mask & (trk_2prong_qual.pfp.trk.len > 50)

    # quality cut
    mu_mask = mu_mask & (np.abs(trk_2prong_qual.pfp.trk.rangeP.p_muon - trk_2prong_qual.pfp.trk.mcsP.fwdP_muon)/trk_2prong_qual.pfp.trk.rangeP.p_muon < 0.5)

    mu_candidates = trk_2prong_qual[mu_mask]
    #choose longest
    mu_candidates = mu_candidates.sort_values(by=('pfp','trk','len'), ascending=False)
    mu_candidate = mu_candidates.groupby(level=[0,1]).nth(0)
    
    # choose proton from remaining pfps
    notmu_trk_idx = trk_2prong_qual.index.difference(mu_candidate.index)
    notmu_trk = trk_2prong_qual.loc[notmu_trk_idx]

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
    slcdf_mux = slcdf_2prong.loc[mu_idx]

    # --- select mu-p slices ----
    p_idx = p_candidates.reset_index(level=[2]).index.unique()
    mu_p_idx = mu_idx.intersection(p_idx)

    slcdf_final = slcdf_mux.loc[mu_p_idx]
    mudf = mu_candidate.reset_index(level=[2]).loc[mu_p_idx].reset_index().set_index(trkdf_index_names)
    pdf = p_candidates.reset_index(level=[2]).loc[mu_p_idx].reset_index().set_index(trkdf_index_names)

    # --- kinematic phase space cuts ----
    mudf.columns = pd.MultiIndex.from_tuples([tuple(["mu"] + list(c)) for c in mudf.columns])
    slcdf_final = multicol_merge(slcdf_final, mudf.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")
    pdf.columns = pd.MultiIndex.from_tuples([tuple(["p"] + list(c)) for c in pdf.columns])
    slcdf_final = multicol_merge(slcdf_final, pdf.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")

    # --- calculate TKI for reco slices ----
    slc_mudf = slcdf_final.mu.pfp.trk
    slc_pdf = slcdf_final.p.pfp.trk
    slc_P_mu_col = pad_column_name(("P", "p_muon"), slc_mudf)
    slc_P_p_col = pad_column_name(("P", "p_proton"), slc_pdf)
    tki_reco = get_cc1p0pi_tki(slc_mudf, slc_pdf, slc_P_mu_col, slc_P_p_col)

    slcdf_final = multicol_add(slcdf_final, tki_reco["del_alpha"].rename("del_alpha"))
    slcdf_final = multicol_add(slcdf_final, tki_reco["del_phi"].rename("del_phi"))
    slcdf_final = multicol_add(slcdf_final, tki_reco["del_Tp"].rename("del_Tp"))
    slcdf_final = multicol_add(slcdf_final, tki_reco["del_p"].rename("del_p"))
    slcdf_final = multicol_add(slcdf_final, tki_reco["del_Tp_x"].rename("del_Tp_x"))
    slcdf_final = multicol_add(slcdf_final, tki_reco["del_Tp_y"].rename("del_Tp_y"))

    # Now, select which slcdf to use based on sel_level
    # sel_level can be: "all", "clearcosmic", "fv", "nu", "2prong", "mux", "final"
    if sel_level == "all":
        slcdf = slcdf_all
        
    elif sel_level == "clearcosmic":
        slcdf = slcdf_clearcosmic

    elif sel_level == "fv":
        slcdf = slcdf_fv

    elif sel_level == "nu":
        slcdf = slcdf_nu

    elif sel_level == "2prong":
        # save longer and shorter tracks
        t1df = track1_2prong.reset_index().set_index(trkdf_index_names)
        t2df = track2_2prong.reset_index().set_index(trkdf_index_names)
        t1df.columns = pd.MultiIndex.from_tuples([tuple(["t1"] + list(c)) for c in t1df.columns])
        slcdf_2prong = multicol_merge(slcdf_2prong, t1df.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")
        t2df.columns = pd.MultiIndex.from_tuples([tuple(["t2"] + list(c)) for c in t2df.columns])
        slcdf_2prong = multicol_merge(slcdf_2prong, t2df.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")

        slcdf = slcdf_2prong

    elif sel_level == "2prong_qual":
        # save longer and shorter tracks
        t1df = track1_2prong_qual.reset_index().set_index(trkdf_index_names)
        t2df = track2_2prong_qual.reset_index().set_index(trkdf_index_names)
        t1df.columns = pd.MultiIndex.from_tuples([tuple(["t1"] + list(c)) for c in t1df.columns])
        slcdf_2prong_qual = multicol_merge(slcdf_2prong_qual, t1df.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")
        t2df.columns = pd.MultiIndex.from_tuples([tuple(["t2"] + list(c)) for c in t2df.columns])
        slcdf_2prong_qual = multicol_merge(slcdf_2prong_qual, t2df.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")

        slcdf = slcdf_2prong_qual

    elif sel_level == "mux":
        slcdf = slcdf_mux

    elif sel_level == "final":
        slcdf = slcdf_final

    else:
        raise ValueError(f"Unknown sel_level: {sel_level}")

    ##### end of reco selection here #####


    # --- calculate TKI for nu MC ---
    mc_mudf = mcdf.mu
    mc_pdf = mcdf.p
    mc_P_mu_col = pad_column_name(("totp",), mc_mudf)
    mc_P_p_col = pad_column_name(("totp",), mc_pdf)
    tki_mc = get_cc1p0pi_tki(mc_mudf, mc_pdf, mc_P_mu_col, mc_P_p_col)

    mcdf = multicol_add(mcdf, tki_mc["del_alpha"].rename("mc_del_alpha"))
    mcdf = multicol_add(mcdf, tki_mc["del_phi"].rename("mc_del_phi"))
    mcdf = multicol_add(mcdf, tki_mc["del_Tp"].rename("mc_del_Tp"))
    mcdf = multicol_add(mcdf, tki_mc["del_p"].rename("mc_del_p"))
    mcdf = multicol_add(mcdf, tki_mc["del_Tp_x"].rename("mc_del_Tp_x"))
    mcdf = multicol_add(mcdf, tki_mc["del_Tp_y"].rename("mc_del_Tp_y"))

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
    # if include_weights:
    #     df, _ = mcstat(df, hdrdf, n_universes=multisim_nuniv)

    df = df.set_index(slcdf.index.names, verify_integrity=True)
    
    return df

