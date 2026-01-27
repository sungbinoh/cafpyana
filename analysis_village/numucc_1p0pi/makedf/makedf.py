from pyanalib.pandas_helpers import *
from pyanalib.variable_calculator import *
from makedf.util import *
import pandas as pd
import numpy as np
from makedf.makedf import *
from makedf.constants import *

from analysis_village.numucc_1p0pi.makedf.util import *

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


def make_pandora_evtdf_wgts(f, include_weights=True, multisim_nuniv=1000, wgt_types=["bnb","genie"], slim=True, 
                       trkScoreCut=False, trkDistCut=10., cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df

def make_pandora_evtdf_mup_wgts(f, sel_level="mup", include_weights=True, multisim_nuniv=200, wgt_types=["bnb","genie","g4"], slim=True, 
                       trkScoreCut=False, trkDistCut=100., cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df

def make_pandora_evtdf_mup(f, sel_level="mup", include_weights=False, multisim_nuniv=0, wgt_types=[], slim=True, 
                       trkScoreCut=False, trkDistCut=100., cutClearCosmic=True, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, cutClearCosmic=cutClearCosmic, **trkArgs)
    return df

# GENIE dfs
qe_genie_systematics = [
'GENIEReWeight_SBN_v1_multisim_RPA_CCQE',
'GENIEReWeight_SBN_v1_multisim_CoulombCCQE',
]
def make_pandora_evtdf_mup_wgts_CCQE(f, sel_level="mup", include_weights=True, multisim_nuniv=200, wgt_types=["genie"], slim=False, 
                       trkScoreCut=False, trkDistCut=100., cutClearCosmic=True, genie_systematics=None, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, cutClearCosmic=cutClearCosmic, genie_systematics=qe_genie_systematics, **trkArgs)
    return df

def make_mcnudf_CCQE(f, include_weights=True, multisim_nuniv=100, wgt_types=["genie"], slim=False, genie_systematics=None):
    df = make_mcnudf(f, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, genie_systematics=qe_genie_systematics)
    return df

mec_genie_systematics = [
'GENIEReWeight_SBN_v1_multisim_NormCCMEC',
'GENIEReWeight_SBN_v1_multisim_NormNCMEC',
"GENIEReWeight_SBN_v1_multisigma_DecayAngMEC",
]

def make_pandora_evtdf_mup_wgts_MEC(f, sel_level="mup", include_weights=True, multisim_nuniv=200, wgt_types=["genie"], slim=False, 
                       trkScoreCut=False, trkDistCut=100., cutClearCosmic=True, genie_systematics=None, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, cutClearCosmic=cutClearCosmic, genie_systematics=mec_genie_systematics, **trkArgs)

def make_mcnudf_MEC(f, include_weights=True, multisim_nuniv=100, wgt_types=["genie"], slim=False, genie_systematics=None):
    df = make_mcnudf(f, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, genie_systematics=mec_genie_systematics)
    return df

res_genie_systematics = [
'GENIEReWeight_SBN_v1_multisim_RDecBR1gamma',
'GENIEReWeight_SBN_v1_multisim_RDecBR1eta',
"GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi",
"GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad",

"GENIEReWeight_SBN_v1_multisigma_MaCCRES",
"GENIEReWeight_SBN_v1_multisigma_MaNCRES",
"GENIEReWeight_SBN_v1_multisigma_MvCCRES",
"GENIEReWeight_SBN_v1_multisigma_MvNCRES",
]
def make_pandora_evtdf_mup_wgts_RES(f, sel_level="mup", include_weights=True, multisim_nuniv=200, wgt_types=["genie"], slim=False, 
                       trkScoreCut=False, trkDistCut=100., cutClearCosmic=True, genie_systematics=None, **trkArgs):
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, cutClearCosmic=cutClearCosmic, genie_systematics=res_genie_systematics, **trkArgs)

def make_mcnudf_RES(f, include_weights=True, multisim_nuniv=100, wgt_types=["genie"], slim=False, genie_systematics=None):
    df = make_mcnudf(f, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, genie_systematics=res_genie_systematics)
    return df

def make_pandora_evtdf_mup_wgts_NonRES(f, sel_level="mup", include_weights=True, multisim_nuniv=200, wgt_types=["genie"], slim=True, 
                       trkScoreCut=False, trkDistCut=100., cutClearCosmic=True, genie_systematics=None, **trkArgs):
    this_genie_systematics = [
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi',
    ]
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, cutClearCosmic=cutClearCosmic, genie_systematics=this_genie_systematics, **trkArgs)

def make_pandora_evtdf_mup_wgts_DIS(f, sel_level="mup", include_weights=True, multisim_nuniv=200, wgt_types=["genie"], slim=False, 
                       trkScoreCut=False, trkDistCut=100., cutClearCosmic=True, genie_systematics=None, **trkArgs):
    this_genie_systematics = [
    'GENIEReWeight_SBN_v1_multisigma_AhtBY',
    'GENIEReWeight_SBN_v1_multisigma_BhtBY',
    'GENIEReWeight_SBN_v1_multisigma_CV1uBY',
    'GENIEReWeight_SBN_v1_multisigma_CV2uBY',
    ]
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, cutClearCosmic=cutClearCosmic, genie_systematics=this_genie_systematics, **trkArgs)

def make_pandora_evtdf_mup_wgts_Other(f, sel_level="mup", include_weights=True, multisim_nuniv=200, wgt_types=["genie"], slim=False, 
                       trkScoreCut=False, trkDistCut=100., cutClearCosmic=True, genie_systematics=None, **trkArgs):
    this_genie_systematics = [
    'GENIEReWeight_SBN_v1_multisigma_MFP_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrCEx_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrInel_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrAbs_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi',
    'GENIEReWeight_SBN_v1_multisigma_MFP_N',
    'GENIEReWeight_SBN_v1_multisigma_FrCEx_N',
    'GENIEReWeight_SBN_v1_multisigma_FrInel_N',
    'GENIEReWeight_SBN_v1_multisigma_FrAbs_N',
    'GENIEReWeight_SBN_v1_multisigma_FrPiProd_N',
    "GENIEReWeight_SBN_v1_multisigma_NormCCCOH", # Handled by re-tuning
    "GENIEReWeight_SBN_v1_multisigma_NormNCCOH",
    'GENIEReWeight_SBN_v1_multisigma_MaNCEL',
    'GENIEReWeight_SBN_v1_multisigma_EtaNCEL',
    ]
    df = make_pandora_evtdf(f, sel_level=sel_level, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, 
                            trkScoreCut=trkScoreCut, trkDistCut=trkDistCut, cutClearCosmic=cutClearCosmic, genie_systematics=this_genie_systematics, **trkArgs)

def make_pandora_evtdf(f, sel_level="all", include_weights=True, multisim_nuniv=1000, wgt_types=[], slim=True, genie_systematics=None,
                       trkScoreCut=False, trkDistCut=100., cutClearCosmic=True, **trkArgs):

    """
    sel_level:
        "all": all slices, no cuts
        "clearcosmic": cosmic rejection
        "fv": vertex in FV
        "nu": n-score cut
        "2prong": 2-prong slices
        "2prong_qual": 2-prong slices with quality cuts on the 2 prongs
        "muX": muon-X cuts
        "mup": final selection
    """

    if sel_level not in ["all", "clearcosmic", "fv", "nu", "2prong", "2prong_contained", "2prong_trackscore", "2prong_vtxdist", "muX", "mup"]:
        raise ValueError("Invalid sel_level: {}".format(sel_level))

    def truth_match(this_evtdf, this_mcdf):
        # ---- truth match ----
        bad_tmatch = np.invert(this_evtdf.slc.tmatch.eff > 0.5) & (this_evtdf.slc.tmatch.idx >= 0)
        # this_evtdf.loc[bad_tmatch, ("slc","tmatch","idx", "", "", "", "")] = np.nan
        this_evtdf.loc[bad_tmatch, pad_column_name(("slc","tmatch","idx"), this_evtdf)] = np.nan

        nlevels = this_evtdf.columns.nlevels

        this_mcdf.columns = pd.MultiIndex.from_tuples([tuple(["mc"] + list(c) +[""] * (nlevels-len(c)-1)) for c in this_mcdf.columns])     # match # of column levels
        df = multicol_merge(this_evtdf.reset_index(), 
                    this_mcdf.reset_index(),
                    left_on=[("entry", "", "",), 
                            ("slc", "tmatch", "idx")], 
                    right_on=[("entry", "", ""), 
                                ("rec.mc.nu..index", "", "")], 
                    how="left"
                    ) 

        df = df.set_index(this_evtdf.index.names, verify_integrity=True) 
        return df

    # SBND Gen 1 analysis
    DETECTOR = "SBND_nohighyz"

    # event selection cuts
    # slice cuts
    nu_score_th = 0.45
    save_ntrks = 2

    # track cuts
    trackscore_th = 0.5
    vtxdist_th = 1.2

    # muon candidate cuts
    mu_chi2mu_th = 30
    mu_chi2p_th = 100
    mu_len_th = 50
    qual_th = 0.2

    # proton candidate cuts
    p_chi2p_th = 90
    p_len_th = 0

    # track kinematics cuts
    mu_Plo_th = 0.22
    mu_Phi_th = 1
    p_Plo_th = 0.3
    p_Phi_th = 1

    mcdf = make_mcnudf(f, include_weights=include_weights, multisim_nuniv=multisim_nuniv, wgt_types=wgt_types, slim=slim, genie_systematics=genie_systematics)
    # calculate TKI for MC 
    tki_var_names = ["del_alpha", "del_phi", "del_Tp", "del_p", "del_Tp_x", "del_Tp_y"]
    mc_mudf = mcdf.mu
    mc_pdf = mcdf.p
    mc_P_mu_col = pad_column_name(("totp",), mc_mudf)
    mc_P_p_col = pad_column_name(("totp",), mc_pdf)
    tki_mc = get_cc1p0pi_tki(mc_mudf, mc_pdf, mc_P_mu_col, mc_P_p_col)
    for var_name in tki_var_names:
        mcdf = multicol_add(mcdf, tki_mc[var_name].rename("{}".format(var_name)))

    slcdf = make_slcdf(f)
    trkdf = make_trkdf(f, det=DETECTOR, scoreCut=trkScoreCut, **trkArgs)

    if sel_level == "all":
        return truth_match(slcdf, mcdf)

    slcdf = cut_clear_cosmic(slcdf)
    if sel_level == "clearcosmic":
        return truth_match(slcdf, mcdf)

    slcdf = cut_vertex_in_fv(slcdf, det=DETECTOR)
    if sel_level == "fv":
        return truth_match(slcdf, mcdf)

    slcdf = cut_nu_score(slcdf, nu_score_th)
    if sel_level == "nu":
        return truth_match(slcdf, mcdf)

    trkdf = get_valid_trks(trkdf)
    trkdf = match_trkdf_to_slcdf(trkdf, slcdf)
    evtdf = get_trk_info(slcdf, trkdf, save_ntrks)

    evtdf = cut_2prong(evtdf)
    if sel_level == "2prong":
        return truth_match(evtdf, mcdf)

    evtdf = cut_2prong_contained(evtdf, det=DETECTOR)
    if sel_level == "2prong_contained":
        return truth_match(evtdf, mcdf)

    evtdf = cut_2prong_trackscore(evtdf, trackscore_th)
    if sel_level == "2prong_trackscore":
        return truth_match(evtdf, mcdf)

    evtdf = cut_2prong_vtxdist(evtdf, vtxdist_th)
    if sel_level == "2prong_vtxdist":
        return truth_match(evtdf, mcdf)

    evtdf = get_mu_p_candidate(evtdf, 
                                mu_chi2mu_th=mu_chi2mu_th, mu_chi2p_th=mu_chi2p_th, mu_len_th=mu_len_th, qual_th=qual_th, 
                                p_chi2mu_th=-1, p_chi2p_th=p_chi2p_th, p_len_th=p_len_th)

    evtdf = cut_has_mu(evtdf)
    evtdf = cut_mu_kinematics(evtdf, mu_Plo_th=mu_Plo_th, mu_Phi_th=mu_Phi_th)
    if sel_level == "muX":
        return truth_match(evtdf, mcdf)

    evtdf = cut_has_p(evtdf)
    evtdf = cut_p_kinematics(evtdf, p_Plo_th=p_Plo_th, p_Phi_th=p_Phi_th)

    # calculate TKI for reco slices
    slc_mudf = evtdf.mu.pfp.trk
    slc_pdf = evtdf.p.pfp.trk
    slc_P_mu_col = pad_column_name(("P", "p_muon"), slc_mudf)
    slc_P_p_col = pad_column_name(("P", "p_proton"), slc_pdf)
    tki_reco = get_cc1p0pi_tki(slc_mudf, slc_pdf, slc_P_mu_col, slc_P_p_col)
    for var_name in tki_var_names:
        evtdf = multicol_add(evtdf, tki_reco[var_name].rename(var_name))

    if sel_level == "mup":
        return truth_match(evtdf, mcdf)