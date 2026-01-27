from pyanalib.pandas_helpers import *
from pyanalib.variable_calculator import *
from makedf.util import *
import pandas as pd
import numpy as np
from makedf.makedf import *
from makedf.constants import *
from makedf import chi2pid
from analysis_village.gump.kinematics import *
from analysis_village.gump.gump_cuts import *

# to do: make_pandora_with_cuts using correct formatting 
# and can be turned into ttree for PROfit

def make_pandora_no_cuts_df(f):
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det
    if det.empty:
        return pd.DataFrame()

    if (1 == det.unique()):
        DETECTOR = "SBND"
    elif (2 == det.unique()):
        DETECTOR = "ICARUS"
    else:
        print("Detector unclear, check rec.hdr.det!")

    slcdf = make_slcdf(f)
    StartingRows = len(slcdf)

    trkdf = make_trkdf(f, False)
    trkdf = multicol_add(trkdf, dmagdf(slcdf.slc.vertex, trkdf.pfp.trk.start).rename(("pfp", "dist_to_vertex")))
    trkdf = trkdf[trkdf.pfp.dist_to_vertex < 10]

    # redo chi2 for ICARUS
    if DETECTOR == "ICARUS":
        trkhitdf = make_trkhitdf(f)

        # systematic variations
        dedx_redo = chi2pid.dedx(trkhitdf, gain="ICARUS", calibrate="ICARUS")
        trkhitdf["dedx_redo"] = dedx_redo

        dedx_hi = chi2pid.dedx(trkhitdf, gain="ICARUS", calibrate="ICARUS", scale=1.01)
        trkhitdf["dedx_hi"] = dedx_hi
        dedx_lo = chi2pid.dedx(trkhitdf, gain="ICARUS", calibrate="ICARUS", scale=0.99)
        trkhitdf["dedx_lo"] = dedx_lo
        dedx_smear = chi2pid.dedx(trkhitdf, gain="ICARUS", calibrate="ICARUS", smear=0.05)
        trkhitdf["dedx_smear"] = dedx_smear

        trkdf["chi2u"] = chi2pid.chi2u(trkhitdf, dedxname="dedx_redo")[0]
        trkdf["chi2p"] = chi2pid.chi2p(trkhitdf, dedxname="dedx_redo")[0]

        trkdf["chi2u_lo"] = chi2pid.chi2u(trkhitdf, dedxname="dedx_lo")[0]
        trkdf["chi2p_lo"] = chi2pid.chi2p(trkhitdf, dedxname="dedx_lo")[0]

        trkdf["chi2u_hi"] = chi2pid.chi2u(trkhitdf, dedxname="dedx_hi")[0]
        trkdf["chi2p_hi"] = chi2pid.chi2p(trkhitdf, dedxname="dedx_hi")[0]

        trkdf["chi2u_smear"] = chi2pid.chi2u(trkhitdf, dedxname="dedx_smear")[0]
        trkdf["chi2p_smear"] = chi2pid.chi2p(trkhitdf, dedxname="dedx_smear")[0]
    else:
        trkdf["chi2u"] = trkdf.pfp.trk.chi2pid.I2.chi2_muon
        trkdf["chi2p"] = trkdf.pfp.trk.chi2pid.I2.chi2_proton

        # TODO: implement
        trkdf["chi2u_lo"] = trkdf.chi2u
        trkdf["chi2u_hi"] = trkdf.chi2u
        trkdf["chi2u_smear"] = trkdf.chi2u
        trkdf["chi2p_lo"] = trkdf.chi2p
        trkdf["chi2p_hi"] = trkdf.chi2p
        trkdf["chi2p_smear"] = trkdf.chi2p

    trkdf[("pfp", "trk", "chi2pid", "I2", "mu_over_p", "")] = trkdf.chi2u / trkdf.chi2p

    # track containment
    trkdf[("pfp", "trk", "is_contained", "", "", "")] = trkfv_cut(trkdf.pfp.trk.start, DETECTOR) & trkfv_cut(trkdf.pfp.trk.end, DETECTOR)

    # reco momentum -- range-only
    trkdf[("pfp", "trk", "P", "p_muon", "", "")] = trkdf[("pfp", "trk", "rangeP", "p_muon", "", "")]
    trkdf[("pfp", "trk", "P", "p_pion", "", "")] = trkdf[("pfp", "trk", "rangeP", "p_pion", "", "")]
    trkdf[("pfp", "trk", "P", "p_proton", "", "")] = trkdf[("pfp", "trk", "rangeP", "p_proton", "", "")]

    # mu candidate is track pfp with smallest chi2_mu/chi2_p
    mudf = trkdf[(trkdf.pfp.trackScore> 0.0)].sort_values(trkdf.pfp.index.names[:-1] + [("pfp", "trk", "chi2pid","I2","mu_over_p", "")]).groupby(level=[0, 1]).head(1)
    # mudf = trkdf[(trkdf.pfp.trackScore> 0.0)].sort_values(trkdf.pfp.index.names[:-1] + [("pfp", "trk", "chi2pid","I2","mu_over_p", "")]).groupby(level=[0, 1]).head(1)
    mudf.columns = pd.MultiIndex.from_tuples([tuple(["mu"] + list(c)) for c in mudf.columns])
    slcdf = multicol_merge(slcdf, mudf.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")
    idx_mu = mudf.index

    # p candidate is track pfp with largest chi2_mu/chi2_p of remaining pfps
    idx_pfps = trkdf.pfp.index
    idx_not_mu = idx_pfps.difference(idx_mu)
    notmudf = trkdf.loc[idx_not_mu]
    pdf = notmudf[(notmudf.pfp.trackScore > 0.0)].sort_values(notmudf.pfp.index.names[:-1] + [("pfp", "trk", "chi2pid", "I2", "mu_over_p", "")]).groupby(level=[0,1]).tail(1)
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

    if slcdf.empty:
        print("found empty slice!")
        empty_index = pd.MultiIndex(
            levels=[[], []],
            codes=[[], []],
            names=['entry','rec.slc..index']
        )
        del_p = pd.Series(dtype='float', name='del_p', index=empty_index)
        del_Tp = pd.Series(dtype='float', name='del_Tp', index=empty_index)
        del_phi = pd.Series(dtype='float', name='del_phi', index=empty_index)
        del_alpha = pd.Series(dtype='float', name='del_alpha', index=empty_index)
        mu_E = pd.Series(dtype='float', name='mu_E', index=empty_index)
        p_E = pd.Series(dtype='float', name='p_E', index=empty_index)
        nu_E_calo = pd.Series(dtype='float', name='nu_E_calo', index=empty_index)
        has_stub = pd.Series(dtype='float', name='has_stub', index=empty_index)
        is_contained = pd.Series(dtype='float', name='is_contained', index=empty_index)
    else:
        tki = transverse_kinematics(slcdf.mu.pfp.trk.P.p_muon, slcdf.mu.pfp.trk.dir, slcdf.p.pfp.trk.P.p_proton, slcdf.p.pfp.trk.dir)
        nu_E_calo = neutrino_energy(slcdf.mu.pfp.trk.P.p_muon, slcdf.mu.pfp.trk.dir, slcdf.p.pfp.trk.P.p_proton, slcdf.p.pfp.trk.dir)
        del_p = tki['del_p']

        del_Tp = tki['del_Tp']
        del_phi = tki['del_phi']
        del_alpha = tki['del_alpha']
        mu_E = tki['mu_E']
        p_E = tki['p_E']
        is_contained = slcdf.p.pfp.trk.is_contained & slcdf.mu.pfp.trk.is_contained

    ######## (9) - c: slc.tmatch.idx for truth matching
    bad_tmatch = np.invert(slcdf.slc.tmatch.eff > 0.5) & (slcdf.slc.tmatch.idx >= 0)
    slcdf.loc[bad_tmatch, ("slc","tmatch","idx", "", "", "", "")] = np.nan
    tmatch_idx_series = slcdf.slc.tmatch.idx

    slc_vtx = slcdf.slc.vertex
    nu_score = slcdf.slc.nu_score
    true_pdg = slcdf.slc.truth.pdg
    crlongtrkdiry = slcdf.slc.nuid.crlongtrkdiry
    is_clear_cosmic = slcdf.slc.is_clear_cosmic
    other_shw_length = slcdf.other_shw_length
    other_trk_length = slcdf.other_trk_length
    mu_chi2_of_mu_cand = slcdf.mu.chi2u
    prot_chi2_of_mu_cand = slcdf.mu.chi2p
    mu_chi2_of_prot_cand = slcdf.p.chi2u
    prot_chi2_of_prot_cand = slcdf.p.chi2p
    mu_len = slcdf.mu.pfp.trk.len
    p_len = slcdf.p.pfp.trk.len

    stubdf = make_stubs(f)
    has_any_stub_series = stubdf.groupby('rec.slc..index')['pass_proton_stub'].transform('any')
    slc_has_stub_series = pd.Series(index=mu_len.index, dtype=bool)

    for k in slc_has_stub_series.keys():
        try:
            threek = k + (0,)
            slc_has_stub_series[k] = has_any_stub_series[threek]
        except KeyError:
            slc_has_stub_series[k] = False

    ## (10) create a slice-based reco df
    slcdf = pd.DataFrame({
        'other_shw_length': other_shw_length,
        'other_trk_length': other_trk_length,
        'slc_vtx_x': slc_vtx.x,
        'slc_vtx_y': slc_vtx.y,
        'slc_vtx_z': slc_vtx.z,
        'is_clear_cosmic': is_clear_cosmic,
        'nu_score': nu_score,
        'true_pdg': true_pdg,
        'is_cosmic': (true_pdg == -1),
        'is_contained': is_contained,
        'crlongtrkdiry': crlongtrkdiry,

        'mu_chi2_of_mu_cand': mu_chi2_of_mu_cand,
        'mu_chi2_of_prot_cand': mu_chi2_of_prot_cand,
        'prot_chi2_of_mu_cand': prot_chi2_of_mu_cand,
        'prot_chi2_of_prot_cand': prot_chi2_of_prot_cand,

        'mu_chi2lo_of_mu_cand': slcdf.mu.chi2u_lo,
        'mu_chi2hi_of_mu_cand': slcdf.mu.chi2u_hi,
        'mu_chi2smear_of_mu_cand': slcdf.mu.chi2u_smear,
        'mu_chi2lo_of_prot_cand': slcdf.p.chi2u_lo,
        'mu_chi2hi_of_prot_cand': slcdf.p.chi2u_hi,
        'mu_chi2smear_of_prot_cand': slcdf.p.chi2u_smear,
        'prot_chi2lo_of_mu_cand': slcdf.mu.chi2p_lo,
        'prot_chi2hi_of_mu_cand': slcdf.mu.chi2p_hi,
        'prot_chi2smear_of_mu_cand': slcdf.mu.chi2p_smear,
        'prot_chi2lo_of_prot_cand': slcdf.p.chi2p_lo,
        'prot_chi2hi_of_prot_cand': slcdf.p.chi2p_hi,
        'prot_chi2smear_of_prot_cand': slcdf.p.chi2p_smear,

        'p_len': p_len,
        'mu_len': mu_len,
        'nu_E_calo': nu_E_calo,
        'mu_E': mu_E,
        'mu_T': mu_E - MUON_MASS,
        'p_E': p_E,
        'p_T': p_E - PROTON_MASS,
        'mu_end_x': slcdf.mu.pfp.trk.end.x,
        'mu_end_y': slcdf.mu.pfp.trk.end.y,
        'mu_end_z': slcdf.mu.pfp.trk.end.z,
        'p_end_x': slcdf.p.pfp.trk.end.x,
        'p_end_y': slcdf.p.pfp.trk.end.y,
        'p_end_z': slcdf.p.pfp.trk.end.z,
        'mu_dir_x': slcdf.mu.pfp.trk.dir.x,
        'mu_dir_y': slcdf.mu.pfp.trk.dir.y,
        'mu_dir_z': slcdf.mu.pfp.trk.dir.z,
        'p_dir_x': slcdf.p.pfp.trk.dir.x,
        'p_dir_y': slcdf.p.pfp.trk.dir.y,
        'p_dir_z': slcdf.p.pfp.trk.dir.z,
        'mu_true_p': magdf(slcdf.mu.pfp.trk.truth.p.genp),
        'mu_true_pdg': slcdf.mu.pfp.trk.truth.p.pdg,
        'p_true_p': magdf(slcdf.p.pfp.trk.truth.p.genp),
        'p_true_pdg': slcdf.p.pfp.trk.truth.p.pdg,
        'baseline': slcdf.slc.truth.baseline, 
        'nu_E_true': slcdf.slc.truth.E,
        'del_p': del_p,
        'del_Tp': del_Tp,
        'del_phi': del_phi,
        # 'tmatch_eff': slcdf.slc.tmatch.eff, 
        # 'tmatch_pur': slcdf.slc.tmatch.pur, 
        'tmatch_idx': tmatch_idx_series,
        'has_stub': slc_has_stub_series
    })

    # include some meta-data
    slcdf['detector'] = DETECTOR

    # Add in crt hit matching for ICARUS
    if DETECTOR == "ICARUS":
        crt = make_crthitdf(f)
        slcdf = slcdf.join(((crt.time > -1) & (crt.time < 1.8) & (crt.plane != 50)).groupby(level=[0]).any().rename("crthit"))
    else:
        slcdf["crthit"] = False

    # add in stub info, per range bin
    stubdf = stubdf[stubdf.plane == 2]

    stub_length_bins = [0, 0.5, 1, 2, 3, 4]
    stub_length_name = ["l0_5cm", "l1cm", "l2cm", "l3cm", "l4cm"]
    tosave = ["dedx", "charge"] # ["dedx_callo", "dedx_calhi", "Q", "length", "inc_charge"]

    for blo, bhi, name in zip(stub_length_bins[:-1], stub_length_bins[1:], stub_length_name):
        stub_tosave = stubdf.dedx[(stubdf.length > blo) & (stubdf.length < bhi)].groupby(level=[0,1]).idxmax()
        for col in tosave:
            s = stubdf.loc[stub_tosave, col]
            s.name = "stub_%s_%s" % (name, col)
            s.index = s.index.droplevel(-1)
            slcdf = slcdf.join(s, how="left", validate="one_to_one")

    return slcdf

gump_genie_systematics = [
    # CCQE
    "GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape",
    'GENIEReWeight_SBN_v1_multisigma_RPA_CCQE',
    'GENIEReWeight_SBN_v1_multisigma_CoulombCCQE',

    'GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE',
    'GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE',
    'GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE',
    'GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE',

    # MEC
    # "GENIEReWeight_SBN_v1_multisigma_NormNCMEC",
    'GENIEReWeight_SBN_v1_multisigma_NormCCMEC',
    'GENIEReWeight_SBN_v1_multisigma_NormNCMEC',
    "GENIEReWeight_SBN_v1_multisigma_DecayAngMEC",

    # RES
    "GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi",
    "GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad",
    "GENIEReWeight_SBN_v1_multisigma_MaCCRES",
    "GENIEReWeight_SBN_v1_multisigma_MaNCRES",
    "GENIEReWeight_SBN_v1_multisigma_MvCCRES",
    "GENIEReWeight_SBN_v1_multisigma_MvNCRES",

    # Non-Res

    # DIS
    # "GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse",
    'GENIEReWeight_SBN_v1_multisigma_AhtBY',
    'GENIEReWeight_SBN_v1_multisigma_BhtBY',
    'GENIEReWeight_SBN_v1_multisigma_CV1uBY',
    'GENIEReWeight_SBN_v1_multisigma_CV2uBY',

    # COH
    "GENIEReWeight_SBN_v1_multisigma_NormCCCOH",
    "GENIEReWeight_SBN_v1_multisigma_NormNCCOH",

    # FSI
    # "GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse",
    # "GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse",
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

    # NCEL
    'GENIEReWeight_SBN_v1_multisigma_MaNCEL',
    'GENIEReWeight_SBN_v1_multisigma_EtaNCEL',
]

# additional (re-weights)
gump_genie_reknob_systematics = gump_genie_systematics + [
    'ZExpPCAWeighter_SBNNuSyst_multisigma_D_ZExp_b1',
    'ZExpPCAWeighter_SBNNuSyst_multisigma_D_ZExp_b3',
    'ZExpPCAWeighter_SBNNuSyst_multisigma_D_ZExp_b2',
    'ZExpPCAWeighter_SBNNuSyst_multisigma_D_ZExp_b4',

    'ZExpPCAWeighter_SBNNuSyst_multisigma_MvA_ZExp_b1',
    'ZExpPCAWeighter_SBNNuSyst_multisigma_MvA_ZExp_b3',
    'ZExpPCAWeighter_SBNNuSyst_multisigma_MvA_ZExp_b2',
    'ZExpPCAWeighter_SBNNuSyst_multisigma_MvA_ZExp_b4',

    'CCQETemplateReweight_SBNNuSyst_multisigma_SF_q0bin1',
    'CCQETemplateReweight_SBNNuSyst_multisigma_SF_q0bin2',
    'CCQETemplateReweight_SBNNuSyst_multisigma_SF_q0bin3',
    'CCQETemplateReweight_SBNNuSyst_multisigma_SF_q0bin4',
    'CCQETemplateReweight_SBNNuSyst_multisigma_SF_q0bin5',

    'CCQETemplateReweight_SBNNuSyst_multisigma_CRPA_q0bin1',
    'CCQETemplateReweight_SBNNuSyst_multisigma_CRPA_q0bin2',
    'CCQETemplateReweight_SBNNuSyst_multisigma_CRPA_q0bin3',
    'CCQETemplateReweight_SBNNuSyst_multisigma_CRPA_q0bin4',
    'CCQETemplateReweight_SBNNuSyst_multisigma_CRPA_q0bin5',

    'QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_0',
    'QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_1',
    'QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_2',
    'QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_3',
    'QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_4',
    'QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_5',

    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrG4_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrINCL_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrG4LoE_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrINCLLoE_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrG4M1E_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrINCLM1E_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrG4M2E_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrINCLM2E_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrG4HiE_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrINCLHiE_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_MFPLoE_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_MFPM1E_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_MFPM2E_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_MFPHiE_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrKin_PiProFix_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrKin_PiProBias_N',

    'MECq0q3InterpWeighting_SBNNuSyst_multisigma_ValQ0Bin0_MECResponse',
    'MECq0q3InterpWeighting_SBNNuSyst_multisigma_ValQ0Bin1_MECResponse',
    'MECq0q3InterpWeighting_SBNNuSyst_multisigma_ValQ0Bin2_MECResponse',
    'MECq0q3InterpWeighting_SBNNuSyst_multisigma_ValQ0Bin3_MECResponse',
    'MECq0q3InterpWeighting_SBNNuSyst_multisigma_ValQ0Bin4_MECResponse',

    'MECq0q3InterpWeighting_SBNNuSyst_multisigma_MarQ0Bin0_MECResponse',
    'MECq0q3InterpWeighting_SBNNuSyst_multisigma_MarQ0Bin1_MECResponse',
    'MECq0q3InterpWeighting_SBNNuSyst_multisigma_MarQ0Bin2_MECResponse',
    'MECq0q3InterpWeighting_SBNNuSyst_multisigma_MarQ0Bin3_MECResponse',
    'MECq0q3InterpWeighting_SBNNuSyst_multisigma_MarQ0Bin4_MECResponse',
]

def make_gump_nuslimwgtdf(f):
    return make_mcnudf(f, include_weights=True, slim=True, genie_systematics=gump_genie_systematics)

def make_gump_nuwgtdf(f):
    return make_mcnudf(f, include_weights=True, slim=False, genie_systematics=gump_genie_systematics)

def make_gump_nurewgtdf(f):
    return make_mcnudf(f, include_weights=True, slim=False, genie_systematics=gump_genie_reknob_systematics)

def make_gump_nudf(f, is_slc=False):
    # note: setting is_slc to false results in pdg for the slice not being used
    # and instead only mcnu pdg info gets saved, but this excludes the -1 pdg for cosmic
    nudf = make_mcdf(f, slc_mcbranches, slc_mcprimbranches) if is_slc else make_mcdf(f)
    nudf["ind"] = nudf.index.get_level_values(1)

    # wgtdf = pd.concat([bnbsyst.bnbsyst(f, nudf.ind), geniesyst.geniesyst_sbnd(f, nudf.ind)], axis=1)
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det

    if det.empty:
        return pd.DataFrame()

    if (1 == det.unique()):
        DETECTOR = "SBND"
    elif (2 == det.unique()):
        DETECTOR = "ICARUS"
    else:
        print("Detector unclear, check rec.hdr.det!")

    is_fv = vtxfv_cut(nudf.position, DETECTOR)
    is_cc = nudf.iscc
    is_nc = (nudf.iscc == 0)
    genie_mode = nudf.genie_mode
    pdg = nudf.pdg
    nmu = nudf.nmu_27MeV
    np = nudf.np_50MeV
    npi = nudf.npi_30MeV
    npi0 = nudf.npi0
    nn = nudf.nn_0MeV
    is_1p0pi = (nudf.nmu_27MeV == 1) & (nudf.np_50MeV == 1) & (nudf.npi_30MeV == 0) & (nudf.npi0 == 0) 
    is_numu = (nudf.pdg == 14)
    is_other_numucc = (is_numu & is_cc & (is_1p0pi == 0) & is_fv)
    is_contained = trkfv_cut(nudf.mu.start, DETECTOR) & trkfv_cut(nudf.p.start, DETECTOR)
    is_sig = is_fv & is_1p0pi & is_numu & is_cc & is_contained

    nudf['nuint_categ'] = genie_mode 

    muon_p_series = magdf(nudf.mu.genp)
    proton_p_series = magdf(nudf.p.genp)
    
    true_tki = transverse_kinematics(muon_p_series, nudf.mu.genp, proton_p_series, nudf.mu.genp)
    true_del_p = true_tki['del_p']

    true_nu_E = nudf.E

    this_nudf = pd.DataFrame({
        'nu_E': true_nu_E,
        'del_p': true_del_p,
        'genie_mode': genie_mode, 
        'is_sig': is_sig, 
        'is_contained': is_contained,
        'is_nc': is_nc, 
        'is_other_numucc': is_other_numucc, 
        'is_fv': is_fv, 
        'pos_x' : nudf.position.x,
        'pos_y' : nudf.position.y,
        'pos_z' : nudf.position.z,
        'pdg': pdg,
        'nmu': nmu,
        'nn': nn,
        'np': np,
        'npi': npi,
        'npi0': npi0 
    })

    this_nudf.columns = pd.MultiIndex.from_tuples([(col, '') for col in this_nudf.columns])

    return this_nudf
