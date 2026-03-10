from . import getsyst
import pandas as pd
import numpy as np

#AR23+ systematic variations
r23p_genie_systematics = [
    'ZExpPCAWeighter_SBNNuSyst_multisigma_MvA_ZExp_b1',
    'ZExpPCAWeighter_SBNNuSyst_multisigma_MvA_ZExp_b2',
    'ZExpPCAWeighter_SBNNuSyst_multisigma_MvA_ZExp_b3',
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
    'QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_0',
    'QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_1',
    'QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_2',
    'QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_3',
    'QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_4',
    'QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_5',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_VecFFCCQEshape',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_CoulombCCQE',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_NormCCMEC',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_NormNCMEC',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_DecayAngMEC',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_MFP_pi',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrCEx_pi',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrInel_pi',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrAbs_pi',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrPiProd_pi',
    'MECq0q3InterpWeighting_SuSAv2ToValenica_q0binned_MECResponse_q0bin0',
    'MECq0q3InterpWeighting_SuSAv2ToValenica_q0binned_MECResponse_q0bin1',
    'MECq0q3InterpWeighting_SuSAv2ToValenica_q0binned_MECResponse_q0bin2',
    'MECq0q3InterpWeighting_SuSAv2ToValenica_q0binned_MECResponse_q0bin3',
    'MECq0q3InterpWeighting_SuSAv2ToMartini_q0binned_MECResponse_q0bin0',
    'MECq0q3InterpWeighting_SuSAv2ToMartini_q0binned_MECResponse_q0bin1',
    'MECq0q3InterpWeighting_SuSAv2ToMartini_q0binned_MECResponse_q0bin2',
    'MECq0q3InterpWeighting_SuSAv2ToMartini_q0binned_MECResponse_q0bin3'
]


# grouped syst knobs
qe_genie_systematics = [
    'GENIEReWeight_SBN_v1_multisim_RPA_CCQE',
    'GENIEReWeight_SBN_v1_multisim_CoulombCCQE',
]

mec_genie_systematics = [
    'GENIEReWeight_SBN_v1_multisim_NormCCMEC',
    'GENIEReWeight_SBN_v1_multisim_NormNCMEC',
    "GENIEReWeight_SBN_v1_multisigma_DecayAngMEC",
]

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

nonres_genie_systematics = [
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

dis_genie_systematics = [
    'GENIEReWeight_SBN_v1_multisigma_AhtBY',
    'GENIEReWeight_SBN_v1_multisigma_BhtBY',
    'GENIEReWeight_SBN_v1_multisigma_CV1uBY',
    'GENIEReWeight_SBN_v1_multisigma_CV2uBY',
]

other_genie_systematics = [
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


def geniesyst(f, nuind, multisim_nuniv=100, slim=False, systematics=None):
    if systematics is None:
        systematics = r23p_genie_systematics + qe_genie_systematics + mec_genie_systematics + res_genie_systematics + nonres_genie_systematics + dis_genie_systematics + other_genie_systematics

    geniewgtdf = getsyst.getsyst(f, systematics, nuind, multisim_nuniv=multisim_nuniv, slim=slim, slimname="GENIE")

    if slim:  # keep only the multiplied "GENIE.univ_" columns
        genie_cols = [c for c in geniewgtdf.columns if c[0] == "GENIE" or "SBNNuSyst" in c or "SuSAv2" in c]
        geniewgtdf = geniewgtdf[genie_cols]
    return geniewgtdf