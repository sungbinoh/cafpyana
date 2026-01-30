from . import getsyst
import pandas as pd
import numpy as np

# Regen systematic variations
# use mulitisim where possible
regen_systematics = [
    # CCQE
    "GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape",
    'GENIEReWeight_SBN_v1_multisigma_RPA_CCQE',
    'GENIEReWeight_SBN_v1_multisigma_CoulombCCQE',

    'GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE',
    'GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE',
    'GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE',
    'GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE',
    # 'ZExpPCAWeighter_SBNnusyst_b1',
    # 'ZExpPCAWeighter_SBNnusyst_b2', 
    # 'ZExpPCAWeighter_SBNnusyst_b3',
    # 'ZExpPCAWeighter_SBNnusyst_b4'

    # "GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse",

    # MEC
    # "GENIEReWeight_SBN_v1_multisigma_NormNCMEC",
    'GENIEReWeight_SBN_v1_multisigma_NormCCMEC',
    'GENIEReWeight_SBN_v1_multisigma_NormNCMEC',
    "GENIEReWeight_SBN_v1_multisigma_DecayAngMEC",

    # RES
    # "GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse",
    # "GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse",
    'GENIEReWeight_SBN_v1_multisim_RDecBR1gamma',
    'GENIEReWeight_SBN_v1_multisim_RDecBR1eta',
    "GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi",
    "GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad",

    "GENIEReWeight_SBN_v1_multisigma_MaCCRES",
    "GENIEReWeight_SBN_v1_multisigma_MaNCRES",
    "GENIEReWeight_SBN_v1_multisigma_MvCCRES",
    "GENIEReWeight_SBN_v1_multisigma_MvNCRES",

    # Non-Res
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

    # DIS
    # "GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse",
    'GENIEReWeight_SBN_v1_multisigma_AhtBY',
    'GENIEReWeight_SBN_v1_multisigma_BhtBY',
    'GENIEReWeight_SBN_v1_multisigma_CV1uBY',
    'GENIEReWeight_SBN_v1_multisigma_CV2uBY',

    # COH
    "GENIEReWeight_SBN_v1_multisigma_NormCCCOH", # Handled by re-tuning
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
    # "GENIEReWeight_SBN_v1_multisim_NCELVariationResponse",
    'GENIEReWeight_SBN_v1_multisigma_MaNCEL',
    'GENIEReWeight_SBN_v1_multisigma_EtaNCEL',
]

def geniesyst(f, nuind, multisim_nuniv=100, slim=False, systematics=None):
    if systematics is None:
        systematics = regen_systematics

    geniewgtdf = getsyst.getsyst(f, systematics, nuind, multisim_nuniv=multisim_nuniv, slim=slim, slimname="GENIE")

    if slim:  # keep only the multiplied "GENIE.univ_" columns
        genie_cols = [c for c in geniewgtdf.columns if c[0] == "GENIE"]
        geniewgtdf = geniewgtdf[genie_cols]
        
    return geniewgtdf
