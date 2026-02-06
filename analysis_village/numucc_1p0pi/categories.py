import numpy as np
import pandas as pd
import sys
sys.path.append('../../')
from makedf.util import *


DETECTOR = "SBND_nohighyz"
# DETECTOR = "SBND"

# ==== definitions for event categories ===-

def IsNu(df):
    return (np.abs(df.mc.pdg) == 14) | (np.abs(df.mc.pdg) == 12)

def IsCosmic(df):
    return ~IsNu(df)

def IsNuOutFV(df):
    return IsNu(df) & ~InFV(df.mc.position, det=DETECTOR)

def IsNuInFV(df):
    return IsNu(df) & InFV(df.mc.position, det=DETECTOR)

def IsNuInFV_NuOther(df):
    return IsNuInFV(df) & (df.mc.iscc == 1) & (df.mc.pdg != 14)

def IsNuInFV_NumuNC(df):
    return IsNuInFV(df) & (df.mc.iscc == 0)

# ---- numu CC in FV, breakdown in topology
def Is_1p0pi(df):
    return (df.mc.nmu_220MeVc == 1) & (df.mc.np_300MeVc == 1) & (df.mc.npi_70MeVc == 0) & (df.mc.npi0 == 0) &\
        (np.sqrt(df.mc.mu.genp.x**2 + df.mc.mu.genp.y**2 + df.mc.mu.genp.z**2) < 1) &\
        (np.sqrt(df.mc.p.genp.x**2 + df.mc.p.genp.y**2 + df.mc.p.genp.z**2) < 1) &\
            InFV(df.mc.mu.start, det=DETECTOR) & InFV(df.mc.p.start, det=DETECTOR) &\
            InFV(df.mc.mu.end, det=DETECTOR) & InFV(df.mc.p.end, det=DETECTOR) 

def Is_Np0pi(df):
    return (df.mc.nmu_220MeVc == 1) & (df.mc.np_300MeVc > 1) & (df.mc.npi_70MeVc == 0) & (df.mc.npi0 == 0) 

def IsNuInFV_NumuCC_Other(df):
    return IsNuInFV(df) & (df.mc.pdg == 14) & (df.mc.iscc == 1) &\
              ~Is_1p0pi(df) & ~Is_Np0pi(df)

def IsNuInFV_NumuCC_Np0pi(df):
    return IsNuInFV(df) & (df.mc.pdg == 14) & (df.mc.iscc == 1) &\
              Is_Np0pi(df)

def IsNuInFV_NumuCC_1p0pi(df):
    return IsNuInFV(df) & (df.mc.pdg == 14) & (df.mc.iscc == 1) &\
              Is_1p0pi(df)

# --- numu CC in FV, breakdown in interaction mode (GENIE)
def IsNuInFV_NumuCC_QE(df):
    return IsNuInFV(df) & (df.mc.pdg == 14) & (df.mc.iscc == 1) &\
              (df.mc.genie_mode == 0)

def IsNuInFV_NumuCC_MEC(df):
    return IsNuInFV(df) & (df.mc.pdg == 14) & (df.mc.iscc == 1) &\
              (df.mc.genie_mode == 10)

def IsNuInFV_NumuCC_RES(df):
    return IsNuInFV(df) & (df.mc.pdg == 14) & (df.mc.iscc == 1) &\
              (df.mc.genie_mode == 1)

def IsNuInFV_NumuCC_DIS(df):
    return IsNuInFV(df) & (df.mc.pdg == 14) & (df.mc.iscc == 1) &\
              (df.mc.genie_mode == 2)

# def IsNuInFV_NumuCC_COH(df):
#     return IsNuInFV(df) & (df.pdg == 14) & (df.iscc == 1) &\
#               (df.genie_mode == 3)

def IsNuInFV_NumuCC_OtherMode(df):
    return IsNuInFV(df) & (df.mc.pdg == 14) & (df.mc.iscc == 1) &\
              ~(df.mc.genie_mode == 0) & ~(df.mc.genie_mode == 1) & ~(df.mc.genie_mode == 2) & ~(df.mc.genie_mode == 10)


# --- numu CC in FV, breakdown in interaction mode (GiBUU)
def IsNuInFV_NumuCC_QE_GiBUU(df):
    return IsNuInFV(df) & (df.mc.pdg == 14) & (df.mc.iscc == 1) &\
              (df.mc.genie_mode == 1)

def IsNuInFV_NumuCC_MEC_GiBUU(df):
    return IsNuInFV(df) & (df.mc.pdg == 14) & (df.mc.iscc == 1) &\
              ((df.mc.genie_mode == 35) | (df.mc.genie_mode == 36))

def IsNuInFV_NumuCC_RES_GiBUU(df):
    return IsNuInFV(df) & (df.mc.pdg == 14) & (df.mc.iscc == 1) &\
              ((df.mc.genie_mode >= 2) & (df.mc.genie_mode <= 31))

def IsNuInFV_NumuCC_DIS_GiBUU(df):
    return IsNuInFV(df) & (df.mc.pdg == 14) & (df.mc.iscc == 1) &\
              ((df.mc.genie_mode == 32) | (df.mc.genie_mode == 33) | (df.mc.genie_mode == 34) | (df.mc.genie_mode == 37))

def IsNuInFV_NumuCC_OtherMode_GiBUU(df):
    return IsNuInFV(df) & (df.mc.pdg == 14) & (df.mc.iscc == 1) &\
              ~IsNuInFV_NumuCC_QE_GiBUU(df) & ~IsNuInFV_NumuCC_MEC_GiBUU(df) & ~IsNuInFV_NumuCC_RES_GiBUU(df) & ~IsNuInFV_NumuCC_DIS_GiBUU(df)


# --- for event topoloy breakdown ---
def IsNu(df):
    return ~df.mc.pdg.isna()


def IsSignal(df): # definition                                                                                                                                                                                                                                                                         
    is_fv = InFV(df.mc.position, det=DETECTOR)
    is_1mu1p0pi = (df.mc.nmu_220MeVc == 1) & (df.mc.npi_70MeVc == 0) & (df.mc.np_300MeVc == 1) & (df.mc.npi0 == 0) & (df.mc.mu.totp < 1) & (df.mc.p.totp < 1) # & (df.np_20MeVc == 1) : add with stubs
    return is_fv & is_1mu1p0pi


def Is1muNp0pi(df): # definition                                                                                                                                                                                                                                                                         
    is_fv = InFV(df.mc.position, det=DETECTOR)
    is_1mu1p0pi = (df.mc.nmu_220MeVc == 1) & (df.mc.npi_70MeVc == 0) & (df.mc.np_300MeVc > 1) & (df.mc.npi0 == 0) & (df.mc.mu.totp < 1) & (df.mc.p.totp < 1) #& (df.mu.genE > 0.25) # & (df.np_20MeVc == 1) : add with stubs
    return is_fv & is_1mu1p0pi


def Is1muNcpi(df): # definition                                                                                                                                                                                                                                                                         
    is_fv = InFV(df.mc.position, det=DETECTOR)
    is_1mu1p0pi = (df.mc.nmu_220MeVc == 1) & (df.mc.npi_70MeVc > 0) & (df.mc.npi0 == 0) #& (df.mu.genE > 0.25) # & (df.np_20MeVc == 1) : add with stubs
    return is_fv & is_1mu1p0pi


# ==== functions to get breadkdowns of event categories ====

def get_int_category(df, ret_cuts=False, print_summary=False):
    cut_cosmic = IsCosmic(df)
    cut_nu_outfv = IsNuOutFV(df)
    cut_nu_infv_nu_other = IsNuInFV_NuOther(df)
    cut_nu_infv_numu_nc = IsNuInFV_NumuNC(df)
    cut_nu_infv_numu_cc_other = IsNuInFV_NumuCC_Other(df)
    cut_nu_infv_numu_cc_np0pi = IsNuInFV_NumuCC_Np0pi(df)
    cut_nu_infv_numu_cc_1p0pi = IsNuInFV_NumuCC_1p0pi(df)

    # assert there's no overlap between the categories, AND that all categories are covered just in case i messed something up...
    assert (cut_cosmic & cut_nu_outfv & cut_nu_infv_nu_other & cut_nu_infv_numu_nc & cut_nu_infv_numu_cc_other & cut_nu_infv_numu_cc_np0pi & cut_nu_infv_numu_cc_1p0pi).sum() == 0
    assert (cut_cosmic | cut_nu_outfv | cut_nu_infv_nu_other | cut_nu_infv_numu_nc | cut_nu_infv_numu_cc_other | cut_nu_infv_numu_cc_np0pi | cut_nu_infv_numu_cc_1p0pi).sum() == len(df)

    # category 1 NEEDS TO BE THE SIGNAL MODE
    nuint_categ = pd.Series(10, index=df.index)
    nuint_categ[cut_cosmic] = -1  # not nu
    nuint_categ[cut_nu_outfv] = 0  # nu out of FV
    nuint_categ[cut_nu_infv_numu_cc_1p0pi] = 1    # nu in FV, signal
    nuint_categ[cut_nu_infv_numu_cc_np0pi] = 2  # 1mu, 0cpi, Np, 0pi0
    nuint_categ[cut_nu_infv_numu_cc_other] = 3  # 1mu, Ncpi, 0pi0
    nuint_categ[cut_nu_infv_numu_nc] = 4  # nu in FV, numu NC
    nuint_categ[cut_nu_infv_nu_other] = 5  # nu in FV, other

    if print_summary:
        print(nuint_categ.value_counts())

    if ret_cuts:
        cuts = [cut_cosmic, cut_nu_outfv, cut_nu_infv_nu_other, cut_nu_infv_numu_nc, 
                cut_nu_infv_numu_cc_other, cut_nu_infv_numu_cc_np0pi, cut_nu_infv_numu_cc_1p0pi]
        return cuts

    return nuint_categ


def get_genie_category(df, ret_cuts=False, print_summary=False):
    cut_cosmic = IsCosmic(df)
    cut_nu_outfv = IsNuOutFV(df)
    cut_nu_infv_nu_other = IsNuInFV_NuOther(df)
    cut_nu_infv_numu_nc = IsNuInFV_NumuNC(df)
    # cut_nu_infv_numu_coh = IsNuInFV_NumuCC_COH(evtdf)
    cut_nu_infv_numu_othermode = IsNuInFV_NumuCC_OtherMode(df)
    cut_nu_infv_numu_cc_dis = IsNuInFV_NumuCC_DIS(df)
    cut_nu_infv_numu_cc_res = IsNuInFV_NumuCC_RES(df)
    cut_nu_infv_numu_cc_me = IsNuInFV_NumuCC_MEC(df)
    cut_nu_infv_numu_cc_qe = IsNuInFV_NumuCC_QE(df)

    assert (cut_cosmic & cut_nu_outfv & cut_nu_infv_nu_other & cut_nu_infv_numu_nc & cut_nu_infv_numu_othermode & cut_nu_infv_numu_cc_dis & cut_nu_infv_numu_cc_res & cut_nu_infv_numu_cc_me & cut_nu_infv_numu_cc_qe).sum() == 0
    assert (cut_cosmic | cut_nu_outfv | cut_nu_infv_nu_other | cut_nu_infv_numu_nc | cut_nu_infv_numu_othermode | cut_nu_infv_numu_cc_dis | cut_nu_infv_numu_cc_res | cut_nu_infv_numu_cc_me | cut_nu_infv_numu_cc_qe).sum() == len(df)

    genie_categ = pd.Series(10, index=df.index)
    genie_categ[cut_cosmic] = -1  # not nu
    genie_categ[cut_nu_outfv] = 0  # nu out of FV
    genie_categ[cut_nu_infv_numu_cc_qe] = 1  # nu in FV, QE
    genie_categ[cut_nu_infv_numu_cc_me] = 2  # nu in FV, MEC
    genie_categ[cut_nu_infv_numu_cc_res] = 3  # nu in FV, RES
    genie_categ[cut_nu_infv_numu_cc_dis] = 4  # nu in FV, DIS
    genie_categ[cut_nu_infv_numu_othermode] = 5  # nu in FV, other mode
    genie_categ[cut_nu_infv_numu_nc] = 6  # nu in FV, numu NC
    genie_categ[cut_nu_infv_nu_other] = 7  # nu in FV, other

    if print_summary:
        print(genie_categ.value_counts())

    if ret_cuts:
        cuts = [cut_cosmic, cut_nu_outfv, cut_nu_infv_nu_other, cut_nu_infv_numu_nc, 
                cut_nu_infv_numu_othermode, cut_nu_infv_numu_cc_dis, cut_nu_infv_numu_cc_res, 
                cut_nu_infv_numu_cc_me, cut_nu_infv_numu_cc_qe]
        return cuts

    return genie_categ

def get_genie_sb_category(df, ret_cuts=False, print_summary=False):
    cut_cosmic = IsCosmic(df)
    cut_nu_outfv = IsNuOutFV(df)
    cut_nu_infv_nu_other = IsNuInFV_NuOther(df)
    cut_nu_infv_numu_nc = IsNuInFV_NumuNC(df)
    # break down numu modes into signal and background
    cut_nu_infv_numu_othermode_s = IsNuInFV_NumuCC_OtherMode(df) & IsNuInFV_NumuCC_1p0pi(df)
    cut_nu_infv_numu_othermode_b = IsNuInFV_NumuCC_OtherMode(df) & ~IsNuInFV_NumuCC_1p0pi(df)
    cut_nu_infv_numu_cc_dis_s = IsNuInFV_NumuCC_DIS(df) & IsNuInFV_NumuCC_1p0pi(df)
    cut_nu_infv_numu_cc_dis_b = IsNuInFV_NumuCC_DIS(df) & ~IsNuInFV_NumuCC_1p0pi(df)
    cut_nu_infv_numu_cc_res_s = IsNuInFV_NumuCC_RES(df) & IsNuInFV_NumuCC_1p0pi(df)
    cut_nu_infv_numu_cc_res_b = IsNuInFV_NumuCC_RES(df) & ~IsNuInFV_NumuCC_1p0pi(df)
    cut_nu_infv_numu_cc_mec_s = IsNuInFV_NumuCC_MEC(df) & IsNuInFV_NumuCC_1p0pi(df)
    cut_nu_infv_numu_cc_mec_b = IsNuInFV_NumuCC_MEC(df) & ~IsNuInFV_NumuCC_1p0pi(df)
    cut_nu_infv_numu_cc_qe_s = IsNuInFV_NumuCC_QE(df) & IsNuInFV_NumuCC_1p0pi(df)
    cut_nu_infv_numu_cc_qe_b = IsNuInFV_NumuCC_QE(df) & ~IsNuInFV_NumuCC_1p0pi(df)

    assert (cut_cosmic & cut_nu_outfv & cut_nu_infv_nu_other & cut_nu_infv_numu_nc & cut_nu_infv_numu_othermode_s & cut_nu_infv_numu_othermode_b & cut_nu_infv_numu_cc_dis_s & cut_nu_infv_numu_cc_dis_b & cut_nu_infv_numu_cc_res_s & cut_nu_infv_numu_cc_res_b & cut_nu_infv_numu_cc_mec_s & cut_nu_infv_numu_cc_mec_b & cut_nu_infv_numu_cc_qe_s & cut_nu_infv_numu_cc_qe_b).sum() == 0
    assert (cut_cosmic | cut_nu_outfv | cut_nu_infv_nu_other | cut_nu_infv_numu_nc | cut_nu_infv_numu_othermode_s | cut_nu_infv_numu_othermode_b | cut_nu_infv_numu_cc_dis_s | cut_nu_infv_numu_cc_dis_b | cut_nu_infv_numu_cc_res_s | cut_nu_infv_numu_cc_res_b | cut_nu_infv_numu_cc_mec_s | cut_nu_infv_numu_cc_mec_b | cut_nu_infv_numu_cc_qe_s | cut_nu_infv_numu_cc_qe_b).sum() == len(df)

    cuts = [cut_cosmic, cut_nu_outfv, cut_nu_infv_nu_other, cut_nu_infv_numu_nc, 
            cut_nu_infv_numu_othermode_b, cut_nu_infv_numu_othermode_s,
            cut_nu_infv_numu_cc_dis_b, cut_nu_infv_numu_cc_dis_s,
            cut_nu_infv_numu_cc_res_b, cut_nu_infv_numu_cc_res_s,
            cut_nu_infv_numu_cc_mec_b, cut_nu_infv_numu_cc_mec_s,
            cut_nu_infv_numu_cc_qe_b, cut_nu_infv_numu_cc_qe_s]

    # if print_summary:
    #     print(genie_sb_categ.value_counts())

    if ret_cuts:
        return cuts

    # TODO: return category series
    return cuts


# --- colors & labelsfor plotting ---
# nu / cosmic breakdown
nu_cosmics_labels = ["Cosmic", r"Out-FV $\nu$", r"FV $\nu$"]
nu_cosmics_colors = ["gray", "C0", "C1"]

# signal / backgroundtopology breakdown
# the signal mode code MUST be the first item in the list for all the code below to work
topology_list = [1, # signal
                 2, 3, 4, 5, 
                 0,
                 -1 # cosmic
                 ]
topology_labels = [r"$\nu_{\mu}$ CC 1p0$\pi$", 
                    r"$\nu_{\mu}$ CC Np0$\pi$", r"$\nu_{\mu}$ CC Other", r"$\nu$ NC", r"Other $\nu$", 
                    r"Out-FV $\nu$", 
                    "Cosmic"]
topology_colors = ["mediumslateblue", 
                    "darkslateblue", "coral", 
                    "darkgreen", 
                    "crimson", 
                    "sienna", 
                    "gray"] 

# --- GENIE interaction mode breakdown ---
genie_mode_list = [1, # CCQE
                   2, 3, 4, 5, 6, 7, 
                   0,
                   -1   # cosmic
                   ]
genie_mode_labels = [r'$\nu_{\mu}$ CC QE', r'$\nu_{\mu}$ CC MEC', r'$\nu_{\mu}$ CC RES', r'$\nu_{\mu}$ CC SIS/DIS', r'$\nu_{\mu}$ CC Other', 
                     r"$\nu$ NC", 
                     r"FV other $\nu$", 
                     r"Out-FV $\nu$", 
                     "Cosmic"]
genie_mode_colors = ["#9b5580", "#390C1E", "#2c7c94", "#D88A3B", "#BFB17C", 
                     "darkgreen", 
                     "crimson", 
                     "sienna",
                     "gray"] 

# --- GENIE SB interaction mode breakdown ---
genie_sb_mode_labels = [r'$\nu_{\mu}$ CC QE', r'$\nu_{\mu}$ CC QE',
                     r'$\nu_{\mu}$ CC MEC', r'$\nu_{\mu}$ CC MEC', 
                     r'$\nu_{\mu}$ CC RES', r'$\nu_{\mu}$ CC RES', 
                     r'$\nu_{\mu}$ CC SIS/DIS', r'$\nu_{\mu}$ CC SIS/DIS',
                     r'$\nu_{\mu}$ CC Other', r'$\nu_{\mu}$ CC Other',
                     r"$\nu$ NC", 
                     r"FV other $\nu$", 
                     r"Out-FV $\nu$", 
                     "Cosmic"]

genie_sb_mode_colors = ["#9b5580", "#9b5580",
                        "#390C1E", "#390C1E",
                        "#2c7c94", "#2c7c94",
                        "#D88A3B", "#D88A3B",
                        "#BFB17C", "#BFB17C",
                        "darkgreen", 
                        "crimson", 
                        "sienna",
                        "gray"] 


# --- GiBUU interaction mode breakdown ---
gibuu_mode_labels = [r'$\nu_{\mu}$ CC QE', r'$\nu_{\mu}$ CC QE (2p2h)', r'$\nu_{\mu}$ CC RES', r'$\nu_{\mu}$ CC DIS', r'$\nu_{\mu}$ CC Other', 
                     r"$\nu$ NC", r"FV other $\nu$", r"Out-FV $\nu$", "Cosmic"]
gibuu_mode_colors = ["#9b5580", "#390C1E", "#2c7c94", "#D88A3B", "#BFB17C", 
                     "darkgreen", "crimson", "sienna","gray"] 