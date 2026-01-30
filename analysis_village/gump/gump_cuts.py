# Standard library imports
import os
import sys

# Third-party imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Add the head direcoty to sys.path
workspace_root = os.getcwd()
sys.path.insert(0, workspace_root + "/../../")

# Local imports
import analysis_village.gump.kinematics
from makedf.util import *

# Fiducial volume cuts for SBND and ICARUS
SBNDFVCuts = {
    "lowYZ": {
        "x": {"min": -200., "max": 200.},
        "y": {"min": -200., "max": 200.},
        "z": {"min": 0., "max": 250.}
    },
    "highYZ": {
        "x": {"min": -200., "max": 200.},
        "y": {"min": -200., "max": 100},
        "z": {"min": 250., "max": 500.}
    }
}

ICARUSFVCuts = {
    "C0": {
        "x": {"min": -358.49, "max": -61.94},
        "y": {"min": -181.86, "max": 134.96},
        "z": {"min": -894.950652270838, "max": 894.950652270838}
    },
    "C1": {
        "x": {"min": 61.94, "max": 358.49},
        "y": {"min": -181.86, "max": 134.96},
        "z": {"min": -894.950652270838, "max": 894.950652270838}
    }
}

def vtxfv_cut(df, det):
    return _fv_cut(df, det, inzback=50)

def slcfv_cut(df, det):
    vtx = pd.DataFrame({'x': df.slc_vtx_x,
                           'y': df.slc_vtx_y,
                           'z': df.slc_vtx_z})
    return vtxfv_cut(vtx, det)

def trkfv_cut(df, det):
    return _fv_cut(df, det, inzback=10)

def mufv_cut(df, det):
    vtx = pd.DataFrame({'x': df.mu_end_x,
                           'y': df.mu_end_y,
                           'z': df.mu_end_z})
    return trkfv_cut(vtx, det)

def pfv_cut(df, det):
    vtx = pd.DataFrame({'x': df.p_end_x,
                           'y': df.p_end_y,
                           'z': df.p_end_z})
    return trkfv_cut(vtx, det)

def _fv_cut(df, det, inx=10, iny=10, inzfront=10, inzback=50):
    if det == "ICARUS":
        return (((df.x < (ICARUSFVCuts['C0']['x']['max'] - inx)) & (df.x > (ICARUSFVCuts['C0']['x']['min'] + inx))) |\
                ((df.x < (ICARUSFVCuts['C1']['x']['max'] - inx)) & (df.x > (ICARUSFVCuts['C1']['x']['min'] + inx)))) &\
                 (df.y < (ICARUSFVCuts['C0']['y']['max'] - iny)) & (df.y > (ICARUSFVCuts['C0']['y']['min'] + iny)) &\
                 (df.z < (ICARUSFVCuts['C0']['z']['max'] - inzback)) & (df.z > (ICARUSFVCuts['C0']['z']['min'] + inzfront))


    elif det == "SBND":

        return ((df.x < SBNDFVCuts['lowYZ']['x']['max'] - inx) & (df.x > SBNDFVCuts['lowYZ']['x']['min'] + inx) &\
                (df.y < SBNDFVCuts['lowYZ']['y']['max'] - iny) & (df.y > SBNDFVCuts['lowYZ']['y']['min'] + iny) &\
                (df.z < SBNDFVCuts['lowYZ']['z']['max'] - inzback) & (df.z > SBNDFVCuts['lowYZ']['z']['min']) + inzfront) |\
               ((df.x < SBNDFVCuts['highYZ']['x']['max'] - inx) & (df.x > SBNDFVCuts['highYZ']['x']['min'] + inx) &\
                (df.y < SBNDFVCuts['highYZ']['y']['max'] - iny) & (df.y > SBNDFVCuts['highYZ']['y']['min'] + iny) &\
                (df.z < SBNDFVCuts['highYZ']['z']['max'] - inzback) & (df.z > SBNDFVCuts['highYZ']['z']['min']) + inzfront)

    else:
        raise NameError("DETECTOR not valid, should be SBND or ICARUS")

def cosmic_cut(df):
    return (df.nu_score > 0.4)

def twoprong_cut(df):
    return (np.isnan(df.other_shw_length) & np.isnan(df.other_trk_length))

def pid_cut_df(df):
    return pid_cut(df.mu_chi2_of_mu_cand, df.mu_chi2_of_prot_cand,
        df.prot_chi2_of_mu_cand, df.prot_chi2_of_prot_cand, df.mu_len)

def pid_cut(mu_chi2_mu_cand, mu_chi2_prot_cand, prot_chi2_mu_cand,
            prot_chi2_prot_cand, mu_len):

    MUSEL_MUSCORE_TH, MUSEL_PSCORE_TH, MUSEL_LEN_TH = 15, 90, 50
    mu_cut = (mu_chi2_mu_cand < MUSEL_MUSCORE_TH) & \
             (prot_chi2_mu_cand > MUSEL_PSCORE_TH) & \
             (mu_len > MUSEL_LEN_TH)

    PSEL_MUSCORE_TH, PSEL_PSCORE_TH = 0, 90
    p_cut = (mu_chi2_prot_cand > PSEL_MUSCORE_TH) & \
            (prot_chi2_prot_cand < PSEL_PSCORE_TH)

    return mu_cut & p_cut

def stub_cut(df):
    cut = (df.has_stub == 0)
    return cut

def clear_cosmic_cut(df):
    cut = (df.is_clear_cosmic == 0)
    return cut

def contained_cut(df):
    cut = (df.is_contained == 1)
    return cut

def crthitveto_cut(df):
    return ~df.crthit

mode_list = [0, 10, 1, 2, 3]
mode_labels = ['QE', 'MEC', 'RES', 'SIS/DIS', 'COH', "other"]

def breakdown_mode(var, df):
    """Break down variable by interaction mode."""
    ret = [var[df.genie_mode == i] for i in mode_list]
    ret.append(var[sum([df.genie_mode == i for i in mode_list]) == 0])
    return ret

top_labels = ["Signal",
              "Other numu CC",
              "NC",
              "Out of FV",
              "Cosmic",
              "Other"]

def breakdown_top(var, df):
    print(df.is_cosmic)
    ret = [var[df.is_sig == True],
           var[df.is_other_numucc == True],
           var[df.is_nc == True],
           var[df.is_fv == False],
           var[df.is_cosmic == True],
           var[(df.is_sig != True) & (df.is_other_numucc != True) & (df.is_nc != True) & (df.is_fv != False) & (df.is_cosmic != True)]
           ]
    return ret
