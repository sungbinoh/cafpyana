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
        "x": {"min": -190., "max": 190.},
        "y": {"min": -190., "max": 190.},
        "z": {"min": 10., "max": 250.}
    },
    "highYZ": {
        "x": {"min": -190., "max": 190.},
        "y": {"min": -190., "max": 100},
        "z": {"min": 250., "max": 450.}
    }
}

ICARUSFVCuts = {
    "C0": {
        "x": {"min": -358.49 + 10, "max": -61.94 - 10},
        "y": {"min": -181.86 + 10, "max": 134.96 - 10},
        "z": {"min": -894.950652270838 + 10, "max": 894.950652270838 - 50}
    },
    "C1": {
        "x": {"min": 61.94, "max": 358.49},
        "y": {"min": -181.86 + 10, "max": 134.96 - 10},
        "z": {"min": -894.950652270838 + 10, "max": 894.950652270838 - 50}
    }
}

def fv_cut(df, det):
    if det == "ICARUS":
        return (((df.x < ICARUSFVCuts['C0']['x']['max']) & (df.x > ICARUSFVCuts['C0']['x']['min'])) |\
                ((df.x < ICARUSFVCuts['C1']['x']['max']) & (df.x > ICARUSFVCuts['C1']['x']['min']))) &\
                 (df.y < ICARUSFVCuts['C0']['y']['max']) & (df.y > ICARUSFVCuts['C0']['y']['min']) &\
                 (df.z < ICARUSFVCuts['C0']['z']['max']) & (df.z > ICARUSFVCuts['C0']['z']['min'])


    elif det == "SBND":

        return ((df.x < SBNDFVCuts['lowYZ']['x']['max']) & (df.x > SBNDFVCuts['lowYZ']['x']['min']) &\
                (df.y < SBNDFVCuts['lowYZ']['y']['max']) & (df.y > SBNDFVCuts['lowYZ']['y']['min']) &\
                (df.z < SBNDFVCuts['lowYZ']['z']['max']) & (df.z > SBNDFVCuts['lowYZ']['z']['min'])) |\
               ((df.x < SBNDFVCuts['highYZ']['x']['max']) & (df.x > SBNDFVCuts['highYZ']['x']['min']) &\
                (df.y < SBNDFVCuts['highYZ']['y']['max']) & (df.y > SBNDFVCuts['highYZ']['y']['min']) &\
                (df.z < SBNDFVCuts['highYZ']['z']['max']) & (df.z > SBNDFVCuts['highYZ']['z']['min']))

    else:
        raise NameError("DETECTOR not valid, should be SBND or ICARUS")

def cosmic_cut(df):
    return (df.nu_score > 0.5)
    #return (df.crlongtrkdiry > 0.5)

def twoprong_cut(df):
    return ((df.p_trackScore > 0.5) & (df.mu_trackScore > 0.5) & np.isnan(df.other_shw_length) & np.isnan(df.other_trk_length))

def pid_cut(mu_chi2_of_mu_cand, mu_chi2_of_prot_cand, prot_chi2_of_mu_cand,
            prot_chi2_of_prot_cand, mu_len):

    MUSEL_MUSCORE_TH, MUSEL_PSCORE_TH, MUSEL_LEN_TH = 25, 100, 50
    mu_cut = (mu_chi2_of_mu_cand < MUSEL_MUSCORE_TH) & \
             (prot_chi2_of_mu_cand > MUSEL_PSCORE_TH) & \
             (mu_len > MUSEL_LEN_TH)

    PSEL_MUSCORE_TH, PSEL_PSCORE_TH = 0, 90
    p_cut = (mu_chi2_of_prot_cand > PSEL_MUSCORE_TH) & \
            (prot_chi2_of_prot_cand < PSEL_PSCORE_TH)

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
    ret = [var[df.is_sig == True],
           var[df.is_other_numucc == True],
           var[df.is_nc == True],
           var[df.is_fv == False],
           var[df.is_cosmic == True],
           var[(df.is_sig != True) & (df.is_other_numucc != True) & (df.is_nc != True) & (df.is_fv != False) & (df.is_cosmic != True)]
           ]
    return ret
 
