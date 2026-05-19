# Standard library imports
import os
import sys

# Third-party imports
import pandas as pd
import numpy as np

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
    "highYZEast": {
        "x": {"min": -200., "max": 0.},
        "y": {"min": -200., "max": 100},
        "z": {"min": 250., "max": 500.}
    },
    "highYZWest": {
        "x": {"min": 0., "max": 200.},
        "y": {"min": -200., "max": 200},
        "z": {"min": 250., "max": 500.}
    },
    "highYZ": {
        "x": {"min": -200., "max": 200.},
        "y": {"min": -200., "max": 100},
        "z": {"min": 250., "max": 500.}
    }
}

ICARUSRun2FVCuts = {
    "C0": {
        "x": {"min": -210.22, "max": -61.94}, # exluce EE in Run 2
        "y": {"min": -181.86, "max": 134.96},
        "z": {"min": -894.950652270838, "max": 894.950652270838}
    },
    "C1": {
        "x": {"min": 61.94, "max": 358.49},
        "y": {"min": -181.86, "max": 134.96},
        "z": {"min": -894.950652270838, "max": 894.950652270838}
    }
}

ICARUSRun4FVCuts = {
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

def true_fv_cut(df, det):
    vtx = pd.DataFrame({
                           'Run': df.Run,
                           'x': df.position.x,
                           'y': df.position.y,
                           'z': df.position.z}, index=df.index)
    return vtxfv_cut(vtx, det)

def slcfv_cut(df, det):
    vtx = pd.DataFrame({
                           'Run': df.Run,
                           'x': df.slc_vtx_x,
                           'y': df.slc_vtx_y,
                           'z': df.slc_vtx_z}, index=df.index)
    return vtxfv_cut(vtx, det)

def trkfv_cut(df, det):
    return _fv_cut(df, det, inzback=10)

def true_trkstartfv_cut(df, det):
    mu_vtx = pd.DataFrame({
                           'Run': df['Run'],
                           'x': df.mu.start.x,
                           'y': df.mu.start.y,
                           'z': df.mu.start.z}, index=df.index)

    p_vtx = pd.DataFrame({
                           'Run': df['Run'],
                           'x': df.mu.start.x,
                           'y': df.mu.start.y,
                           'z': df.mu.start.z}, index=df.index)

    return trkfv_cut(mu_vtx, det) & trkfv_cut(p_vtx, det)

def true_trkendfv_cut(df, det):
    mu_vtx = pd.DataFrame({
                           'Run': df['Run'],
                           'x': df.mu.end.x,
                           'y': df.mu.end.y,
                           'z': df.mu.end.z}, index=df.index)

    p_vtx = pd.DataFrame({
                           'Run': df['Run'],
                           'x': df.mu.end.x,
                           'y': df.mu.end.y,
                           'z': df.mu.end.z}, index=df.index)

    return trkfv_cut(mu_vtx, det) & trkfv_cut(p_vtx, det)

def trkstartfv_cut(df, det):
    vtx = pd.DataFrame({
                           'Run': df['Run'],
                           'x': df.pfp.trk.start.x,
                           'y': df.pfp.trk.start.y,
                           'z': df.pfp.trk.start.z}, index=df.index)
    return trkfv_cut(vtx, det)

def trkendfv_cut(df, det):
    vtx = pd.DataFrame({
                           'Run': df['Run'],
                           'x': df.pfp.trk.end.x,
                           'y': df.pfp.trk.end.y,
                           'z': df.pfp.trk.end.z}, index=df.index)
    return trkfv_cut(vtx, det)

def mufv_cut(df, det):
    vtx = pd.DataFrame({
                           'Run': df.Run,
                           'x': df.mu_end_x,
                           'y': df.mu_end_y,
                           'z': df.mu_end_z}, index=df.index)
    return trkfv_cut(vtx, det)

def pfv_cut(df, det):
    vtx = pd.DataFrame({
                           'Run': df.Run,
                           'x': df.p_end_x,
                           'y': df.p_end_y,
                           'z': df.p_end_z}, index=df.index)
    return trkfv_cut(vtx, det)

def _fv_cut(df, det, inx=10, iny=10, inzfront=10, inzback=50):
    if "ICARUS" in det:
        FVRun2 = (((df.x < (ICARUSRun2FVCuts['C0']['x']['max'] - inx)) & (df.x > (ICARUSRun2FVCuts['C0']['x']['min'] + inx))) |\
                ((df.x < (ICARUSRun2FVCuts['C1']['x']['max'] - inx)) & (df.x > (ICARUSRun2FVCuts['C1']['x']['min'] + inx)))) &\
                 (df.y < (ICARUSRun2FVCuts['C0']['y']['max'] - iny)) & (df.y > (ICARUSRun2FVCuts['C0']['y']['min'] + iny)) &\
                 (df.z < (ICARUSRun2FVCuts['C0']['z']['max'] - inzback)) & (df.z > (ICARUSRun2FVCuts['C0']['z']['min'] + inzfront))
        FVRun4 = (((df.x < (ICARUSRun4FVCuts['C0']['x']['max'] - inx)) & (df.x > (ICARUSRun4FVCuts['C0']['x']['min'] + inx))) |\
                ((df.x < (ICARUSRun4FVCuts['C1']['x']['max'] - inx)) & (df.x > (ICARUSRun4FVCuts['C1']['x']['min'] + inx)))) &\
                 (df.y < (ICARUSRun4FVCuts['C0']['y']['max'] - iny)) & (df.y > (ICARUSRun4FVCuts['C0']['y']['min'] + iny)) &\
                 (df.z < (ICARUSRun4FVCuts['C0']['z']['max'] - inzback)) & (df.z > (ICARUSRun4FVCuts['C0']['z']['min'] + inzfront))
        if det == "ICARUS":
            ret = FVRun2
            ret[df.Run == 4] = FVRun4[df.Run == 4]
            return ret
        elif det == "ICARUS Run2":
            return FVRun2
        elif det == "ICARUS Run4":
            return FVRun4
        else:
            raise NameError("DETECTOR not valid, should be SBND or ICARUS Run2 or ICARUS Run4")
    elif det == "SBND":
        return ((df.x < SBNDFVCuts['lowYZ']['x']['max'] - inx) & (df.x > SBNDFVCuts['lowYZ']['x']['min'] + inx) &\
                (df.y < SBNDFVCuts['lowYZ']['y']['max'] - iny) & (df.y > SBNDFVCuts['lowYZ']['y']['min'] + iny) &\
                (df.z < SBNDFVCuts['lowYZ']['z']['max']) & (df.z > SBNDFVCuts['lowYZ']['z']['min'] + inzfront)) |\
               ((df.x < SBNDFVCuts['highYZ']['x']['max'] - inx) & (df.x > SBNDFVCuts['highYZ']['x']['min'] + inx) &\
                (df.y < SBNDFVCuts['highYZ']['y']['max'] - iny) & (df.y > SBNDFVCuts['highYZ']['y']['min'] + iny) &\
                (df.z < SBNDFVCuts['highYZ']['z']['max'] - inzback) & (df.z > SBNDFVCuts['highYZ']['z']['min']))

        # should switch to this soon, but right now need to keep using older cuts
        #return ((df.x < SBNDFVCuts['lowYZ']['x']['max'] - inx) & (df.x > SBNDFVCuts['lowYZ']['x']['min'] + inx) &\
        #        (df.y < SBNDFVCuts['lowYZ']['y']['max'] - iny) & (df.y > SBNDFVCuts['lowYZ']['y']['min'] + iny) &\
        #        (df.z < SBNDFVCuts['lowYZ']['z']['max']) & (df.z > SBNDFVCuts['lowYZ']['z']['min'] + inzfront)) |\
        #       ((df.x < SBNDFVCuts['highYZWest']['x']['max'] - inx) & (df.x > SBNDFVCuts['highYZWest']['x']['min']) &\
        #        (df.y < SBNDFVCuts['highYZWest']['y']['max'] - iny) & (df.y > SBNDFVCuts['highYZWest']['y']['min'] + iny) &\
        #        (df.z < SBNDFVCuts['highYZWest']['z']['max'] - inzback) & (df.z > SBNDFVCuts['highYZWest']['z']['min'])) |\
        #       ((df.x < SBNDFVCuts['highYZEast']['x']['max']) & (df.x > SBNDFVCuts['highYZEast']['x']['min'] + inx) &\
        #        (df.y < SBNDFVCuts['highYZEast']['y']['max'] - iny) & (df.y > SBNDFVCuts['highYZEast']['y']['min'] + iny) &\
        #        (df.z < SBNDFVCuts['highYZEast']['z']['max'] - inzback) & (df.z > SBNDFVCuts['highYZEast']['z']['min']))

    else:
        raise NameError("DETECTOR not valid, should be SBND or ICARUS Run2 or ICARUS Run4")

def flash_cut(df, det):
    if det == "SBND":
        return df.flash_maxpe > 2000.
    elif det == "ICARUS Run2":
        return df.flash_maxpe > 6000.
    elif det == "ICARUS Run4":
        return df.flash_maxpe > 1000.
    elif det == "ICARUS":
        return ((df.flash_maxpe > 3000) & (df.Run == 2)) | ((df.flash_maxpe > 3000) & (df.Run == 4))


def cosmic_cut(df, is_old=False):
    if is_old:
        return (df.nu_score > 0.4)
    else:
        df = add_opening_angle_mu_p(df)
        return (df.nu_score > 0.4) & (df["mu_p_opening_angle_deg"] < 155)

def add_opening_angle_mu_p(df, out_col="mu_p_opening_angle_deg", degrees=True):
    mu = df[["mu_dir_x", "mu_dir_y", "mu_dir_z"]].to_numpy(dtype=float)
    p  = df[["p_dir_x",  "p_dir_y",  "p_dir_z"]].to_numpy(dtype=float)

    mu_norm = np.linalg.norm(mu, axis=1)
    p_norm  = np.linalg.norm(p,  axis=1)

    valid = np.isfinite(mu_norm) & np.isfinite(p_norm) & (mu_norm > 0) & (p_norm > 0)

    cosang = np.full(len(df), np.nan, dtype=float)
    dot = np.einsum("ij,ij->i", mu, p)
    cosang[valid] = dot[valid] / (mu_norm[valid] * p_norm[valid])
    cosang = np.clip(cosang, -1.0, 1.0)

    ang = np.arccos(cosang)
    if degrees:
        ang = np.degrees(ang)

    df[out_col] = ang
    return df

def del_p_cut(df):
    return (df.del_p <= 0.6)

def twoprong_cut(df):
    return (np.isnan(df.other_shw_length) & np.isnan(df.other_trk_length))

def pid_cut_df(df):
    return pid_cut(df.mu_chi2_of_mu_cand, df.mu_chi2_of_prot_cand,
        df.prot_chi2_of_mu_cand, df.prot_chi2_of_prot_cand, df.mu_len)

def pid_cut(mu_chi2_mu_cand, mu_chi2_prot_cand, prot_chi2_mu_cand,
            prot_chi2_prot_cand, mu_len, is_old=False):
    if is_old:
        MUSEL_MUSCORE_TH, MUSEL_PSCORE_TH, MUSEL_LEN_TH = 15, 90, 50
    else:
        MUSEL_MUSCORE_TH, MUSEL_PSCORE_TH, MUSEL_LEN_TH = 30, 80, 25

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
    ret = [var[df.is_sig == True],
           var[df.is_other_numucc == True],
           var[df.is_nc == True],
           var[df.is_fv == False],
           var[df.is_cosmic == True],
           var[(df.is_sig != True) & (df.is_other_numucc != True) & (df.is_nc != True) & (df.is_fv != False) & (df.is_cosmic != True)]
           ]
    return ret

def old_cuts(recodf, DETECTOR, det_run=False):
    if det_run:
        recodf['Run'] = det_run

    ## fv cut
    recodf = recodf[slcfv_cut(recodf, DETECTOR)]

    ### NuScore cut
    recodf = recodf[cosmic_cut(recodf, is_old=True)]

    ### Two prong cut
    recodf = recodf[twoprong_cut(recodf)]

    ### containment cut
    recodf = recodf[mufv_cut(recodf, DETECTOR)]
    recodf = recodf[pfv_cut(recodf, DETECTOR)]

    ### PID cut
    recodf = recodf[pid_cut(recodf.mu_chi2_of_mu_cand, recodf.mu_chi2_of_prot_cand,
                            recodf.prot_chi2_of_mu_cand, recodf.prot_chi2_of_prot_cand,
                            recodf.mu_len, is_old=True)]

    ### crthitveto cut
    if DETECTOR == "ICARUS":
        recodf = recodf[crthitveto_cut(recodf)]

    return recodf

def cathode_cut(df):
    p_start = df[['slc_vtx_x', 'slc_vtx_y', 'slc_vtx_z']].values
    p_mu = df[['mu_end_x', 'mu_end_y', 'mu_end_z']].values
    df = df[~intersects_prism_vectorized(p_start, p_mu, (-5., -200., 0.), (5., 200., 500.))]

    p_start = df[['slc_vtx_x', 'slc_vtx_y', 'slc_vtx_z']].values
    p_prot = df[['p_end_x', 'p_end_y', 'p_end_z']].values
    df = df[~intersects_prism_vectorized(p_start, p_prot, (-5., -200., 0.), (5., 200., 500.))]

    return df

def intersects_prism_vectorized(p1_array, p2_array, prism_min=(-200., 100., 250.), prism_max=(200., 200., 500.), solid=True):
    """
    Determines intersection for multiple line segments simultaneously.
    
    p1_array, p2_array: NumPy arrays of shape (N, 3)
    prism_min, prism_max: Tuples or arrays of (x, y, z)
    """
    p1 = np.asarray(p1_array)
    p2 = np.asarray(p2_array)
    p_min = np.asarray(prism_min)
    p_max = np.asarray(prism_max)

    # Initialize t_min and t_max for each segment
    t_min = np.zeros(len(p1))
    t_max = np.ones(len(p1))

    direction = p2 - p1

    inside_bool = np.zeros(len(p1))

    p_mins = np.array([p_min]*len(p1))
    p_maxs = np.array([p_max]*len(p1))

    if solid:
        inside_bool = ((p_mins < p1) & (p1 < p_maxs)).all(axis=1) | ((p_mins < p2) & (p2 < p_maxs)).all(axis=1)

    for i in range(3): # Iterate over X, Y, Z dimensions
        # Use a small epsilon to avoid true division by zero
        # or handle it via numpy's error handling
        inv_dir = 1.0 / np.where(direction[:, i] == 0, 1e-9, direction[:, i])

        tmin = (p_min[i] - p1[:, i]) * inv_dir
        tmax = (p_max[i] - p1[:, i]) * inv_dir

        t_near = np.minimum(tmin, tmax)
        t_far = np.maximum(tmin, tmax)

        # Update entry/exit for the entire batch
        t_min = np.maximum(t_min, t_near)
        t_max = np.minimum(t_max, t_far)
        # If the line is parallel to the axis (dir=0), manually check bounds
        parallel_mask = (direction[:, i] == 0)
        outside_bounds = (p1[:, i] < p_min[i]) | (p1[:, i] > p_max[i])

        # If parallel and outside, invalidate the t range so it returns False
        t_min = np.where(parallel_mask & outside_bounds, 1.0, t_min)
        t_max = np.where(parallel_mask & outside_bounds, 0.0, t_max)

    return (t_min <= t_max) | inside_bool

def all_cuts(recodf, DETECTOR, det_run=False):
    if det_run:
        recodf['Run'] = det_run

    ## fv cut
    recodf = recodf[slcfv_cut(recodf, DETECTOR)]

    ### cosmic cut
    recodf = recodf[cosmic_cut(recodf)]

    ### flash cut
    #recodf = recodf[flash_cut(recodf, DETECTOR)]

    ### Two prong cut
    recodf = recodf[twoprong_cut(recodf)]

    ### containment cut
    recodf = recodf[mufv_cut(recodf, DETECTOR)]
    recodf = recodf[pfv_cut(recodf, DETECTOR)]

    ### PID cut
    recodf = recodf[pid_cut(recodf.mu_chi2_of_mu_cand, recodf.mu_chi2_of_prot_cand,
                            recodf.prot_chi2_of_mu_cand, recodf.prot_chi2_of_prot_cand,
                            recodf.mu_len)]

    ### Cathode cut
    #recodf = cathode_cut(recodf)

    return recodf

