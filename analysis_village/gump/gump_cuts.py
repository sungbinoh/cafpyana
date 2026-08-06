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

def vtxfv_cut(df):
    return _fv_cut(df, inzback=50)

def true_fv_cut(df):
    vtx = pd.DataFrame({
                           'detector': df.detector,
                           'Run': df.Run,
                           'x': df.position.x,
                           'y': df.position.y,
                           'z': df.position.z}, index=df.index)
    return vtxfv_cut(vtx)

def slcfv_cut(df):
    vtx = pd.DataFrame({
                           'detector': df.detector,
                           'Run': df.Run,
                           'x': df.slc_vtx_x,
                           'y': df.slc_vtx_y,
                           'z': df.slc_vtx_z}, index=df.index)
    return vtxfv_cut(vtx)

def trkfv_cut(df):
    return _fv_cut(df, inzback=10)

def true_trkstartfv_cut(df):
    mu_vtx = pd.DataFrame({
                           'detector': df.detector,
                           'Run': df.Run,
                           'x': df.mu.start.x,
                           'y': df.mu.start.y,
                           'z': df.mu.start.z}, index=df.index)

    p_vtx = pd.DataFrame({
                           'detector': df.detector,
                           'Run': df.Run,
                           'x': df.mu.start.x,
                           'y': df.mu.start.y,
                           'z': df.mu.start.z}, index=df.index)

    return trkfv_cut(mu_vtx) & trkfv_cut(p_vtx)

def true_trkendfv_cut(df):
    mu_vtx = pd.DataFrame({
                           'detector': df.detector,
                           'Run': df.Run,
                           'x': df.mu.end.x,
                           'y': df.mu.end.y,
                           'z': df.mu.end.z}, index=df.index)

    p_vtx = pd.DataFrame({
                           'detector': df.detector,
                           'Run': df.Run,
                           'x': df.mu.end.x,
                           'y': df.mu.end.y,
                           'z': df.mu.end.z}, index=df.index)

    return trkfv_cut(mu_vtx) & trkfv_cut(p_vtx)

def trkstartfv_cut(df):
    vtx = pd.DataFrame({
                           'detector': df.detector,
                           'Run': df.Run,
                           'x': df.pfp.trk.start.x,
                           'y': df.pfp.trk.start.y,
                           'z': df.pfp.trk.start.z}, index=df.index)
    return trkfv_cut(vtx)

def trkendfv_cut(df):
    vtx = pd.DataFrame({
                           'detector': df.detector,
                           'Run': df.Run,
                           'x': df.pfp.trk.end.x,
                           'y': df.pfp.trk.end.y,
                           'z': df.pfp.trk.end.z}, index=df.index)
    return trkfv_cut(vtx)

def mufv_cut(df):
    vtx = pd.DataFrame({
                           'detector': df.detector,
                           'Run': df.Run,
                           'x': df.mu_end_x,
                           'y': df.mu_end_y,
                           'z': df.mu_end_z}, index=df.index)
    return trkfv_cut(vtx)

def pfv_cut(df):
    vtx = pd.DataFrame({
                           'detector': df.detector,
                           'Run': df.Run,
                           'x': df.p_end_x,
                           'y': df.p_end_y,
                           'z': df.p_end_z}, index=df.index)
    return trkfv_cut(vtx)

def _fv_cut(df, inx=10, iny=10, inzfront=10, inzback=50):

    det_col = "detector"
    if det_col not in df.columns:
        raise KeyError(
            f"Could not find a detector column ('det' or 'detector') in the DataFrame."
        )

    valid_detectors = {"ICARUS Run2", "ICARUS Run4", "ICARUS", "SBND"}
    
    present_detectors = set(df[det_col].unique())
    invalid_detectors = present_detectors - valid_detectors

    if invalid_detectors:
        bad_rows = df[df[det_col].isin(invalid_detectors)].index.tolist()[:5]
        raise ValueError(
            f"DETECTOR type not valid! Found unrecognized labels: {invalid_detectors} "
            f"at dataframe rows {bad_rows}. Must be strictly 'SBND', 'ICARUS', 'ICARUS Run2', or 'ICARUS Run4'."
        )

    FVRun2 = (((df.x < (ICARUSRun2FVCuts['C0']['x']['max'] - inx)) & (df.x > (ICARUSRun2FVCuts['C0']['x']['min'] + inx))) |\
            ((df.x < (ICARUSRun2FVCuts['C1']['x']['max'] - inx)) & (df.x > (ICARUSRun2FVCuts['C1']['x']['min'] + inx)))) &\
             (df.y < (ICARUSRun2FVCuts['C0']['y']['max'] - iny)) & (df.y > (ICARUSRun2FVCuts['C0']['y']['min'] + iny)) &\
             (df.z < (ICARUSRun2FVCuts['C0']['z']['max'] - inzback)) & (df.z > (ICARUSRun2FVCuts['C0']['z']['min'] + inzfront))

    FVRun4 = (((df.x < (ICARUSRun4FVCuts['C0']['x']['max'] - inx)) & (df.x > (ICARUSRun4FVCuts['C0']['x']['min'] + inx))) |\
            ((df.x < (ICARUSRun4FVCuts['C1']['x']['max'] - inx)) & (df.x > (ICARUSRun4FVCuts['C1']['x']['min'] + inx)))) &\
             (df.y < (ICARUSRun4FVCuts['C0']['y']['max'] - iny)) & (df.y > (ICARUSRun4FVCuts['C0']['y']['min'] + iny)) &\
             (df.z < (ICARUSRun4FVCuts['C0']['z']['max'] - inzback)) & (df.z > (ICARUSRun4FVCuts['C0']['z']['min'] + inzfront))

    FVSBND = ((df.x < SBNDFVCuts['lowYZ']['x']['max'] - inx) & (df.x > SBNDFVCuts['lowYZ']['x']['min'] + inx) &\
            (df.y < SBNDFVCuts['lowYZ']['y']['max'] - iny) & (df.y > SBNDFVCuts['lowYZ']['y']['min'] + iny) &\
            (df.z < SBNDFVCuts['lowYZ']['z']['max']) & (df.z > SBNDFVCuts['lowYZ']['z']['min'] + inzfront)) |\
           ((df.x < SBNDFVCuts['highYZ']['x']['max'] - inx) & (df.x > SBNDFVCuts['highYZ']['x']['min'] + inx) &\
            (df.y < SBNDFVCuts['highYZ']['y']['max'] - iny) & (df.y > SBNDFVCuts['highYZ']['y']['min'] + iny) &\
            (df.z < SBNDFVCuts['highYZ']['z']['max'] - inzback) & (df.z > SBNDFVCuts['highYZ']['z']['min']))

    conditions = [
        (df[det_col] == "ICARUS Run2")
        | ((df[det_col] == "ICARUS") & (df.Run == 2)),
        (df[det_col] == "ICARUS Run4")
        | ((df[det_col] == "ICARUS") & (df.Run == 4)),
        (df[det_col] == "SBND"),
    ]

    choices = [
        FVRun2,
        FVRun4,
        FVSBND,
    ]

    np_mask = np.select(conditions, choices, default=False)

    return pd.Series(np_mask, index=df.index) 

def flash_cut(df):
    det_col = "detector"
    if det_col not in df.columns:
        raise KeyError(
            f"Could not find a detector column ('det' or 'detector') in the DataFrame."
        )

    valid_detectors = {"ICARUS Run2", "ICARUS Run4", "ICARUS", "SBND"}
    present_detectors = set(df[det_col].unique())
    invalid_detectors = present_detectors - valid_detectors

    if invalid_detectors:
        bad_rows = df[df[det_col].isin(invalid_detectors)].index.tolist()[:5]
        raise ValueError(
            f"DETECTOR type not valid in flash_cut! Found: {invalid_detectors} "
            f"at dataframe rows {bad_rows}."
        )

    conditions = [
        df[det_col] == "SBND",
        df[det_col] == "ICARUS Run2",
        df[det_col] == "ICARUS Run4",
        df[det_col] == "ICARUS",
    ]

    choices = [
        df.flash_maxpe > 2000.0,  # SBND
        df.flash_maxpe > 5000.0,  # ICARUS Run2
        df.flash_maxpe > 1000.0,  # ICARUS Run4
        ((df.flash_maxpe > 5000.0) & (df.Run == 2))
        | ((df.flash_maxpe > 1000.0) & (df.Run == 4)),  # ICARUS (Generic)
    ]

    np_mask = np.select(conditions, choices, default=False)

    return pd.Series(np_mask, index=df.index) 

# Final optimized selection
# opt to ratio, min mu len 40, shared track score, per-detector PID
# opt used the most updated calo treatment (no EMB IC, no sqmear15, double SBND 13)

SBND_CUTS = {
    "nu_score_th": 0.35,
    "max_opening_angle": 160,
    "musel_track_score_min": 0.5,
    "musel_muscore_th": 38,
    "musel_pscore_th": 82,
    "musel_len_th_min": 40,
    "musel_len_th_max": 400,
    "psel_muscore_th": 0,
    "psel_pscore_th": 141,
}

ICARUS_CUTS = {
    "nu_score_th": 0.35,
    "max_opening_angle": 160,
    "musel_track_score_min": 0.5,
    "musel_muscore_th": 111,
    "musel_pscore_th": 74,
    "musel_len_th_min": 40,
    "musel_len_th_max": 400,
    "psel_muscore_th": 0,
    "psel_pscore_th": 92,
}

def _det_cut_th(detector, key):
    return np.where(detector == "SBND", SBND_CUTS[key], ICARUS_CUTS[key])

def cosmic_cut(df):
    df = add_opening_angle_mu_p(df)
    return (df.nu_score > _det_cut_th(df.detector, "nu_score_th")) & \
           (df["mu_p_opening_angle_deg"] < _det_cut_th(df.detector, "max_opening_angle"))

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

def pid_cut(df):
    return pid_cut_df(df.mu_chi2_of_mu_cand, df.mu_chi2_of_prot_cand,
        df.prot_chi2_of_mu_cand, df.prot_chi2_of_prot_cand, df.mu_len,
        df.mu_track_score, df.detector)

def pid_cut_df(mu_chi2_mu_cand, mu_chi2_prot_cand, prot_chi2_mu_cand,
            prot_chi2_prot_cand, mu_len, mu_track_score, detector):
    mu_cut = (mu_track_score > _det_cut_th(detector, "musel_track_score_min")) & \
             (mu_chi2_mu_cand < _det_cut_th(detector, "musel_muscore_th")) & \
             (prot_chi2_mu_cand > _det_cut_th(detector, "musel_pscore_th")) & \
             (mu_len > _det_cut_th(detector, "musel_len_th_min")) & \
             (mu_len < _det_cut_th(detector, "musel_len_th_max"))

    p_cut = (mu_chi2_prot_cand > _det_cut_th(detector, "psel_muscore_th")) & \
            (prot_chi2_prot_cand < _det_cut_th(detector, "psel_pscore_th"))

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

def cathode_cut(df):
    if df.detector.iloc[0] == "SBND":
        p_start = df[['slc_vtx_x', 'slc_vtx_y', 'slc_vtx_z']].values
        p_mu = df[['mu_end_x', 'mu_end_y', 'mu_end_z']].values
        

        p_start = df[['slc_vtx_x', 'slc_vtx_y', 'slc_vtx_z']].values
        p_prot = df[['p_end_x', 'p_end_y', 'p_end_z']].values

        return ~intersects_prism_vectorized(p_start, p_prot, (-5., -200., 0.), (5., 200., 500.)) & ~intersects_prism_vectorized(p_start, p_mu, (-5., -200., 0.), (5., 200., 500.))
    else:
        return pd.Series(True, df.index)

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

def containment_cut(df):
    return mufv_cut(df) & pfv_cut(df)

def presel_cut(df):
    return slcfv_cut(df) & containment_cut(df) & cathode_cut(df)

def all_cuts(recodf, DETECTOR=None, det_run=None):
    if DETECTOR:
        print(f"manual detector: {DETECTOR}")
        recodf['detector'] = DETECTOR
    if det_run:
        print(f"manual run: {det_run}")
        recodf['Run'] = det_run

    ## presel cut
    presel_mask = presel_cut(recodf)

    ### cosmic cut
    cosmic_mask = cosmic_cut(recodf)

    ### flash cut
    flash_mask = flash_cut(recodf)

    ### Two prong cut
    two_prong_mask = twoprong_cut(recodf)

    ### PID cut
    pid_mask = pid_cut(recodf)

    return presel_mask & cosmic_mask & flash_mask & two_prong_mask & pid_mask

