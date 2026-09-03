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

# Final optimized selection
# opt to ratio, min mu len 40, shared track score, per-detector PID
# opt used the most updated calo treatment (no EMB IC, no sqmear15, double SBND 13)

SBND_CUTS = {
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
    "max_opening_angle": 160,
    "musel_track_score_min": 0.5,
    # The chi2-under-muon-hypothesis (chi2u) cut on the ICARUS muon candidate is
    # deliberately disabled: the ICARUS muon candidate is PID'd by the
    # proton-hypothesis chi2 (musel_pscore_th) alone. The key is kept -- rather
    # than deleted -- because _det_cut_th() indexes SBND_CUTS and ICARUS_CUTS
    # with the same key. NB inf, not "no cut": a NaN chi2u still fails.
    "musel_muscore_th": np.inf,
    "musel_pscore_th": 74,
    "musel_len_th_min": 40,
    "musel_len_th_max": 400,
    "psel_muscore_th": 0,
    "psel_pscore_th": 92,
}

def maple_v_calovariation(df, variation):
    """Copy of df with the nominal chi2 candidate columns replaced by the
    {variation}-suffixed ones (GUMP nb v_calovariation semantics)."""
    df = df.copy()
    for c in CHI2_CAND_COLS:
        df[c % ""] = df[c % variation]
    return df

# function for picking cut val based on detector
def _det_cut_th(detector, key):
    return np.where(detector == "SBND", SBND_CUTS[key], ICARUS_CUTS[key])

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

# Dangling cable in TPC WW (West cryostat, West TPC -- the high-x TPC of C1,
# bounded below by the west cathode at 210.22). The high-y, downstream-z corner
# is unusable and is removed from the FV in both Run 2 and Run 4.
ICARUSCableCuts = {
    "WW": {
        "x": {"min": 210.22, "max": 358.49},
        "y": {"min": 70.},
        "z": {"min": 0.}
    }
}

def prefix_fv_cut(df, prefix, is_trk=False):
    """Generic helper to extract coordinates via prefix and apply FV cuts."""
    vtx = pd.DataFrame({
        'detector': df.detector,
        'Run': df.Run,
        'x': df[f'{prefix}_x'],
        'y': df[f'{prefix}_y'],
        'z': df[f'{prefix}_z']
    }, index=df.index)
    inz = 10 if is_trk else 50
    return _fv_cut(vtx, inzback=inz)

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

    # Dangling-cable region in TPC WW -- removed in both runs. The boundaries are
    # absolute, so the inx/iny/inz insets are deliberately not applied here.
    inCable = (df.x > ICARUSCableCuts['WW']['x']['min']) & (df.x < ICARUSCableCuts['WW']['x']['max']) &\
              (df.y > ICARUSCableCuts['WW']['y']['min']) &\
              (df.z > ICARUSCableCuts['WW']['z']['min'])

    FVRun2 = (((df.x < (ICARUSRun2FVCuts['C0']['x']['max'] - inx)) & (df.x > (ICARUSRun2FVCuts['C0']['x']['min'] + inx))) |\
            ((df.x < (ICARUSRun2FVCuts['C1']['x']['max'] - inx)) & (df.x > (ICARUSRun2FVCuts['C1']['x']['min'] + inx)))) &\
             (df.y < (ICARUSRun2FVCuts['C0']['y']['max'] - iny)) & (df.y > (ICARUSRun2FVCuts['C0']['y']['min'] + iny)) &\
             (df.z < (ICARUSRun2FVCuts['C0']['z']['max'] - inzback)) & (df.z > (ICARUSRun2FVCuts['C0']['z']['min'] + inzfront)) &\
             ~inCable

    FVRun4 = (((df.x < (ICARUSRun4FVCuts['C0']['x']['max'] - inx)) & (df.x > (ICARUSRun4FVCuts['C0']['x']['min'] + inx))) |\
            ((df.x < (ICARUSRun4FVCuts['C1']['x']['max'] - inx)) & (df.x > (ICARUSRun4FVCuts['C1']['x']['min'] + inx)))) &\
             (df.y < (ICARUSRun4FVCuts['C0']['y']['max'] - iny)) & (df.y > (ICARUSRun4FVCuts['C0']['y']['min'] + iny)) &\
             (df.z < (ICARUSRun4FVCuts['C0']['z']['max'] - inzback)) & (df.z > (ICARUSRun4FVCuts['C0']['z']['min'] + inzfront)) &\
             ~inCable

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

def cosmic_cut(df):
    return (df["mu_p_opening_angle_deg"] < _det_cut_th(df.detector, "max_opening_angle")) & (df["crthit"] == 0)

def slcfv_cut(df):
    return prefix_fv_cut(df, "slc_vtx")

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

def sbnd_cathode_crossing(vtx_x, vtx_y, vtx_z, end_x, end_y, end_z):
    """Per-track SBND cathode-crossing flag (GUMP cathode_cut semantics).

    True where the segment (slice vertex -> track end) intersects the cathode
    prism x in (-5, 5). NaN coordinates yield False.
    """

    SBND_CATHODE_PRISM_MIN = (-5.0, -200.0, 0.0)
    SBND_CATHODE_PRISM_MAX = (5.0, 200.0, 500.0)

    p1 = np.stack([np.asarray(vtx_x, dtype=float),
                   np.asarray(vtx_y, dtype=float),
                   np.asarray(vtx_z, dtype=float)], axis=1)
    p2 = np.stack([np.asarray(end_x, dtype=float),
                   np.asarray(end_y, dtype=float),
                   np.asarray(end_z, dtype=float)], axis=1)

    return intersects_prism_vectorized(p1, p2, SBND_CATHODE_PRISM_MIN, SBND_CATHODE_PRISM_MAX)

def sanity_cut(df):
    return df.slc_vtx_x.notna() & df.slc_vtx_y.notna() & df.slc_vtx_z.notna()

# Requires at least two pfp candidates (a muon and a candidate proton), as the
# legacy GUMP preselection did implicitly through the mu/p endpoint FV cuts.
def presel_cut(df):
    return sanity_cut(df) & slcfv_cut(df) & df.cut_contained & df.cut_cathode & \
        df.has_muon & df.cut_np

def trk_cut(df):
    # The shower/other veto (cut_0shwother, from the MAPLE id_pfp
    # classification with VTX_MAX_DIST=50) is NOT applied: multiplicity is
    # enforced by the candidate count alone (n_pfp == 2 for gump_sel,
    # n_pfp > 2 for maple_sel). n_shower/n_other/cut_0shwother remain in the
    # evt df as diagnostics.
    return df.cut_np & df.has_muon

def get_base_muon_mask(df, level="slc"):
    if level == "trk":
        pref = ""
    elif level == "slc":
        pref = "mu_"

    base_mask = (
            df[pref+"end_x"].notna()
            & df[pref+"len"].notna()
            & (df[pref+"trackScore"] >= _det_cut_th(df.detector, "musel_track_score_min"))
            & (df[pref+"dist_start"] <= 10.0)
            & df[pref+"len"].between(_det_cut_th(df.detector, "musel_len_th_min"), _det_cut_th(df.detector, "musel_len_th_max"))
            & df[pref+"prim_pfp"]
            & df[pref+"contained10"]
            & (df[pref+"end_x"] * df["slc_vtx_x"] > 0)
        )
    return base_mask

def get_muon_mask(df, variation=None):
    """Generates a boolean mask filtering candidate muon tracks.

    Applies the variation suffix 'v' strictly to chi2 variables.
    """

    # Format variation suffix (e.g., None -> "", "_sys" -> "_sys")
    v = "" if variation is None else variation

    # Standard kinematic/topological selection
    base_mask = get_base_muon_mask(df, level="slc")

    # Apply muon/proton chi2 cuts using the variation suffix on chi2 columns.
    # Column names match the makedf outputs (mu_chi2{v}_of_mu_cand /
    # prot_chi2{v}_of_mu_cand); a missing column raises rather than silently
    # skipping the cut (which is how the muon chi2 PID was lost in the
    # sbn-rewgted-16/17/18 productions).
    chi2_muon_col = f"mu_chi2{v}_of_mu_cand"
    chi2_prot_col = f"prot_chi2{v}_of_mu_cand"

    chi2_mask = (df[chi2_muon_col] < _det_cut_th(df.detector, "musel_muscore_th")) & \
                (df[chi2_prot_col] > _det_cut_th(df.detector, "musel_pscore_th"))

    return base_mask & chi2_mask

def get_proton_mask(df, variation=None):
    """Generates a boolean mask filtering candidate proton tracks on chi2 score."""

    # Format variation suffix (e.g., None -> "", "_sys" -> "_sys")
    v = "" if variation is None else variation

    # Dynamic column name targeting only prot_chi2
    chi2_prot_col = f"prot_chi2{v}_of_prot_cand"

    return df[chi2_prot_col] < _det_cut_th(df.detector, "psel_pscore_th")

def pid_cut(df, variation=None):
    cut_muon = get_muon_mask(df, variation=variation)
    cut_protons = get_proton_mask(df, variation=variation)
    return cut_muon & cut_protons

def gumple_base_cuts(recodf, DETECTOR=None, det_run=None, variation=None, do_cosmic_cut=True):
    """The shared (multiplicity-inclusive) cut chain: presel & cosmic & flash
    & trk & PID.  The GUMP/MAPLE selections split this on n_pfp (the number
    of candidate pfps in the slice, muon included)."""
    if DETECTOR:
        print(f"manual detector: {DETECTOR}")
        recodf['detector'] = DETECTOR
    if det_run:
        print(f"manual run: {det_run}")
        recodf['Run'] = det_run

    ## presel cut
    presel_mask = presel_cut(recodf)

    ### cosmic cut
    if do_cosmic_cut:
        cosmic_mask = cosmic_cut(recodf)
    else:
        cosmic_mask = pd.Series(True, index=recodf.index)

    ### flash cut
    flash_mask = flash_cut(recodf)

    ### trk cut
    trk_mask = trk_cut(recodf)

    ### PID cut
    pid_mask = pid_cut(recodf, variation=variation)

    return presel_mask & cosmic_mask & flash_mask & trk_mask & pid_mask

def all_gump_cuts(recodf, DETECTOR=None, det_run=None, variation=None):
    """GUMP (1u1p): base chain + exactly two candidate pfps (muon + proton)."""
    return gumple_base_cuts(recodf, DETECTOR=DETECTOR, det_run=det_run, variation=variation, do_cosmic_cut=True) & \
        (recodf.n_pfp == 2)

def all_maple_cuts(recodf, DETECTOR=None, det_run=None, variation=None):
    """MAPLE (1uNp): base chain + more than two candidate pfps.

    NB semantic change: this is the multiplicity-INCLUSIVE chain.
    """
    cut_far_shw = np.isnan(recodf.max_far_shw_len)
    return gumple_base_cuts(recodf, DETECTOR=DETECTOR, det_run=det_run, variation=variation, do_cosmic_cut=True) & \
        cut_far_shw

def all_maplemp_cuts(recodf, DETECTOR=None, det_run=None, variation=None):
    """MAPLE (1uN>1p): base chain + more than two candidate pfps.

    NB semantic change: this used to be the multiplicity-INCLUSIVE chain
    (now gumple_base_cuts); it is now exclusive of the GUMP (n_pfp == 2)
    selection."""

    cut_far_shw = np.isnan(recodf.max_far_shw_len)
    return maple_base_cuts(recodf, DETECTOR=DETECTOR, det_run=det_run, variation=variation, do_cosmic_cut=False) & \
        cut_far_shw & (recodf.n_pfp > 2)

def maple_cut_chain(recodf, DETECTOR=None, det_run=None, variation=None):
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

    ### trk cut
    trk_mask = trk_cut(recodf)

    ### PID cut
    cut_muon = get_muon_mask(recodf, variation=variation)
    cut_protons = get_proton_mask(recodf, variation=variation)

    ### far-shower veto (MAPLE side only): reject slices with a displaced primary
    ### shower (trackScore<0.5, 10-50 cm from the vertex). Passes iff none exists.
    cut_far_shw = np.isnan(recodf.max_far_shw_len)

    ### base chain, split on candidate-pfp multiplicity into the exclusive
    ### gump (1u1p, n_pfp==2) and maple (1uN>1p, n_pfp>2) selections
    base_sel = presel_mask & cosmic_mask & flash_mask & trk_mask & cut_muon & cut_protons

    return pd.DataFrame({
        "cut_presel": presel_mask,
        "cut_cosmic": cosmic_mask,
        "cut_flash": flash_mask,
        "cut_trk": trk_mask,
        "cut_muon": cut_muon,
        "cut_protons": cut_protons,
        "cut_far_shw": cut_far_shw,
        "gump_sel": base_sel & (recodf.n_pfp == 2),
        "maple_sel": base_sel & cut_far_shw & (recodf.n_pfp > 2),
    }, index=recodf.index)
