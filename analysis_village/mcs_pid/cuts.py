# Standard library imports
import os
import sys

# Third-party imports
import pandas as pd
import numpy as np

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

def containment_cut(df):
    interaction_levels = ['__ntuple', 'entry', 'rec.slc..index']

    # 3. Track score cut: evaluated ONLY across real particles
    # (Dummy rows are assigned True so they don't break the group .all() condition)
    contained_track = prefix_fv_cut(df, "end") | (df['trackScore'] != -5.0)
    slc_contained_mask = contained_track.groupby(level=interaction_levels).transform('all')

    return slc_contained_mask

def twoprong_cut(df):
    interaction_levels = ['__ntuple', 'entry', 'rec.slc..index']

    # 1. Define real particle PFPs
    is_real_particle = df['trackScore'] != -5.0

    # 2. Count ONLY real particles in each slice
    real_pfp_count = is_real_particle.groupby(level=interaction_levels).transform('sum')
    twopfp_mask = (real_pfp_count == 2)

    # 3. Track score cut
    valid_track_score = (df['trackScore'] > 0.5) | (~is_real_particle)
    slc_track_mask = valid_track_score.groupby(level=interaction_levels).transform('all')

    # 4. Connection distance cut
    valid_dist = (df['dist_start'] < 5.0) | (~is_real_particle)
    slc_connected_mask = valid_dist.groupby(level=interaction_levels).transform('all')

    # 5. Track length cuts
    in_twopfp_slice = is_real_particle & twopfp_mask
    real_lengths = df['len'].where(in_twopfp_slice)   # NaN for dummies / other slices

    slc_max_len = real_lengths.groupby(level=interaction_levels).transform('max')
    slc_min_len = real_lengths.groupby(level=interaction_levels).transform('min')

    is_muon_candidate   = in_twopfp_slice & (real_lengths == slc_max_len)
    is_proton_candidate = in_twopfp_slice & (real_lengths == slc_min_len)

    muon_length_ok   = (~is_muon_candidate)   | (df['len'] >= 140.0)
    proton_length_ok = (~is_proton_candidate) | ((df['len'] >= 10.0) & (df['len'] <= 50.0))

    slc_muon_len_mask   = muon_length_ok.groupby(level=interaction_levels).transform('all')
    slc_proton_len_mask = proton_length_ok.groupby(level=interaction_levels).transform('all')

    slice_passes = (twopfp_mask & slc_track_mask & slc_connected_mask
                    & slc_muon_len_mask & slc_proton_len_mask)

    return slice_passes & is_real_particle

def nuscore_cut(df):
    return (df.nu_score > 0.5)

def baryscore_cut(df):
    return (df.baryscore > 1e-6)

def geom1u1p_cuts(df):

    twoprong_mask = twoprong_cut(df)

    nuscore_mask = nuscore_cut(df)

    baryscore_mask = baryscore_cut(df)

    fv_mask = slcfv_cut(df)

    containment_mask = containment_cut(df)

    return twoprong_mask & fv_mask & nuscore_mask & containment_mask
