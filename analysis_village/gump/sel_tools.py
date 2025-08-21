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
import kinematics
from makedf.util import *


# Fiducial volume cuts for SBND and ICARUS
SBNDFVCuts = {
    "x": {"min": -199.15 + 10, "max": 199.15 - 10},
    "y": {"min": -200. + 10, "max": 200. - 10},
    "z": {"min": 0.0 + 10, "max": 500. - 50}
}

ICARUSFVCuts = {
    "C0": {
        "x": {"min": -358.49 + 10, "max": -61.94 - 10},
        "y": {"min": -181.86 + 10, "max": 134.96 - 10},
        "z": {"min": -894.950652270838 + 10, "max": 894.950652270838 - 50}
    },
    "C1": {
        "x": {"min": 358.49, "max": 61.94},
        "y": {"min": -181.86 + 10, "max": 134.96 - 10},
        "z": {"min": -894.950652270838 + 10, "max": 894.950652270838 - 50}
    }
}

def tmatchdf(df: pd.DataFrame, mcdf: pd.DataFrame, DETECTOR: str) -> pd.DataFrame:
    """
    Match reconstructed and MC truth DataFrames, add baseline info.
    Args:
        df: Reconstructed DataFrame.
        mcdf: MC truth DataFrame.
        DETECTOR: 'SBND' or 'ICARUS'.
    Returns:
        Merged DataFrame with baseline and truth columns.
    """
    # filter the columns from the mcdf we want
    # map old name to new name
    tosave = {
      "E": "trueE",
      "Q2": "trueQ2",
      "pdg": "truepdg",
      "iscc": "iscc",
      "genie_mode": "genie_mode",
      "np_50MeV": "np_50MeV",
      "np_20MeV": "np_20MeV",
    }

    def savecol(c):
        return c[0] in list(tosave.keys()) or is_weightcol(c[0])

    mcdf_cols = [c for c in mcdf.columns if savecol(c)]
    mcdf = mcdf[mcdf_cols]

    # Get the column depth matching on both sides
    mcdf.columns = pd.MultiIndex.from_tuples([(tosave[c[0]], c[1]) if c[0] in tosave else (c[0], c[1]) for c in mcdf.columns])
    df.columns = pd.MultiIndex.from_tuples([(c, "") for c in df.columns])

    tmatch_df = pd.merge(df, mcdf, how="left", left_on=["__ntuple", "entry", "tmatch"], right_index=True)

    # Fill nan's for not-matching columns
    for c in tosave.values():
        tmatch_df[(c, "")] = tmatch_df[(c, "")].fillna(0.)

    # Systematic weights are all 1 for no truth match
    tmatch_df.fillna(1, inplace=True)

    # TODO: get actual baseline
    if DETECTOR == "SBND":
        tmatch_df[("baseline", "")] = 110.
        print("Using SBND Baseline of 110")
    elif DETECTOR == "ICARUS":
        tmatch_df[("baseline", "")] = 600.
        print("Using ICARUS Baseline of 600")
    else:
        print("No detector configured!")
        sys.exit()

    return tmatch_df

def is_weightcol(c: str) -> bool:
    """Return True if column is a weight column."""
    return c.startswith("GENIEReWeight") or c.endswith("Flux")

# Reconstructed variables we want to save
def recodf(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extract reconstructed variables to save.
    Args:
        df: DataFrame with reconstructed event info.
    Returns:
        DataFrame with selected variables.
    """
    dPt = df.del_p.rename("dPt")
    caloE = kinematics.neutrino_energy(df.mu.pfp.trk.P.p_muon, df.mu.pfp.trk.dir, df.p.pfp.trk.P.p_proton, df.p.pfp.trk.dir).rename("caloE")
    muonE = kinematics.muon_energy(df.mu.pfp.trk.P.p_muon, df.mu.pfp.trk.dir).rename("muonE")
    # Lookup truth match
    # TODO: make this better, include cut on purity
    tmatch = df['rec.mc.nu..index'].fillna(-1).astype(int).rename("tmatch")

    return pd.concat([dPt, caloE, muonE, tmatch], axis=1) # .droplevel(0, 0)

def SelFV(df: pd.DataFrame, det: str, inzback: int = 50) -> pd.Series:
    """
    Fiducial volume cut for SBND or ICARUS.
    Args:
        df: DataFrame with x, y, z columns.
        det: Detector name ('SBND' or 'ICARUS').
        inzback: Z-back cut (unused).
    Returns:
        Boolean mask for events in FV.
    """
    if det == "ICARUS":
        return (((df.x < ICARUSFVCuts['C0']['x']['max']) & (df.x > ICARUSFVCuts['C0']['x']['min'])) |\
                ((df.x < ICARUSFVCuts['C1']['x']['max']) & (df.x > ICARUSFVCuts['C1']['x']['min']))) &\
                 (df.y < ICARUSFVCuts['C0']['y']['max']) & (df.y > ICARUSFVCuts['C0']['y']['min']) &\
                 (df.z < ICARUSFVCuts['C0']['z']['max']) & (df.z > ICARUSFVCuts['C0']['z']['min'])

    elif det == "SBND":
        return (df.x < SBNDFVCuts['x']['max']) & (df.x > SBNDFVCuts['x']['min']) &\
               (df.y < SBNDFVCuts['y']['max']) & (df.y > SBNDFVCuts['y']['min']) &\
               (df.z < SBNDFVCuts['z']['max']) & (df.z > SBNDFVCuts['z']['min'])

    else:
        raise NameError("DETECTOR not valid, should be SBND or ICARUS")

def InBeam(t):
    return (t > 0.) & (t < 1.800)

def is_cosmic(df: pd.DataFrame) -> pd.Series:
    """Return mask for cosmic events."""
    return (df.slc.truth.pdg == -1)

def is_FV(df: pd.DataFrame, det: str) -> pd.Series:
    """Return mask for events in fiducial volume."""
    return SelFV(df.slc.vertex, det)

def is_numu(df: pd.DataFrame) -> pd.Series:
    """Return mask for numu events."""
    return (np.abs(df.slc.truth.pdg) == 14)

def is_CC(df: pd.DataFrame) -> pd.Series:
    """Return mask for CC events."""
    return (df.slc.truth.iscc == 1)

def is_NC(df: pd.DataFrame) -> pd.Series:
    """Return mask for NC events."""
    return (df.slc.truth.iscc == 0)

def is_1p0pi(df: pd.DataFrame) -> pd.Series:
    """Return mask for 1mu, 1p, 0pi events."""
    return (df.slc.truth.nmu_27MeV == 1) & (df.slc.truth.np_50MeV == 1) & (df.slc.truth.npi_30MeV == 0) & (df.slc.truth.npi0 == 0)

def is_signal(df: pd.DataFrame, det: str) -> pd.Series:
    """Return mask for signal events."""
    return is_numu(df) & is_CC(df) & is_1p0pi(df) & is_FV(df, det)

def is_outFV(df: pd.DataFrame, det: str) -> pd.Series:
    """Return mask for signal events outside FV."""
    return is_numu(df) & is_CC(df) & is_1p0pi(df) & np.invert(is_FV(df, det))

def is_othernumuCC(df: pd.DataFrame, det: str) -> pd.Series:
    """Return mask for other numu CC events in FV."""
    return is_numu(df) & is_CC(df) & np.invert(is_1p0pi(df)) & is_FV(df, det)

mode_list = [0, 10, 1, 2, 3]
mode_labels = ['QE', 'MEC', 'RES', 'SIS/DIS', 'COH', "other"]

def breakdown_mode(var: pd.Series, df: pd.DataFrame) -> list:
    """Break down variable by interaction mode."""
    ret = [var[df.genie_mode == i] for i in mode_list]
    ret.append(var[sum([df.genie_mode == i for i in mode_list]) == 0])
    return ret

def breakdown_top(var: pd.Series, df: pd.DataFrame, det: str) -> list:
    """Break down variable by topological category."""
    ret = [var[is_signal(df, det)],
            var[is_othernumuCC(df, det)],
            var[is_NC(df)],
            var[is_outFV(df, det)],
            var[is_cosmic(df)],
            var[np.invert(is_signal(df, det) | is_othernumuCC(df, det) | is_NC(df) | is_outFV(df, det) | is_cosmic(df))]
            ]
    return ret

def multicol_merge(lhs: pd.DataFrame, rhs: pd.DataFrame, **panda_kwargs) -> pd.DataFrame:
    """
    Merge two DataFrames with possibly different column levels.
    """
    lhs_col = lhs.columns
    rhs_col = rhs.columns

    nlevel = max(lhs_col.nlevels, rhs_col.nlevels)

    def pad(c):
        nc = 1 if isinstance(c, str) else len(c)
        c0 = [c] if isinstance(c, str) else list(c)
        return tuple(c0 + [""]*(nlevel - nc))

    lhs.columns = pd.MultiIndex.from_tuples([pad(c) for c in lhs_col])
    rhs.columns = pd.MultiIndex.from_tuples([pad(c) for c in rhs_col])

    return lhs.merge(rhs, **panda_kwargs)

def add_weights(df: pd.DataFrame, wgtdf: pd.DataFrame, tmatch_col=("slc", "tmatch", "idx")) -> pd.DataFrame:
    """
    Add weights from wgtdf to df, filling non-matched events with 1.
    """
    na_set = {}
    nlevel = max(df.columns.nlevels, wgtdf.columns.nlevels)
    for w in wgtdf.columns:
        nc = len(w)
        na_set[tuple(list(w) + [""]*(nlevel-nc))] = 1 # don't weight non-truth matched events
    return multicol_merge(df, wgtdf, how="left", left_on=wgtdf.index.names[:2] + [tmatch_col], right_index=True).fillna(na_set)
