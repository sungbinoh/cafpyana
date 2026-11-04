"""Best-plane shower energy for SBND, with the per-plane scale correction.

Kept out of cuts.py because it is reconstruction, not selection. cuts.py
imports it for the ENERGY_WINDOW cut; plot_composition.py uses it directly.
"""

import numpy as np
import pandas as pd

# 1 / median(E_plane / E_true_photon) on selected signal; mildly circular.
SHW_ENERGY_SCALE = (1.20, 1.20, 1.15)   # (I0 induction-1, I1 induction-2, I2 collection)


def shw_best_plane(df):
    """Index of the best plane for energy: the one with the most hits."""
    nh = np.vstack([np.nan_to_num(df.shw_nhits_I0.to_numpy(dtype=float), nan=-1.),
                    np.nan_to_num(df.shw_nhits_I1.to_numpy(dtype=float), nan=-1.),
                    np.nan_to_num(df.shw_nhits_I2.to_numpy(dtype=float), nan=-1.)])
    return np.argmax(nh, axis=0)


def shw_energy_corrected(df, scale=SHW_ENERGY_SCALE):
    """Best-plane shower energy [GeV] with the per-plane scale applied.

    Reproduces the CAF's shw_energy_best to 99.97% of showers before scaling,
    so it can be used on dataframes that predate this correction.
    """
    ep = np.vstack([df.shw_energy_I0.to_numpy(dtype=float),
                    df.shw_energy_I1.to_numpy(dtype=float),
                    df.shw_energy_I2.to_numpy(dtype=float)])
    bp = shw_best_plane(df)
    idx = np.arange(len(df))
    return pd.Series(np.asarray(scale)[bp] * ep[bp, idx], index=df.index)


def visible_energy(df):
    """Shower energy plus proton KE from range [GeV]; shower only with no track."""
    return shw_energy_corrected(df) + df.trk_ke_p.where(df.n_trk > 0, 0.).fillna(0.)
