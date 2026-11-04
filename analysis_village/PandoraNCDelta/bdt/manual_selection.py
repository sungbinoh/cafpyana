"""A loose cut-based filter to run before the BDT.

Its job is only to remove background whose "shower" is not a shower: cosmic
muons, protons and pions filed as shower-like PFPs. It is deliberately high
efficiency -- the BDT does the real work on what survives.

Operating point chosen at ~92% signal efficiency, where 1g1p has a clear knee
(background rejection 57% -> 84% for three points of efficiency, then flat).
Measured at topology + pre-selection + energy window:

                 sig eff   bkg rej   non-EM rej
  1g1p            92.2%     84.1%      82.0%
  1g0p            92.3%     52.7%      67.1%

shw_density_max is implemented here rather than in cuts.py, so the manual
selection stays self-contained and cuts.py is untouched.
"""

import numpy as np

MANUAL = {
    "1g1p": {"nu_score_min": 0.40, "shw_trackScore_max": 0.50,
             "shw_len_min": 14.0, "shw_density_max": 10.0},
    "1g0p": {"nu_score_min": 0.40, "shw_trackScore_max": 0.44,
             "shw_len_min": 8.0},
}


def manual_cuts(df, topo, cuts):
    """Mask for the manual cuts alone (no topology, pre-selection or window)."""
    c = MANUAL[topo]
    m = np.ones(len(df), dtype=bool)
    if c.get("nu_score_min") is not None:
        m &= (df.nu_score > c["nu_score_min"]).to_numpy()
    if c.get("shw_trackScore_max") is not None:
        m &= (df.shw_trackScore < c["shw_trackScore_max"]).to_numpy()
    if c.get("shw_len_min") is not None:
        m &= (df.shw_len > c["shw_len_min"]).to_numpy()
    if c.get("shw_density_max") is not None:
        m &= (cuts.sv(df, "shw_density") < c["shw_density_max"]).to_numpy()
    return m


def stage(df, topo, cuts):
    """Topology + pre-selection + energy window + the manual cuts."""
    pre = (cuts.PRESELECTION[topo](df) & cuts.energy_window_cut(df, topo)).to_numpy()
    return pre & manual_cuts(df, topo, cuts)
