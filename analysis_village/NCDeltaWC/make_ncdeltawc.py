"""Dataframe makers for the Wire-Cell NC Delta -> N gamma (1g1p / 1g0p) selection.

The Wire-Cell counterpart of analysis_village/NCDelta: one row per Wire-Cell
neutrino candidate, index (entry, nu_index), instead of one per Pandora slice.
Columns:
  * the T_tagger summary (vertex, flash, flags, numu/nue scores) and every
    shw_sp_* leaf (the single-photon tagger), raw from the CAF. shw_sp_* are
    only computed when shw_sp_filled == 1; otherwise they hold 0 / -1 defaults.
  * T_kine energies (kine_*, MeV as in the CAF)
  * topology from the T_kine particles counted in the energy (GeV, like
    NCDelta): n_shw, n_trk, per-species counts, leading shower and proton
  * on MC, the matched true neutrino (tmatch_idx -> rec.mc.nu) with the same
    true_* columns as make_ncdelta_nudf, so NCDelta's cuts.truth_category and
    cuts.is_signal apply unchanged.

Only trees that data has too are used for reco columns: T_pf is MC-only.
"""

import numpy as np
import pandas as pd

from makedf import wirecell as wc

from analysis_village.NCDelta import cuts
from analysis_village.NCDelta.make_ncdelta import fetch_metadata, make_ncdelta_nudf

MEV = 1e-3               # MeV -> GeV
PROTON_KE_THRESHOLD = 0.035  # GeV; the MicroBooNE Wire-Cell proton threshold
SHW_E_THRESHOLD = 0.020      # GeV; as shw_sp_n_20mev_showers
NU_LEVELS = wc.WC_NU_LEVELS

SHW_PDG = 11   # Wire-Cell labels every shower 11, photons included
TRK_PDGS = {13: "muon", 211: "pion", 2212: "proton"}


def _topology(f, index):
    """Counts and leading energies from the T_kine particles in the energy sum.

    kine_energy_included == 1 are the particles summed into kine_reco_Enu.
    On MC every such track and ~3/4 of such showers are also T_pf particles;
    the rest (== 3) are mostly few-MeV fragments.
    """
    k = wc.make_wc_kineparticledf(f)
    k = k[k.kine_energy_included == 1]
    E = k.kine_energy_particle * MEV
    pdg = k.kine_particle_type.abs()
    is_shw = pdg == SHW_PDG

    g = lambda m: m.groupby(level=NU_LEVELS).sum()
    T = pd.DataFrame(index=index)
    T["n_shw"] = g(is_shw)
    T["n_shw_20MeV"] = g(is_shw & (E > SHW_E_THRESHOLD))
    T["n_trk"] = g(pdg.isin(TRK_PDGS.keys()))
    for p, name in TRK_PDGS.items():
        T["n_" + name] = g(pdg == p)
    T["n_proton_35MeV"] = g((pdg == 2212) & (E > PROTON_KE_THRESHOLD))
    T = T.fillna(0).astype(int)

    T["shw_energy"] = E[is_shw].groupby(level=NU_LEVELS).max()
    T["shw_energy_2nd"] = E[is_shw].groupby(level=NU_LEVELS).apply(
        lambda s: s.nlargest(2).iloc[-1] if len(s) > 1 else np.nan)
    T["trk_ke_p"] = E[pdg == 2212].groupby(level=NU_LEVELS).max()
    return T


def make_ncdeltawc_df(f, include_truth=None):
    DETECTOR, RUN, ismc = fetch_metadata(f)
    if include_truth is None:
        include_truth = ismc

    leaves = wc.wctaggerbranches + wc.wc_tagger_leaves(f, ["shw_sp_"])
    D = wc.make_wc_nudf(f, tagger_leaves=leaves, include_tmatch=include_truth)
    D["detector"] = DETECTOR
    D["Run"] = RUN
    D["ismc"] = ismc
    D["reco_Enu"] = D.kine_reco_Enu * MEV
    D["vtx_dist_to_wall"] = cuts.dist_to_active_boundary(D.nu_x, D.nu_y, D.nu_z, DETECTOR)

    D = D.join(_topology(f, D.index))

    if include_truth:
        N = make_ncdelta_nudf(f).drop(columns=["detector", "Run"])
        N.index.names = ["entry", "tmatch_idx"]
        D = D.join(N, on=["entry", "tmatch_idx"])
    return D


def make_ncdeltawc_df_mc(f):
    return make_ncdeltawc_df(f, include_truth=True)


def make_ncdeltawc_df_data(f):
    return make_ncdeltawc_df(f, include_truth=False)


def make_ncdeltawc_nudf(f):
    """make_ncdelta_nudf + wc_valid: drop wc_valid == 0 events from efficiency
    denominators, Wire-Cell never ran on them."""
    N = make_ncdelta_nudf(f)
    valid = wc.make_wc_evtdf(f).valid
    N["wc_valid"] = valid.reindex(N.index.get_level_values(0)).values
    return N
