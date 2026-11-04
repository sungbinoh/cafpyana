"""Dataframe makers for the Wire-Cell NC Delta -> N gamma (1g1p / 1g0p) selection.

One row per Wire-Cell neutrino candidate, index (entry, nu_index).
Self-contained: uses makedf (wirecell, make_mcdf) and this directory only.
Columns:
  * the T_tagger summary (vertex, flash, flags, numu/nue scores) and every
    shw_sp_* leaf (the single-photon tagger), raw from the CAF. shw_sp_* are
    only computed when shw_sp_filled == 1; otherwise they hold 0 / -1 defaults.
  * T_kine energies (kine_*, MeV as in the CAF)
  * topology from the T_kine particles counted in the energy (GeV): n_shw,
    n_trk, per-species counts, leading shower and proton
  * on MC, the matched true neutrino (tmatch_idx -> rec.mc.nu) with the
    true_* columns of make_ncdeltawc_truthdf, for cuts.truth_category and
    cuts.is_signal.

Only trees that data has too are used for reco columns: T_pf is MC-only.
"""

import numpy as np
import pandas as pd

from makedf import wirecell as wc
from makedf.makedf import make_mcdf
from pyanalib.pandas_helpers import loadbranches

from analysis_village.NCDeltaWC import cuts

MEV = 1e-3               # MeV -> GeV
PROTON_KE_THRESHOLD = 0.035  # GeV; the MicroBooNE Wire-Cell proton threshold
SHW_E_THRESHOLD = 0.020      # GeV; as shw_sp_n_20mev_showers
NU_LEVELS = wc.WC_NU_LEVELS

# kine_particle_type is the start segment's PID label. EM showers (photons
# included) are 11; a shower object can also carry 13/211/2212, and those
# count as tracks here, as Wire-Cell adds their rest mass like a track's.
SHW_PDG = 11
TRK_PDGS = {13: "muon", 211: "pion", 2212: "proton"}


# GeV. Photon thresholds for the truth counts: GENIE primaries include ~3 MeV
# nuclear de-excitation gammas, which these separate from a Delta -> N gamma
# photon (O(100) MeV).
GAMMA_THRESHOLDS = (0.020, 0.100)
GAMMA_COUNT_COLS = ["true_ngamma_%dMeV" % (1000. * t) for t in GAMMA_THRESHOLDS]
NU_TRUTH_LEVELS = [0, 1]   # (entry, rec.mc.nu index)


def fetch_metadata(f):
    """(DETECTOR, RUN, ismc) from rec.hdr."""
    hdr = lambda leaf: loadbranches(f["recTree"], ["rec.hdr." + leaf]).rec.hdr[leaf]
    dets = pd.unique(hdr("det"))
    if len(dets) != 1:
        raise ValueError("fetch_metadata: rec.hdr.det = %s; need exactly one detector." % dets)
    if dets[0] == 1:
        DETECTOR = "SBND"
    elif dets[0] == 2:
        DETECTOR = "ICARUS"
    else:
        raise ValueError("rec.hdr.det must be 1 (SBND) or 2 (ICARUS); got %s" % dets[0])
    RUN = 1 if DETECTOR == "SBND" else (2 if hdr("run").iloc[0] < 12960 else 4)
    ismc = bool(hdr("ismc").iloc[0])
    return DETECTOR, RUN, ismc


def make_ncdeltawc_truthdf(f):
    """One row per true neutrino (rec.mc.nu), index (entry, nu index): the
    GENIE truth columns used by cuts.py. T_truth_nu has the same order, so
    candidates join on (entry, tmatch_idx)."""
    DETECTOR, RUN, ismc = fetch_metadata(f)
    t = make_mcdf(f)
    N = pd.DataFrame({
        "true_E": t.E,                 # GeV
        "true_pdg": t.pdg,
        "true_iscc": t.iscc,
        "true_isnc": t.isnc,
        "true_genie_mode": t.genie_mode,
        "true_w": t.w,                 # GeV
        "true_vtx_x": t.position.x,
        "true_vtx_y": t.position.y,
        "true_vtx_z": t.position.z,
        "true_npi0": t.npi0,
        "true_ngamma": t.ng,
        "true_np": t.np,
        "true_np_20MeV": t.np_20MeV,
        "true_np_50MeV": t.np_50MeV,
        "true_npi_30MeV": t.npi_30MeV,
        "true_nmu_27MeV": t.nmu_27MeV,
    })

    prim = loadbranches(f["recTree"], ["rec.mc.nu.prim.pdg", "rec.mc.nu.prim.genE"]).rec.mc.nu.prim
    genE = prim.genE[prim.pdg == 22]
    N["true_leading_gamma_E"] = genE.groupby(level=NU_TRUTH_LEVELS).max().reindex(N.index)
    for thr, col in zip(GAMMA_THRESHOLDS, GAMMA_COUNT_COLS):
        N[col] = (genE > thr).groupby(level=NU_TRUTH_LEVELS).sum().reindex(N.index).fillna(0)
    N["true_npi0"] = N.true_npi0.fillna(0)
    # GENIE resonance number; 0 is the Delta(1232)
    N["true_resnum"] = loadbranches(f["recTree"], ["rec.mc.nu.resnum"]).rec.mc.nu.resnum.reindex(N.index)
    N["detector"] = DETECTOR
    N["Run"] = RUN
    return N


def _topology(f, index):
    """Counts and leading energies from the T_kine particles in the energy sum.

    Every T_kine row is summed into kine_reco_Enu, whatever this flag says.
    1 = connected to the candidate; 3 = a leftover shower whose connection
    is unclear (PRShower association type 3), mostly few-MeV fragments. Only
    the 1 rows are counted here. On MC every such track and ~3/4 of such
    showers are also T_pf particles.
    """
    k = wc.make_wc_kineparticledf(f)
    k = k[k.kine_energy_included == 1]
    E = k.kine_energy_particle * MEV
    pdg = k.kine_particle_type.abs()
    is_shw = pdg == SHW_PDG   # EM showers only, see SHW_PDG

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
        N = make_ncdeltawc_truthdf(f).drop(columns=["detector", "Run"])
        N.index.names = ["entry", "tmatch_idx"]
        D = D.join(N, on=["entry", "tmatch_idx"])
    return D


def make_ncdeltawc_df_mc(f):
    return make_ncdeltawc_df(f, include_truth=True)


def make_ncdeltawc_df_data(f):
    return make_ncdeltawc_df(f, include_truth=False)


def make_ncdeltawc_nudf(f):
    """make_ncdeltawc_truthdf + wc_valid: drop wc_valid == 0 events from
    efficiency denominators, they have no Wire-Cell output."""
    N = make_ncdeltawc_truthdf(f)
    valid = wc.make_wc_evtdf(f).valid
    N["wc_valid"] = valid.reindex(N.index.get_level_values(0)).values
    return N
