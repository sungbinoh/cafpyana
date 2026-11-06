"""Signal definition and truth categories for the inclusive single-photon
search with Wire-Cell (SBND).

Signal, following the MicroBooNE inclusive photon-like LEE search
(PRL, doi:10.1103/89qs-4lcp), which is modelled on MiniBooNE's single
Cherenkov ring:

  * any true neutrino vertex position: a single photon from a vertex
    outside the active volume is signal, as in MicroBooNE (category
    "1g: vertex outside AV")
  * exactly one reconstructable photon shower. A reconstructable photon has
    true E >= 20 MeV and starts its shower (pair conversion) >= 3 cm inside
    the TPC boundary. Photons from the same vertex closer than 20 deg are
    one shower.
  * no primary muon with KE >= 100 MeV (below MiniBooNE's threshold)
  * no primary electron/positron (photon-like only)
  * any number of protons, charged pions, neutrons; any number of photons
    that are not reconstructable (soft, escaping, converting near a wall)

So NC pi0 events where one photon is lost or the two merge are signal, as
are NC Delta -> N gamma and CC pi0 with a soft muon. Photons are the prompt
ones from the interaction: GENIE primaries and the decay photons of primary
hadrons (pi0, eta, ...); photons from secondary interactions in the argon
do not count.

Truth columns are those of make_truth.make_sp_truthdf, joined to candidates
on (entry, tmatch_idx); a candidate with no match is "cosmic".

Reco selections (one row per Wire-Cell candidate, make_wc_candidates.py):
PRESELECTION and the MicroBooNE reference UB_SELECTION, cut chains applied
in order by passes().
"""

import numpy as np
import pandas as pd

# --- Geometry (SBND only)

# Active volume: the TPC box, no insets [cm]
ACTIVE_VOLUME = {"x": (-200.0, 200.0), "y": (-200.0, 200.0), "z": (0.0, 500.0)}

# --- Signal definition

PHOTON_E_MIN = 0.020       # GeV, true photon energy
SHOWER_INSET = 3.0         # cm, shower start from the TPC boundary (as MicroBooNE)
MERGE_ANGLE = 20.0         # deg, photons closer than this are one shower
MUON_KE_MAX = 0.100        # GeV, primary muons below this are allowed


def dist_to_active_boundary(x, y, z):
    """Signed distance [cm] to the active-volume boundary, positive inside.
    NaN coordinates give NaN."""
    x, y, z = (np.asarray(v, dtype=float) for v in (x, y, z))
    b = ACTIVE_VOLUME
    d = np.minimum.reduce([
        x - b["x"][0], b["x"][1] - x,
        y - b["y"][0], b["y"][1] - y,
        z - b["z"][0], b["z"][1] - z,
    ])
    d[~(np.isfinite(x) & np.isfinite(y) & np.isfinite(z))] = np.nan
    return d


def true_in_active(df):
    return pd.Series(
        dist_to_active_boundary(df.true_vtx_x, df.true_vtx_y, df.true_vtx_z) > 0,
        index=df.index)


def is_cosmic(df):
    """No true neutrino behind the row (never true in the truth table)."""
    m = df.true_E.isna()
    if "tmatch_idx" in df.columns:
        m |= df.tmatch_idx < 0
    return m


def is_signal(df):
    """MicroBooNE-style: one reconstructable photon shower in the TPC, any
    neutrino vertex position (out-of-AV 1g is signal, as in MicroBooNE)."""
    return (df.true_nshower == 1) & (df.true_nmu_100MeV == 0) & (df.true_ne == 0)


# --- Truth categories

RESNUM_DELTA_1232 = 0   # genie::EResonance kP33_1232
ORIGIN_PRIMARY = 22     # true_shw1_origin of a GENIE-primary photon

CATEGORY_LABELS = {
    0: "1g: NC Delta->Ngamma",
    1: "1g: NC pi0",
    2: "1g: CC pi0, mu < 100 MeV",
    3: "1g: other",
    10: "1g: vertex outside AV",
    4: ">=2 showers",
    5: "0 showers",
    6: "CC numu, mu >= 100 MeV",
    7: "CC nue (primary e)",
    8: "Outside AV, not 1g",
    9: "Cosmic / unmatched",
}
SIGNAL_CATEGORIES = [0, 1, 2, 3, 10]


def is_signal_cat(cat):
    return cat.isin(SIGNAL_CATEGORIES)


def truth_category(df):
    """Integer category per row (CATEGORY_LABELS); first match wins."""
    o = df.true_shw1_origin
    sig = is_signal(df)
    ncdelta = (df.true_isnc == 1) & (o == ORIGIN_PRIMARY) & (df.true_resnum == RESNUM_DELTA_1232)
    conds = [
        (is_cosmic(df), 9),
        (sig & ~true_in_active(df), 10),
        (~true_in_active(df), 8),
        (df.true_ne > 0, 7),
        (df.true_nmu_100MeV > 0, 6),
        (sig & ncdelta, 0),
        (sig & (df.true_isnc == 1) & (o == 111), 1),
        (sig & (df.true_iscc == 1) & (o == 111), 2),
        (sig, 3),
        (df.true_nshower >= 2, 4),
        (df.true_nshower == 0, 5),
    ]
    cat = pd.Series(-1, index=df.index)
    for m, c in reversed(conds):
        cat[m.fillna(False).astype(bool)] = c
    return cat


# --- Reco preselection (Wire-Cell candidate columns)
#
#   wc_valid    Wire-Cell output exists for the event (wcevt valid = 1; 0 is
#               a production failure, not a reco outcome). An event cut:
#               valid = 0 events have no candidates, so cutflow applies it
#               to the true neutrinos (EVENT_CUTS)
#   vertex      has_vertex == 1: a vertex was found, so the taggers ran
#   generic     MicroBooNE Wire-Cell generic neutrino selection: main
#               cluster not a through-going muon, stopping muon or light
#               mismatch, and > 15 cm. Every candidate already has a matched
#               flash; tgm/stm/lm fire on < 0.1% of SBND candidates (mostly
#               applied upstream), so in practice this is the length cut.
#               The cosmic tagger (cosmict_*) is not cut on: as in MicroBooNE
#               it feeds the "other" BDT.
#   fv          reco vertex >= 5 cm inside the active volume
#   em_shower   >= 1 connected EM shower (T_kine PID 11)

FV_INSET = 5.0                 # cm
MIN_CLUSTER_LENGTH = 15.0      # cm


def reco_in_fv(df, inset=FV_INSET):
    return pd.Series(dist_to_active_boundary(df.nu_x, df.nu_y, df.nu_z) >= inset, index=df.index)


EVENT_CUTS = {"wc_valid"}

PRESELECTION = {
    "wc_valid": lambda df: df.wc_valid == 1,
    "vertex": lambda df: df.has_vertex == 1,
    "generic": lambda df: ((df.main_cluster_tgm == 0) & (df.main_cluster_stm == 0) & (df.main_cluster_lm == 0)
                           & (df.main_cluster_length_cm > MIN_CLUSTER_LENGTH)),
    "fv": reco_in_fv,
    "em_shower": lambda df: df.n_em >= 1,
}


# --- MicroBooNE reference selection (microboone_single_photon_selection.yaml)
#
# On the MicroBooNE numu, other and NC pi0 BDT scores, columns ub_<name>
# (microboone_bdt.score; NaN = not evaluable), plus the SBND Wire-Cell
# tagger scores: numu_score (removes CC numu the MicroBooNE numu BDT keeps)
# and nue_score instead of the MicroBooNE nue BDT (~2x the signal at the
# same purity). The MicroBooNE FV is not used.
#   vertex_shower   >= 1 photon-like shower > 20 MeV at the vertex
#   bdts_evaluable  shw_sp_filled and every input finite, the three BDTs used
#   wc_numu         not numu CC-like: numu_score < 0.5
#   ub_numu         MicroBooNE numu BDT > 1 (MicroBooNE: 0.4; tighter on SBND,
#                   tuned out-of-fold, as wc_numu already cut)
#   ub_other/ncpi0  score above the MicroBooNE threshold (log10 score scale,
#                   see microboone_bdt). Re-tuning on SBND MC (3 cm inset)
#                   found these already near-optimal.
#   wc_nue          not nue-like: nue_score < 7 (MicroBooNE's nue selection
#                   is nue_score > 7)
#   one_shower      exactly one reconstructable photon-like shower

# Keep the order: the MicroBooNE BDTs are a cascade, each trained on events
# passing the previous ones (the NC pi0 BDT has no separation before "other").
UB_SCORES = ["ub_numu", "ub_other", "ub_ncpi0"]
UB_NUMU_MIN = 1.0
UB_OTHER_MIN = 0.2
UB_NCPI0_MIN = -0.05
WC_NUMU_SCORE_MAX = 0.5
WC_NUE_SCORE_MAX = 7.0

UB_SELECTION = {
    "vertex_shower": lambda df: df.shw_sp_n_20mev_showers > 0,
    "bdts_evaluable": lambda df: df[UB_SCORES].notna().all(axis=1),
    "wc_numu": lambda df: df.numu_score < WC_NUMU_SCORE_MAX,
    "ub_numu": lambda df: df.ub_numu > UB_NUMU_MIN,
    "ub_other": lambda df: df.ub_other > UB_OTHER_MIN,
    "ub_ncpi0": lambda df: df.ub_ncpi0 > UB_NCPI0_MIN,
    "wc_nue": lambda df: df.nue_score < WC_NUE_SCORE_MAX,
    "one_shower": lambda df: df.shw_sp_n_20br1_showers == 1,
}

# Cut chains for cutflow.py
CHAINS = {
    "presel": PRESELECTION,
    "uboone": UB_SELECTION,
    "presel+uboone": {**PRESELECTION, **UB_SELECTION},
}


def passes(df, chain=PRESELECTION, upto=None):
    """Candidates passing the cuts of `chain`, through `upto` (inclusive)."""
    m = pd.Series(True, index=df.index)
    for name, cut in chain.items():
        m &= cut(df)
        if name == upto:
            break
    return m
