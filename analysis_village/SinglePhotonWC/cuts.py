"""Signal definition, truth categories and reco cuts for the SBND inclusive
single-photon search with Wire-Cell.

Signal, adapted from the MicroBooNE inclusive photon-like LEE search (PRL,
doi:10.1103/89qs-4lcp):
  * exactly one reconstructable photon: true E >= 20 MeV, shower start >=
    SHOWER_INSET inside the TPC. No other such photon anywhere in the TPC;
    a pair < 20 deg apart is two photons (MicroBooNE counts it as one
    shower: MERGED_PAIRS_ARE_SIGNAL).
  * true vertex and shower start >= SIGNAL_FV_INSET = 10 cm inside the TPC
  * no primary muon or e+- (MicroBooNE allows muons < 100 MeV KE); any
    number of protons, charged pions and neutrons
Photons are the prompt ones: GENIE primaries and decay photons of primary
hadrons. Unlike MicroBooNE, single photons from outside the AV (dirt) are
background (category 12), as are those from 0-10 cm inside (category 11).

Truth columns come from make_truth.make_sp_truthdf; candidates join them on
(entry, tmatch_idx) and an unmatched candidate is cosmic. Reco cut chains
(PRESELECTION, UB_SELECTION) are applied in order by passes().
"""

import numpy as np
import pandas as pd

# --- Geometry and signal definition

ACTIVE_VOLUME = {"x": (-200.0, 200.0), "y": (-200.0, 200.0), "z": (0.0, 500.0)}  # cm, SBND TPC

PHOTON_E_MIN = 0.020         # GeV, true photon energy
SHOWER_INSET = 3.0           # cm, shower start from the TPC boundary for a visible photon
SIGNAL_FV_INSET = 10.0       # cm, FV for the signal vertex and photon conversion
MERGE_ANGLE = 20.0           # deg, photons closer than this are one shower (true_nshower)
MERGED_PAIRS_ARE_SIGNAL = False  # True: count showers (MicroBooNE), False: count photons

RESNUM_DELTA_1232 = 0        # genie::EResonance kP33_1232
ORIGIN_PRIMARY = 22          # true_shw1_origin of a GENIE-primary photon


def dist_to_active_boundary(x, y, z):
    """Signed distance [cm] to the active-volume boundary, positive inside
    (NaN for NaN coordinates)."""
    x, y, z = (np.asarray(v, dtype=float) for v in (x, y, z))
    b = ACTIVE_VOLUME
    d = np.minimum.reduce([
        x - b["x"][0], b["x"][1] - x,
        y - b["y"][0], b["y"][1] - y,
        z - b["z"][0], b["z"][1] - z,
    ])
    d[~(np.isfinite(x) & np.isfinite(y) & np.isfinite(z))] = np.nan
    return d


def _true_vtx_dist(df):
    return dist_to_active_boundary(df.true_vtx_x, df.true_vtx_y, df.true_vtx_z)


def true_in_active(df):
    return pd.Series(_true_vtx_dist(df) > 0, index=df.index)


def true_in_fv(df):
    return pd.Series(_true_vtx_dist(df) >= SIGNAL_FV_INSET, index=df.index)


def is_cosmic(df):
    """No true neutrino behind the row (never true in the truth table)."""
    m = df.true_E.isna()
    if "tmatch_idx" in df.columns:
        m |= df.tmatch_idx < 0
    return m


def n_visible(df):
    """Reconstructable showers (merged pairs as one) or photons."""
    return df.true_nshower if MERGED_PAIRS_ARE_SIGNAL else df.true_ngamma_rec


def one_shower_in_fv(df):
    """Exactly one visible photon in the TPC, its shower starting in the FV.
    true_ngamma_edge (converting 0.5-3 cm in) also counts as visible."""
    conv = dist_to_active_boundary(df.true_shw1_conv_x, df.true_shw1_conv_y, df.true_shw1_conv_z)
    return (n_visible(df) == 1) & (df.true_ngamma_edge == 0) & (conv >= SIGNAL_FV_INSET)


def is_signal(df):
    return one_shower_in_fv(df) & true_in_fv(df) & (df.true_nmu == 0) & (df.true_ne == 0)


# --- Truth categories

CATEGORY_LABELS = {
    0: "1g: NC Delta->Ngamma",
    1: "1g: NC pi0",
    3: "1g: other",
    4: ">=2 visible photons",
    5: "0 showers",
    6: "Primary muon (CC numu)",
    7: "CC nue (primary e)",
    8: "Outside AV, not 1g",
    9: "Cosmic / unmatched",
    11: "1g, not signal (vertex or conversion < 10 cm in)",
    12: "1g from outside the AV (dirt)",
}
SIGNAL_CATEGORIES = [0, 1, 3]


def is_signal_cat(cat):
    return cat.isin(SIGNAL_CATEGORIES)


def truth_category(df):
    """Integer category per row (CATEGORY_LABELS); first match wins."""
    o = df.true_shw1_origin
    sig = is_signal(df)
    nvis = n_visible(df)
    ncdelta = (df.true_isnc == 1) & (o == ORIGIN_PRIMARY) & (df.true_resnum == RESNUM_DELTA_1232)
    conds = [
        (is_cosmic(df), 9),
        (~true_in_active(df) & one_shower_in_fv(df), 12),
        (~true_in_active(df), 8),
        (df.true_ne > 0, 7),
        (df.true_nmu > 0, 6),
        (sig & ncdelta, 0),
        (sig & (df.true_isnc == 1) & (o == 111), 1),
        (sig, 3),
        (nvis + df.true_ngamma_edge >= 2, 4),
        (nvis == 1, 11),
        (nvis == 0, 5),
    ]
    cat = pd.Series(-1, index=df.index)
    for m, c in reversed(conds):
        cat[m.fillna(False).astype(bool)] = c
    return cat


# --- Reco preselection (Wire-Cell candidate columns)
#
# Not cut, unlike MicroBooNE: the TGM/STM/light-mismatch flags (Wire-Cell
# makes no candidate from a flagged cluster; 0.04% effect), the > 15 cm
# main-cluster length (-20% signal, no purity gain) and >= 1 EM shower
# (implied by one_shower).

FV_INSET = SIGNAL_FV_INSET   # cm, from every AV face
CATHODE_INSET = 5.0          # cm, from the cathode (x = 0); purity 27% -> 32%


def reco_in_fv(df, inset=FV_INSET):
    d = dist_to_active_boundary(df.nu_x, df.nu_y, df.nu_z)
    return pd.Series(d >= inset, index=df.index) & (df.nu_x.abs() >= CATHODE_INSET)


# Cuts on true neutrinos (cutflow applies them to mcnu): wc_valid = 0 is a
# production failure, such events have no candidates.
EVENT_CUTS = {"wc_valid"}

PRESELECTION = {
    "wc_valid": lambda df: df.wc_valid == 1,
    "vertex": lambda df: df.has_vertex == 1,   # the taggers ran
    "fv": reco_in_fv,
}


# --- MicroBooNE reference selection (microboone_single_photon_selection.yaml)
#
# MicroBooNE numu / other / NC pi0 BDTs (ub_*, microboone_bdt.score, NaN =
# not evaluable, which fails the cut) plus the SBND Wire-Cell numu_score and
# nue_score (the latter replaces the MicroBooNE nue BDT: ~2x the signal at
# the same purity). The cuts are ANDed; the order follows the MicroBooNE BDT
# cascade, each BDT trained on events passing the previous ones.

WC_NUMU_SCORE_MAX = 0.4  # 0.2-3 within 10% of the best eff x purity (Nov14 scan)
UB_NUMU_MIN = 1.0      # MicroBooNE: 0.4; best of the Nov14 scan (0.8-1.2 equivalent)
UB_OTHER_MIN = 0.2     # as MicroBooNE; best of the Nov14 scan
UB_NCPI0_MIN = -0.1    # highest purity on both file halves (Nov14; MicroBooNE: -0.05)
WC_NUE_SCORE_MAX = 7.0  # MicroBooNE's nue selection is nue_score > 7
BACK_DIST_MIN = 80.0   # cm; removes 24/25 dirt single photons for ~9% of signal (Nov8)

UB_SELECTION = {
    "wc_numu": lambda df: df.numu_score < WC_NUMU_SCORE_MAX,
    "ub_numu": lambda df: df.ub_numu > UB_NUMU_MIN,
    "ub_other": lambda df: df.ub_other > UB_OTHER_MIN,
    "ub_ncpi0": lambda df: df.ub_ncpi0 > UB_NCPI0_MIN,
    "wc_nue": lambda df: df.nue_score < WC_NUE_SCORE_MAX,
    # one EM shower > 20 MeV at the vertex passing the SP tagger's first check
    "one_shower": lambda df: df.shw_sp_n_20br1_showers == 1,
    # no shower labelled mu/pi/p: a second photon often gets a mu/pi start
    # segment, and one_shower does not count it
    "em_showers": lambda df: df.n_shower <= df.n_em,
    # active volume behind the shower start along its axis; short for dirt
    "back_traced": lambda df: df.shw1_back_dist > BACK_DIST_MIN,
}

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
