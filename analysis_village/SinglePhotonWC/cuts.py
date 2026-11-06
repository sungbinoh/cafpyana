"""Signal definition and truth categories for the inclusive single-photon
search with Wire-Cell (SBND).

Signal, following the MicroBooNE inclusive photon-like LEE search
(PRL, doi:10.1103/89qs-4lcp), which is modelled on MiniBooNE's single
Cherenkov ring:

  * exactly one reconstructable photon shower. A reconstructable photon has
    true E >= 20 MeV and starts its shower (pair conversion) >= 3 cm inside
    the TPC boundary. Photons from the same vertex closer than 20 deg are
    one shower.
  * vertex in the active volume: no primary muon with KE >= 100 MeV (below
    MiniBooNE's threshold) and no primary electron/positron
  * vertex outside the active volume: any interaction type (MicroBooNE's
    "out of FV 1g": a photon entering and showering in the TPC). The
    definition is what MiniBooNE could see, and a lepton from outside the
    detector would not have been
  * any number of protons, charged pions, neutrons; any number of photons
    that are not reconstructable (soft, escaping, converting near a wall)

The signal categories are MicroBooNE's: NC pi0 1g (one photon lost or the
two merged), NC Delta 1g, NC other 1g, numu CC 1g (muon < 100 MeV, any
photon origin) and out of FV 1g. Photons are the prompt
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
    """MicroBooNE-style: one reconstructable photon shower in the TPC; the
    muon and electron vetoes apply only to vertices in the active volume."""
    lepton_ok = (df.true_nmu_100MeV == 0) & (df.true_ne == 0)
    return (df.true_nshower == 1) & (lepton_ok | ~true_in_active(df))


# --- Truth categories

RESNUM_DELTA_1232 = 0   # genie::EResonance kP33_1232
ORIGIN_PRIMARY = 22     # true_shw1_origin of a GENIE-primary photon

CATEGORY_LABELS = {
    0: "1g: NC Delta->Ngamma",
    1: "1g: NC pi0",
    2: "1g: numu CC, mu < 100 MeV",
    3: "1g: NC other",
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
        (sig & (df.true_iscc == 1), 2),
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
#   fv          reco vertex >= 5 cm from the anodes and the cathode (drift
#               direction only, as MicroBooNE)
#
# Not cut here, unlike MicroBooNE:
#   * the generic neutrino selection's TGM/STM/light-mismatch flags: SBND
#     Wire-Cell makes no candidate from a flagged cluster (nu_skip_cosmic),
#     so the cut removed 0.04% of candidates. The cosmic tagger (cosmict_*)
#     feeds the "other" BDT, as in MicroBooNE.
#   * the > 15 cm main-cluster length: it removed ~20% of signal (soft
#     photons) at no purity gain. Wire-Cell's own 15 cm floor exempts the
#     event's longest cluster, so it does not do this.
#   * >= 1 EM shower: nearly implied by vertex_shower below (removing it
#     changed the final efficiency by +0.05%).

FV_INSET = 5.0                 # cm


def reco_in_fv(df, inset=FV_INSET):
    ax = df.nu_x.abs()
    return (ax >= inset) & (ax <= ACTIVE_VOLUME["x"][1] - inset)


EVENT_CUTS = {"wc_valid"}

PRESELECTION = {
    "wc_valid": lambda df: df.wc_valid == 1,
    "vertex": lambda df: df.has_vertex == 1,
    "fv": reco_in_fv,
}


# --- MicroBooNE reference selection (microboone_single_photon_selection.yaml)
#
# On the MicroBooNE numu, other and NC pi0 BDT scores, columns ub_<name>
# (microboone_bdt.score; NaN = not evaluable), plus the SBND Wire-Cell
# tagger scores: numu_score (removes CC numu the MicroBooNE numu BDT keeps)
# and nue_score instead of the MicroBooNE nue BDT (~2x the signal at the
# same purity). The yaml's MicroBooNE-geometry FV is replaced by
# the SBND fv cut of the preselection.
#   wc_numu         not numu CC-like: numu_score < 0.5
#   ub_numu         MicroBooNE numu BDT > 1 (MicroBooNE: 0.4; tighter on SBND,
#                   tuned out-of-fold, as wc_numu already cut)
#   ub_other        MicroBooNE other BDT > -0.25 (MicroBooNE: 0.2; it barely
#                   separates on SBND, so looser keeps ~30% more signal at
#                   the same purity; tuned two-fold on the Nov6 MC)
#   ub_ncpi0        MicroBooNE NC pi0 BDT > -0.05, as MicroBooNE (log10
#                   score scale, see microboone_bdt)
#   wc_nue          not nue-like: nue_score < 7 (MicroBooNE's nue selection
#                   is nue_score > 7)
#   one_shower      exactly one EM shower > 20 MeV at the vertex passing the
#                   single-photon tagger's first reconstruction check
#
# MicroBooNE's ">= 1 shower > 20 MeV at the vertex" and "BDTs evaluable"
# steps are not cut: one_shower and the BDT thresholds (NaN = not
# evaluable fails them) imply both.
# The cuts are ANDed, so the order only shapes the cut flow. It follows the
# MicroBooNE BDT cascade (each trained on events passing the previous ones;
# the NC pi0 BDT has no separation before "other"), so tune them in order.
UB_NUMU_MIN = 1.0
UB_OTHER_MIN = -0.25
UB_NCPI0_MIN = -0.05
WC_NUMU_SCORE_MAX = 0.5
WC_NUE_SCORE_MAX = 7.0

UB_SELECTION = {
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
