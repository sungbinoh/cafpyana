"""Geometry and truth categories for the Wire-Cell NC Delta -> N gamma analysis.

Self-contained: nothing here imports another analysis_village directory.
Truth columns are those of make_ncdeltawc.make_ncdeltawc_truthdf (GENIE truth
from rec.mc.nu), joined to candidates through the Wire-Cell truth match
(tmatch_idx); a candidate with no match is "cosmic".
"""

from collections import OrderedDict

import numpy as np
import pandas as pd

# --- Geometry (SBND only)

# Active volume: the TPC box, no insets [cm]
ACTIVE_VOLUME = {"x": (-200.0, 200.0), "y": (-200.0, 200.0), "z": (0.0, 500.0)}


def _check_sbnd(detector):
    det = np.unique(np.asarray(detector))
    if not all(str(d) == "SBND" for d in det):
        raise ValueError("NCDeltaWC is SBND-only; got detector = %s" % det)


def dist_to_active_boundary(x, y, z, detector="SBND"):
    """Signed distance [cm] to the active-volume boundary, positive inside.

    NaN coordinates give NaN.
    """
    _check_sbnd(detector)
    x, y, z = (np.asarray(v, dtype=float) for v in (x, y, z))
    b = ACTIVE_VOLUME
    d = np.minimum.reduce([
        x - b["x"][0], b["x"][1] - x,
        y - b["y"][0], b["y"][1] - y,
        z - b["z"][0], b["z"][1] - z,
    ])
    d[~(np.isfinite(x) & np.isfinite(y) & np.isfinite(z))] = np.nan
    return d


# --- Truth

GENIE_MODE_COH = 3
RESNUM_DELTA_1232 = 0   # genie::EResonance kP33_1232; genie_mode == kRes is ANY resonance

CATEGORY_LABELS = OrderedDict([
    (0, "NC Delta->Ngamma"),
    (1, "NC 1pi0 coherent"),
    (2, "NC 1pi0 non-coherent"),
    (3, "CC numu 1pi0"),
    (4, "CC nue/nuebar"),
    (5, "BNB other"),
    (6, "Dirt (outside TPC)"),
    (7, "Cosmic / unmatched"),
])


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


def is_ncdelta_rad(df):
    """NC Delta(1232) -> N gamma with a single photon, from GENIE truth.

    Exactly one primary photon above 20 MeV: GENIE primaries also include
    ~3-9 MeV nuclear de-excitation gammas, which do not shower. pi0 decay
    photons are not GENIE primaries, hence npi0 == 0.
    """
    return ((df.true_isnc == 1)
            & (df.true_resnum == RESNUM_DELTA_1232)
            & (df.true_ngamma_20MeV == 1)
            & (df.true_npi0 == 0))


def truth_category(df):
    """Integer category per row; see CATEGORY_LABELS."""
    cosmic = is_cosmic(df)
    nu = ~cosmic
    av = true_in_active(df) & nu
    nc = df.true_isnc == 1
    cc = df.true_iscc == 1
    npi0 = df.true_npi0
    coh = df.true_genie_mode == GENIE_MODE_COH
    # Mutually exclusive in this order; anything else is 5, "BNB other"
    conds = [
        cosmic,                                          # 7
        nu & ~av,                                        # 6 dirt
        av & is_ncdelta_rad(df),                         # 0 signal
        av & nc & (npi0 == 1) & coh,                     # 1
        av & nc & (npi0 == 1) & ~coh,                    # 2
        av & cc & (df.true_pdg == 14) & (npi0 == 1),     # 3; numubar -> other
        av & cc & (df.true_pdg.abs() == 12),             # 4
    ]
    return pd.Series(np.select(conds, [7, 6, 0, 1, 2, 3, 4], default=5), index=df.index)


def is_signal(df):
    """True NC Delta -> N gamma with vertex in the active volume."""
    return is_ncdelta_rad(df) & true_in_active(df)


def is_signal_1g0p(df):
    """Signal whose final state is the photon and nucleons only, with no
    proton above 35 MeV KE (neutrons and softer protons allowed).
    Needs the make_truth.py columns."""
    return is_signal(df) & (df.true_nother == 0) & (df.true_np_35MeV == 0)


def is_signal_1g1p(df):
    """As is_signal_1g0p, with exactly one proton above 35 MeV KE."""
    return is_signal(df) & (df.true_nother == 0) & (df.true_np_35MeV == 1)
