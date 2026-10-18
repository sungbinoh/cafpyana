"""Pandora NC Delta -> N gamma selection, adapted from MicroBooNE.

References
----------
[1] MicroBooNE, arXiv:2502.05750. Its Pandora channels are unchanged from [2].
[2] MicroBooNE, PRL 128, 111801 (2022), arXiv:2110.00409 (+ supplemental).
[3] Numeric pre-cuts: final Runs 1-3 Hive configs,
    github.com/markrosslonergan/hellstroms_hive, hive/xml/
      1g1p_xmls/hive2.5/1g1p_CinnabarWorldOrder_v33.3_15_04_30_Combined_Run123_v1_updatedBDTfiles.xml
      1g0p_xmls/hive2.5/1g0p_CinnabarWorldOrder_v33.3_4_29_20_CombinedRuns13_updatedBDTfiles.xml

Stages in [2]:
  1. Topological: one shower and one (1g1p) or zero (1g0p) tracks among the
     slice's primary PFPs.
  2. Pre-selection (PRESEL_*). Stages 1+2 are applied at dataframe creation.
  3. Background rejection. Here: rectangular cuts on CAF features
     (CUTBASED_*), one group per background class.
  4. Signal window in reconstructed shower energy (ENERGY_WINDOW).

Energies are MeV, lengths cm.
"""

from collections import OrderedDict

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Thresholds
# ---------------------------------------------------------------------------

# Track/shower classification. PFPs with NaN or negative score are counted
# as "unclassified", not as tracks or showers.
TRK_SCORE_CUT = 0.5

# MicroBooNE cut on distance to the space-charge boundary (SCB); here it is
# the distance to the active-volume boundary (dist_to_active_boundary()).
PRESEL_1G1P = {
    "shw_energy_min":      40.0,  # #GE# > 40 MeV (collection-plane energy)
    "vtx_dist_min":         2.0,  # reco_vertex_dist_to_SCB > 2
    "trk_dist_min":         2.0,  # min(track start, end dist to SCB) > 2
    "trk_dedx_min":         2.0,  # truncated-mean dE/dx > 2 MeV/cm
    "abs_cos_trk_shw_max":  0.99, # |cos(track, implied shower dir)| < 0.99
    "shw_start_dist_min":   7.0,  # reco_shower_start_dist_to_SCB > 7
    "trk_ke_max":         500.0,  # proton KE from range < 0.5 GeV
    "trk_calo_E_max":     400.0,  # track calorimetric energy < 0.4 GeV
}

PRESEL_1G0P = {
    "shw_energy_max_min": 30.0,   # reco_shower_energy_max > 30 MeV
    "shw_start_dist_min":  2.0,   # reco_shower_start_dist_to_SCB > 2
}

ENERGY_WINDOW = {
    "1g1p": (0.0, 600.0),
    "1g0p": (100.0, 700.0),
}

TOPOLOGIES = ("1g1p", "1g0p")

# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

# Outer active-volume boxes, as in analysis_village/mcs_pid/cuts.py without
# the FV insets. ICARUS has one box per cryostat.
ACTIVE_VOLUMES = {
    "SBND": [
        {"x": (-200.0, 200.0), "y": (-200.0, 200.0), "z": (0.0, 500.0)},
    ],
    "ICARUS": [
        {"x": (-358.49, -61.94), "y": (-181.86, 134.96), "z": (-894.950652270838, 894.950652270838)},
        {"x": (61.94, 358.49),   "y": (-181.86, 134.96), "z": (-894.950652270838, 894.950652270838)},
    ],
}


def _detector_key(detector):
    det = str(detector)
    if det.startswith("ICARUS"):
        return "ICARUS"
    if det == "SBND":
        return "SBND"
    raise ValueError("Unknown detector %r; expected 'SBND' or 'ICARUS[ RunN]'." % det)


def _dist_to_boxes(x, y, z, boxes):
    """Signed distance to the nearest face of the containing box.

    Positive inside (exact). Negative outside every box, with magnitude equal
    to the largest face violation (a lower bound). NaN coordinates give NaN.
    """
    best = np.full(len(x), -np.inf)
    for b in boxes:
        d = np.minimum.reduce([
            x - b["x"][0], b["x"][1] - x,
            y - b["y"][0], b["y"][1] - y,
            z - b["z"][0], b["z"][1] - z,
        ])
        best = np.maximum(best, d)
    best[~(np.isfinite(x) & np.isfinite(y) & np.isfinite(z))] = np.nan
    return best


def dist_to_active_boundary(x, y, z, detector):
    """Signed distance [cm] to the active-volume boundary.

    `detector` is a scalar or array-like aligned with the coordinates.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)

    if np.ndim(detector) == 0:
        return _dist_to_boxes(x, y, z, ACTIVE_VOLUMES[_detector_key(detector)])

    det = np.asarray(detector)
    out = np.full(len(x), np.nan)
    for d in pd.unique(det):
        m = det == d
        out[m] = _dist_to_boxes(x[m], y[m], z[m], ACTIVE_VOLUMES[_detector_key(d)])
    return out


# ---------------------------------------------------------------------------
# Topology
# ---------------------------------------------------------------------------

def topology_1g1p(df):
    return (df.n_shw == 1) & (df.n_trk == 1)


def topology_1g0p(df):
    return (df.n_shw == 1) & (df.n_trk == 0)


def topology_cut(df):
    return topology_1g1p(df) | topology_1g0p(df)


TOPOLOGY_CUTS = {"1g1p": topology_1g1p, "1g0p": topology_1g0p}

# ---------------------------------------------------------------------------
# Pre-selection (NaN comparisons are False, so missing quantities fail)
# ---------------------------------------------------------------------------

def presel_1g1p_cuts(df, cfg=PRESEL_1G1P):
    """Ordered dict of the individual 1g1p pre-cuts (topology excluded)."""
    c = cfg
    return OrderedDict([
        ("shw_energy",    df.shw_energy > c["shw_energy_min"]),
        ("vtx_contained", df.vtx_dist_to_wall > c["vtx_dist_min"]),
        ("trk_contained", np.minimum(df.trk_start_dist_to_wall, df.trk_end_dist_to_wall) > c["trk_dist_min"]),
        ("trk_dedx",      df.trk_dedx_trunc > c["trk_dedx_min"]),
        ("trk_shw_angle", df.cos_trk_shw_implied.abs() < c["abs_cos_trk_shw_max"]),
        ("shw_contained", df.shw_start_dist_to_wall > c["shw_start_dist_min"]),
        ("trk_ke",        df.trk_ke_p < c["trk_ke_max"]),
        ("trk_calo_E",    df.trk_calo_E < c["trk_calo_E_max"]),
    ])


def presel_1g0p_cuts(df, cfg=PRESEL_1G0P):
    """Ordered dict of the individual 1g0p pre-cuts (topology excluded)."""
    c = cfg
    return OrderedDict([
        ("shw_energy",    df.shw_energy_max > c["shw_energy_max_min"]),
        ("shw_contained", df.shw_start_dist_to_wall > c["shw_start_dist_min"]),
    ])


PRESEL_CUTS = {"1g1p": presel_1g1p_cuts, "1g0p": presel_1g0p_cuts}


def _all(masks, index):
    out = pd.Series(True, index=index)
    for m in masks.values():
        out &= m
    return out


def preselection_1g1p(df):
    return topology_1g1p(df) & _all(presel_1g1p_cuts(df), df.index)


def preselection_1g0p(df):
    return topology_1g0p(df) & _all(presel_1g0p_cuts(df), df.index)


def preselection(df):
    """Topology + pre-selection for either channel. Applied at df creation."""
    return preselection_1g1p(df) | preselection_1g0p(df)


PRESELECTION = {"1g1p": preselection_1g1p, "1g0p": preselection_1g0p}

# ---------------------------------------------------------------------------
# Cut-based background rejection
# ---------------------------------------------------------------------------
# Cuts are grouped by the background they target. Unless a source is given,
# values are untuned physics-motivated starting points: tune them on SBND MC
# with nminus1_masks(). Set a value to None to disable.

CUTBASED_1G1P = {
    # Cosmics
    "nu_score_min":            0.5,
    "bary_score_min":          1e-6,   # as mcs_pid baryscore_cut; SBND only
    "vtx_fv":                  True,   # mcs_pid FV on the slice vertex
    # CC nue: a photon shower starts at ~2 MIP (~4 MeV/cm), an electron at ~1
    "shw_dedx_min":            3.0,
    "shw_conv_dist_min":       1.0,    # MicroBooNE 2g1p pre-cut [2, suppl.]
    # NC pi0: veto a second (non-primary) shower in the slice
    "n_nonprim_shw_max":       0,
    # Proton + photon invariant mass window [MeV]. Off by default: it depends
    # on the shower energy scale, so locate the signal peak before using it.
    "delta_mass_range":        None,   # e.g. (1000., 1450.)
    # CC numu: track more proton-like than muon-like, and attached to vertex
    "trk_proton_like":         True,   # chi2_proton < chi2_muon
    "trk_start_dist_to_vtx_max": 10.0, # MicroBooNE 2g1p pre-cut [2, suppl.]
}

CUTBASED_1G0P = {
    "nu_score_min":      0.5,
    "bary_score_min":    1e-6,
    "shw_fv":            True,   # FV on the shower start (no track defines the vertex)
    "shw_dedx_min":      3.0,
    "n_nonprim_shw_max": 0,
}

CUTBASED = {"1g1p": CUTBASED_1G1P, "1g0p": CUTBASED_1G0P}


# Fiducial volume. A verbatim copy of the FV in analysis_village/mcs_pid/cuts.py
# (SBND high-yz exclusion, ICARUS Run 2 East-East exclusion and WW cable
# region), vendored so this file has no dependencies outside numpy/pandas and
# can be run outside cafpyana. Keep the two in sync if either changes.

SBND_FV = {
    "lowYZ":  {"x": (-200., 200.), "y": (-200., 200.), "z": (0., 250.)},
    "highYZ": {"x": (-200., 200.), "y": (-200., 100.), "z": (250., 500.)},
}
ICARUS_RUN2_FV = {
    "C0": {"x": (-210.22, -61.94), "y": (-181.86, 134.96), "z": (-894.950652270838, 894.950652270838)},
    "C1": {"x": (61.94, 358.49),   "y": (-181.86, 134.96), "z": (-894.950652270838, 894.950652270838)},
}
ICARUS_RUN4_FV = {
    "C0": {"x": (-358.49, -61.94), "y": (-181.86, 134.96), "z": (-894.950652270838, 894.950652270838)},
    "C1": {"x": (61.94, 358.49),   "y": (-181.86, 134.96), "z": (-894.950652270838, 894.950652270838)},
}
ICARUS_WW_CABLE = {"x": (210.22, 358.49), "y_min": 70., "z_min": 0.}


def fv_cut(x, y, z, detector, run, inx=10., iny=10., inzfront=10., inzback=50.):
    """True inside the fiducial volume; insets in cm from the FV boundaries."""
    x, y, z = (np.asarray(v, dtype=float) for v in (x, y, z))
    det = np.broadcast_to(np.asarray(detector), x.shape)
    run = np.broadcast_to(np.asarray(run), x.shape)

    lo, hi = SBND_FV["lowYZ"], SBND_FV["highYZ"]
    sbnd = (((x < lo["x"][1] - inx) & (x > lo["x"][0] + inx) &
             (y < lo["y"][1] - iny) & (y > lo["y"][0] + iny) &
             (z < lo["z"][1]) & (z > lo["z"][0] + inzfront)) |
            ((x < hi["x"][1] - inx) & (x > hi["x"][0] + inx) &
             (y < hi["y"][1] - iny) & (y > hi["y"][0] + iny) &
             (z < hi["z"][1] - inzback) & (z > hi["z"][0])))

    cable = ((x > ICARUS_WW_CABLE["x"][0]) & (x < ICARUS_WW_CABLE["x"][1]) &
             (y > ICARUS_WW_CABLE["y_min"]) & (z > ICARUS_WW_CABLE["z_min"]))

    def icarus(fv):
        c0, c1 = fv["C0"], fv["C1"]
        in_x = (((x < c0["x"][1] - inx) & (x > c0["x"][0] + inx)) |
                ((x < c1["x"][1] - inx) & (x > c1["x"][0] + inx)))
        return (in_x & (y < c0["y"][1] - iny) & (y > c0["y"][0] + iny) &
                (z < c0["z"][1] - inzback) & (z > c0["z"][0] + inzfront) & ~cable)

    is_run2 = (det == "ICARUS Run2") | ((det == "ICARUS") & (run == 2))
    is_run4 = (det == "ICARUS Run4") | ((det == "ICARUS") & (run == 4))
    return np.select([is_run2, is_run4, det == "SBND"],
                     [icarus(ICARUS_RUN2_FV), icarus(ICARUS_RUN4_FV), sbnd], default=False)


def _fv(df, prefix):
    return pd.Series(fv_cut(df[prefix + "_x"], df[prefix + "_y"], df[prefix + "_z"],
                            df.detector, df.Run), index=df.index)


def cutbased_cuts(df, topology, cfg=None):
    """Ordered dict of the cut-based stage masks for one channel.

    The flash-match cut is skipped when `bary_score` is absent (ICARUS CAFs).
    """
    c = dict(CUTBASED[topology]) if cfg is None else cfg
    on = lambda k: c.get(k) is not None and c.get(k) is not False
    out = OrderedDict()

    # cosmics
    if on("nu_score_min"):
        out["nu_score"] = df.nu_score > c["nu_score_min"]
    if on("bary_score_min") and "bary_score" in df.columns:
        out["flash_match"] = df.bary_score > c["bary_score_min"]
    if on("vtx_fv"):
        out["vtx_fv"] = _fv(df, "slc_vtx")
    if on("shw_fv"):
        out["shw_fv"] = _fv(df, "shw_start")

    # CC nue
    if on("shw_dedx_min"):
        out["shw_dedx"] = df.shw_dedx_best > c["shw_dedx_min"]
    if on("shw_conv_dist_min"):
        out["shw_conv_dist"] = df.shw_conv_dist > c["shw_conv_dist_min"]

    # NC pi0
    if on("n_nonprim_shw_max"):
        out["no_second_shw"] = df.n_nonprim_shw <= c["n_nonprim_shw_max"]
    if on("delta_mass_range"):
        lo, hi = c["delta_mass_range"]
        out["delta_mass"] = (df.delta_mass > lo) & (df.delta_mass < hi)

    # CC numu
    if on("trk_proton_like"):
        out["trk_proton_like"] = df.trk_chi2_proton < df.trk_chi2_muon
    if on("trk_start_dist_to_vtx_max"):
        out["trk_attached"] = df.trk_start_dist_to_vtx < c["trk_start_dist_to_vtx_max"]

    return out


def energy_window_cut(df, topology):
    lo, hi = ENERGY_WINDOW[topology]
    return (df.shw_energy >= lo) & (df.shw_energy < hi)


def stage_cuts(df, topology, cutbased=True, energy_window=True):
    """Ordered dict of every selection stage mask for one channel.

    cutbased: True for the defaults, a dict to override CUTBASED[topology]
    (copy and edit it), or False to skip the stage.
    """
    stages = OrderedDict([("topology", TOPOLOGY_CUTS[topology](df))])
    stages.update(PRESEL_CUTS[topology](df))
    if cutbased is not False:
        stages.update(cutbased_cuts(df, topology, None if cutbased is True else cutbased))
    if energy_window:
        stages["energy_window"] = energy_window_cut(df, topology)
    return stages


def selection(df, topology, **kwargs):
    """Full selection mask for one channel ("1g1p" or "1g0p")."""
    return _all(stage_cuts(df, topology, **kwargs), df.index)


def nminus1_masks(df, topology, **kwargs):
    """{stage: mask of every other stage}, for tuning one cut at a time."""
    stages = stage_cuts(df, topology, **kwargs)
    return OrderedDict(
        (name, _all(OrderedDict((k, m) for k, m in stages.items() if k != name), df.index))
        for name in stages
    )


def _weights(df, weight):
    if weight is None:
        return pd.Series(1.0, index=df.index)
    if isinstance(weight, str):
        return df[weight]
    return weight


def cutflow(df, topology, weight=None, **kwargs):
    """Cumulative yield per stage, by truth category if truth is present.

    On preselected dataframes the topology and pre-cut rows are all equal;
    use the *_nopresel dataframes to see those stages.
    """
    w = _weights(df, weight)
    has_truth = "true_pdg" in df.columns
    cat = truth_category(df) if has_truth else None

    rows = []
    running = pd.Series(True, index=df.index)
    for name, m in stage_cuts(df, topology, **kwargs).items():
        running = running & m
        row = {"stage": name, "total": w[running].sum()}
        if has_truth:
            for code, label in CATEGORY_LABELS.items():
                row[label] = w[running & (cat == code)].sum()
        rows.append(row)
    return pd.DataFrame(rows).set_index("stage")


# ---------------------------------------------------------------------------
# Truth categories ([2], Table II). Works on the slice df and the mcnu df.
# ---------------------------------------------------------------------------

CATEGORY_LABELS = OrderedDict([
    (0, "NC Delta->Ngamma"),
    (1, "NC 1pi0 coherent"),
    (2, "NC 1pi0 non-coherent"),
    (3, "CC numu 1pi0"),
    (4, "CC nue/nuebar"),
    (5, "BNB other"),
    (6, "Dirt (outside TPC)"),
    (7, "Cosmic"),
])

GENIE_MODE_RES = 1
GENIE_MODE_COH = 3


def true_in_active(df):
    return pd.Series(
        dist_to_active_boundary(df.true_vtx_x, df.true_vtx_y, df.true_vtx_z, df.detector) > 0,
        index=df.index)


def is_cosmic(df):
    """Slice not matched to a neutrino (never true in the mcnu df)."""
    m = df.true_pdg == -1
    if "tmatch_idx" in df.columns:
        m |= df.tmatch_idx < 0
    return m


def is_ncdelta_rad(df):
    """NC Delta -> N gamma from GENIE truth: NC, RES, >= 1 primary photon.

    CAF truth has no resonance PDG, so heavier-resonance radiative decays can
    pass. Photons from pi0 decay are not primaries.
    """
    return ((df.true_isnc == 1) & (df.true_genie_mode == GENIE_MODE_RES)
            & (df.true_ngamma >= 1))


def truth_category(df):
    """Integer category per row; see CATEGORY_LABELS."""
    cosmic = is_cosmic(df)
    nu = ~cosmic
    av = true_in_active(df) & nu
    nc = df.true_isnc == 1
    cc = df.true_iscc == 1
    npi0 = df.true_npi0

    conds = [
        cosmic,
        nu & ~av,
        av & is_ncdelta_rad(df),
        av & nc & (npi0 == 1) & (df.true_genie_mode == GENIE_MODE_COH),
        av & nc & (npi0 == 1),
        av & cc & (df.true_pdg == 14) & (npi0 == 1),  # numubar -> other, as in [3]
        av & cc & (df.true_pdg.abs() == 12),
    ]
    codes = [7, 6, 0, 1, 2, 3, 4]
    return pd.Series(np.select(conds, codes, default=5), index=df.index)


def is_signal(df):
    """True NC Delta -> N gamma with vertex in the active volume."""
    return is_ncdelta_rad(df) & true_in_active(df)
