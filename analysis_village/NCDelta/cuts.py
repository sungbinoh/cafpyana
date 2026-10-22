"""Pandora NC Delta(1232) -> N gamma selection for SBND, adapted from [1].

Energies GeV, lengths cm; dE/dx and PIDA keep the conventional MeV/cm.
Cut values are provisional: not holdout-validated, on a small signal sample.

[1] MicroBooNE, arXiv:2502.05750 (Pandora channels unchanged from [2]).
[2] MicroBooNE, PRL 128, 111801 (2022), arXiv:2110.00409 (+ supplemental).
[3] Numeric pre-cuts: Runs 1-3 Hive configs, github.com/markrosslonergan/
    hellstroms_hive, hive/xml/{1g1p,1g0p}_xmls/hive2.5
"""

from collections import OrderedDict

import numpy as np
import pandas as pd

# cuts.py is imported both as analysis_village.NCDelta.cuts and, by the
# standalone scripts, as plain "cuts" with this directory on sys.path.
try:
    from analysis_village.NCDelta.shower_energy import (      # noqa: F401
        SHW_ENERGY_SCALE, shw_best_plane, shw_energy_corrected)
except ImportError:
    from shower_energy import (                               # noqa: F401
        SHW_ENERGY_SCALE, shw_best_plane, shw_energy_corrected)

# --- Thresholds

# PFPs with NaN or negative score count as neither track nor shower.
TRK_SCORE_CUT = 0.5

# Applied at dataframe creation, so kept minimal: anything here is baked in.
PRESEL_1G1P = {
    "not_clear_cosmic": True,   # Pandora's neutrino candidate only, as in [3]
    "vtx_fv": True,
    "shw_fv": True,
}

PRESEL_1G0P = {
    "not_clear_cosmic": True,
    "vtx_fv": True,
    "shw_fv": True,
}

# GeV. The measurement's signal region, so common to both working points.
ENERGY_WINDOW = {
    "1g1p": (0.020, 0.700),
    "1g0p": (0.120, 0.600),
}

TOPOLOGIES = ("1g1p", "1g0p")

# --- Geometry (SBND only)

# Active volume: the TPC box, no insets.
ACTIVE_VOLUME = {"x": (-200.0, 200.0), "y": (-200.0, 200.0), "z": (0.0, 500.0)}

# Fiducial volume: a copy of the SBND one in analysis_village/mcs_pid/cuts.py.
SBND_FV = {
    "lowYZ":  {"x": (-200., 200.), "y": (-200., 200.), "z": (0., 250.)},
    "highYZ": {"x": (-200., 200.), "y": (-200., 100.), "z": (250., 500.)},
}


def _check_sbnd(detector):
    """This analysis is SBND-only; fail loudly on anything else."""
    det = np.unique(np.asarray(detector))
    if not all(str(d) == "SBND" for d in det):
        raise ValueError("NCDelta is SBND-only; got detector = %s" % det)


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


def fv_cut(x, y, z, detector="SBND", inx=10., iny=10., inzfront=10., inzback=50.):
    """True inside the fiducial volume; insets in cm from the FV boundaries."""
    _check_sbnd(detector)
    x, y, z = (np.asarray(v, dtype=float) for v in (x, y, z))
    lo, hi = SBND_FV["lowYZ"], SBND_FV["highYZ"]
    return (((x < lo["x"][1] - inx) & (x > lo["x"][0] + inx) &
             (y < lo["y"][1] - iny) & (y > lo["y"][0] + iny) &
             (z < lo["z"][1]) & (z > lo["z"][0] + inzfront)) |
            ((x < hi["x"][1] - inx) & (x > hi["x"][0] + inx) &
             (y < hi["y"][1] - iny) & (y > hi["y"][0] + iny) &
             (z < hi["z"][1] - inzback) & (z > hi["z"][0])))


def _fv(df, prefix):
    return pd.Series(fv_cut(df[prefix + "_x"], df[prefix + "_y"], df[prefix + "_z"],
                            df.detector), index=df.index)


def extrap_det_path(df):
    """Distance back along the shower direction to the AV edge [cm].

    How far the photon could have travelled inside the detector before
    converting. Small for dirt, which enters through a nearby wall.
    """
    b = ACTIVE_VOLUME
    start = (sv(df, "shw_start_x").to_numpy(dtype=float),
             sv(df, "shw_start_y").to_numpy(dtype=float),
             sv(df, "shw_start_z").to_numpy(dtype=float))
    back = (-sv(df, "shw_dir_x").to_numpy(dtype=float),
            -sv(df, "shw_dir_y").to_numpy(dtype=float),
            -sv(df, "shw_dir_z").to_numpy(dtype=float))
    t = np.full(len(df), np.inf)
    for c, dc, (lo, hi) in zip(start, back, (b["x"], b["y"], b["z"])):
        with np.errstate(divide="ignore", invalid="ignore"):
            for edge in (lo, hi):
                ti = (edge - c) / dc
                t = np.minimum(t, np.where(ti > 0, ti, np.inf))
    t[~np.isfinite(t)] = np.nan        # direction parallel to every face
    return pd.Series(t, index=df.index)


# --- Topology

def topology_1g1p(df):
    return (df.n_shw == 1) & (df.n_trk == 1)


def topology_1g0p(df):
    return (df.n_shw == 1) & (df.n_trk == 0)


def topology_cut(df):
    return topology_1g1p(df) | topology_1g0p(df)


TOPOLOGY_CUTS = {"1g1p": topology_1g1p, "1g0p": topology_1g0p}

# --- Pre-selection (NaN comparisons are False, so missing quantities fail)

def presel_1g1p_cuts(df, cfg=PRESEL_1G1P):
    """Ordered dict of the individual 1g1p pre-cuts (topology excluded)."""
    out = OrderedDict()
    if cfg.get("not_clear_cosmic"):
        out["not_clear_cosmic"] = ~df.is_clear_cosmic.astype(bool)
    if cfg.get("vtx_fv"):
        out["vtx_fv"] = _fv(df, "slc_vtx")
    if cfg.get("shw_fv"):
        out["shw_fv"] = _fv(df, "shw_start")
    return out


def presel_1g0p_cuts(df, cfg=PRESEL_1G0P):
    """Ordered dict of the individual 1g0p pre-cuts (topology excluded)."""
    out = OrderedDict()
    if cfg.get("not_clear_cosmic"):
        out["not_clear_cosmic"] = ~df.is_clear_cosmic.astype(bool)
    if cfg.get("vtx_fv"):
        out["vtx_fv"] = _fv(df, "slc_vtx")
    if cfg.get("shw_fv"):
        out["shw_fv"] = _fv(df, "shw_start")
    return out


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

# --- Cut-based background rejection
# Every cut a working point can set. A config lists only the ones it uses;
# anything absent is off. Thresholds are repeated per working point on purpose
# -- each config is readable on its own.
CUT_KEYS = (
    "nu_score_min", "bary_score_min", "opt0_score_min",
    "shw_dedx_min", "shw_dedx_max", "shw_nhits_min",
    "shw_len_min", "shw_open_angle_max", "shw_trackScore_max",
    "shw_conv_dist_min", "shw_conv_dist_max", "shw_start_dist_min",
    "extrap_det_path_min", "vtx_dist_min", "trk_dist_min",
    "trk_ke_max", "trk_calo_E_max", "trk_pida_max", "trk_proton_like",
    "trk_start_dist_to_vtx_max", "abs_cos_trk_shw_max", "n_nonprim_shw_max",
    "delta_mass_range", "energy_window",
)

# Maximizes S/sqrt(S+B).
BALANCED_1G1P = {
    "nu_score_min":        0.45,
    "bary_score_min":      0.05,
    "opt0_score_min":      100.,   # higher = better match; poorly validated
    "shw_dedx_min":        3.5,    # MeV/cm; photon ~2 MIP, electron ~1
    "shw_dedx_max":        6.0,    # shw_dedx_min only bounds below; protons pass
    "shw_conv_dist_max":   60.,    # X0 = 14 cm, so 60 cm is misreco
    "shw_trackScore_max":  0.44,   # track-like fakes pile up under TRK_SCORE_CUT
    "trk_proton_like":     True,   # chi2_proton < chi2_muon
}

BALANCED_1G0P = {
    "bary_score_min":      0.05,   # nu_score_min is inert here, so omitted
    "opt0_score_min":      700.,
    "shw_dedx_min":        3.0,    # looser: dedx is measured worse with no track
    "shw_dedx_max":        10.0,   # looser: two-proton events are rarer
    "shw_trackScore_max":  0.42,
    "extrap_det_path_min": 60.,    # dirt; the FV cannot do this on a shower start
}

# Maximizes purity.
HIGH_PURITY_1G1P = {
    "nu_score_min":        0.50,
    "bary_score_min":      0.05,
    "opt0_score_min":      100.,
    "shw_dedx_min":        3.5,
    "shw_conv_dist_min":   4.0,    # shw_dedx_max off: len/open_angle already
    "shw_conv_dist_max":   50.,    # drop the stubby protons it would catch
    "shw_len_min":         30.,
    "shw_open_angle_max":  0.15,
    "shw_trackScore_max":  0.42,
    "trk_calo_E_max":      0.300,
    "trk_proton_like":     True,
}

HIGH_PURITY_1G0P = {
    "nu_score_min":        0.50,
    "bary_score_min":      0.10,
    "opt0_score_min":      700.,
    "shw_dedx_min":        2.5,
    "shw_dedx_max":        10.0,
    "shw_open_angle_max":  0.20,
    "shw_trackScore_max":  0.44,
    "extrap_det_path_min": 80.,
}

SELECTIONS = {
    "balanced":    {"1g1p": BALANCED_1G1P,    "1g0p": BALANCED_1G0P},
    "high_purity": {"1g1p": HIGH_PURITY_1G1P, "1g0p": HIGH_PURITY_1G0P},
}
DEFAULT_SELECTION = "balanced"

# Back-compat names for evaluate.py and plot_composition.py.
CUTBASED_1G1P = BALANCED_1G1P
CUTBASED_1G0P = BALANCED_1G0P
CUTBASED = SELECTIONS[DEFAULT_SELECTION]


def selection_config(name, topology):
    """The named working point's cuts for one channel."""
    if name not in SELECTIONS:
        raise KeyError("Unknown selection %r; have %s" % (name, sorted(SELECTIONS)))
    return dict(SELECTIONS[name][topology])


def cutbased_cuts(df, topology, cfg=None):
    """Ordered dict of the cut-based stage masks for one channel.

    The flash-match cut is skipped when `bary_score` is absent.
    """
    c = dict(CUTBASED[topology]) if cfg is None else cfg
    on = lambda k: c.get(k) is not None and c.get(k) is not False
    out = OrderedDict()

    # cosmics
    if on("nu_score_min"):
        out["nu_score"] = df.nu_score > c["nu_score_min"]
    if on("bary_score_min") and "bary_score" in df.columns:
        out["flash_match"] = sv(df, "bary_score") > c["bary_score_min"]
    if on("opt0_score_min"):
        out["opt0_score"] = df.opt0_score > c["opt0_score_min"]

    # containment (distances to the active-volume boundary)
    if on("vtx_dist_min"):
        out["vtx_contained"] = df.vtx_dist_to_wall > c["vtx_dist_min"]
    if on("trk_dist_min"):
        out["trk_contained"] = np.minimum(df.trk_start_dist_to_wall,
                                          df.trk_end_dist_to_wall) > c["trk_dist_min"]
    if on("shw_start_dist_min"):
        out["shw_contained"] = df.shw_start_dist_to_wall > c["shw_start_dist_min"]

    # track energy / PID thresholds
    if on("abs_cos_trk_shw_max"):
        out["trk_shw_angle"] = df.cos_trk_shw_implied.abs() < c["abs_cos_trk_shw_max"]
    if on("trk_ke_max"):
        out["trk_ke"] = df.trk_ke_p < c["trk_ke_max"]
    if on("trk_calo_E_max"):
        out["trk_calo_E"] = df.trk_calo_E < c["trk_calo_E_max"]

    # CC nue
    if on("shw_dedx_min"):
        out["shw_dedx"] = sv(df, "shw_dedx_best") > c["shw_dedx_min"]
    if on("shw_dedx_max"):
        out["shw_dedx_max"] = sv(df, "shw_dedx_best") < c["shw_dedx_max"]
    if on("shw_conv_dist_min"):
        out["shw_conv_dist"] = df.shw_conv_dist > c["shw_conv_dist_min"]

    # shower quality
    if on("shw_conv_dist_max"):
        out["shw_conv_dist_max"] = df.shw_conv_dist < c["shw_conv_dist_max"]
    if on("shw_len_min"):
        out["shw_len"] = df.shw_len > c["shw_len_min"]
    if on("shw_trackScore_max"):
        out["shw_trackScore"] = df.shw_trackScore < c["shw_trackScore_max"]
    if on("extrap_det_path_min"):
        out["extrap_det_path"] = extrap_det_path(df) > c["extrap_det_path_min"]
    if on("shw_nhits_min"):
        out["shw_nhits"] = (df.shw_nhits_I0 + df.shw_nhits_I1
                            + df.shw_nhits_I2) > c["shw_nhits_min"]
    if on("shw_open_angle_max"):
        out["shw_open_angle"] = sv(df, "shw_open_angle") < c["shw_open_angle_max"]

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
    if on("trk_pida_max"):
        out["trk_pida"] = df.trk_pida < c["trk_pida_max"]

    return out


def energy_window(topology, cfg=None):
    """(lo, hi) for a channel: the CUTBASED override if set, else the default."""
    if cfg is not None and cfg.get("energy_window") is not None:
        return tuple(cfg["energy_window"])
    return ENERGY_WINDOW[topology]


def energy_window_cut(df, topology, cfg=None):
    lo, hi = energy_window(topology, cfg)
    e = shw_energy_corrected(df)
    return (e >= lo) & (e < hi)


def stage_cuts(df, topology, cutbased=True, energy_window=True):
    """Ordered dict of every selection stage mask for one channel.

    cutbased: True for the defaults, a dict to override CUTBASED[topology]
    (copy and edit it), or False to skip the stage. An "energy_window" key in
    that dict overrides ENERGY_WINDOW[topology], so a whole configuration --
    signal window included -- can be given as one dict or one --config JSON.
    """
    cfg = None if cutbased in (True, False) else cutbased
    stages = OrderedDict([("topology", TOPOLOGY_CUTS[topology](df))])
    stages.update(PRESEL_CUTS[topology](df))
    if cutbased is not False:
        stages.update(cutbased_cuts(df, topology, cfg))
    if energy_window:
        stages["energy_window"] = energy_window_cut(df, topology, cfg)
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


# --- Truth categories ([2], Table II). Works on the slice df and the mcnu df.

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

# CAF "no value" markers. tmatch_idx/true_pdg excluded: is_cosmic() reads them.
SENTINELS = {
    "shw_dedx_best": (-5.,), "shw_dedx_I2": (-999.,),
    "shw_density": (-5.,), "shw_open_angle": (-5.,),
    "shw_dir_x": (-999.,), "shw_dir_y": (-999.,), "shw_dir_z": (-999.,),
    "shw_start_x": (-999.,), "shw_start_y": (-999.,), "shw_start_z": (-999.,),
    "bary_score": (-10000.,), "bary_deltaZ": (-10000.,),
    "crlongtrkdiry": (-10000.,),
}


def sv(df, col):
    """Column with its CAF sentinels masked to NaN. See SENTINELS."""
    v = df[col]
    bad = SENTINELS.get(col)
    return v.mask(v.isin(bad)) if bad else v


GENIE_MODE_RES = 1
GENIE_MODE_COH = 3

# genie::EResonance (BaryonResonance.h). genie_mode == kRes is ANY resonance.
RESNUM_NO_RESONANCE = -1
RESNUM_DELTA_1232 = 0      # kP33_1232
RESNUM = {
    -999: "no res (unmatched slice)",   # CAF sentinel: no truth interaction
    -1: "none",                         # kNoResonance: not resonant production
    0:  "P33(1232)",   1:  "S11(1535)",   2:  "D13(1520)",   3:  "S11(1650)",
    4:  "D13(1700)",   5:  "D15(1675)",   6:  "S31(1620)",   7:  "D33(1700)",
    8:  "P11(1440)",   9:  "P33(1600)",   10: "P13(1720)",   11: "F15(1680)",
    12: "P31(1910)",   13: "P33(1920)",   14: "F35(1905)",   15: "F37(1950)",
    16: "P11(1710)",   17: "F17(1970)",
}


def true_in_active(df):
    return pd.Series(
        dist_to_active_boundary(df.true_vtx_x, df.true_vtx_y, df.true_vtx_z) > 0,
        index=df.index)


def is_cosmic(df):
    """Slice not matched to a neutrino (never true in the mcnu df)."""
    m = df.true_pdg == -1
    if "tmatch_idx" in df.columns:
        m |= df.tmatch_idx < 0
    return m

def is_ncdelta_rad(df):
    """NC Delta(1232) -> N gamma from GENIE truth.

    resnum pins the Delta(1232); genie_mode == kRes is any resonance. GENIE
    primaries include ~3 MeV de-excitation gammas, hence the hard photon.
    """
    return ((df.true_isnc == 1)
            & (df.true_resnum == RESNUM_DELTA_1232)
            & (df.true_ngamma_100MeV >= 1)
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
    # Mutually exclusive by construction; anything unmatched is 5, "BNB other".
    conds = [
        cosmic,                                                     # 7
        nu & ~av,                                                   # 6 dirt
        av & is_ncdelta_rad(df),                                    # 0 signal
        av & nc & (npi0 == 1) & coh,                                # 1
        av & nc & (npi0 == 1) & ~coh,                               # 2
        av & cc & (df.true_pdg == 14) & (npi0 == 1),  # 3; numubar -> other, as in [3]
        av & cc & (df.true_pdg.abs() == 12),                        # 4
    ]
    codes = [7, 6, 0, 1, 2, 3, 4]
    return pd.Series(np.select(conds, codes, default=5), index=df.index)


def is_signal(df):
    """True NC Delta -> N gamma with vertex in the active volume."""
    return is_ncdelta_rad(df) & true_in_active(df)
