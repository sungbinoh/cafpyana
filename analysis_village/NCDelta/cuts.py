"""Pandora NC Delta -> N gamma selection for SBND, adapted from MicroBooNE.

References
----------
[1] MicroBooNE, arXiv:2502.05750. Its Pandora channels are unchanged from [2].
[2] MicroBooNE, PRL 128, 111801 (2022), arXiv:2110.00409 (+ supplemental).
[3] Numeric pre-cuts: final Runs 1-3 Hive configs,
    github.com/markrosslonergan/hellstroms_hive, hive/xml/
      1g1p_xmls/hive2.5/1g1p_CinnabarWorldOrder_v33.3_15_04_30_Combined_Run123_v1_updatedBDTfiles.xml
      1g0p_xmls/hive2.5/1g0p_CinnabarWorldOrder_v33.3_4_29_20_CombinedRuns13_updatedBDTfiles.xml

Stages:
  1. Topology: one shower and one (1g1p) or zero (1g0p) primary tracks.
  2. Pre-selection (PRESEL_*): the fiducial volume only. Stages 1+2 are what
     the dataframe maker applies, so keeping stage 2 minimal means the
     dataframes do not bake in a tuning choice.
  3. Background rejection: rectangular cuts on CAF features (CUTBASED_*).
  4. Signal window in reconstructed shower energy (ENERGY_WINDOW).

Energies are GeV and lengths cm, except dE/dx and PIDA which keep the
conventional MeV/cm. Shower energy is shw_energy_corrected(): the
best plane's energy (best = most hits, as sbncode FillReco.cxx defines it for
SBND) with the per-plane scale SHW_ENERGY_SCALE applied. The MicroBooNE
"#GE#" convention -- collection plane with a max-plane fallback -- is not
used: SBND's wire planes are not MicroBooNE's, and post-selection it resolves
the photon energy slightly worse.

Tuning
------
CUTBASED_* was tuned on 500 nopresel dataframes (1.8e20 POT). Every change had
to hold up on the even and the odd file split separately and on a true
holdout. Signal is is_signal(): true vertex in the *active* volume, so
efficiencies here and from evaluate.py are AV-based; [1] uses the FV, which is
a factor ~1.33 larger (15519 true signal neutrinos instead of 21450).

Two working points, quoted separately rather than summed (see SELECTIONS).
FV-based efficiency / purity, against [1] Table I:

                this file              [1] Pandora
  1g1p   2.27% / 11.1% (balanced)   4.24% / 14.8%
         0.69% / 20.9% (high_purity)
  1g0p   3.24% /  8.2% (balanced)   5.52% /  4.0%
         2.07% / 10.5% (high_purity)

Kept despite the figures of merit preferring otherwise, because an MC-only
figure of merit cannot see what they cost: the fiducial volume (keeps badly
modelled regions out) and trk_proton_like (dropping it inflates the 1g1p
CC numu 1pi0 background 8x). Conversely opt0_score_min was dropped despite
being worth ~0.2 points of purity, for the same reason in reverse -- the score
is not well validated and costs real MC/data comparison work.

Where the ceiling is
--------------------
Topology plus the FV caps efficiency at 10.7% (1g1p) and 6.9% (1g0p) at ~0.5%
purity, so everything happens in stage 3 -- and that is where the limit bites.
It is the input variables, not the cut shape:

  - the best single variable rejects only ~20% of background at 90% signal
    efficiency;
  - a 17-variable naive-Bayes likelihood ratio, trained on one split and
    evaluated on the other, lands on the same efficiency-purity frontier;
  - against NC pi0 (80% of the 1g1p background, 39% of the 1g0p one) the best
    variable rejects 18% -- correctly, since a pi0 that loses one photon
    really is a single photon in the detector.

Closing that needs a second-shower / unassociated-hit variable, which [2] used
and these CAFs lack (n_nonprim_shw is identically 0 for 1g0p and counts only
the leading track's daughters for 1g1p). That is a make_ncdelta.py change plus
a dataframe remake.
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

# Pre-selection. Applied at dataframe creation, so it is deliberately limited
# to the fiducial volume and Pandora's clear-cosmic flag: anything tighter
# bakes a tuning choice into the dataframes and cannot be loosened afterwards.
# Every other cut lives in CUTBASED_* below, re-tunable without a remake.
PRESEL_1G1P = {
    "not_clear_cosmic": True,   # Pandora's neutrino candidate only, as in [3]
    "vtx_fv": True,             # FV on the slice vertex (the track fixes it)
}

PRESEL_1G0P = {
    "not_clear_cosmic": True,
    "shw_fv": True,             # FV on the shower start (no track defines the vertex)
}

# Signal window in shw_energy_corrected(). Common to both working points in
# SELECTIONS: the window defines the signal region of the measurement, so
# narrowing it per selection would inflate purity by discarding signal phase
# space rather than by rejecting background.
# GeV, on shw_energy_corrected().
ENERGY_WINDOW = {
    "1g1p": (0.020, 0.700),
    "1g0p": (0.120, 0.950),
}

TOPOLOGIES = ("1g1p", "1g0p")

# ---------------------------------------------------------------------------
# Geometry (SBND only)
# ---------------------------------------------------------------------------

# Active volume: the TPC box, no insets.
ACTIVE_VOLUME = {"x": (-200.0, 200.0), "y": (-200.0, 200.0), "z": (0.0, 500.0)}

# Fiducial volume: the SBND part of the FV in analysis_village/mcs_pid/cuts.py,
# i.e. the active volume with insets and the high-yz region excluded.
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


# ---------------------------------------------------------------------------
# Shower energy
# ---------------------------------------------------------------------------

# Per-plane energy scale, applied to the best plane's energy. SBND's best
# plane for energy is the one with the most hits (sbncode FillReco.cxx), and
# which plane that is varies event to event, so all three factors are in play.
#
# Derived from selected signal showers truth-matched to a photon, as
# 1 / median(E_plane / E_true_photon) over the events where that plane is the
# best plane. Stable to +-0.01 across both working points and both file
# splits. Applying them centres the energy scale within 0.5% of truth
# out-of-sample; the fractional resolution (~13%) is unchanged, since a
# multiplicative factor cannot change it.
#
# Caveat: the factors are derived on selected events, which is also where they
# are applied, so there is mild circularity -- the inclusive bias is larger
# (median 0.74 rather than 0.85), the selection keeping better-reconstructed
# showers. Re-derive if the selection changes substantially.
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


def extrap_det_path(df):
    """Extrapolated detector path length [cm].

    Distance from the shower start, projected *backwards* along the shower
    direction, to where it leaves the active volume: how far the photon could
    have travelled inside the detector before converting.

    Small for dirt -- a neutrino interacting in the surrounding rock sends its
    photon in through a nearby wall, so there is little detector behind the
    conversion point. Large for a genuine in-detector interaction. This is the
    handle the fiducial-volume cut cannot provide, because with no track the
    FV is applied to the shower start, which for dirt is legitimately inside
    (median 23 cm in from the wall, while the true vertex is ~29 cm outside).
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
    out = OrderedDict()
    if cfg.get("not_clear_cosmic"):
        out["not_clear_cosmic"] = ~df.is_clear_cosmic.astype(bool)
    if cfg.get("vtx_fv"):
        out["vtx_fv"] = _fv(df, "slc_vtx")
    return out


def presel_1g0p_cuts(df, cfg=PRESEL_1G0P):
    """Ordered dict of the individual 1g0p pre-cuts (topology excluded)."""
    out = OrderedDict()
    if cfg.get("not_clear_cosmic"):
        out["not_clear_cosmic"] = ~df.is_clear_cosmic.astype(bool)
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

# ---------------------------------------------------------------------------
# Cut-based background rejection
# ---------------------------------------------------------------------------
# Cuts are grouped by the background they target. Unless a source is given,
# values are untuned physics-motivated starting points: tune them on SBND MC
# with nminus1_masks(). Set a value to None to disable.

CUTBASED_1G1P = {
    "nu_score_min":            0.35,
    "bary_score_min":          1e-6,   # mcs_pid baryscore_cut
    # OpT0 flash match. Dropped from every working point despite being worth
    # ~0.2 points of purity on both splits: the OpT0 score is not well
    # validated and carries a real MC/data comparison burden, which an MC-only
    # figure of merit cannot see. It was also doing different things per
    # channel -- in 1g0p it removed cosmics and dirt (75 of 124 slices, no
    # signal), in 1g1p only 4 NC pi0 slices and no cosmics at all. For the
    # record, higher score = better match: at 1g0p topology level cosmics have
    # median score 10.6 with 35% having no match, signal median 656 with 3%.
    "opt0_score_min":          None,
    "shw_dedx_min":            3.5,    # MeV/cm; photon ~2 MIP (peak 4.2), electron ~1
    "shw_conv_dist_max":       80.0,   # X0 = 14 cm, so 80 cm is misreco
    "trk_proton_like":         True,   # chi2_proton < chi2_muon
    "energy_window":           None,   # None = ENERGY_WINDOW[topology]

    # Switched off: a leave-one-out test showed each was worth <0.3 points of
    # purity once the cuts above are applied (several removed nothing at all).
    # Kept as None to document what is available. Signal/background slices
    # each removed, out of 350/2734 selected:
    "vtx_dist_min":            None,   # 2.0:    0 / 0
    "trk_dist_min":            None,   # 2.0:    0 / 2
    "shw_start_dist_min":      None,   # 2.0:    0 / 6
    "trk_ke_max":              None,   # 0.7 GeV: 0 / 6
    "n_nonprim_shw_max":       None,   # 1:      0 / 3
    "shw_nhits_min":           None,   # 150:    0 / 1
    "trk_start_dist_to_vtx_max": None, # 10.0:   0 / 0
    "abs_cos_trk_shw_max":     None,   # 0.99:   1 / 7
    # On in high_purity only:
    "trk_calo_E_max":          None,   # GeV
    "trk_pida_max":            None,
    "shw_len_min":             None,
    "shw_open_angle_max":      None,
    "shw_conv_dist_min":       None,
    # Pandora's track/shower score for the leading shower. TRK_SCORE_CUT = 0.5
    # already bounds it from above by construction, but genuine photon showers
    # sit near 0.39 while CC numu events whose muon was filed as shower-like
    # pile up at ~0.47 -- just under the boundary, which is how they reach
    # 1g0p with no track to cut on. Tightening it targets that background.
    "shw_trackScore_max":      None,
    "extrap_det_path_min":     None,   # 1g0p only: dirt is just 1-3% of this
                                       # channel, the track already vetoes it
    # p+gamma mass window [GeV]: the peak is at ~1.143, below the PDG
    # Delta(1232), and NC pi0 peaks at 1.108, so a window buys ~0.2% purity.
    "delta_mass_range":        None,
}

CUTBASED_1G0P = {
    "nu_score_min":       0.35,
    "bary_score_min":     0.02,
    # OpT0 flash match. Dropped from every working point despite being worth
    # ~0.2 points of purity on both splits: the OpT0 score is not well
    # validated and carries a real MC/data comparison burden, which an MC-only
    # figure of merit cannot see. It was also doing different things per
    # channel -- in 1g0p it removed cosmics and dirt (75 of 124 slices, no
    # signal), in 1g1p only 4 NC pi0 slices and no cosmics at all. For the
    # record, higher score = better match: at 1g0p topology level cosmics have
    # median score 10.6 with 35% having no match, signal median 656 with 3%.
    "opt0_score_min":     None,
    "shw_dedx_min":        2.0,  # looser than 1g1p: shw_dedx_best is measured
                                 # badly enough that tighter costs more signal
                                 # than the CC nue it removes
    "energy_window":      None,  # None = ENERGY_WINDOW[topology]

    # Switched off (see the 1g1p note); removed out of 525/5941 selected:
    "shw_start_dist_min": None,  # 2.0:  0 / 0
    "n_nonprim_shw_max":  None,  # 0: structurally dead, n_nonprim_shw is
                                 # identically 0 for every 1g0p slice
    "shw_nhits_min":      None,  # 100:  5 / 131
    "shw_conv_dist_max":  None,  # 5.0:  1 / 65
    "shw_len_min":        None,  # on in high_purity only
    "shw_open_angle_max": None,  # on in high_purity only
    # Dirt rejection, 1g0p only. With no track the FV is applied to the shower
    # start, which for a dirt event is legitimately inside the detector (median
    # 23 cm in) even though the true vertex is ~29 cm outside -- the photon
    # travels in and converts. This cut asks how much detector lies behind the
    # conversion point instead, and roughly halves the dirt background.
    # Pandora's track/shower score for the leading shower. TRK_SCORE_CUT = 0.5
    # already bounds it from above by construction, but genuine photon showers
    # sit near 0.39 while CC numu events whose muon was filed as shower-like
    # pile up at ~0.47 -- just under the boundary, which is how they reach
    # 1g0p with no track to cut on. Tightening it targets that background.
    "shw_trackScore_max": 0.48,
    "extrap_det_path_min": 40.0,
}

CUTBASED = {"1g1p": CUTBASED_1G1P, "1g0p": CUTBASED_1G0P}


# The defaults above are "balanced": they maximize S/sqrt(S+B), which is the
# same as maximizing purity x efficiency, since (S/sqrt(S+B))^2 = S * purity
# and purity x efficiency = S * purity / N_true with N_true fixed. Verified by
# optimizing each separately and getting identical cuts in both channels.
# "high_purity" maximizes purity, equivalently S/B, since
# purity = (S/B) / (1 + S/B).
#
# Use selection_config(name, topo), dump one to JSON for evaluate.py --config,
# or pass --selection to plot_composition.py.
SELECTIONS = {
    "balanced": {"1g1p": {}, "1g0p": {}},
    "high_purity": {
        # Same variables as the original tuning -- the re-convergence after the
        # energy change had added shw_nhits_min / n_nonprim_shw_max /
        # trk_start_dist_to_vtx_max here and shw_conv_dist_max / shw_len_min /
        # shw_nhits_min / shw_start_dist_min to 1g0p, which is a change of
        # variables rather than of values. Reverted: it costs 0.12 points of
        # purity in 1g1p, and in 1g0p trades 2.3 points for 29% more signal.
        "1g1p": {"nu_score_min": 0.50, "bary_score_min": 0.02,
                 "shw_conv_dist_min": 4.0, "shw_conv_dist_max": 50.0,
                 "shw_len_min": 30.0, "shw_open_angle_max": 0.25,
                 "shw_trackScore_max": 0.42, "trk_calo_E_max": 0.300},
        "1g0p": {"nu_score_min": 0.50, "shw_dedx_min": 2.5,
                 "shw_open_angle_max": 0.20, "shw_trackScore_max": 0.44,
                 "extrap_det_path_min": 80.0},
    },
}


def selection_config(name, topology):
    """CUTBASED[topology] with the named working point's overrides applied."""
    if name not in SELECTIONS:
        raise KeyError("Unknown selection %r; have %s" % (name, sorted(SELECTIONS)))
    cfg = dict(CUTBASED[topology])
    cfg.update(SELECTIONS[name][topology])
    return cfg


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

# CAF "no value" markers, left in the dataframe and masked here instead. A cut
# on NaN is False either way, which is what a missing measurement should do.
# tmatch_idx and true_pdg are excluded: is_cosmic() reads their sentinels.
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

# genie::EResonance, from rec.mc.nu.resnum. genie_mode == GENIE_MODE_RES is ANY
# resonance, so resnum is what picks out the Delta(1232).
# Ref: GENIE Framework/ParticleData/BaryonResonance.h, sbnanaobj
# StandardRecord/SRTrueInteraction.h.
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

    genie_mode == kRes is any resonance, so resnum pins the Delta(1232).
    GENIE primaries include ~3 MeV de-excitation gammas, hence the hard-photon
    requirement. npi0 == 0 drops Delta -> N pi0 with a photon from elsewhere,
    which is indistinguishable from signal and overlaps the NC 1pi0 category.
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
