"""Collection-plane end-point topology ("coltopo") score of the muon candidate.

A per-detector LightGBM (models in coltopo_models/) on six features of the
muon candidate's collection-plane calo points within its last 25 cm of
residual range (the make_maple_evt_df muon end hits, see _mu_endhits):
  step_med   median 3D distance between consecutive points
  ang2_med   median angle between the 2 cm chords before and after each
             point of a 1 cm arc-length (rr) grid
  wire_span  |signed wire displacement| over the segment (pitch-projected)
  rr_span    rr range covered by the points
  ang1_med   as ang2_med with 1 cm chords
  cos_pitch  wire_span * pitch / rr_span
The output is P(gamma/e shower) vs stopping muon: high = shower-like.
Training (trk-points/mcd_bdt_pfail.py --sep --loo, set l_a4rb = "wd6"):
contained mu / gamma-e candidates failing the proton chi2 cut, even entries,
SBND and ICARUS (Run2+Run4) separately. Features are a port of
trk-points/r2_features.feats.

coltopo_trkscore_a2shift is the same score with ang2_med scaled by
(1 + A2_SHIFT[det]) -- the data/MC ang2_med scale difference fitted at the
PID stage (with the track-score cut removed) -- for a MC shift systematic.
It is NaN on data.
"""
import os

import numpy as np
import pandas as pd

COLTOPO_FEATURES = ["step_med", "ang2_med", "wire_span", "rr_span", "ang1_med", "cos_pitch"]
A2_SHIFT = {"SBND": 0.028, "ICARUS": 0.047}
MODEL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "coltopo_models")
PITCH = 0.3  # cm, collection wire pitch (SBND and ICARUS)
NGRID = 25  # 1 cm arc-length grid, s = 0.5 .. 24.5
_BOOSTERS = {}


def _booster(DETECTOR):
    if DETECTOR not in _BOOSTERS:
        import lightgbm as lgb
        _BOOSTERS[DETECTOR] = lgb.Booster(model_file=os.path.join(MODEL_DIR, "coltopo_wd6_%s.txt" % DETECTOR))
    return _BOOSTERS[DETECTOR]


def _pitch_dir(DETECTOR, x):
    # unit pitch direction (y, z components) of the collection plane. SBND:
    # vertical wires -> along z. ICARUS: +-60 deg wires -> -+30 deg from z,
    # alternating by TPC.
    if DETECTOR == "SBND":
        a = np.zeros_like(x)
    else:
        a = np.where((x < -210) | ((x > 0) & (x < 210)), -30.0, 30.0)
    r = np.radians(a)
    return np.sin(r), np.cos(r)


def _angle(a, b):
    na = np.linalg.norm(a, axis=-1); nb = np.linalg.norm(b, axis=-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = (a * b).sum(-1) / (na * nb)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def _nanmedian(a):
    with np.errstate(invalid="ignore"), np.testing.suppress_warnings() as sup:
        sup.filter(RuntimeWarning)
        return np.nanmedian(a, axis=1)


def coltopo_features(hits, DETECTOR):
    """The six features per slice from the muon end hits (columns x, y, z,
    rr; index (entry, slc, hit)). Returns a DataFrame indexed (entry, slc)."""
    lvl = list(hits.index.names[:2])
    mh = hits[["x", "y", "z", "rr"]].reset_index(level=lvl).reset_index(drop=True)
    mh = mh.sort_values(lvl + ["rr"], kind="mergesort").reset_index(drop=True)
    gcode = mh.groupby(lvl, sort=False).ngroup().to_numpy()
    keys = mh.loc[np.r_[True, gcode[1:] != gcode[:-1]], lvl]
    g = gcode
    ng = len(keys)
    x, y, z, rr = (mh[c].to_numpy(float) for c in ["x", "y", "z", "rr"])
    P = np.stack([x, y, z], 1)
    N = np.bincount(g, minlength=ng).astype(float)
    starts = np.r_[0, np.cumsum(N[:-1]).astype(int)]
    same = np.r_[False, g[1:] == g[:-1]]
    out = {}
    span = np.fmax.reduceat(rr, starts) - np.fmin.reduceat(rr, starts)
    out["rr_span"] = span

    # ---- pitch-projected steps ----
    sy, sz = _pitch_dir(DETECTOR, x)
    proj = y * sy + z * sz
    w = np.where(same, np.r_[np.nan, np.diff(proj)], np.nan) / PITCH
    W = np.bincount(g, np.nan_to_num(w), minlength=ng)  # total signed wire displacement
    with np.errstate(invalid="ignore", divide="ignore"):
        out["wire_span"] = np.abs(W)
        out["cos_pitch"] = np.abs(W) * PITCH / span

    # ---- 3D steps ----
    st = np.where(same, np.r_[np.nan, np.linalg.norm(np.diff(P, axis=0), axis=1)], np.nan)
    out["step_med"] = pd.Series(st).groupby(g).quantile(0.5).reindex(range(ng)).to_numpy()

    # ---- arc-length resampling on a 1 cm grid ----
    sg = np.arange(NGRID) + 0.5
    key = g * 100.0 + rr
    qk = (np.arange(ng)[:, None] * 100.0 + sg[None, :]).ravel()
    idx = np.searchsorted(key, qk)
    gq = np.repeat(np.arange(ng), NGRID)
    valid = (idx > 0) & (idx < len(key))
    idx_c = np.clip(idx, 1, len(key) - 1)
    valid &= (g[idx_c] == gq) & (g[idx_c - 1] == gq)
    r0, r1 = rr[idx_c - 1], rr[idx_c]
    with np.errstate(invalid="ignore", divide="ignore"):
        t = (qk - gq * 100.0 - r0) / (r1 - r0)
        Pg = P[idx_c - 1] + t[:, None] * (P[idx_c] - P[idx_c - 1])
    Pg[~valid] = np.nan
    Pg = Pg.reshape(ng, NGRID, 3)
    for L in [1, 2]:
        a = Pg[:, L:NGRID - L] - Pg[:, :NGRID - 2 * L]
        b = Pg[:, 2 * L:] - Pg[:, L:NGRID - L]
        out["ang%i_med" % L] = _nanmedian(_angle(a, b))
    return pd.DataFrame(out, index=pd.MultiIndex.from_frame(keys))[COLTOPO_FEATURES]


def coltopo_scores(hits, DETECTOR, ismc):
    """coltopo_trkscore (and the ang2_med-shifted coltopo_trkscore_a2shift,
    MC only) per slice, indexed (entry, slc), from the muon end hits."""
    cols = ["coltopo_trkscore", "coltopo_trkscore_a2shift"]
    if len(hits) == 0:
        return pd.DataFrame({c: pd.Series(dtype=float) for c in cols}, index=hits.index.droplevel(-1)[:0])
    F = coltopo_features(hits, DETECTOR)
    X = F.to_numpy(np.float64)
    allnan = np.isnan(X).all(axis=1)
    bst = _booster(DETECTOR)
    s = bst.predict(X, num_threads=1)
    s[allnan] = np.nan
    s_sh = np.full(len(X), np.nan)
    if ismc:
        X[:, COLTOPO_FEATURES.index("ang2_med")] *= 1 + A2_SHIFT[DETECTOR]
        s_sh = bst.predict(X, num_threads=1)
        s_sh[allnan] = np.nan
    return pd.DataFrame({cols[0]: s, cols[1]: s_sh}, index=F.index)
