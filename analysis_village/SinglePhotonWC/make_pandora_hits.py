"""Pandora track hits per plane, summarised per track, for the dQ/dx
calibration checks (plot_dqdx_data_mc.py). Pandora's hits come from the same
signal-processing output (simtpc2d / sptpc2d : dnnsp) as Wire-Cell's
imaging, so their widths, multiplicities and integral / ADC-sum ratios probe
the signal processing independently of the Wire-Cell reconstruction.

One row per (entry, slice, pfp, plane) for Pandora tracks >= TRK_MIN_LEN:
  slice / pfp: is_clear_cosmic, nu_score, track_score, len, dir_x/y/z,
               start_x/y/z, end_x/y/z, chi2_mu, chi2_p (collection plane),
               true_pdg (MC; NaN on data)
  hits:        n_hits; medians over all hits of dqdx (the CAF's calibrated
               dQ/dx), integral, sumadc, int_over_sum (integral / sumadc per
               hit), width (ticks), pitch (cm), phi (deg, the hit's angle to
               the drift field), x; mean mult (hits fitted together);
               dqdx_mip, width_mip, mult_mip, int_over_sum_mip over the hits
               with residual range > RR_MIP_MIN (no Bragg rise)

Reco only from recTree via pyanalib.loadbranches, so it runs on MC and data
(no CVMFS imports).
"""

import numpy as np
import pandas as pd

from pyanalib.pandas_helpers import loadbranches

TRK = "rec.slc.reco.pfp.trk."
TRK_MIN_LEN = 30.0   # cm
RR_MIP_MIN = 30.0    # cm
HIT_LEAVES = ["dqdx", "integral", "sumadc", "width", "mult", "pitch", "phi", "rr", "x"]
TRK_LEAVES = {"len": "len", "dir.x": "dir_x", "dir.y": "dir_y", "dir.z": "dir_z",
              "start.x": "start_x", "start.y": "start_y", "start.z": "start_z",
              "end.x": "end_x", "end.y": "end_y", "end.z": "end_z",
              "chi2pid.2.chi2_muon": "chi2_mu", "chi2pid.2.chi2_proton": "chi2_p"}
LEVELS = ["entry", "rec.slc..index", "rec.slc.reco.pfp..index"]


def _tracks(tree):
    keys = set(tree.keys())
    t = loadbranches(tree, [TRK + k for k in TRK_LEAVES] + ["rec.slc.reco.pfp.trackScore"])
    t.columns = [".".join(c for c in col if c) for col in t.columns]
    # loadbranches names an array index n as "In" (chi2pid.2 -> chi2pid.I2)
    t = t.rename(columns={TRK + k.replace("chi2pid.2", "chi2pid.I2"): v for k, v in TRK_LEAVES.items()}
                 | {"rec.slc.reco.pfp.trackScore": "track_score"})
    ismc = bool(loadbranches(tree, ["rec.hdr.ismc"]).rec.hdr.ismc.iloc[0]) if tree.num_entries else False
    if ismc and TRK + "truth.p.pdg" in keys:
        t = t.join(loadbranches(tree, [TRK + "truth.p.pdg"]).rec.slc.reco.pfp.trk.truth.p.pdg.rename("true_pdg"))
    s = loadbranches(tree, ["rec.slc.is_clear_cosmic", "rec.slc.nu_score"]).rec.slc
    t = t.join(s, on=["entry", "rec.slc..index"])
    return t[t.len >= TRK_MIN_LEN]


def _plane_summary(tree, plane, keep):
    leaves = [TRK + "calo.%d.points.%s" % (plane, h) for h in HIT_LEAVES]
    h = loadbranches(tree, leaves)
    h.columns = [col[-1] for col in h.columns]
    h = h[h.index.droplevel(-1).isin(keep)]
    if not len(h):
        return pd.DataFrame()
    h["int_over_sum"] = h.integral / h.sumadc.replace(0, np.nan)
    g = h.groupby(level=LEVELS)
    out = pd.DataFrame({
        "n_hits": g.size(),
        "dqdx": g.dqdx.median(), "integral": g.integral.median(), "sumadc": g.sumadc.median(),
        "int_over_sum": g.int_over_sum.median(), "width": g.width.median(), "mult": g.mult.mean(),
        "pitch": g.pitch.median(), "phi": np.degrees(g.phi.median()), "x": g.x.median(),
    })
    m = h[h.rr > RR_MIP_MIN].groupby(level=LEVELS)
    out = out.join(pd.DataFrame({"dqdx_mip": m.dqdx.median(), "width_mip": m.width.median(),
                                 "mult_mip": m.mult.mean(), "int_over_sum_mip": m.int_over_sum.median()}))
    out["plane"] = plane
    return out.set_index("plane", append=True)


def _empty():
    """No rows, but the index names (an unnamed empty frame strips them from the whole split)."""
    return pd.DataFrame(index=pd.MultiIndex.from_arrays([[]] * 4, names=LEVELS + ["plane"]))


def make_pandora_trkhit_df(f):
    """Never raises: a CAF it cannot read gives an empty table (with a
    warning) rather than failing the whole dataframe job."""
    try:
        return _make(f)
    except Exception as e:
        print("make_pandora_trkhit_df: skipped (%s: %s)" % (type(e).__name__, e))
        return _empty()


def _make(f):
    if "recTree" not in f or not f["recTree"].num_entries:
        return _empty()
    tree = f["recTree"]
    t = _tracks(tree)
    if not len(t):
        return _empty()
    planes = [p for p in (_plane_summary(tree, pl, t.index) for pl in range(3)) if len(p)]
    if not planes:
        return _empty()
    hits = pd.concat(planes)
    out = hits.join(t, on=LEVELS)
    return out.astype(np.float32)
