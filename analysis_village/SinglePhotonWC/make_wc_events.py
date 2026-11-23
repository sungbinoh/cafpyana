"""One row per CAF event, index entry: Wire-Cell bookkeeping, so events
without a candidate (no row in the wc table) can be accounted for.

  valid        1: Wire-Cell output exists for the event. 0: none (a step-2
               job that failed every try, a partial merge, or an event lost
               upstream); a production failure, not a reconstruction outcome
  n_nu         candidates (T_tagger rows)
  nu_n_flashes, nu_n_flashes_in_window, nu_n_bundles   from Trun
  n_bundle_reason_<r>   T_bundle rows (flash-cluster bundles) per reason:
               0 candidate (main), 1 candidate (demoted main), 2 all
               cosmic-vetoed, 3 below the 15 cm length floor, 4 STM-only,
               5 nothing eligible, 6 flash-group dedup

valid = 1 with n_nu = 0 is a processed event Wire-Cell found no candidate in;
the reason counts say which step dropped its bundles (none: no bundle).
Joins to mcnu on (file, __ntuple, entry).

make_hdrdf: the rec.hdr header per event (pot, noffbeambnb, ...) for the
data POT / off-beam gate normalisation.

make_wc_chanrun_df: long single-channel signals from wc.T_proj_data (the
deconvolved charge per channel and 2 us time slice that Wire-Cell's imaging
used, with the charge its 3D solution predicts there), for the AC-coupling
/ signal-processing checks: one row per run of >= CHANRUN_MIN_SLICES
consecutive slices of one cluster on one channel. Columns: cluster_id,
nu_index (the candidate whose main cluster it is, -1 none), channel, tpc,
plane (SBND channel map, CHANNEL_BLOCKS), ts_start, n_slices, q / qpred
(summed measured / predicted), and two profiles: q_rel_<i> / qp_rel_<i>,
the run cut into CHANRUN_NREL equal parts, and q_abs_<i> / qp_abs_<i>, the
charge in CHANRUN_ABS_SLICES-slice blocks from the start of the run (NaN past
its end). Measured / predicted vs time within a run is the in-signal sag;
the stored charge is >= 0 and ROI-limited, so the after-signal undershoot
is not visible here.
"""

import numpy as np
import pandas as pd

from makedf import wirecell as wc
from pyanalib.pandas_helpers import loadbranches

BUNDLE_REASONS = range(7)
RUN_COLUMNS = ["nu_n_flashes", "nu_n_flashes_in_window", "nu_n_bundles"]


def make_wc_event_df(f):
    ev = wc.make_wc_evtdf(f)[["run", "subrun", "evt", "valid", "n_nu"]]
    run = wc.make_wc_rundf(f).reset_index("row", drop=True)[RUN_COLUMNS]
    reason = wc.make_wc_bundledf(f, ["reason"]).reason
    n = (pd.crosstab(reason.index.get_level_values("entry"), reason)
         .reindex(index=ev.index, columns=BUNDLE_REASONS, fill_value=0)
         .add_prefix("n_bundle_reason_"))
    n.columns.name = None
    # valid = 0 events have no Trun row: their counts are 0
    return ev.join(run.reindex(ev.index).fillna(0).astype(int)).join(n)


def make_hdrdf(f):
    """rec.hdr pot and noffbeambnb (as makedf.make_hdrdf, without its CVMFS imports)."""
    return loadbranches(f["recTree"], ["rec.hdr.run", "rec.hdr.subrun", "rec.hdr.evt", "rec.hdr.pot",
                                       "rec.hdr.noffbeambnb", "rec.hdr.nbnbinfo", "rec.hdr.first_in_subrun"]).rec.hdr


# (first channel, plane, tpc) per block, from the Pandora calo points' channel vs (wire, plane, tpc): per TPC
# 1984 U + 1984 V + 1670 Y channels
CHANNEL_BLOCKS = [(0, 0, 0), (1984, 1, 0), (3968, 2, 0), (5638, 0, 1), (7622, 1, 1), (9606, 2, 1), (11276, -1, -1)]
CHANRUN_MIN_SLICES = 10      # 2 us slices: >= 20 us of signal on one channel
CHANRUN_NREL = 10
CHANRUN_ABS_SLICES = 5       # 10 us blocks
CHANRUN_NABS = 20            # up to 200 us


def channel_plane_tpc(ch):
    ch = np.asarray(ch)
    starts = np.array([b[0] for b in CHANNEL_BLOCKS])
    i = np.clip(np.searchsorted(starts, ch, side="right") - 1, 0, len(CHANNEL_BLOCKS) - 1)
    return np.array([CHANNEL_BLOCKS[k][1] for k in i]), np.array([CHANNEL_BLOCKS[k][2] for k in i])


def _chanrun_cols():
    return (["cluster_id", "nu_index", "channel", "tpc", "plane", "ts_start", "n_slices", "q", "qpred"]
            + ["q_rel_%d" % i for i in range(CHANRUN_NREL)] + ["qp_rel_%d" % i for i in range(CHANRUN_NREL)]
            + ["q_abs_%d" % i for i in range(CHANRUN_NABS)] + ["qp_abs_%d" % i for i in range(CHANRUN_NABS)])


def _empty_chanruns():
    return pd.DataFrame(columns=_chanrun_cols(), index=pd.MultiIndex.from_arrays([[]] * 2, names=["entry", "run"]),
                        dtype=np.float32)


def make_wc_chanrun_df(f):
    """Never raises: an unreadable CAF gives an empty table (with a warning)."""
    try:
        return _chanruns(f)
    except Exception as e:
        print("make_wc_chanrun_df: skipped (%s: %s)" % (type(e).__name__, e))
        return _empty_chanruns()


def _chanruns(f):
    if not wc.has_wc(f):
        return _empty_chanruns()
    t = f[wc.WC_TREE]
    pre = "wc.T_proj_data."
    if pre + "charge.v" not in t.keys():
        return _empty_chanruns()
    a = t.arrays([pre + "cluster_id", pre + "channel.v..length", pre + "channel.v", pre + "time_slice.v",
                  pre + "charge.v", pre + "charge_pred.v"], library="np")
    nu = wc.make_wc_nudf(f, tagger_leaves=["cluster_id"], include_tmatch=False) if t.num_entries else None
    main = {}
    if nu is not None and len(nu) and "cluster_id" in nu:
        for (e, n), c in nu.cluster_id.items():
            main[(e, int(c))] = n
    rows = []
    for e in range(t.num_entries):
        L = np.asarray(a[pre + "channel.v..length"][e], int)
        if not L.sum():
            continue
        cl = np.repeat(np.asarray(a[pre + "cluster_id"][e], int), L)
        ch = np.asarray(a[pre + "channel.v"][e], int)
        ts = np.asarray(a[pre + "time_slice.v"][e], int)
        q = np.asarray(a[pre + "charge.v"][e], float)
        qp = np.asarray(a[pre + "charge_pred.v"][e], float)
        o = np.lexsort((ts, ch, cl))
        cl, ch, ts, q, qp = cl[o], ch[o], ts[o], q[o], qp[o]
        new = np.r_[True, (cl[1:] != cl[:-1]) | (ch[1:] != ch[:-1]) | (ts[1:] != ts[:-1] + 1)]
        starts = np.flatnonzero(new)
        ends = np.r_[starts[1:], len(ts)]
        for s0, s1 in zip(starts, ends):
            n = s1 - s0
            if n < CHANRUN_MIN_SLICES:
                continue
            qq, pp = q[s0:s1], qp[s0:s1]
            rel = np.array_split(np.arange(n), CHANRUN_NREL)
            qa, pa = np.full(CHANRUN_NABS, np.nan), np.full(CHANRUN_NABS, np.nan)
            for k in range(min(CHANRUN_NABS, n // CHANRUN_ABS_SLICES)):
                sl = slice(k * CHANRUN_ABS_SLICES, (k + 1) * CHANRUN_ABS_SLICES)
                qa[k], pa[k] = qq[sl].sum(), pp[sl].sum()
            plane, tpc = channel_plane_tpc([ch[s0]])
            rows.append([e, cl[s0], main.get((e, int(cl[s0])), -1), ch[s0], tpc[0], plane[0], ts[s0], n, qq.sum(), pp.sum()]
                        + [qq[r].sum() for r in rel] + [pp[r].sum() for r in rel] + list(qa) + list(pa))
    if not rows:
        return _empty_chanruns()
    df = pd.DataFrame(rows, columns=["entry"] + _chanrun_cols())
    df["run"] = df.groupby("entry").cumcount()
    return df.set_index(["entry", "run"]).astype(np.float32)
