"""One row per Wire-Cell neutrino candidate, index (entry, nu_index).

Everything Wire-Cell writes about a candidate in one row, for selection work:

  * event: run/subrun/evt, candidates in the event (n_nu), Trun flash and
    bundle counts
  * every scalar T_tagger leaf (~1030: the summary flags and scores plus
    the inputs of Wire-Cell's BDTs, shw_sp_*, pio_*, numu_cc_*, cosmict_*,
    the nue-tagger blocks, ssm_*), under their T_tagger names; T_kine
    scalars (energies in MeV, as in the CAF); on MC the T_match level-3
    truth match
  * main cluster (T_cluster row of T_tagger.cluster_id): main_cluster_*
  * flash bundle (T_bundle row with this nu_index): bundle_*; bundle_reason
    0 main candidate, 1 demoted main (codes in make_wc_events.py)
  * particle list (T_kine vectors), summarised over the particles connected
    to the candidate (kine_energy_included == 1): counts per PID label and
    leading energies. The unclear-link leftover showers (== 3) are counted
    apart, n_unclear / E_unclear; kine_reco_Enu includes both.
  * segments (T_segment) and charge points (T_rec_charge): counts, lengths
  * leading EM shower geometry from its charge points (T_rec_charge points
    mapped to particles through T_segment): shw1_start_* (point nearest the
    reco vertex), shw1_dir_* (principal axis of the points within
    SHW_AXIS_RADIUS of the start, pointing away from it), shw1_back_dist
    (from the start back along the axis to the active-volume boundary; short
    for a photon entering from outside), shw1_start_vtx_dist, shw1_npts;
    vtx_other_q (candidate charge within VTX_RADIUS of the reco vertex not
    in that shower) and other_q_frac (candidate charge not in it)
  * MC only: the T_match truth match, tmatch_idx (row in rec.mc.nu, -1: no
    match) and its charge / point IoU. Absent on data.

Reco only: truth is the separate mcnu table (make_truth.py), joined on
(file, __ntuple, entry, tmatch_idx) by cutflow.load. Events without a
candidate have no row here; mcnu gives the efficiency denominators.
"""

import numpy as np
import pandas as pd

from makedf import wirecell as wc
from analysis_village.SinglePhotonWC import cuts

LEVELS = wc.WC_NU_LEVELS
PID = {"em": 11, "muon": 13, "pion": 211, "proton": 2212, "nopid": 0}
# T_tagger identifiers not carried (run/subrun/evt come from the event table)
TAGGER_SKIP = {"run", "subrun", "event", "nu_index", "sel_cluster_id"}
SHW_AXIS_RADIUS = 30.0   # cm, points used for the shower axis
VTX_RADIUS = 5.0         # cm, vertex-activity sphere


def _by_nu(df):
    return df.groupby(level=LEVELS)


def _nth(E, n):
    """n-th largest value per candidate (n = 1 is the maximum)."""
    s = E.sort_values(ascending=False)
    return s[s.groupby(level=LEVELS).cumcount() == n - 1]


def _particles(f, index):
    k = wc.make_wc_kineparticledf(f)
    k = k.reset_index("kine_particle", drop=True)
    E = k.kine_energy_particle
    pdg = k.kine_particle_type.abs()
    conn = k.kine_energy_included == 1

    T = pd.DataFrame(index=index)
    T["n_particles"] = _by_nu(conn).sum()
    T["n_shower"] = _by_nu(conn & (k.kine_particle_is_shower == 1)).sum()
    for name, code in PID.items():
        T["n_" + name] = _by_nu(conn & (pdg == code)).sum()
    T["n_unclear"] = _by_nu(~conn).sum()
    T = T.fillna(0).astype(int)
    T["E_unclear"] = _by_nu(E.where(~conn, 0)).sum().reindex(index).fillna(0)

    for name in ["em", "muon", "pion", "proton"]:
        Ei = E[conn & (pdg == PID[name])]
        T["E_%s_1" % name] = _nth(Ei, 1)
    T["E_em_2"] = _nth(E[conn & (pdg == 11)], 2)
    return T


def _segments(f, index):
    s = wc.make_wc_segmentdf(f)
    s = s[s.nu_index >= 0].reset_index("row", drop=True).set_index("nu_index", append=True)
    trk = s.flag_shower == 0
    T = pd.DataFrame(index=index)
    T["n_segments"] = _by_nu(s).size()
    T["n_segments_in_enu"] = _by_nu(s.in_enu == 1).sum()
    T = T.fillna(0).astype(int)
    T["trk_len_max"] = _by_nu(s.length_cm[trk]).max()
    T["trk_len_sum"] = _by_nu(s.length_cm[trk]).sum().reindex(index).fillna(0)
    T["shw_len_max"] = _by_nu(s.length_cm[~trk]).max()
    return T


def _charge(f, index):
    c = wc.make_wc_chargedf(f, ["nu_index"])
    n = c.reset_index("row", drop=True).set_index("nu_index", append=True).groupby(level=LEVELS).size()
    return pd.DataFrame({"n_charge_points": n}, index=index).fillna(0).astype(int)


def _back_dist(p, u):
    """Distance from p along -u to the active-volume boundary (p inside)."""
    b = cuts.ACTIVE_VOLUME
    t = np.inf
    for i, ax in enumerate("xyz"):
        d = -u[i]
        if d > 1e-9: t = min(t, (b[ax][1] - p[i]) / d)
        elif d < -1e-9: t = min(t, (b[ax][0] - p[i]) / d)
    return max(t, 0.0)


def _shower_geom(f, D):
    cols = ["shw1_start_x", "shw1_start_y", "shw1_start_z", "shw1_dir_x", "shw1_dir_y", "shw1_dir_z",
            "shw1_back_dist", "shw1_start_vtx_dist", "shw1_npts", "vtx_other_q", "other_q_frac"]
    out = pd.DataFrame(np.nan, index=D.index, columns=cols)
    k = wc.make_wc_kineparticledf(f)
    em = k[(k.kine_particle_type.abs() == 11) & (k.kine_particle_is_shower == 1)]
    if not len(em):
        return out
    lead = em.kine_energy_particle.groupby(level=LEVELS).idxmax().map(lambda i: i[-1])
    seg = wc.make_wc_segmentdf(f, ["nu_index", "segment_id", "kine_index"])
    seg = seg.reset_index("row", drop=True).set_index(["nu_index", "segment_id"], append=True).kine_index
    ch = wc.make_wc_chargedf(f, ["x", "y", "z", "q", "nu_index", "real_cluster_id"]).reset_index("row", drop=True)
    ch = ch[ch.nu_index >= 0]
    ch["kine_index"] = seg.reindex(pd.MultiIndex.from_arrays([ch.index, ch.nu_index, ch.real_cluster_id])).values
    for (entry, nu), g in ch.groupby([ch.index, "nu_index"]):
        key = (entry, nu)
        if key not in lead.index or key not in D.index:
            continue
        xyz = g[["x", "y", "z"]].to_numpy(float); q = g.q.to_numpy(float)
        vtx = D.loc[key, ["nu_x", "nu_y", "nu_z"]].to_numpy(float)
        sh = (g.kine_index == lead[key]).to_numpy()
        tot = q.sum()
        near = np.linalg.norm(xyz - vtx, axis=1) < VTX_RADIUS
        out.loc[key, "vtx_other_q"] = q[near & ~sh].sum()
        out.loc[key, "other_q_frac"] = q[~sh].sum() / tot if tot > 0 else np.nan
        if sh.sum() < 3:
            continue
        p = xyz[sh]; w = np.clip(q[sh], 0, None)
        start = p[np.argmin(np.linalg.norm(p - vtx, axis=1))]
        m = np.linalg.norm(p - start, axis=1) < SHW_AXIS_RADIUS
        if m.sum() < 3:
            continue
        if w[m].sum() <= 0:
            w = np.ones_like(w)
        c = np.average(p[m], axis=0, weights=w[m])
        u = np.linalg.svd((p[m] - c) * np.sqrt(w[m])[:, None], full_matrices=False)[2][0]
        if np.dot(c - start, u) < 0:
            u = -u
        out.loc[key, ["shw1_start_x", "shw1_start_y", "shw1_start_z"]] = start
        out.loc[key, ["shw1_dir_x", "shw1_dir_y", "shw1_dir_z"]] = u
        out.loc[key, "shw1_back_dist"] = _back_dist(start, u)
        out.loc[key, "shw1_start_vtx_dist"] = np.linalg.norm(start - vtx)
        out.loc[key, "shw1_npts"] = int(sh.sum())
    return out


def make_wc_candidate_df(f):
    extra = [l for l in wc.wc_tagger_leaves(f) if l not in TAGGER_SKIP and l not in wc.wctaggerbranches]
    nu = wc.make_wc_nudf(f, tagger_leaves=wc.wctaggerbranches + extra)
    if not len(nu):
        return nu
    entry = nu.index.get_level_values("entry")

    ev = wc.make_wc_evtdf(f)[["run", "subrun", "evt", "n_nu"]]
    run = wc.make_wc_rundf(f).reset_index("row", drop=True)[
        ["nu_n_flashes", "nu_n_flashes_in_window", "nu_n_bundles"]]
    front = ev.join(run).reindex(entry)
    front.index = nu.index
    D = front.join(nu)

    cl = wc.make_wc_clusterdf(f).reset_index("row", drop=True).set_index("cluster_id", append=True)
    cl = cl[~cl.index.duplicated()][["npoints", "length_cm", "fc", "tgm", "stm", "lm"]]
    cl.columns = ["main_cluster_" + c for c in cl.columns]
    D = D.join(cl, on=["entry", "cluster_id"])

    b = wc.make_wc_bundledf(f)
    b = b[b.nu_index >= 0].reset_index("row", drop=True).set_index("nu_index", append=True)
    b = b[["reason", "n_main"]].add_prefix("bundle_")
    D = D.join(b)

    D = D.join(_particles(f, D.index)).join(_segments(f, D.index)).join(_charge(f, D.index))
    return D.join(_shower_geom(f, D))
