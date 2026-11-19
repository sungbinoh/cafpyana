"""One row per Wire-Cell neutrino candidate, index (entry, nu_index).

Columns:
  * event: run/subrun/evt, n_nu, Trun flash and bundle counts
  * every scalar T_tagger leaf (~1030: flags, scores and BDT inputs) and the
    T_kine scalars (MeV), under their Wire-Cell names; on MC the T_match
    truth match (tmatch_idx: row in rec.mc.nu, -1 none; tmatch_iou_q/_n)
  * main_cluster_*: the T_cluster row of cluster_id; bundle_*: the T_bundle
    row (reason 0 main, 1 demoted main; codes in make_wc_events.py)
  * T_kine particles connected to the candidate (kine_energy_included ==
    1): counts per PID and leading energies; the unclear-link leftovers
    apart (n_unclear, E_unclear, E_em_unclear_1); E_em_all_2 is the second
    EM energy over both
  * pf_*: e/gamma counts from the particle-flow tree (MC only, see _pf)
  * T_segment / T_rec_charge counts and lengths
  * track particles by length and dQ/dx (see _tracks): the longest track,
    the longest muon-labelled one, and the longest MIP-like ones
  * leading EM shower from its charge points (see _shower_geom): start,
    axis, shw1_back_dist (back along the axis to the AV boundary; short for
    a photon entering from outside), vtx_other_q, other_q_frac, and shape:
    shw1_a68 (68% charge opening angle), shw1_width, shw1_lobe_sep (two
    photons merged into one shower show up as two lobes)

Reco only: truth is the mcnu table (make_truth.py), joined by cutflow.load.
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
SHAPE_MIN_R = 3.0        # cm, shower points closer to the start do not enter the shape
SHAPE_MIN_POINTS = 10
LOBE_MIN_FRAC = 0.15     # charge share each lobe needs for shw1_lobe_sep
DQDX_MIP = 3350.0        # T_rec_charge q/nq of a MIP: median of true muon tracks > 10 cm (Nov16)
MIP_LIKE_MAX = 1.6       # MIPs; protons sit at ~3 (16% quantile ~2.1), muons at 1 (84% ~1.2)
DQDX_START_CM = 4.0      # cm, start of a track for mu1_dqdx_start
DQDX_END_CM = 3.0        # cm, end of a track for mu1_end_ratio (Bragg peak)


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
    # every EM particle, also the unclear-link leftovers (not attached to the candidate's particle flow)
    T["E_em_all_2"] = _nth(E[pdg == 11], 2)
    T["E_em_unclear_1"] = _nth(E[~conn & (pdg == 11)], 1)
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


def _pf_ancestors(t, max_depth=50):
    """Per T_pf node: real (non-pseudo) nodes above it, and whether a pseudo
    neutron is above it."""
    entry, nu = t.index.get_level_values(0), t.index.get_level_values(1)
    cur = t.parent_id.to_numpy()
    real_above = np.zeros(len(t), int)
    neutron_above = np.zeros(len(t), bool)
    for _ in range(max_depth):
        up = cur >= 0
        if not up.any():
            break
        p = t.reindex(pd.MultiIndex.from_arrays([entry, nu, cur]))
        ok = up & p.parent_id.notna().to_numpy()
        pseudo = p.is_pseudo.to_numpy() == 1
        real_above += ok & ~pseudo
        neutron_above |= ok & pseudo & (p.pdg.abs().to_numpy() == 2112)
        cur = np.where(ok, p.parent_id.fillna(-1).to_numpy(), -1).astype(int)
    return real_above, neutron_above


def _pf(f, index):
    """e/gamma counts from the T_pf particle-flow tree. A shower is primary
    when no real particle and no pseudo neutron is above it: attached to the
    vertex (depth 0) or displaced through pseudo photons / pi0s. NaN when
    the file has no T_pf (data)."""
    cols = ["pf_n_em", "pf_n_em_primary", "pf_n_em_primary_20MeV", "pf_n_em_attached",
            "pf_n_em_displaced", "pf_E_em_primary_1", "pf_E_em_primary_2"]
    t = wc.make_wc_pfdf(f)
    if not len(t):
        return pd.DataFrame(np.nan, index=index, columns=cols)
    t = t.reset_index("row", drop=True).set_index(["nu_index", "node_id"], append=True)
    t = t[~t.index.duplicated()]
    real_above, neutron_above = _pf_ancestors(t)
    em = (t.pdg.abs() == 11) & (t.is_pseudo == 0)
    prim = em & (real_above == 0) & ~neutron_above
    by = lambda m: m.groupby(level=LEVELS).sum()
    T = pd.DataFrame(index=index)
    T["pf_n_em"] = by(em)
    T["pf_n_em_primary"] = by(prim)
    T["pf_n_em_primary_20MeV"] = by(prim & (t.ke_mev > 20))
    T["pf_n_em_attached"] = by(prim & (t.depth == 0))
    T["pf_n_em_displaced"] = by(prim & (t.depth > 0))
    T = T.fillna(0).astype(int)
    E = t.ke_mev[prim].droplevel("node_id")
    T["pf_E_em_primary_1"] = _nth(E, 1)
    T["pf_E_em_primary_2"] = _nth(E, 2)
    return T


def _back_dist(p, u):
    """Distance from p along -u to the active-volume boundary (p inside)."""
    b = cuts.ACTIVE_VOLUME
    t = np.inf
    for i, ax in enumerate("xyz"):
        d = -u[i]
        if d > 1e-9: t = min(t, (b[ax][1] - p[i]) / d)
        elif d < -1e-9: t = min(t, (b[ax][0] - p[i]) / d)
    return max(t, 0.0)


def _shower_axis(p, q, vtx):
    """(start, unit axis) of a shower's charge points p with charges q, or
    None with too few points: the start is the point nearest the vertex, the
    axis the charge-weighted principal axis of the points within
    SHW_AXIS_RADIUS of it, pointing away from it."""
    w = np.clip(q, 0, None)
    start = p[np.argmin(np.linalg.norm(p - vtx, axis=1))]
    m = np.linalg.norm(p - start, axis=1) < SHW_AXIS_RADIUS
    if m.sum() < 3:
        return None
    if w[m].sum() <= 0:
        w = np.ones_like(w)
    c = np.average(p[m], axis=0, weights=w[m])
    u = np.linalg.svd((p[m] - c) * np.sqrt(w[m])[:, None], full_matrices=False)[2][0]
    return start, (-u if np.dot(c - start, u) < 0 else u)


def _two_lobes(v, w, min_frac=LOBE_MIN_FRAC, seeds=5):
    """Largest angle [deg] between two charge-weighted direction clusters
    (2-means on the unit vectors v) that each hold >= min_frac of the
    charge; 0 when there is no such second lobe."""
    best = 0.0
    for seed in range(seeds):
        cen = v[np.random.RandomState(seed).choice(len(v), 2, replace=False)]
        for _ in range(20):
            lab = np.argmax(v @ cen.T, axis=1)
            new = np.array([np.average(v[lab == k], axis=0, weights=w[lab == k] + 1e-9) if (lab == k).any() else cen[k]
                            for k in (0, 1)])
            new /= np.linalg.norm(new, axis=1)[:, None]
            if np.allclose(new, cen):
                break
            cen = new
        frac = np.array([w[lab == k].sum() for k in (0, 1)]) / w.sum()
        if frac.min() >= min_frac:
            best = max(best, np.degrees(np.arccos(np.clip(cen[0] @ cen[1], -1, 1))))
    return best


def _shower_shape(p, q, vtx):
    """(a68, width, lobe_sep) of a shower's charge points, or None with fewer
    than SHAPE_MIN_POINTS beyond SHAPE_MIN_R of the start: the angle around
    the charge-weighted principal axis (from the start, the point nearest
    the vertex) holding 68% of the charge [deg], the transverse /
    longitudinal spread, and the angle between two charge lobes [deg]."""
    q = np.clip(q, 0, None) + 1e-9
    start = p[np.argmin(np.linalg.norm(p - vtx, axis=1))]
    d = p - start
    r = np.linalg.norm(d, axis=1)
    m = r > SHAPE_MIN_R
    if m.sum() < SHAPE_MIN_POINTS:
        return None
    c = np.average(p, axis=0, weights=q)
    _, sv, vt = np.linalg.svd((p - c) * np.sqrt(q)[:, None], full_matrices=False)
    ax = vt[0] if np.dot(c - start, vt[0]) > 0 else -vt[0]
    ang = np.degrees(np.arccos(np.clip((d[m] @ ax) / r[m], -1, 1)))
    o = np.argsort(ang)
    a68 = ang[o][np.searchsorted(np.cumsum(q[m][o]) / q[m].sum(), 0.68)]
    return a68, sv[1] / sv[0], _two_lobes(d[m] / r[m][:, None], q[m])


def _shower_geom(f, D):
    """Leading EM shower geometry from the T_rec_charge points mapped to it
    through T_segment, plus the candidate charge outside it: vtx_other_q
    (within VTX_RADIUS of the vertex) and other_q_frac (all), and its shape
    (_shower_shape): shw1_a68, shw1_width, shw1_lobe_sep."""
    cols = ["shw1_start_x", "shw1_start_y", "shw1_start_z", "shw1_dir_x", "shw1_dir_y", "shw1_dir_z",
            "shw1_back_dist", "shw1_start_vtx_dist", "shw1_npts", "vtx_other_q", "other_q_frac",
            "shw1_a68", "shw1_width", "shw1_lobe_sep"]
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
        shape = _shower_shape(xyz[sh], q[sh], vtx) if sh.sum() >= SHAPE_MIN_POINTS else None
        if shape is not None:
            out.loc[key, ["shw1_a68", "shw1_width", "shw1_lobe_sep"]] = shape
        axis = _shower_axis(xyz[sh], q[sh], vtx) if sh.sum() >= 3 else None
        if axis is None:
            continue
        start, u = axis
        out.loc[key, ["shw1_start_x", "shw1_start_y", "shw1_start_z"]] = start
        out.loc[key, ["shw1_dir_x", "shw1_dir_y", "shw1_dir_z"]] = u
        out.loc[key, "shw1_back_dist"] = _back_dist(start, u)
        out.loc[key, "shw1_start_vtx_dist"] = np.linalg.norm(start - vtx)
        out.loc[key, "shw1_npts"] = int(sh.sum())
    return out


def _track_profile(xyz, dqdx, vtx):
    """(median, first-DQDX_START_CM median, last-DQDX_END_CM median / median)
    of a track's dQ/dx, ordered from the end nearer the vertex."""
    c = xyz.mean(0)
    u = np.linalg.svd(xyz - c, full_matrices=False)[2][0]
    t = (xyz - c) @ u
    if abs(t.min() - (vtx - c) @ u) > abs(t.max() - (vtx - c) @ u):
        t = -t
    o = np.argsort(t)
    d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(xyz[o], axis=0), axis=1))]
    q = dqdx[o]
    med = np.median(q)
    return med, np.median(q[d <= DQDX_START_CM]), np.median(q[d >= d[-1] - DQDX_END_CM]) / med


def _tracks(f, D):
    """Track particles (connected, not showers) by length and dQ/dx, from the
    T_rec_charge fit points (q / nq, file charge units per cm; DQDX_MIP is
    the median of true muon tracks):
      trk1_len, trk1_dqdx, trk1_pdg       longest track of any label
      mu1_len, mu1_dqdx, mu1_dqdx_start, mu1_end_ratio
                                          longest muon-labelled track
      mu_mip1_len                         longest muon-labelled track with a
                                          MIP-like median (< MIP_LIKE_MAX MIPs)
      mip1_len                            the same for any label"""
    cols = ["trk1_len", "trk1_dqdx", "trk1_pdg", "mu1_len", "mu1_dqdx", "mu1_dqdx_start", "mu1_end_ratio",
            "mu_mip1_len", "mip1_len"]
    out = pd.DataFrame(np.nan, index=D.index, columns=cols)
    k = wc.make_wc_kineparticledf(f)
    k = k[(k.kine_energy_included == 1) & (k.kine_particle_is_shower == 0)]
    if not len(k):
        return out
    seg = wc.make_wc_segmentdf(f, ["nu_index", "segment_id", "kine_index", "length_cm"]).reset_index("row", drop=True)
    seg = seg[seg.nu_index >= 0]
    length = seg.groupby([seg.index, "nu_index", "kine_index"]).length_cm.sum()
    sid = seg.set_index(["nu_index", "segment_id"], append=True).kine_index
    ch = wc.make_wc_chargedf(f, ["x", "y", "z", "q", "nq", "nu_index", "real_cluster_id"]).reset_index("row", drop=True)
    ch = ch[(ch.nu_index >= 0) & (ch.nq > 0)]
    ch["kine_index"] = sid.reindex(pd.MultiIndex.from_arrays([ch.index, ch.nu_index, ch.real_cluster_id])).values
    ch = ch[ch.kine_index.notna()]
    rows = []
    for (entry, nu, ki), g in ch.groupby([ch.index, "nu_index", "kine_index"]):
        key = (entry, nu, int(ki))
        if key not in k.index or len(g) < 3 or (entry, nu) not in D.index:
            continue
        vtx = D.loc[(entry, nu), ["nu_x", "nu_y", "nu_z"]].to_numpy(float)
        med, start, end = _track_profile(g[["x", "y", "z"]].to_numpy(float), (g.q / g.nq).to_numpy(float), vtx)
        rows.append((entry, nu, abs(int(k.loc[key, "kine_particle_type"])), length.get(key, 0.0), med, start, end))
    if not rows:
        return out
    t = pd.DataFrame(rows, columns=["entry", "nu_index", "pdg", "len", "dqdx", "dqdx_start", "end_ratio"])
    t = t.sort_values("len", ascending=False)
    mip = t.dqdx < MIP_LIKE_MAX * DQDX_MIP
    first = lambda m: t[m].groupby(["entry", "nu_index"]).first()
    a, mu = first(np.ones(len(t), bool)), first(t.pdg == PID["muon"])
    out["trk1_len"], out["trk1_dqdx"], out["trk1_pdg"] = a.len, a.dqdx, a.pdg
    out["mu1_len"], out["mu1_dqdx"], out["mu1_dqdx_start"], out["mu1_end_ratio"] = mu.len, mu.dqdx, mu.dqdx_start, mu.end_ratio
    out["mu_mip1_len"] = first((t.pdg == PID["muon"]) & mip).len
    out["mip1_len"] = first(mip).len
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

    D = D.join(_particles(f, D.index)).join(_pf(f, D.index)).join(_segments(f, D.index)).join(_charge(f, D.index))
    return D.join(_shower_geom(f, D)).join(_tracks(f, D))
