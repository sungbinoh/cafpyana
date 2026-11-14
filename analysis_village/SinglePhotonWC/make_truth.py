"""Truth table for the inclusive single-photon search (signal definition:
cuts.py). One row per true neutrino (rec.mc.nu), index (entry, nu), so
candidates join on (entry, tmatch_idx). Columns, all prefixed true_:

  * rec.mc.nu: E, pdg, iscc, isnc, genie_mode, w, resnum, vtx_x/y/z
  * primary counts: npi0, neta, npipm, np_35MeV, nn, nmu, nmu_100MeV, ne
    (e+-), mu_KE (leading muon, GeV)
  * prompt photons (GENIE primaries and decay photons of primary hadrons):
    ngamma_20MeV (E >= 20 MeV); of those ngamma_rec (shower starts >=
    cuts.SHOWER_INSET inside the TPC) and ngamma_out (the rest), of which
    ngamma_edge convert inside the TPC; ndecayg_soft (decay photons < 20
    MeV); nshower (ngamma_rec after merging photons < 20 deg apart) and
    nmerged = ngamma_rec - nshower
  * wc_q_dep / wc_q_img: charge the neutrino deposited / that was imaged;
    wc_f_cand / _sameflash / _otherflash / _noflash: fractions of its
    cluster-matched charge in its candidate's particle flow, other clusters
    on that candidate's flash, clusters on another flash, unmatched clusters
    (MC truth-match trees; NaN without them)
  * shw1_* / shw2_*, the two most energetic showers: E [GeV], gap (vertex to
    nearest conversion, cm), costh (to the beam), conv_x/y/z, nphot, origin
    (22: GENIE primary, else the parent hadron pdg)

The shower start is the photon trajectory end: SBND's Geant4 keeps no EM
shower daughters and the CAF stores only in-TPC points, so an escaping
photon ends on the boundary and fails the inset.
"""

import numpy as np
import pandas as pd

from pyanalib.pandas_helpers import loadbranches

from makedf import wirecell as wc

from analysis_village.SinglePhotonWC import cuts

PROTON_MASS = 0.938272   # GeV
MUON_MASS = 0.105658     # GeV
MUON_KE_100MEV = 0.100   # GeV, nmu_100MeV (MicroBooNE's muon threshold)
G4_DECAY = 3             # caf::g4_process_ kG4Decay
LEPTONS = [11, 12, 13, 14, 15, 16]
MAX_DECAY_DEPTH = 4      # eta -> 3pi0 -> 6gamma needs 2
EDGE_EPS = 0.5           # cm; trajectory ends closer to the boundary are escapes


def _load(f, prefix, leaves):
    df = loadbranches(f["recTree"], [prefix + l for l in leaves])
    for level in prefix.strip(".").split("."):
        df = df[level]
    df.columns = [("_".join(c for c in col if c) if isinstance(col, tuple) else col).replace(".", "_")
                  for col in df.columns]
    return df


def _prompt_photons(prim, tp):
    """Photons from the neutrino interaction: primaries, and the decay chains
    of primary hadrons. Columns entry, nu, origin plus the tp columns."""
    p = prim.reset_index()
    t = tp.reset_index()
    t.columns = ["entry", "itp"] + list(t.columns[2:])
    key = ["entry", "G4ID"]

    out = [t.merge(p.loc[p.pdg == 22, key + ["nu"]], on=key).assign(origin=cuts.ORIGIN_PRIMARY)]
    front = p.loc[~p.pdg.abs().isin(LEPTONS + [22]), key + ["nu", "pdg"]].rename(columns={"pdg": "origin"})
    dec = t[t.start_process == G4_DECAY]
    for _ in range(MAX_DECAY_DEPTH):
        if not len(front):
            break
        kids = dec.merge(front.rename(columns={"G4ID": "parent"}), on=["entry", "parent"])
        out.append(kids[kids.pdg == 22])
        front = kids.loc[~kids.pdg.abs().isin(LEPTONS + [22]), key + ["nu", "origin"]]
    g = pd.concat(out, ignore_index=True)
    p_ = np.sqrt(g.genp_x**2 + g.genp_y**2 + g.genp_z**2)
    for c in "xyz":
        g["d" + c] = g["genp_" + c] / p_
    return g


def _merge_showers(r):
    """Shower label per reconstructable photon: photons closer than
    MERGE_ANGLE to each other (single linkage) form one shower."""
    r = r.sort_values(["entry", "nu", "genE"], ascending=[True, True, False]).reset_index(drop=True)
    r["shower"] = r.groupby(["entry", "nu"]).cumcount()
    multi = r[r.groupby(["entry", "nu"]).genE.transform("size") > 1]
    for _, g in multi.groupby(["entry", "nu"]):
        d = g[["dx", "dy", "dz"]].to_numpy()
        ang = np.degrees(np.arccos(np.clip(d @ d.T, -1, 1)))
        lab = np.arange(len(g))
        for i in range(len(g)):
            for j in range(i):
                if ang[i, j] < cuts.MERGE_ANGLE:
                    lab[lab == lab[i]] = lab[j]
        r.loc[g.index, "shower"] = lab
    return r


def _nested(t, base):
    """vector<vector> leaf -> per entry, list of per-row arrays."""
    a = t.arrays([base, base + "..length", base + "..idx"], library="np")
    return [[v[i:i + l] for i, l in zip(I, L)] for v, L, I in zip(a[base], a[base + "..length"], a[base + "..idx"])]


def _cluster_charge_by_nu(l0, match_obj, match_q):
    """{true nu: [(cluster_id, shared charge)]} from one entry's level-0
    T_match rows (side 0 clusters, side 1 interactions)."""
    truth = l0[l0.side == 1]
    truth_nu = dict(zip(truth.obj, truth.truth_nu_idx))
    acc = {}
    for row, r in l0[l0.side == 0].iterrows():
        for o, sq in zip(match_obj[row], match_q[row]):
            nu = truth_nu.get(o, -1)
            if nu >= 0:
                acc.setdefault(nu, []).append((int(r.cluster_id), float(sq)))
    return acc


def _charge_fractions(clusters, cand_clusters, flash, cand_flash):
    """Fractions of the charge in `clusters` on the candidate, its flash,
    another flash and no flash."""
    tot = sum(x for _, x in clusters) or np.nan
    fc = fs = fo = fn = 0.0
    for c, x in clusters:
        if c in cand_clusters: fc += x
        elif flash.get(c, -1) < 0: fn += x
        elif cand_flash is not None and flash.get(c) == cand_flash: fs += x
        else: fo += x
    return [fc / tot, fs / tot, fo / tot, fn / tot]


def _wc_charge_destination(f, N):
    """Where each true neutrino's charge went (true_wc_* columns)."""
    cols = ["wc_q_dep", "wc_q_img", "wc_f_cand", "wc_f_sameflash", "wc_f_otherflash", "wc_f_noflash"]
    out = pd.DataFrame(np.nan, index=N.index, columns=cols)
    if not wc.has_wc(f) or "wc.T_match.level" not in f[wc.WC_TREE].keys():
        return out
    t = f[wc.WC_TREE]
    tt = wc.make_wc_tabledf(f, "T_truth_track", ["nu_idx", "q_readout", "q_in_blobs"])
    tt = tt[tt.nu_idx >= 0].reset_index("row", drop=True).set_index("nu_idx", append=True)
    q = tt.groupby(level=["entry", "nu_idx"]).sum()
    q.index.names = N.index.names
    out["wc_q_dep"] = q.q_readout.reindex(N.index)
    out["wc_q_img"] = q.q_in_blobs.reindex(N.index)

    m = wc.make_wc_tabledf(f, "T_match", ["level", "side", "obj", "cluster_id", "truth_nu_idx"])
    mo, ms = _nested(t, "wc.T_match.match_obj"), _nested(t, "wc.T_match.match_shared_q")
    cl = wc.make_wc_clusterdf(f, ["cluster_id", "matched_flash_gid"])
    tag = wc.make_wc_taggerdf(f, ["matched_flash_gid"])
    tm = wc.make_wc_tmatchdf(f)
    ov = wc.make_wc_tabledf(f, "T_overlap", ["nu_index", "cluster_id", "in_pf"])
    has = lambda df, entry: entry in df.index.get_level_values(0)
    for entry in np.unique(m.index.get_level_values("entry")):
        me = m.loc[entry]
        acc = _cluster_charge_by_nu(me[me.level == 0], mo[entry], ms[entry])
        flash = dict(zip(cl.loc[entry].cluster_id, cl.loc[entry].matched_flash_gid)) if has(cl, entry) else {}
        cand = {}   # true nu -> matched candidates
        if has(tm, entry):
            for nu_index, r in tm.loc[entry].iterrows():
                if r.tmatch_idx >= 0:
                    cand.setdefault(int(r.tmatch_idx), []).append(nu_index)
        for nu, clusters in acc.items():
            key = (entry, nu)
            if key not in out.index:
                continue
            cset, cflash = set(), None
            if nu in cand:
                ci = cand[nu][0]
                if has(ov, entry):
                    oe = ov.loc[entry]
                    cset = set(oe[(oe.nu_index == ci) & (oe.in_pf == 1)].cluster_id)
                if (entry, ci) in tag.index:
                    cflash = tag.loc[(entry, ci), "matched_flash_gid"]
            out.loc[key, cols[2:]] = _charge_fractions(clusters, cset, flash, cflash)
    return out


def make_sp_truthdf(f):
    N = _load(f, "rec.mc.nu.", ["E", "pdg", "iscc", "isnc", "genie_mode", "w", "resnum",
                                "position.x", "position.y", "position.z"])
    N = N.rename(columns={"position_x": "vtx_x", "position_y": "vtx_y", "position_z": "vtx_z"})
    if not len(N):
        return N.add_prefix("true_")
    names = N.index.names
    N.index.names = ["entry", "nu"]
    cnt = lambda m, by: m.groupby(by).sum().reindex(N.index).fillna(0).astype(int)

    # primaries
    prim = _load(f, "rec.mc.nu.prim.", ["pdg", "G4ID", "genE"])
    prim.index.names = ["entry", "nu", "iprim"]
    pdg, nu = prim.pdg.abs(), ["entry", "nu"]
    muKE = (prim.genE - MUON_MASS).where(pdg == 13)
    N["npi0"] = cnt(pdg == 111, nu)
    N["neta"] = cnt(pdg == 221, nu)
    N["npipm"] = cnt(pdg == 211, nu)
    N["np_35MeV"] = cnt((pdg == 2212) & (prim.genE - PROTON_MASS > 0.035), nu)
    N["nn"] = cnt(pdg == 2112, nu)
    N["nmu"] = cnt(pdg == 13, nu)
    N["nmu_100MeV"] = cnt(muKE >= MUON_KE_100MEV, nu)
    N["mu_KE"] = muKE.groupby(nu).max().reindex(N.index)
    N["ne"] = cnt(pdg == 11, nu)

    # prompt photons
    tp = _load(f, "rec.true_particles.", ["pdg", "G4ID", "parent", "start_process", "genE",
                                          "genp.x", "genp.y", "genp.z", "end.x", "end.y", "end.z"])
    g = _prompt_photons(prim, tp).join(N[["vtx_x", "vtx_y", "vtx_z"]], on=nu)
    g["gap"] = np.sqrt((g.end_x - g.vtx_x)**2 + (g.end_y - g.vtx_y)**2 + (g.end_z - g.vtx_z)**2)
    hard = g.genE >= cuts.PHOTON_E_MIN
    g["rec"] = hard & (cuts.dist_to_active_boundary(g.end_x, g.end_y, g.end_z) >= cuts.SHOWER_INSET)
    by = [g.entry, g.nu]
    N["ngamma_20MeV"] = cnt(hard, by)
    N["ngamma_rec"] = cnt(g.rec, by)
    N["ngamma_out"] = cnt(hard & ~g.rec, by)
    # escaping photons end on the boundary (distance ~0), so require > EDGE_EPS inside
    N["ngamma_edge"] = cnt(hard & ~g.rec & (cuts.dist_to_active_boundary(g.end_x, g.end_y, g.end_z) > EDGE_EPS), by)
    N["ndecayg_soft"] = cnt(~hard & (g.origin != cuts.ORIGIN_PRIMARY), by)

    # reconstructable showers, most energetic first
    s = _merge_showers(g[g.rec]).groupby(["entry", "nu", "shower"]).agg(
        E=("genE", "sum"), gap=("gap", "min"), costh=("dz", "first"),
        conv_x=("end_x", "first"), conv_y=("end_y", "first"), conv_z=("end_z", "first"),
        nphot=("genE", "size"), origin=("origin", "first")).reset_index()
    s = s.sort_values(["entry", "nu", "E"], ascending=[True, True, False])
    s["rank"] = s.groupby(nu).cumcount()
    s = s.set_index(nu)
    N["nshower"] = cnt(s["rank"] >= 0, nu)
    N["nmerged"] = N.ngamma_rec - N.nshower
    conv = ["conv_x", "conv_y", "conv_z"]
    N = N.join(s.loc[s["rank"] == 0, ["E", "gap", "costh", "nphot", "origin"] + conv].add_prefix("shw1_"))
    N = N.join(s.loc[s["rank"] == 1, ["E", "gap"] + conv].add_prefix("shw2_"))
    N = N.join(_wc_charge_destination(f, N))
    N.index.names = names
    return N.add_prefix("true_")
