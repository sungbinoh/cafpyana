"""Dataframe makers for the SBND Wire-Cell reconstruction (recTreeWireCell).

The Wire-Cell output is a separate tree, entry-aligned with recTree (same
run/subrun/evt per entry), with branches "wc.<T_tree>.<leaf>". Each T_* tree is
a vector per event; T_tagger and T_kine have one row per Wire-Cell neutrino
candidate, keyed by nu_index. All makers here return frames indexed by
(entry, nu_index) or (entry, nu_index, <row>) so they join directly.
Format and production: github.com/HaiwangYu/sbnd-wirecell-production (aurora/).

- wc.valid = 0 marks a CAF event Wire-Cell failed on (all its trees empty);
  make_wc_evtdf keeps it, the per-candidate makers simply have no rows.
- T_pf, T_match, T_overlap, T_truth_track exist only on MC made with the
  truth match. On data their branches are present but always have 0 rows, so
  reco quantities must not be built from T_pf.

Units are as in the CAF: Wire-Cell energies are MeV (kine_*, ke_mev), the
truth tables (T_truth_nu.E) are GeV, lengths cm.

Only pyanalib is imported, so this module does not need CVMFS.
"""

import numpy as np
import pandas as pd

from pyanalib.pandas_helpers import loadbranches

WC_TREE = "recTreeWireCell"
WC_NU_LEVELS = ["entry", "nu_index"]

# T_tagger leaves loaded by default: vertex, flash, and the summary flags and
# BDT scores. Use wc_tagger_leaves(f, prefixes) for the BDT input variables.
wctaggerbranches = [
    "nu_x", "nu_y", "nu_z", "has_vertex", "neutrino_type",
    "cluster_id", "sel_cluster_id", "vertex_moved_cluster", "matched_flash_gid",
    "flash_time_us", "flash_pe", "flash_tpc", "flash_group",
    "match_isFC", "cosmic_flag", "numu_cc_flag", "photon_flag",
    "numu_score", "nue_score",
]

wckinebranches = [
    "kine_reco_Enu", "kine_reco_add_energy",
    "kine_nu_x_corr", "kine_nu_y_corr", "kine_nu_z_corr",
    "kine_pio_mass", "kine_pio_flag", "kine_pio_vtx_dis", "kine_pio_angle",
    "kine_pio_energy_1", "kine_pio_theta_1", "kine_pio_phi_1", "kine_pio_dis_1",
    "kine_pio_energy_2", "kine_pio_theta_2", "kine_pio_phi_2", "kine_pio_dis_2",
    "kine_energy_excluded", "kine_n_excluded", "kine_main_vertex_id",
]

# Per-particle vectors inside each T_kine row
wckineparticlebranches = [
    "kine_energy_particle", "kine_particle_type", "kine_energy_info",
    "kine_energy_included", "kine_particle_id", "kine_particle_is_shower",
    "kine_particle_cluster_id", "kine_particle_nseg", "kine_particle_pool",
]

# T_pf: the particle-flow hierarchy, MC with truth match only. is_pseudo marks
# inferred particles (a gamma for a gap-converted shower, a neutron for a
# detached proton).
wcpfbranches = [
    "nu_index", "node_id", "parent_id", "depth", "n_children", "pdg",
    "is_pseudo", "kine_index", "ke_mev",
    "start_x", "start_y", "start_z", "end_x", "end_y", "end_z",
]


def has_wc(f):
    return WC_TREE in f


def _wc_tree(f):
    if not has_wc(f):
        raise KeyError("%s not in file; it was not made with the Wire-Cell CAF maker." % WC_TREE)
    return f[WC_TREE]


def _load(f, group, leaves):
    """Flat per-row leaves of wc.<group>, columns = leaf names."""
    return loadbranches(_wc_tree(f), ["wc.%s.%s" % (group, l) for l in leaves]).wc[group]


def _index_by_nu(df):
    """(entry, row) -> (entry, nu_index)."""
    df = df.reset_index(level=1, drop=True).set_index("nu_index", append=True)
    df.index.names = WC_NU_LEVELS
    return df


def make_wc_evtdf(f):
    """One row per CAF event: run, subrun, evt, valid, n_nu (candidates)."""
    names = {"wc.run": "run", "wc.subrun": "subrun", "wc.evt": "evt", "wc.valid": "valid",
             "wc.T_tagger..length": "n_nu"}
    df = _wc_tree(f).arrays(list(names), library="pd").rename(columns=names)
    df.index.name = "entry"
    return df


def wc_tagger_leaves(f, prefixes=None):
    """Scalar-per-candidate T_tagger leaves, optionally filtered by prefix.

    The ~200 nested (vector-per-candidate) leaves are excluded.
    """
    keys = _wc_tree(f).keys()
    head = "wc.T_tagger."
    nested = {k[len(head):].split("..")[0] for k in keys if k.startswith(head) and k.endswith("..idx")}
    leaves = [k[len(head):] for k in keys
              if k.startswith(head) and ".." not in k and k[len(head):] not in nested]
    if prefixes is not None:
        leaves = [l for l in leaves if l.startswith(tuple(prefixes))]
    return leaves


def make_wc_taggerdf(f, leaves=wctaggerbranches):
    leaves = ["nu_index"] + [l for l in leaves if l != "nu_index"]
    return _index_by_nu(_load(f, "T_tagger", leaves))


def make_wc_kinedf(f, leaves=wckinebranches):
    leaves = ["nu_index"] + [l for l in leaves if l != "nu_index"]
    return _index_by_nu(_load(f, "T_kine", leaves))


def make_wc_kineparticledf(f, leaves=wckineparticlebranches):
    """One row per T_kine particle, index (entry, nu_index, kine_particle)."""
    t = _wc_tree(f)
    # loadbranches needs one vector structure per call, so one leaf at a time
    cols = []
    for l in leaves:
        c = loadbranches(t, ["wc.T_kine." + l]).wc.T_kine[l]
        c.index.names = ["entry", "kine_row", "kine_particle"]
        cols.append(c)
    df = pd.concat(cols, axis=1)
    nu = loadbranches(t, ["wc.T_kine.nu_index"]).wc.T_kine.nu_index
    nu.index.names = ["entry", "kine_row"]
    df = df.join(nu, on=["entry", "kine_row"])
    df = df.reset_index(level="kine_row", drop=True).set_index("nu_index", append=True)
    return df.reorder_levels(["entry", "nu_index", "kine_particle"])


def make_wc_pfdf(f, leaves=wcpfbranches):
    """One row per T_pf particle, index (entry, nu_index, pf)."""
    df = _load(f, "T_pf", leaves)
    df.index.names = ["entry", "pf"]
    df = df.set_index("nu_index", append=True)
    return df.reorder_levels(["entry", "nu_index", "pf"])


def make_wc_tmatchdf(f, level=3):
    """Truth neutrino matched to each Wire-Cell candidate, from T_match.

    T_match rows are reco (side 0) or truth (side 1) objects at a level:
    0 cluster, 1 particle (T_kine row vs prong), 2 candidate (T_kine particles
    vs interaction), 3 candidate (particle-flow tree vs interaction). The
    candidate is the side-0 row with nu_index >= 0; its match_obj lists up to
    3 side-1 'obj' ids at the same level, best point IoU first.

    tmatch_idx is the index into rec.mc.nu (T_truth_nu is in the same order);
    -1 when nothing matched. tmatch_iou_q / _n are the charge / point IoU.
    """
    t = _wc_tree(f)
    if "wc.T_match.level" not in t.keys():
        raise KeyError("no T_match: this file was made without the Wire-Cell truth match (WC_TRUTH_MATCH=0 or schema v1).")
    m = _load(f, "T_match", ["level", "side", "obj", "nu_index", "truth_nu_idx"])
    m.index.names = ["entry", "row"]
    nested = {}
    for l in ["match_obj", "match_iou_q", "match_iou_n"]:
        c = loadbranches(t, ["wc.T_match." + l]).wc.T_match[l]
        c.index.names = ["entry", "row", "imatch"]
        nested[l] = c
    nested = pd.DataFrame(nested)

    m = m[m.level == level]
    reco = m[(m.side == 0) & (m.nu_index >= 0)]
    truth = m[m.side == 1].reset_index()[["entry", "obj", "truth_nu_idx"]]

    # partners are ranked by point IoU, so the first is the best match
    cand = nested[nested.index.droplevel(-1).isin(reco.index)]
    best = cand.xs(0, level="imatch")

    out = reco[["nu_index"]].join(best)
    out = out.reset_index().merge(truth.rename(columns={"obj": "match_obj"}),
                                  on=["entry", "match_obj"], how="left")
    out = out.set_index(WC_NU_LEVELS)
    return pd.DataFrame({
        "tmatch_idx": out.truth_nu_idx.fillna(-1).astype(int),
        "tmatch_iou_q": out.match_iou_q,
        "tmatch_iou_n": out.match_iou_n,
    })


def make_wc_nudf(f, tagger_leaves=wctaggerbranches, kine_leaves=wckinebranches, include_tmatch=None):
    """T_tagger + T_kine (+ T_match truth match on MC), one row per candidate."""
    tag = make_wc_taggerdf(f, tagger_leaves)
    kine = make_wc_kinedf(f, kine_leaves)
    df = tag.join(kine.drop(columns=[c for c in kine.columns if c in tag.columns]), how="left")
    if include_tmatch is None:
        ismc = bool(loadbranches(f["recTree"], ["rec.hdr.ismc"]).rec.hdr.ismc.iloc[0])
        include_tmatch = ismc and "wc.T_match.level" in _wc_tree(f).keys()
    if include_tmatch:
        df = df.join(make_wc_tmatchdf(f), how="left")
        df["tmatch_idx"] = df.tmatch_idx.fillna(-1).astype(int)
    return df
