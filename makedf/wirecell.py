"""Dataframe makers for the SBND Wire-Cell reconstruction (recTreeWireCell).

The Wire-Cell tree is entry-aligned with recTree; its branches are
"wc.<T_tree>.<leaf>", each T_* tree a vector per event. T_tagger and T_kine
have one row per neutrino candidate (nu_index), so the per-candidate makers
index by (entry, nu_index) and the table makers by (entry, row).
Format: github.com/HaiwangYu/sbnd-wirecell-production (aurora/); source:
github.com/HaiwangYu/wire-cell-toolkit, branch sbnd-wc-flatcaf.

- wc.valid = 0: no Wire-Cell output for the event (failed job or lost
  upstream), all trees empty. valid = 1 without T_tagger rows is a processed
  event with no candidate; T_bundle.reason says why.
- T_pf, T_match, T_overlap and T_truth_track are written only with the
  truth match (MC). T_pf is reco but gated the same way: empty on data.

Units as in the CAF: Wire-Cell energies MeV, truth tables GeV, lengths cm.
Only pyanalib is imported, so this module does not need CVMFS.
"""

import numpy as np
import pandas as pd

from pyanalib.pandas_helpers import loadbranches

WC_TREE = "recTreeWireCell"
WC_NU_LEVELS = ["entry", "nu_index"]

# T_tagger defaults (wc_tagger_leaves(f) lists every scalar leaf).
# cosmict_flag is the overall cosmic verdict; neutrino_type bits: 2 cosmic,
# 4 numu CC, 8 NC, 32 nue (0: no vertex). cluster_id is the main cluster
# after the vertex finder.
wctaggerbranches = [
    "nu_x", "nu_y", "nu_z", "has_vertex", "neutrino_type",
    "cluster_id", "matched_flash_gid",
    "flash_time_us", "flash_pe", "flash_tpc",
    "match_isFC", "cosmict_flag", "numu_cc_flag", "photon_flag",
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
wckineparticlebranches = [   # per-particle vectors inside each T_kine row
    "kine_energy_particle", "kine_particle_type", "kine_energy_info",
    "kine_energy_included", "kine_particle_id", "kine_particle_is_shower",
    "kine_particle_cluster_id", "kine_particle_nseg", "kine_particle_pool",
]
wcrunbranches = [
    "nu_n_flashes", "nu_n_flashes_in_window", "nu_n_bundles", "nu_n_candidates",
    "beam_window_low_us", "beam_window_high_us",
]
# tgm / stm: through-going / stopping muon; lm: -1 unset, 1 low energy, 2 light mismatch
wcclusterbranches = ["cluster_id", "is_main", "beam_flash", "fc", "tgm", "stm", "lm",
                     "npoints", "length_cm", "flash_time_us", "flash_pe"]
wcbundlebranches = ["gid", "reason", "nu_index", "n_main", "sel_cluster_id", "sel_length_cm",
                    "flash_time_us", "flash_pe"]
wcsegmentbranches = ["nu_index", "segment_id", "kine_index", "pdg", "flag_shower", "dirsign",
                     "length_cm", "ke_mev", "in_enu"]
wcchargebranches = ["x", "y", "z", "q", "rr", "nu_index", "sub_cluster_id"]
# T_pf nodes: depth 0 hangs off the vertex; is_pseudo nodes (pdg 22, 111,
# 2112) are invisible links to displaced daughters
wcpfbranches = ["nu_index", "node_id", "parent_id", "depth", "pdg", "is_pseudo", "kine_index", "ke_mev"]


def has_wc(f):
    return WC_TREE in f


def _wc_tree(f):
    if not has_wc(f):
        raise KeyError("%s not in file; it was not made with the Wire-Cell CAF maker." % WC_TREE)
    return f[WC_TREE]


def _index_by_nu(df):
    """(entry, row) -> (entry, nu_index)."""
    df = df.reset_index(level=1, drop=True).set_index("nu_index", append=True)
    df.index.names = WC_NU_LEVELS
    return df


def _nested_leaf(t, group, leaf, names):
    c = loadbranches(t, ["wc.%s.%s" % (group, leaf)]).wc[group][leaf]
    c.index.names = names
    return c


# --- Tables: one row per row of a wc.T_* tree, index (entry, row)

def make_wc_tabledf(f, group, leaves):
    """Flat per-row leaves of wc.<group>. A tree with no rows (MC-only trees
    on data) gives an empty frame with the right columns."""
    t = _wc_tree(f)
    br = ["wc.%s.%s" % (group, l) for l in leaves]
    a = t.arrays(br, library="np")
    n = np.array([len(x) for x in a[br[0]]], dtype=int)
    entry = np.repeat(np.arange(len(n)), n)
    row = np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n)
    cols = {l: (np.concatenate(a[b]) if n.sum() else np.array([])) for l, b in zip(leaves, br)}
    return pd.DataFrame(cols, index=pd.MultiIndex.from_arrays([entry, row], names=["entry", "row"]))


def make_wc_rundf(f, leaves=wcrunbranches):
    return make_wc_tabledf(f, "Trun", leaves)


def make_wc_clusterdf(f, leaves=wcclusterbranches):
    return make_wc_tabledf(f, "T_cluster", leaves)


def make_wc_bundledf(f, leaves=wcbundlebranches):
    return make_wc_tabledf(f, "T_bundle", leaves)


def make_wc_segmentdf(f, leaves=wcsegmentbranches):
    return make_wc_tabledf(f, "T_segment", leaves)


def make_wc_chargedf(f, leaves=wcchargebranches):
    return make_wc_tabledf(f, "T_rec_charge", leaves)


def make_wc_pfdf(f, leaves=wcpfbranches):
    """Particle-flow nodes. MC only (truth-match gate)."""
    return make_wc_tabledf(f, "T_pf", leaves)


# --- Per event and per candidate

def make_wc_evtdf(f):
    """One row per CAF event: run, subrun, evt, valid, n_nu (candidates)."""
    names = {"wc.run": "run", "wc.subrun": "subrun", "wc.evt": "evt", "wc.valid": "valid",
             "wc.T_tagger..length": "n_nu"}
    df = _wc_tree(f).arrays(list(names), library="pd").rename(columns=names)
    df.index.name = "entry"
    return df


def wc_tagger_leaves(f):
    """Scalar-per-candidate T_tagger leaves (the nested ones are excluded)."""
    keys = _wc_tree(f).keys()
    head = "wc.T_tagger."
    nested = {k[len(head):].split("..")[0] for k in keys if k.startswith(head) and k.endswith("..idx")}
    return [k[len(head):] for k in keys
            if k.startswith(head) and ".." not in k and k[len(head):] not in nested]


def make_wc_taggerdf(f, leaves=wctaggerbranches):
    leaves = ["nu_index"] + [l for l in leaves if l != "nu_index"]
    return _index_by_nu(make_wc_tabledf(f, "T_tagger", leaves))


def make_wc_kinedf(f, leaves=wckinebranches):
    leaves = ["nu_index"] + [l for l in leaves if l != "nu_index"]
    return _index_by_nu(make_wc_tabledf(f, "T_kine", leaves))


def make_wc_kineparticledf(f, leaves=wckineparticlebranches):
    """One row per T_kine particle, index (entry, nu_index, kine_particle)."""
    t = _wc_tree(f)
    # loadbranches takes one vector structure per call, so one leaf at a time
    df = pd.concat([_nested_leaf(t, "T_kine", l, ["entry", "kine_row", "kine_particle"]) for l in leaves], axis=1)
    nu = loadbranches(t, ["wc.T_kine.nu_index"]).wc.T_kine.nu_index
    nu.index.names = ["entry", "kine_row"]
    df = df.join(nu, on=["entry", "kine_row"])
    df = df.reset_index(level="kine_row", drop=True).set_index("nu_index", append=True)
    return df.reorder_levels(["entry", "nu_index", "kine_particle"])


def make_wc_tmatchdf(f, level=3):
    """Truth neutrino matched to each candidate, from T_match.

    T_match rows are reco (side 0) or truth (side 1) objects at a level: 0
    cluster, 1 particle, 2 candidate (T_kine particles), 3 candidate
    (particle-flow tree). A candidate's match_obj lists up to 3 truth ids by
    decreasing point IoU (ties: charge IoU), so the first is the best.
    tmatch_idx indexes rec.mc.nu (-1: no match); tmatch_iou_q / _n are the
    charge / point IoU.
    """
    t = _wc_tree(f)
    if "wc.T_match.level" not in t.keys():
        raise KeyError("no T_match: this file was made without the Wire-Cell truth match (WC_TRUTH_MATCH=0 or schema v1).")
    m = make_wc_tabledf(f, "T_match", ["level", "side", "obj", "nu_index", "truth_nu_idx"])
    nested = pd.DataFrame({l: _nested_leaf(t, "T_match", l, ["entry", "row", "imatch"])
                           for l in ["match_obj", "match_iou_q", "match_iou_n"]})

    m = m[m.level == level]
    reco = m[(m.side == 0) & (m.nu_index >= 0)]
    truth = m[m.side == 1].reset_index()[["entry", "obj", "truth_nu_idx"]]
    best = nested[nested.index.droplevel(-1).isin(reco.index)].xs(0, level="imatch")

    out = reco[["nu_index"]].join(best)
    out = out.reset_index().merge(truth.rename(columns={"obj": "match_obj"}),
                                  on=["entry", "match_obj"], how="left")
    out = out.set_index(WC_NU_LEVELS)
    return pd.DataFrame({
        "tmatch_idx": out.truth_nu_idx.fillna(-1).astype(int),
        "tmatch_iou_q": out.match_iou_q,
        "tmatch_iou_n": out.match_iou_n,
    })


def make_wc_flashbothtpcdf(f):
    """flash_both_tpc per candidate: 1 if its matched flash's group (flashes
    on different TPCs within Trun.flash_pair_dt_us) has a flash in each TPC."""
    fl = make_wc_tabledf(f, "T_flash", ["gid", "tpc", "flash_group"]).reset_index("row", drop=True)
    ntpc = fl.groupby(["entry", "flash_group"]).tpc.nunique().rename("ntpc")
    fl = fl.join(ntpc, on=["entry", "flash_group"]).set_index("gid", append=True)
    tag = make_wc_taggerdf(f, ["matched_flash_gid"])
    both = tag.join(fl.ntpc, on=["entry", "matched_flash_gid"]).ntpc
    return (both >= 2).astype(int).rename("flash_both_tpc").to_frame()


def make_wc_nudf(f, tagger_leaves=wctaggerbranches, kine_leaves=wckinebranches, include_tmatch=None):
    """T_tagger + T_kine + flash_both_tpc (+ the T_match truth match on MC),
    one row per candidate."""
    tag = make_wc_taggerdf(f, tagger_leaves).join(make_wc_flashbothtpcdf(f))
    kine = make_wc_kinedf(f, kine_leaves)
    df = tag.join(kine.drop(columns=[c for c in kine.columns if c in tag.columns]), how="left")
    if include_tmatch is None:
        ismc = bool(loadbranches(f["recTree"], ["rec.hdr.ismc"]).rec.hdr.ismc.iloc[0])
        include_tmatch = ismc and "wc.T_match.level" in _wc_tree(f).keys()
    if include_tmatch:
        df = df.join(make_wc_tmatchdf(f), how="left")
        df["tmatch_idx"] = df.tmatch_idx.fillna(-1).astype(int)
    return df
