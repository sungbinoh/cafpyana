"""Every scalar T_tagger leaf per Wire-Cell candidate (~1040, index
(entry, nu_index)): the inputs of Wire-Cell's own BDTs, for the BDT here.

  shw_sp_*   single-photon tagger (vertex activity, shower, proton/MIP tracks)
  pio_*      pi0 tagger
  numu_cc_*, numu_*_score       numu CC tagger
  cosmict_*  cosmic sub-taggers
  mip_*, gap_*, br1-4_*, stem_*, vis_*, spt_*, hol_*, lol_*, mgo_*, mgt_*,
  stw_*, anc_*, lem_*, tro_*, cme_*, brm_*, sig_*:  nue tagger inputs
  ssm_*, ssmsp_*                single short-track (kDAR / muon) tagger

The nested (vector) leaves are left out. Identifiers are dropped here.
"""

from makedf import wirecell as wc

DROP = ["run", "subrun", "event", "nu_index", "cluster_id", "matched_flash_gid"]


def make_wc_tagger_df(f):
    leaves = wc.wc_tagger_leaves(f)
    df = wc.make_wc_taggerdf(f, leaves)
    return df.drop(columns=[c for c in DROP if c in df.columns]).add_prefix("tag_")
