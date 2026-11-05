from analysis_village.SinglePhotonWC.make_wc_candidates import make_wc_candidate_df
from analysis_village.SinglePhotonWC.make_truth import make_sp_truthdf

# wc: one row per Wire-Cell candidate, reco only. mcnu: every true neutrino
# with the single-photon truth. POT from histpotdf (dfio.load totals).
DFS = [make_wc_candidate_df, make_sp_truthdf]
NAMES = ["wc", "mcnu"]
