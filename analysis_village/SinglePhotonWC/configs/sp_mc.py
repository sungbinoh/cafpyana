from analysis_village.SinglePhotonWC.make_wc_candidates import make_wc_candidate_df
from analysis_village.SinglePhotonWC.make_truth import make_sp_truthdf
from analysis_village.SinglePhotonWC.make_wc_events import make_wc_event_df

# wc: one row per Wire-Cell candidate, reco only. mcnu: every true neutrino
# with the single-photon truth. wcevt: every CAF event, Wire-Cell valid flag
# and bundle reasons. POT from histpotdf (dfio.load totals).
DFS = [make_wc_candidate_df, make_sp_truthdf, make_wc_event_df]
NAMES = ["wc", "mcnu", "wcevt"]
