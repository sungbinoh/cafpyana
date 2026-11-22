from analysis_village.SinglePhotonWC.make_wc_candidates import make_wc_candidate_df, make_wc_trkseg_df
from analysis_village.SinglePhotonWC.make_truth import make_sp_truthdf
from analysis_village.SinglePhotonWC.make_wc_events import make_wc_event_df
from analysis_village.SinglePhotonWC.make_weights import make_sp_wgtdf

# wc: one row per Wire-Cell candidate, reco only. mcnu: every true neutrino
# with the single-photon truth. wcevt: every CAF event, Wire-Cell valid flag
# and bundle reasons. wgt: systematic weights per true neutrino
# (make_weights.py; as gump's gump_no_cuts_wgt.py), index as mcnu.
# trkseg: dQ/dx along the tracks, 4 cm pieces (make_wc_trkseg_df; calibration checks).
# POT from histpotdf (dfio.load totals).
DFS = [make_wc_candidate_df, make_sp_truthdf, make_wc_event_df, make_sp_wgtdf, make_wc_trkseg_df]
NAMES = ["wc", "mcnu", "wcevt", "wgt", "trkseg"]
