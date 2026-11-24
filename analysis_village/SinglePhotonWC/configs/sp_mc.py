from analysis_village.SinglePhotonWC.make_wc_candidates import make_wc_candidate_df, make_wc_trkseg_df
from analysis_village.SinglePhotonWC.make_truth import make_sp_truthdf
from analysis_village.SinglePhotonWC.make_pandora_hits import make_pandora_trkhit_df, make_pandora_hittrain_df, make_pandora_wire_df
from analysis_village.SinglePhotonWC.make_wc_events import make_wc_event_df, make_wc_chanrun_df
from analysis_village.SinglePhotonWC.make_weights import make_sp_wgtdf

# wc: one row per Wire-Cell candidate, reco only. mcnu: every true neutrino
# with the single-photon truth. wcevt: every CAF event, Wire-Cell valid flag
# and bundle reasons. wgt: systematic weights per true neutrino
# (make_weights.py; as gump's gump_no_cuts_wgt.py), index as mcnu.
# trkseg: dQ/dx along the tracks, 4 cm pieces (make_wc_trkseg_df; calibration checks).
# pandtrk: Pandora track hits per plane, summarised per track (make_pandora_hits.py).
# pandtrain: Pandora hit trains, consecutive hits of a track on one wire (make_pandora_hits.py).
# pandwire: the charge of one Pandora track on one wire, its calo points (ROI snippets) summed (make_pandora_hits.py).
# wcchanrun: long single-channel signals from wc.T_proj_data, measured vs predicted (make_wc_events.py).
# POT from histpotdf (dfio.load totals).
DFS = [make_wc_candidate_df, make_sp_truthdf, make_wc_event_df, make_sp_wgtdf, make_wc_trkseg_df, make_pandora_trkhit_df, make_pandora_hittrain_df, make_pandora_wire_df, make_wc_chanrun_df]
NAMES = ["wc", "mcnu", "wcevt", "wgt", "trkseg", "pandtrk", "pandtrain", "pandwire", "wcchanrun"]
