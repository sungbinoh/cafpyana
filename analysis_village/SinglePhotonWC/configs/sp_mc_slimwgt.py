from analysis_village.SinglePhotonWC.make_wc_candidates import make_wc_candidate_df, make_wc_trkseg_df
from analysis_village.SinglePhotonWC.make_truth import make_sp_truthdf
from analysis_village.SinglePhotonWC.make_pandora_hits import make_pandora_trkhit_df, make_pandora_hittrain_df
from analysis_village.SinglePhotonWC.make_wc_events import make_wc_event_df, make_wc_chanrun_df
from analysis_village.SinglePhotonWC.make_weights import make_sp_slimwgtdf

# sp_mc.py with the slim systematic weights ("slimwgt": products per group,
# make_weights.py; as gump's gump_no_cuts_slimwgt.py) in place of the full ones.
DFS = [make_wc_candidate_df, make_sp_truthdf, make_wc_event_df, make_sp_slimwgtdf, make_wc_trkseg_df, make_pandora_trkhit_df, make_pandora_hittrain_df, make_wc_chanrun_df]
NAMES = ["wc", "mcnu", "wcevt", "slimwgt", "trkseg", "pandtrk", "pandtrain", "wcchanrun"]
