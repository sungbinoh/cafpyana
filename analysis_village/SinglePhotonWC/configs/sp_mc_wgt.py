from analysis_village.SinglePhotonWC.make_wc_candidates import make_wc_candidate_df
from analysis_village.SinglePhotonWC.make_truth import make_sp_truthdf
from analysis_village.SinglePhotonWC.make_wc_events import make_wc_event_df
from analysis_village.SinglePhotonWC.make_weights import make_sp_wgtdf

# sp_mc.py plus the systematic weights per true neutrino ("wgt", make_weights.py;
# as gump's gump_no_cuts_wgt.py), index as mcnu.
DFS = [make_wc_candidate_df, make_sp_truthdf, make_wc_event_df, make_sp_wgtdf]
NAMES = ["wc", "mcnu", "wcevt", "wgt"]
