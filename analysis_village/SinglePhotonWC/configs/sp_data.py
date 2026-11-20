from analysis_village.SinglePhotonWC.make_wc_candidates import make_wc_candidate_df
from analysis_village.SinglePhotonWC.make_wc_events import make_wc_event_df

# Data (on- or off-beam): the sp_mc.py tables without the truth table.
DFS = [make_wc_candidate_df, make_wc_event_df]
NAMES = ["wc", "wcevt"]
