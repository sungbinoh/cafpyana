from analysis_village.SinglePhotonWC.make_wc_candidates import make_wc_candidate_df
from analysis_village.SinglePhotonWC.make_wc_events import make_wc_event_df, make_hdrdf

# Data (on- or off-beam): the sp_mc.py tables without the truth table, plus the
# header (pot, noffbeambnb) for the POT / off-beam gate normalisation.
DFS = [make_wc_candidate_df, make_wc_event_df, make_hdrdf]
NAMES = ["wc", "wcevt", "hdr"]
