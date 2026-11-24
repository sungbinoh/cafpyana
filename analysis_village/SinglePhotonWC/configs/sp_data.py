from analysis_village.SinglePhotonWC.make_wc_candidates import make_wc_candidate_df
from analysis_village.SinglePhotonWC.numu_bdt import make_wc_numuobj_df
from analysis_village.SinglePhotonWC.make_wc_events import make_wc_event_df, make_hdrdf

# Data (on- or off-beam): the sp_mc.py tables without the truth table, plus the
# header (pot, noffbeambnb) for the POT / off-beam gate normalisation.
# numuobj as in sp_mc.py.
DFS = [make_wc_candidate_df, make_wc_event_df, make_hdrdf, make_wc_numuobj_df]
NAMES = ["wc", "wcevt", "hdr", "numuobj"]
