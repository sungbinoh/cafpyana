from analysis_village.SinglePhotonWC.make_wc_candidates import make_wc_candidate_df, make_wc_trkseg_df
from analysis_village.SinglePhotonWC.make_pandora_hits import make_pandora_trkhit_df
from analysis_village.SinglePhotonWC.make_wc_events import make_wc_event_df, make_hdrdf

# Data (on- or off-beam): the sp_mc.py tables without the truth table, plus the
# header (pot, noffbeambnb) for the POT / off-beam gate normalisation, and the
# trkseg table (dQ/dx along the tracks) and the pandtrk table (Pandora track hits).
DFS = [make_wc_candidate_df, make_wc_event_df, make_hdrdf, make_wc_trkseg_df, make_pandora_trkhit_df]
NAMES = ["wc", "wcevt", "hdr", "trkseg", "pandtrk"]
