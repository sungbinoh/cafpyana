from analysis_village.SinglePhotonWC.make_wc_candidates import make_wc_candidate_df, make_wc_trkseg_df
from analysis_village.SinglePhotonWC.make_pandora_hits import make_pandora_trkhit_df, make_pandora_hittrain_df
from analysis_village.SinglePhotonWC.make_wc_events import make_wc_event_df, make_hdrdf, make_wc_chanrun_df

# Data (on- or off-beam): the sp_mc.py tables without the truth table, plus the
# header (pot, noffbeambnb) for the POT / off-beam gate normalisation, and the
# trkseg table (dQ/dx along the tracks), the pandtrk table (Pandora track hits) and
# the pandtrain table (Pandora hit trains on one wire) and the wcchanrun table
# (long single-channel Wire-Cell signals).
DFS = [make_wc_candidate_df, make_wc_event_df, make_hdrdf, make_wc_trkseg_df, make_pandora_trkhit_df, make_pandora_hittrain_df, make_wc_chanrun_df]
NAMES = ["wc", "wcevt", "hdr", "trkseg", "pandtrk", "pandtrain", "wcchanrun"]
