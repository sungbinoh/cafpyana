# MAPLE evt production, data: PID-free MAPLE PRESELECTION applied at
# dataframe-build time (not the full selection), no mcnu table, no calo
# systematic variations. The full selection is applied post-hoc.
from analysis_village.maple.makedf import *

DFS = [make_maple_evt_presel_data_df, make_hdrdf, make_triggerdf, make_potdf_bnb]
NAMES = ["evt", "hdr", "trig", "bnb"]
