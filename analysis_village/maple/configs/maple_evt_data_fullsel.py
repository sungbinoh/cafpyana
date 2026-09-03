# MAPLE evt production, data: FULL selection applied at dataframe-build
# time, no mcnu table, no calo systematic variations.
from analysis_village.maple.makedf import *

DFS = [make_maple_evt_fullsel_data_df, make_hdrdf, make_triggerdf, make_potdf_bnb]
NAMES = ["evt", "hdr", "trig", "bnb"]
