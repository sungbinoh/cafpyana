# MAPLE evt production, data: no mcnu table.
from analysis_village.maple.makedf import *

DFS = [make_maple_evt_nosel_df, make_hdrdf, make_triggerdf, make_potdf_bnb]
NAMES = ["evt", "hdr", "trig", "bnb"]
