# MAPLE evt production, MC: no selection applied (all slices kept, with
# per-cut booleans). PID candidates from the cafpyana-recomputed (gump-style) chi2.
from analysis_village.maple.makedf import *

DFS = [make_maple_evt_nosel_df, make_maple_nudf, make_hdrdf, make_triggerdf, make_potdf_bnb]
NAMES = ["evt", "mcnu", "hdr", "trig", "bnb"]
