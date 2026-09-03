# MAPLE evt production, MC: apply the FULL MAPLE selection during dataframe
# building (1mu + >1p, no pions/showers/other, cafpyana-PID (gump-style) candidates).
from analysis_village.maple.makedf import *

DFS = [make_maple_evt_fullsel_df, make_maple_nudf, make_hdrdf, make_triggerdf, make_potdf_bnb]
NAMES = ["evt", "mcnu", "hdr", "trig", "bnb"]
