from analysis_village.gump.makedf import *
from preprocess.preprocess import Script

DFS = [make_pandora_no_cuts_df, make_gump_nudf, make_stubs, make_crthitdf, make_hdrdf, make_triggerdf, make_potdf_bnb, make_gump_nuwgtdf]
NAMES = ["evt", "mcnu", "stub", "crt", "hdr", "trig", "bnb", "wgt"]
