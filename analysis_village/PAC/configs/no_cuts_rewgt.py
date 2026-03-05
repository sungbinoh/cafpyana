from analysis_village.PAC.makedf import *
from preprocess.preprocess import Script

DFS = [make_pandora_no_cuts_df, make_PAC_nudf, make_stubs, make_crthitdf, make_hdrdf, make_triggerdf, make_potdf_bnb, make_PAC_nuwgtdf]
NAMES = ["evt", "mcnu", "stub", "crt", "hdr", "trig", "bnb", "wgt"]
