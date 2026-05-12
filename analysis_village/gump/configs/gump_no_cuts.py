from analysis_village.gump.makedf import *

DFS = [make_pandora_no_cuts_df, make_gump_nudf, make_stubs, 
    # make_crthitdf, 
    make_hdrdf, make_triggerdf, make_potdf_bnb, make_opflashdf]

NAMES = ["evt", "mcnu", "stub", 
    # "crt", 
    "hdr", "trig", "bnb", "flash"]
