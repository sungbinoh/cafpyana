from analysis_village.gump.makedf import *
from preprocess.preprocess import Script

PREPROCESS = [Script("/exp/sbnd/app/users/gputnam/Ar23-knobs/update_reweight_anywhere.sh")]
PREPROCESS = []

DFS = [make_pandora_no_cuts_df, make_gump_nudf, make_stubs, make_crthitdf, make_hdrdf, make_triggerdf, make_potdf_bnb, make_gump_nurewgtdf, make_opflashdf]
NAMES = ["evt", "mcnu", "stub", "crt", "hdr", "trig", "bnb", "wgt", "flash"]
