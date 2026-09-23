# MAPLE MC CV production: FULL selection applied at dataframe-build time,
# plus the GUMP CV reweight knob set (wgt) and the raw GENIE event record
# (evtrec), as used by the gump MC-CV reweighting workflow.
# evtrec rows link to mcnu via mcnu.genie_evtrec_idx.
from analysis_village.gump.makedf import *

DFS = [make_gump_evt_fullsel_df, make_gump_nudf, make_gump_rewgtdf,
       make_gump_evtrec_df, make_hdrdf, make_triggerdf, make_potdf_bnb]
NAMES = ["evt", "mcnu", "wgt", "evtrec", "hdr", "trig", "bnb"]
