# MAPLE MC CV truth-only production with systematic weights: NO reco (no evt
# table) and NO preselection -- mcnu, the GUMP CV reweight knob set (wgt) and the
# raw GENIE event record (evtrec) are kept for EVERY neutrino in the file.
# evtrec rows link to mcnu via mcnu.genie_evtrec_idx.
from analysis_village.maple.makedf import *

DFS = [make_maple_nudf, make_maple_rewgtdf, make_maple_evtrec_df, make_hdrdf, make_triggerdf]
NAMES = ["mcnu", "wgt", "evtrec", "hdr", "trig"]
