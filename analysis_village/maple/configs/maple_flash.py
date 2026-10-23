# MAPLE flash production, MC (non-CV / detector-variation): per-event in-time
# flash summary (flashevt: flash_sumpe/flash_maxpe/flash_maxpe_cryo0/1 in
# [-5, 5] us) for EVERY event, independent of the GUMPLE selection -- no slice
# reco. Keeps mcnu.
from analysis_village.maple.makedf import *

DFS = [make_maple_flashdf, make_maple_nudf, make_hdrdf, make_triggerdf, make_maple_bnbdf]
NAMES = ["flashevt", "mcnu", "hdr", "trig", "bnb"]
