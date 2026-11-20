# MAPLE flash production, data: per-event in-time flash summary (flashevt:
# flash_sumpe/flash_maxpe/flash_maxpe_cryo0/1 in [-5, 5] us) for EVERY event,
# independent of the GUMPLE selection -- no slice reco, no truth.
from analysis_village.maple.makedf import *

DFS = [make_maple_flashdf, make_hdrdf, make_triggerdf, make_maple_bnbdf]
NAMES = ["flashevt", "hdr", "trig", "bnb"]
