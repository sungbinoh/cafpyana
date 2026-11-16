# MAPLE MC truth-only production: NO reco (no evt table) and NO preselection --
# the mcnu table is kept for EVERY neutrino in the file (+ hdr/trig for POT).
from analysis_village.maple.makedf import *

DFS = [make_maple_nudf, make_hdrdf, make_triggerdf]
NAMES = ["mcnu", "hdr", "trig"]
