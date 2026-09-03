# GUMP evt production, MC: apply the PID-free MAPLE preselection during
# dataframe building (sanity + FV + all-track containment + cathode + a muon
# and >=1 candidate proton; no CRT veto or cryo-light).
from analysis_village.gump.makedf import *

DFS = [make_gump_evt_presel_df, make_gump_nudf, make_hdrdf, make_triggerdf, make_potdf_bnb]
NAMES = ["evt", "mcnu", "hdr", "trig", "bnb"]
