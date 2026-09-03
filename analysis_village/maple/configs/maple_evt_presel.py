# MAPLE evt production, MC: apply the PID-free MAPLE preselection during
# dataframe building (sanity + FV + all-track containment + cathode + a muon
# and >=1 candidate proton; no CRT veto or cryo-light).
# Calorimetric variations are NOT computed here -- they are only needed for the
# CV reweighting workflow (see maple_evtrec_presel.py) and substantially slow
# down processing, so the non-CV / detector-variation MC samples skip them.
from analysis_village.maple.makedf import *

DFS = [make_maple_evt_presel_nocalo_df, make_maple_nudf, make_hdrdf, make_triggerdf, make_potdf_bnb]
NAMES = ["evt", "mcnu", "hdr", "trig", "bnb"]
