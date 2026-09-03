# MAPLE MC CV production: PID-free MAPLE PRESELECTION applied at dataframe-build
# time (sanity + FV + all-track containment + cathode + a muon and >=1
# candidate proton; no CRT veto or cryo-light), NOT the full
# selection -- the full/calo-varied selection is applied post-hoc from the
# preselected slices (see maple_sel.py). Also keeps the GUMP CV reweight knob
# set (wgt) and the raw GENIE event record (evtrec), as used by the gump
# MC-CV reweighting workflow. evtrec rows link to mcnu via mcnu.genie_evtrec_idx.
from analysis_village.maple.makedf import *

DFS = [make_maple_evt_presel_df, make_maple_nudf, make_maple_rewgtdf,
       make_maple_evtrec_df, make_hdrdf, make_triggerdf, make_potdf_bnb]
NAMES = ["evt", "mcnu", "wgt", "evtrec", "hdr", "trig", "bnb"]
