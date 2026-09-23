# MAPLE production with systematic weights: adds the wgt dataframe
# (mcnu + BNB flux + GENIE multisim/multisigma weights).
from analysis_village.gump.makedf import *

DFS = [make_gump_evt_nosel_df, make_gump_nudf, make_gump_wgtdf, make_hdrdf, make_triggerdf, make_potdf_bnb]
NAMES = ["evt", "mcnu", "wgt", "hdr", "trig", "bnb"]
