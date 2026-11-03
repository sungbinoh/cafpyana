# MAPLE production with systematic weights: adds the wgt dataframe
# (mcnu + BNB flux + GENIE multisim/multisigma weights).
from analysis_village.maple.makedf import *

DFS = [make_maple_evt_nosel_df, make_maple_nudf, make_maple_wgtdf, make_hdrdf, make_triggerdf, make_maple_bnbdf]
NAMES = ["evt", "mcnu", "wgt", "hdr", "trig", "bnb"]
