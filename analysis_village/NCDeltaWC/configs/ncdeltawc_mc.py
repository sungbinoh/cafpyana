from analysis_village.NCDeltaWC.make_ncdeltawc import *
from makedf.makedf import make_hdrdf
from makedf.wirecell import make_wc_evtdf

DFS = [make_ncdeltawc_df_mc, make_hdrdf, make_ncdeltawc_nudf, make_wc_evtdf]
NAMES = ["wc", "hdr", "mcnu", "wcevt"]
