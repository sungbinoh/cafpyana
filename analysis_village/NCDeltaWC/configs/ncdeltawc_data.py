from analysis_village.NCDeltaWC.make_ncdeltawc import *
from makedf.makedf import make_hdrdf, make_potdf_bnb
from makedf.wirecell import make_wc_evtdf

DFS = [make_ncdeltawc_df_data, make_hdrdf, make_potdf_bnb, make_wc_evtdf]
NAMES = ["wc", "hdr", "pot", "wcevt"]
