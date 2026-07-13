from analysis_village.hnl_nuee_nupi0.makedf.make_hnldf import *

GRID_PARAMS = {
    "memory":   "2GB",
    "cpu":      5,
    "disk":     "100GB",
    "lifetime": "2h",
}

DFS =   [make_hnldf_mcnu, make_hdrdf, make_mcnulite_df_hnl]
NAMES = ["rec", "hdr", "nulite"]
