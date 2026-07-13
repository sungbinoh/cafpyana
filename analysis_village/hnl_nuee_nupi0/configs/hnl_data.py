from analysis_village.hnl_nuee_nupi0.makedf.make_hnldf import *

GRID_PARAMS = {
    "memory":   "2GB",
    "cpu":      5,
    "disk":     "100GB",
    "lifetime": "2h",
}

DFS =   [make_hnldf_data, make_hdrdf, make_potdf_bnb]
NAMES = ["rec", "hdr", "pot"]
