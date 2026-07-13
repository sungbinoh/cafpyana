from analysis_village.hnl_nuee_nupi0.makedf.make_hnldf import * 

GRID_PARAMS = {
    "memory":   "4GB",
    "cpu":      10,
    "disk":     "100GB",
    "lifetime": "2h",
}

DFS =   [make_hnldf_mcnu_wgt, make_hdrdf, make_potdf_bnb, make_mcnudf]
NAMES = ["rec", "hdr", "pot", "mcnu"]


