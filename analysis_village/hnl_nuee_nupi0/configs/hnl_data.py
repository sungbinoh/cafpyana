from analysis_village.hnl_nuee_nupi0.makedf.make_hnldf import *

# Lightest of the four HNL configs: no truth merge, no weight decode at all.
# MEASURED 2026-07-07 (dtbnb_fix real grid job): actual peak ~1.2GB at ~151
# files/job (60% of the 2GB request) -> ~7.9MB/file (linear-through-origin,
# one data point). Applied to dtbnb_roll/dtoff too (same config, unverified
# for those specific lists). submit_jobs.sh's -ngrid retuned to target ~227
# files/job (90% of 2GB); re-measure once the retuned jobs run.
GRID_PARAMS = {
    "memory":   "2GB",
    "cpu":      5,
    "disk":     "100GB",
    "lifetime": "1h",
}

DFS =   [make_hnldf_data, make_hdrdf, make_potdf_bnb]
NAMES = ["rec", "hdr", "pot"]
