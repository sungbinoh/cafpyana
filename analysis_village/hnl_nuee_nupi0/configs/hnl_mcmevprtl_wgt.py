from analysis_village.hnl_nuee_nupi0.makedf.make_hnldf import *

# Explicit GRID_PARAMS: previously absent, which fell back to run_df_maker.py's
# DEFAULT_GRID_PARAMS (cpu=7, lifetime=1h). Measured ~5.7s/file single-threaded;
# at cpu=7 a single ~4700-file HNL mass-point sample runs right at/over the 1h
# wall, so lifetime is bumped here for margin.
# MEASURED 2026-07-07 (mchnl_nupi0_m260 real grid job): actual peak ~1.6GB at
# ~315 files/job (40% of the 4GB request) -- notably lower than hnl_mcnu_wgt.py's
# 11.2MB/file rate (mcbnb_cv decodes the full 52 GENIE + 13 Flux + 6 G4 knobs;
# HNL isn't GENIE-generated, so this config only decodes Flux+G4, ~5.1MB/file).
# submit_jobs.sh's -ngrid retuned off this rate (~709 files/job, 90% of 4GB);
# only one data point exists so linear-through-origin is assumed, unverified --
# re-measure once the retuned jobs run (detvar showed this assumption can be
# wrong, coming out flat/overhead-dominated instead).
GRID_PARAMS = {
    "memory":   "4GB",
    "cpu":      10,
    "disk":     "100GB",
    "lifetime": "2h",
}

DFS =   [make_hnldf_mevprtl_wgt, make_hdrdf, make_potdf_bnb, make_mevprtldf]
NAMES = ["rec", "hdr", "pot", "mevprtl"]
