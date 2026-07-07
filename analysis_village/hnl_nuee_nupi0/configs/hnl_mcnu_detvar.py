from analysis_village.hnl_nuee_nupi0.makedf.make_hnldf import *

# Detector-variation CV/DV samples: no systematic weights needed here, just
# nominal reco compared across variations, plus the lightweight "nulite"
# table used for cross-file (CV<->DV) event matching in detvar/store.py.
# cpu/memory set lower than hnl_mcnu_wgt.py: this config never runs
# bnbsyst/geniesyst/g4syst (no weight-branch decode, no per-knob columns),
# which is what drove the wgt config's footprint.
# MEASURED 2026-07-07 (detvar_cv_0/detvar_0p94xly_0/detvar_1p19xly_0 real grid
# jobs): actual usage ~1.1GB, FLAT across 161-482 files/job (3x range, no
# measurable growth) -- memory here is dominated by fixed overhead, not
# per-file cost, unlike hnl_mcnu_wgt.py. 2GB request left unchanged (still
# ample margin even at the tested max); -ngrid in submit_jobs.sh was bumped
# to target ~800 files/job -- this extrapolates ~1.7x past the tested max
# (482), so re-check actual memory again once run to see if growth kicks in
# above the tested range before pushing further.
GRID_PARAMS = {
    "memory":   "2GB",
    "cpu":      5,
    "disk":     "100GB",
    "lifetime": "1h",
}

DFS =   [make_hnldf_mcnu, make_hdrdf, make_mcnulite_df_hnl]
NAMES = ["rec", "hdr", "nulite"]
