#!/usr/bin/env bash
# Grid job-submit commands for all sample lists ready in data/sample_lists/sbnd/.
# Run from the cafpyana repo root.
#
# -ngrid values below are re-tuned from a REAL measurement: the mcbnb_cv
# nominal run (2026_07_06_163827__mcbnb_cv, 4000 jobs / 749262 files, i.e.
# ~187.3 files/job) used only 2.1GB actual against its 4GB GRID_PARAMS
# request -- 52% utilization. That gives a marginal cost of ~11.2 MB/file.
# Every -ngrid below is recomputed to target 90% of each config's requested
# memory (10% safety margin against OOM-kill from file-size variance),
# i.e. files/job = 0.9 * mem_request / 11.2MB, ngrid = ceil(list_len / files/job):
#   hnl_mcnu_wgt.py       cpu=10 mem=4GB lifetime=2h  -> target ~321 files/job (MEASURED basis)
#   hnl_mcmevprtl_wgt.py  cpu=10 mem=4GB lifetime=2h  -> target ~709 files/job (MEASURED basis,
#                         2026-07-07: real mchnl_nupi0_m260 job showed ~1.6GB actual at ~315
#                         files/job, i.e. ~5.1MB/file -- lower than hnl_mcnu_wgt.py's 11.2MB/file
#                         since HNL isn't GENIE-generated, so no 52-knob GENIE decode here.)
#   hnl_mcnu_detvar.py    cpu=5  mem=2GB lifetime=2h  -> target ~800 files/job (MEASURED basis,
#                         2026-07-07: real detvar_cv_0/0p94xly/1p19xly jobs showed ~1.1GB actual,
#                         FLAT across 161-482 files/job -- overhead-dominated, not per-file-scaled.
#                         800 extrapolates ~1.7x past the tested max of 482; re-measure once run.)
#   hnl_data.py           cpu=5  mem=2GB lifetime=2h  -> target ~227 files/job (MEASURED basis,
#                         2026-07-07: real dtbnb_fix job showed ~1.2GB actual at ~151 files/job,
#                         i.e. ~7.9MB/file -- applied to dtbnb_roll/dtoff too, unverified for those.)
# All of the above are single-data-point (linear-through-origin) extrapolations except
# hnl_mcnu_detvar.py's flat/overhead-dominated result -- re-check actual memory again
# once each retuned batch runs for real, since detvar showed the linear assumption can fail.
# 2026-07-07: lifetime bumped 1h->2h across the board (hnl_mcmevprtl_wgt.py was already
# 2h) for extra wall-clock margin given the -ngrid retunes above increase files/job.

# --- MC BNB nominal (full weighted dataframe, main analysis) ---
# ALREADY RAN (see 2026_07_06_163827__mcbnb_cv, 4000/4000 shards complete) --
# ngrid below is the re-tuned value for a future rerun, not yet re-submitted.
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcnu_wgt.py -l ./data/sample_lists/sbnd/mcbnb_cv.list -ngrid 2334 -o mcbnb_cv

# --- MC BNB DetVar (CV + light-yield DV variations), hnl_mcnu_detvar.py ---
# Prior -ngrid 1554 run FAILED (not successful) -- unrelated to the memory
# retune, still unexplained; watch this rerun for the same failure.
#
# -o KEEPS the manual "_0" suffix (detvar_cv_0, not detvar_cv) -- this looks
# redundant with run_df_maker.py's own per-shard "_<jobindex>" indexing
# (shards land as detvar_cv_0_<jobindex>.df) but it's load-bearing:
# concat_hdf.py strips exactly one trailing "_<digits>" group before ".df"
# (the shard index), which only leaves a valid "detvar_cv_0.df" name for
# process_detvars.py's `detvar_(.+)_(\d+)\.df` regex (group 2 = the "_0"
# CV/DV pairing index) if that "_0" was there to begin with. Dropping it
# (tried briefly, reverted) makes concat_hdf.py's output collapse to
# "detvar_cv.df", which doesn't match the regex and process_detvars.py
# hard-exits with "no CV files found". Traced 2026-07-07; do not drop the
# "_0" again without also fixing concat_hdf.py/process_detvars.py to match.
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcnu_detvar.py -l ./data/sample_lists/sbnd/mcbnb_cv.list -ngrid 937 -o detvar_cv_0
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcnu_detvar.py -l ./data/sample_lists/sbnd/mcbnb_0p94xLY.list -ngrid 312 -o detvar_0p94xly_0
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcnu_detvar.py -l ./data/sample_lists/sbnd/mcbnb_1p19xLY.list -ngrid 312 -o detvar_1p19xly_0

# --- MC HNL: nu_ee channel ---
# ngrid retuned to ~709 files/job (see header) -- none of these five submitted yet.
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nuee_m35.list -ngrid 7 -o mchnl_nuee_m35
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nuee_m50.list -ngrid 7 -o mchnl_nuee_m50
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nuee_m75.list -ngrid 7 -o mchnl_nuee_m75
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nuee_m100.list -ngrid 7 -o mchnl_nuee_m100
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nuee_m125.list -ngrid 7 -o mchnl_nuee_m125
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nuee_m140.list -ngrid 7 -o mchnl_nuee_m140

# --- MC HNL: nu_pi0 channel ---
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nupi0_m140.list -ngrid 7 -o mchnl_nupi0_m140
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nupi0_m165.list -ngrid 7 -o mchnl_nupi0_m165
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nupi0_m190.list -ngrid 7 -o mchnl_nupi0_m190
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nupi0_m215.list -ngrid 7 -o mchnl_nupi0_m215
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nupi0_m240.list -ngrid 7 -o mchnl_nupi0_m240
# ALREADY SUBMITTED at -ngrid 15 (running as of 2026-07-07, this is the job that gave
# the ~1.6GB measurement above). ngrid below is the re-tuned value for a future rerun.
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nupi0_m260.list -ngrid 7 -o mchnl_nupi0_m260

# --- Data (BNB + offbeam) ---
# dtbnb_fix and dtoff ALREADY SUBMITTED at their old -ngrid (12 and 306, running as of
# 2026-07-07 -- dtbnb_fix gave the ~1.2GB measurement above). ngrid below is the
# re-tuned value (~227 files/job) for a future rerun. dtbnb_roll was never submitted.
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_data.py -l ./data/sample_lists/sbnd/dtbnb_fix.list -ngrid 9 -o dtbnb_fix
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_data.py -l ./data/sample_lists/sbnd/dtbnb_roll.list -ngrid 7 -o dtbnb_roll
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_data.py -l ./data/sample_lists/sbnd/dtoff.list -ngrid 216 -o dtoff
