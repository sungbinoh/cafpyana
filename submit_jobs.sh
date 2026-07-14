#!/usr/bin/env bash
# Grid job-submit commands for all sample lists ready in data/sample_lists/sbnd/.
# Run from the cafpyana repo root.

# --- MC BNB nominal (full weighted dataframe, main analysis) ---
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcnu_wgt.py -l ./data/sample_lists/sbnd/mcbnb_cv.list -ngrid 2334 -o mcbnb_cv

# --- MC BNB DetVar (CV + light-yield DV variations), hnl_mcnu_detvar.py ---
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
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nupi0_m140.list -ngrid 7 -o mchnl_nupi0_m140
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nupi0_m165.list -ngrid 7 -o mchnl_nupi0_m165
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nupi0_m190.list -ngrid 7 -o mchnl_nupi0_m190
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nupi0_m215.list -ngrid 7 -o mchnl_nupi0_m215
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nupi0_m240.list -ngrid 7 -o mchnl_nupi0_m240
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nupi0_m260.list -ngrid 7 -o mchnl_nupi0_m260

# --- Data (BNB + offbeam) ---
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_data.py -l ./data/sample_lists/sbnd/dtbnb_fix.list -ngrid 9 -o dtbnb_fix
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_data.py -l ./data/sample_lists/sbnd/dtbnb_roll.list -ngrid 7 -o dtbnb_roll
#python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_data.py -l ./data/sample_lists/sbnd/dtoff.list -ngrid 216 -o dtoff
