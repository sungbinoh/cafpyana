python3 run_df_maker.py -c analysis_village/mcs_pid/configs/calo_study_mc.py -l data/sample_lists/gen2_cv.list -ncpu 30 -nfile 8000 -o /exp/sbnd/data/users/nrowe/calo_mc_dbg.df
python3 run_df_maker.py -c analysis_village/mcs_pid/configs/calo_study.py -l data/sample_lists/gen2_dev.list -ncpu 30 -o /exp/sbnd/data/users/nrowe/calo_data_dbg.df
