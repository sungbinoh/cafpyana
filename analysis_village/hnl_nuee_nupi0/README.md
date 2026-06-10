#EXAMPLE COMMANDS:

#Make dataframes of MC MeVPrtl without any weights
#python run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl.py -l input.list -o output.df

#Make dataframes of MC MeVPrtl with BNB flux and G4 weights
#python run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l input.list -o output.df

#Make dataframes of MC SM neutrinos without any weights
#python run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcnu.py -l input.list -o output.df

#Make dataframes of MC SM neutrinos with GENIE/BNB flux/G4 weights
#python run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcnu_wgt.py -l input.list -o output.df

#Make dataframes of data
#python run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_data.py -l input.list -o output.df
