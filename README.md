# cafpyana

A set of scripts for analyzing SBN CAF files with the python.
Goal is to provide an easy starting point for SBN analysis to everyone.
For more details and instructions, please check [wiki](https://github.com/sungbinoh/cafpyana/wiki).

## Releases

Please check tags for releases!

## Version compatibility

Current `main` branch is based on python v3.9.15.
It is for running the repository without any issue at gpvm servers with the `spack`.
There is no need to open an SL7 image.
If python version is updated in gpvm servers, compatibility issue should be revisited.

# for the dune plotting style, follow the instructions here
https://github.com/DUNE/dune_plot_style.git

(brief description, working in the virtual python environment)
# dependencies first
python3 -m pip install matplotlib numpy scipy

# this is one way to obtain the tarball, but use any way you like
cd /path/to/install/area
export DUNE_PLOT_STYLE_LATEST_TAG=`curl --silent "https://api.github.com/repos/DUNE/dune_plot_style/releases" | jq -r 'map(select(.prerelease == false)) | first | .tag_name'`
wget --no-check-certificate https://github.com/DUNE/dune_plot_style/archive/refs/tags/${DUNE_PLOT_STYLE_LATEST_TAG}.tar.gz -O dune_plot_style.tar.gz
tar -xvzf dune_plot_style.tar.gz

# obviously adjust the directory name for whatever came out of the tarball
cd /path/to/install/area/dune_plot_style
python3 -m pip install .

################################

# to run the pandora df maker and plot 
# distributions using a jupyter notebook
# using all the pandora events

# bnb data
python run_df_maker.py -c configs/pandora.py -l data/sample_lists/sbnd/prod_2025B/data_MCP2025B_02_DevSample_bnblight_v10_06_00_02_flatcaf_sbnd_35ms.txt -o /exp/sbnd/data/users/apapadop/dfs/v10_06_00_02/data_MCP2025B_02_DevSample_bnblight_v10_06_00_02_flatcaf_sbnd_35ms

# mc bnb
python run_df_maker.py -c configs/pandora.py -l data/sample_lists/sbnd/prod_2025B/mc_MCP2025B_5e18_02_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd.txt -o /exp/sbnd/data/users/apapadop/dfs/v10_06_00_02/mc_MCP2025B_5e18_02_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd

# mc intime 
python run_df_maker.py -c configs/pandora.py -l data/sample_lists/sbnd/prod_2025B/mc_MCP2025B_5e18_02_prodcorsika_proton_intime_sbnd_CV_caf_flat_caf_sbnd.txt -o /exp/sbnd/data/users/apapadop/dfs/v10_06_00_02/mc_MCP2025B_5e18_02_prodcorsika_proton_intime_sbnd_CV_caf_flat_caf_sbnd

################################

# to run the cc2p df maker and apply the cc2p selection

#bnb data
python run_df_maker.py -c configs/pandora_cc2p.py -l data/sample_lists/sbnd/prod_2025B/data_MCP2025B_02_DevSample_bnblight_v10_06_00_02_flatcaf_sbnd_35ms.txt -o /exp/sbnd/data/users/apapadop/dfs/v10_06_00_02/data_MCP2025B_02_DevSample_bnblight_v10_06_00_02_flatcaf_sbnd_35ms_cc2p [-nfile 100]

#mc bnb
python run_df_maker.py -c configs/pandora_cc2p.py -l data/sample_lists/sbnd/prod_2025B/mc_MCP2025B_5e18_02_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd.txt -o /exp/sbnd/data/users/apapadop/dfs/v10_06_00_02/mc_MCP2025B_5e18_02_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_cc2p

#mc intime cosmics
python run_df_maker.py -c configs/pandora_cc2p.py -l data/sample_lists/sbnd/prod_2025B/mc_MCP2025B_5e18_02_prodcorsika_proton_intime_sbnd_CV_caf_flat_caf_sbnd.txt -o /exp/sbnd/data/users/apapadop/dfs/v10_06_00_02/mc_MCP2025B_5e18_02_prodcorsika_proton_intime_sbnd_CV_caf_flat_caf_sbnd_100files_cc2p -nfile 100

################################

# to run the cc2p ttree maker

# bnb data
python run_ttree_maker.py -c configs/cc2p_ttree_data.py -i /exp/sbnd/data/users/apapadop/dfs/v10_06_00_02/data_MCP2025B_02_DevSample_bnblight_v10_06_00_02_flatcaf_sbnd_35ms_cc2p.df -o /exp/sbnd/data/users/apapadop/dfs/v10_06_00_02/data_MCP2025B_02_DevSample_bnblight_v10_06_00_02_flatcaf_sbnd_35ms_cc2p.root

# mc bnb
python run_ttree_maker.py -c configs/cc2p_ttree_mc.py -i /exp/sbnd/data/users/apapadop/dfs/v10_06_00_02/mc_MCP2025B_5e18_02_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_cc2p.df -o /exp/sbnd/data/users/apapadop/dfs/v10_06_00_02/mc_MCP2025B_5e18_02_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd_cc2p.root

# mc intime cosmics
python run_ttree_maker.py -c configs/cc2p_ttree_data.py -i /exp/sbnd/data/users/apapadop/dfs/v10_06_00_02/mc_MCP2025B_5e18_02_prodcorsika_proton_intime_sbnd_CV_caf_flat_caf_sbnd_100files_cc2p.df -o /exp/sbnd/data/users/apapadop/dfs/v10_06_00_02/mc_MCP2025B_5e18_02_prodcorsika_proton_intime_sbnd_CV_caf_flat_caf_sbnd_100files_cc2p.root