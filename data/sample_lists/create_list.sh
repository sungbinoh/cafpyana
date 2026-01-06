#!/bin/bash

# needs to be run outside the virtual environment
# in an sl7 container
# source /exp/$(id -ng)/data/users/vito/podman/start_SL7dev_jsl.sh
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
setup sbndcode v10_10_03_01 -q e26:prof
htgettoken -a htvaultprod.fnal.gov -i sbnd

# example usage
# ./create_list.sh mc_MCP2025B_5e18_02_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd sbnd/prod_2025B/FallValidationII sbnd

# SAMDEFs for 2025B production
# MC: mc_MCP2025B_5e18_02_prodgenie_corsika_proton_rockbox_sbnd_CV_caf_flat_caf_sbnd
# in time cosmics: mc_MCP2025B_5e18_02_prodcorsika_proton_intime_sbnd_CV_caf_flat_caf_sbnd
# data: data_MCP2025B_02_DevSample_bnblight_v10_06_00_02_flatcaf_sbnd
# in time data cosmics: data_MCP2025B_02_InTimeCosmics_offbeamlight_v10_06_00_02_flatcaf_sbnd

SAMDEF=${1}
OUTDIR=${2}
EXP=${3}

( samweb list-definition-files -e ${EXP} ${SAMDEF} | while read FILENAME;\
do PREFIX='enstore:'; SUFFIX='(.*)'; DIR_TO_FILE=$(samweb -e ${EXP} locate-file ${FILENAME} | sed -e "s@${PREFIX}@@" -e "s@${SUFFIX}@@" -e "s@dcache:@@");\
pnfsToXRootD ${DIR_TO_FILE}/${FILENAME}; done ) 2>&1 | tee ${OUTDIR}/${SAMDEF}.txt
