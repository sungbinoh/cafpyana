# Setup grid submission
echo "@@ initial test:"
xrdcp root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/sbn/data_add/sbn_nd/poms_production/mc/MCP2025C_1e20_v10_06_00_09/v10_06_00_09/prodgenie_corsika_proton_rockbox_sbnd/CV/caf/d2/64/caf.flat.caf-3750ae6a-8d07-44db-b813-1a16c0407cf4.root .

outDir=$1
echo "@@ outDir : ${outDir}"
DFPREFIX=$2

nProcess=$PROCESS
echo "@@ nProcess : "${nProcess}

source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh

echo "@@ ls -alh"
ls -alh
echo "@@ git clone cafpyana"
git clone https://github.com/gputnam/cafpyana.git
echo "@@ cd to cafpyana dir"
cd cafpyana
git checkout remotes/origin/N8Dev
echo "@@ git branch -a"
echo "@@ ls -alh"
ls -alh
echo "@@ second test:"
xrdcp root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/sbn/data_add/sbn_nd/poms_production/mc/MCP2025C_1e20_v10_06_00_09/v10_06_00_09/prodgenie_corsika_proton_rockbox_sbnd/CV/caf/d2/64/caf.flat.caf-3750ae6a-8d07-44db-b813-1a16c0407cf4.root second_test.root
ls -alh

# spack load scitokens-cpp@1.0.1
spack load cmake@3.27.7
# spack load hdf5@1.14.3
# spack load xrootd@5.6.1
spack load ifdhc@2.7.2
# spack find --loaded

echo "@@ third test:"
xrdcp root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/sbn/data_add/sbn_nd/poms_production/mc/MCP2025C_1e20_v10_06_00_09/v10_06_00_09/prodgenie_corsika_proton_rockbox_sbnd/CV/caf/d2/64/caf.flat.caf-3750ae6a-8d07-44db-b813-1a16c0407cf4.root third_test.root
ls -alh

thisOutputCreationDir=`pwd`
filesFromSender=${CONDOR_DIR_INPUT}/bin_dir/

#echo "@@ first attempt!"
#echo "@@ source ${filesFromSender}/run_"${nProcess}".sh "
#source run_${nProcess}.sh  &> log_${nProcess}_first.log

echo "@@ run init_grid.sh"
source ./bin/init_grid.sh
echo "@@ ls -alh"
ls -alh

echo "@@ mkdir output"
mkdir output
echo "@@ Done!"
echo "@@ check filesFromSender dir"
ls -alh ${filesFromSender}

#echo "@@ Setup xrootd"
cp -r ${filesFromSender}/XRootD $VIRTUAL_ENV/lib/python3.9/site-packages/
cp -r ${filesFromSender}/pyxrootd $VIRTUAL_ENV/lib/python3.9/site-packages/

echo "@@ new fourth test:"
xrdcp root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/sbn/data_add/sbn_nd/poms_production/mc/MCP2025C_1e20_v10_06_00_09/v10_06_00_09/prodgenie_corsika_proton_rockbox_sbnd/CV/caf/d2/64/caf.flat.caf-3750ae6a-8d07-44db-b813-1a16c0407cf4.root new_fourth_test.root
ls -alh

#export LD_LIBRARY_PATH=$VIRTUAL_ENV/lib/python3.9/site-packages/pyxrootd:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=$VIRTUAL_ENV/lib/python3.9/site-packages/xrootd-5.8.3-py3.10-linux-x86_64.egg/pyxrootd:$LD_LIBRARY_PATH
#
#echo "@@ checking venv"
#ls $VIRTUAL_ENV
#ls $VIRTUAL_ENV/lib
#ls $VIRTUAL_ENV/lib/python3.10
#ls $VIRTUAL_ENV/lib/python3.10/site-packages
#ls $VIRTUAL_ENV/lib/python3.10/site-packages/xrootd-5.8.3-py3.10-linux-x86_64.egg
#ls $VIRTUAL_ENV/lib/python3.10/site-packages/xrootd-5.8.3-py3.10-linux-x86_64.egg/pyxrootd

echo "@@ fifth test:"
xrdcp root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/sbn/data_add/sbn_nd/poms_production/mc/MCP2025C_1e20_v10_06_00_09/v10_06_00_09/prodgenie_corsika_proton_rockbox_sbnd/CV/caf/d2/64/caf.flat.caf-3750ae6a-8d07-44db-b813-1a16c0407cf4.root fifth_test.root
ls -alh

export IFDH_CP_MAXRETRIES=2

echo "@@ outDir : "${outDir}
echo "@@ ifdh  mkdir_p "${outDir}
ifdh  mkdir_p ${outDir} || true

echo "@@ source ${filesFromSender}/run_"${nProcess}".sh "
ls -ltr

cp ${filesFromSender}/run_${nProcess}.sh ./
source run_${nProcess}.sh  &> log_${nProcess}.log
ls -alh
echo "@@ second attempt!"
echo "@@ Check output : ${DFPREFIX}_${nProcess}.df"
ls -alh ${DFPREFIX}_${nProcess}.df
pip list
outFILE=${thisOutputCreationDir}/${DFPREFIX}_${nProcess}.df
if [ -f "$outFILE" ]; then
  echo "ifdh cp ${thisOutputCreationDir}/${DFPREFIX}_${nProcess}.df ${outDir}/${DFPREFIX}_${nProcess}.df"
  ifdh cp ${thisOutputCreationDir}/${DFPREFIX}_${nProcess}.df ${outDir}/${DFPREFIX}_${nProcess}.df
  echo "ifdh cp ${thisOutputCreationDir}/log_${nProcess}.log ${outDir}/log_${nProcess}.log"
  ifdh cp ${thisOutputCreationDir}/log_${nProcess}.log ${outDir}/log_${nProcess}.log
  echo "@@ Done!"
else
  ifdh cp ${thisOutputCreationDir}/log_${nProcess}.log ${outDir}/log_${nProcess}.log
  echo "File not exist"
fi
