# Setup grid submission
outDir=$1
echo "@@ outDir : ${outDir}"
DFPREFIX=$2

nProcess=$PROCESS
echo "@@ nProcess : "${nProcess}

source /cvmfs/larsoft.opensciencegrid.org/spack-fnal-v1.0.0/setup-env.sh

echo "@@ ls -alh"
ls -alh

spack load cmake@3.31.8
spack load ifdhc@2.8.0

spack load hdf5@1.12.2%gcc@12.5.0 arch=linux-almalinux9-x86_64_v2
spack load xrootd@5.6.9%gcc@12.5.0 arch=linux-almalinux9-x86_64_v2
# this sets Python to 3.10.16
spack load root arch=linux-almalinux9-x86_64_v2

thisOutputCreationDir=`pwd`
filesFromSender=${CONDOR_DIR_INPUT}/bin_dir/

echo "@@ copy cafpyana.tar.gz and extract"
cp ${filesFromSender}/cafpyana.tar.gz ./
mkdir cafpyana
mv cafpyana.tar.gz cafpyana
cd cafpyana
tar -xvf cafpyana.tar.gz
echo "@@ ls -alh"
ls -alh

echo "@@ run init_grid.sh"
source ./bin/init_grid.sh
echo "@@ ls -alh"
ls -alh

echo "@@ mkdir output"
mkdir output
echo "@@ Done!"
echo "@@ check filesFromSender dir"
ls -alh ${filesFromSender}

echo "@@ Setup xrootd"
echo $VIRTUAL_ENV
ls envs/
ls envs/venv_py310_cafpyana/
ls envs/venv_py310_cafpyana/lib/
ls envs/venv_py310_cafpyana/lib/python3.10
ls envs/venv_py310_cafpyana/lib/python3.10/site-packages/
cp -r ${filesFromSender}/XRootD $VIRTUAL_ENV/lib/python3.10/site-packages/
cp -r ${filesFromSender}/pyxrootd $VIRTUAL_ENV/lib/python3.10/site-packages/

export IFDH_CP_MAXRETRIES=2

echo "@@ outDir : "${outDir}
echo "@@ ifdh  mkdir_p "${outDir}
ifdh mkdir_p ${outDir} || true

echo "@@ source ${filesFromSender}/run_"${nProcess}".sh "
ls -ltr

cp ${filesFromSender}/run_${nProcess}.sh ./
source run_${nProcess}.sh  &> log_${nProcess}.log
ls -ltr

echo "@@ print log"
more log_${nProcess}.log
echo "@@ Check output : ${DFPREFIX}_${nProcess}.df"
ls -alh ${DFPREFIX}_${nProcess}.df

outFILE=${thisOutputCreationDir}/cafpyana/${DFPREFIX}_${nProcess}.df
if [ -f "$outFILE" ]; then
  echo "ifdh cp ${thisOutputCreationDir}/cafpyana/${DFPREFIX}_${nProcess}.df ${outDir}/${DFPREFIX}_${nProcess}.df"
  ifdh cp ${thisOutputCreationDir}/cafpyana/${DFPREFIX}_${nProcess}.df ${outDir}/${DFPREFIX}_${nProcess}.df
  echo "ifdh cp ${thisOutputCreationDir}/cafpyana/log_${nProcess}.log ${outDir}/log_${nProcess}.log"
  ifdh cp ${thisOutputCreationDir}/cafpyana/log_${nProcess}.log ${outDir}/log_${nProcess}.log
  echo "@@ Done!"
else
  ifdh cp ${thisOutputCreationDir}/log_${nProcess}.log ${outDir}/log_${nProcess}.log
  echo "File not exist"
fi
