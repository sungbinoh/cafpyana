# Setup grid submission
echo "BEARER_TOKEN_FILE is set to: $BEARER_TOKEN_FILE"
outDir=$1
echo "@@ outDir : ${outDir}"
DFPREFIX=$2

nProcess=$PROCESS
echo "@@ nProcess : "${nProcess}
source /cvmfs/fermilab.opensciencegrid.org/packages/common/setup-env.sh
spack load sam-web-client@3.6%gcc@11.4.1

echo "@@ pwd"
pwd
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
echo "@@ check if there is cmake"
spack find cmake
spack load cmake@3.27.7
which cmake
echo "@@ check if other spack packages"

source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
spack load hdf5@1.12.2%gcc@12.5.0 arch=linux-almalinux9-x86_64_v2
spack load xrootd@5.6.9%gcc@12.5.0 arch=linux-almalinux9-x86_64_v2
spack load ifdhc@2.7.2
echo "@@ run init_grid.sh"
source ./bin/init_grid.sh
echo "@@ ls -alh"
ls -alh
echo "@@ mkdir output"
mkdir output
echo "@@ Done!"
thisOutputCreationDir=`pwd`
filesFromSender=${CONDOR_DIR_INPUT}/bin_dir/
echo "@@ check filesFromSender dir"
ls -alh ${filesFromSender}

echo "@@ Setup xrootd"
cp -r ${filesFromSender}/XRootD $VIRTUAL_ENV/lib/python3.9/site-packages/
cp -r ${filesFromSender}/pyxrootd $VIRTUAL_ENV/lib/python3.9/site-packages/
export LD_LIBRARY_PATH=$VIRTUAL_ENV/lib/python3.9/site-packages/pyxrootd:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=$VIRTUAL_ENV/lib/python3.9/site-packages/xrootd-5.6.1-py3.9-linux-x86_64.egg/pyxrootd:$LD_LIBRARY_PATH

export IFDH_CP_MAXRETRIES=2

echo "@@ outDir : "${outDir}
echo "@@ ifdh  mkdir_p "${outDir}
ifdh  mkdir_p ${outDir}

echo "@@ source ${filesFromSender}/run_"${nProcess}".sh "
ls -alh
pwd

#htgettoken -a htvaultprod.fnal.gov -i sbnd

cp ${filesFromSender}/run_${nProcess}.sh ./
source run_${nProcess}.sh  &> log_${nProcess}.log
ls -alh
echo "@@ Check output : ${DFPREFIX}_${nProcess}.df"
ls -alh ${DFPREFIX}_${nProcess}.df

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

echo "BEARER_TOKEN_FILE is set to: $BEARER_TOKEN_FILE"
