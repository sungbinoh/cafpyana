# Setup grid submission
outDir=$1
echo "@@ outDir : ${outDir}"
DFPREFIX=$2

nProcess=$PROCESS
echo "@@ nProcess : "${nProcess}

source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh

echo "@@ ls -alh"
ls -alh
echo "@@ git clone cafpyana"
git clone https://github.com/sungbinoh/cafpyana.git
echo "@@ cd to cafpyana dir"
cd cafpyana
git checkout remotes/origin/grid
echo "@@ git branch -a"
echo "@@ ls -alh"
ls -alh

# spack load scitokens-cpp@1.0.1
spack load cmake@3.27.7
# spack load hdf5@1.14.3
# spack load xrootd@5.6.1
spack load ifdhc@2.7.2
# spack find --loaded

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

echo "@@ Setup xrootd"
cp -r ${filesFromSender}/XRootD $VIRTUAL_ENV/lib/python3.9/site-packages/
cp -r ${filesFromSender}/pyxrootd $VIRTUAL_ENV/lib/python3.9/site-packages/

export IFDH_CP_MAXRETRIES=2

echo "@@ outDir : "${outDir}
echo "@@ ifdh  mkdir_p "${outDir}
ifdh  mkdir_p ${outDir} || true

echo "@@ source ${filesFromSender}/run_"${nProcess}".sh "
ls -ltr

cp ${filesFromSender}/run_${nProcess}.sh ./
source run_${nProcess}.sh  &> log_${nProcess}.log

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
