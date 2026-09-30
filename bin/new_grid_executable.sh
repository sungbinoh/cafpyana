# Setup grid submission
outDir=$1
DFPREFIX=$2
nProcess=$PROCESS
git clone https://github.com/gputnam/cafpyana.git
cd cafpyana
git checkout remotes/origin/N8Dev

source setup.sh
mkdir output
thisOutputCreationDir=`pwd`
filesFromSender=${CONDOR_DIR_INPUT}/bin_dir/

export IFDH_CP_MAXRETRIES=2
ifdh  mkdir_p ${outDir}

echo "@@ source ${filesFromSender}/run_"${nProcess}".sh "
ls -alh
pwd

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
