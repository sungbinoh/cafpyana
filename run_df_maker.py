#!/usr/bin/env python3 
import os,sys,time
import datetime
import pathlib
#from TimeTools import *
import argparse
import tables
from pyanalib.ntuple_glob import NTupleGlob
import pandas as pd
import warnings
import re

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=tables.exceptions.NaturalNameWarning)
pd.set_option('future.no_silent_downcasting', True)

## Arguments
parser = argparse.ArgumentParser(
    description="Data frame maker command: process input flatcaf files and generate output dataframes.",
    epilog="""\
Examples:

  -- Use Pool
  $ python run_df_maker.py -c ./configs/cohpi_slcdf.py -o test_cohpi_slcdf -i input_0.root,input_1.root,...

  -- Use Grid (adding -ngrid to an integer > 0 will automatically submit grid jobs)
  $ python run_df_maker.py -ngrid 2 -c ./configs/cohpi_slcdf.py -o test_cohpi_slcdf -i input_0.root,input_1.root,...

  -- Note!!
  Output df files are sent to /pnfs/<exp>/scratch/users/<User>/cafpyana_out in Grid mode
""",
    formatter_class=argparse.RawTextHelpFormatter  # Ensures line breaks are preserved
)
parser.add_argument('-c', dest='config', default="", help="Path to the data frame maker configuration file in ./configs, i.e.) -c ./configs/mcnu.py.")
parser.add_argument('-o', dest='output', default="", help="output data frame name prefix")
parser.add_argument('-i', dest='inputfiles', default="", help="input root file path, you can submit multiple files using comma, i.e.) -i input_0.root,input_1.root")
parser.add_argument('-l', dest='inputfilelist', default="", help="a file of list for input root files")
parser.add_argument('-ncpu', dest='NCPU', default=-1, type=int, help="Number of CPUs to run on. Default is set to number on server.")
parser.add_argument('-ngrid', dest='NGridJobs', default=0, type=int, help="Number of grid jobs. Default = 0, no grid submission.")
parser.add_argument('-nfile', dest='NFiles', default=0, type=int, help="Number of files to run. Default = 0, run all input files.")
parser.add_argument('-split', dest='SplitSize', default=1.0, type=float, help="Split size in GB before writing to HDF5. Default = 1.0 GB.")

args = parser.parse_args()

def run_pool(output, inputs, nproc):
    os.nice(10)
    ntuples = NTupleGlob(inputs, None)

    # if PREPROCESS doesn't exist, set it to None
    global PREPROCESS
    try:
        PREPROCESS
    except:
        PREPROCESS = []

    dfss = ntuples.dataframes(nproc=nproc, fs=DFS, preprocess=PREPROCESS)
    output = pathlib.Path(output).with_suffix('.df')
    split_margin = args.SplitSize
    with pd.HDFStore(output) as hdf_pd:
        NAMES.append("histpotdf")
        NAMES.append("histgenevtdf")
        split_idx = 0
        split_bytes = 0.0
        split_buffers = {k: [] for k in NAMES}

        def _flush_split():
            nonlocal split_idx, split_bytes, split_buffers
            wrote_any = False
            for k, buffer in split_buffers.items():
                if not buffer:
                    continue
                concat_df = pd.concat(buffer, ignore_index=False)
                this_key = f"{k}_{split_idx}"
                hdf_pd.put(key=this_key, value=concat_df, format="fixed")
                print(f"Saved {this_key}: {concat_df.memory_usage(deep=True).sum() / (1024**3):.4f} GB")
                del concat_df
                wrote_any = True

            if wrote_any:
                split_idx += 1
            split_bytes = 0.0
            split_buffers = {k: [] for k in NAMES}

        for dfs in dfss:
            this_NAMES = NAMES
            if len(dfs) == 2: ## no or empty recTree
                this_NAMES = ["histpotdf", "histgenevtdf"]

            for k, df in zip(reversed(this_NAMES), reversed(dfs)):
                if df is None:
                    continue

                size_bytes = df.memory_usage(deep=True).sum()
                size_gb = size_bytes / (1024**3)
                split_buffers[k].append(df)
                split_bytes += size_gb
                del df

            if split_bytes >= split_margin:
                _flush_split()

        _flush_split()

        # Save the split count metadata
        split_df = pd.DataFrame({"n_split": [split_idx]})
        hdf_pd.put(key="split", value=split_df, format="fixed")
        print(f"Saved split info: {split_df.iloc[0]['n_split']} total splits")

def run_grid(inputfiles):
    # 1) dir/file name style
    JobStartTime = datetime.datetime.now()
    timestamp =  JobStartTime.strftime('%Y_%m_%d_%H%M%S')

    # 2) Define MasterJobDir -- produce grid job submission scripts in $CAFPYANA_GRID_OUT_DIR
    CAFPYANA_GRID_OUT_DIR = os.environ['CAFPYANA_GRID_OUT_DIR']
    MasterJobDir = CAFPYANA_GRID_OUT_DIR + "/logs/" + timestamp + '__' + args.output + "_log"
    OutputDir = CAFPYANA_GRID_OUT_DIR + "/dfs/" + timestamp + '__' + args.output
    os.system('mkdir -p ' + MasterJobDir)

    # 3) grid job is based on number of files
    ngrid = args.NGridJobs
    if(len(inputfiles) <= ngrid):
        ngrid = len(inputfiles)

    NInputfiles = len(inputfiles)
    print("Number of Grid Jobs: %d, number of input caf files: %d" % (ngrid, NInputfiles))

    # 4) prepare bash scripts for each job and make tarball
    flistForEachJob = []
    for i in range(0,ngrid):
        flistForEachJob.append( [] )

    for i_line in range(0,len(inputfiles)):
        flistForEachJob[i_line%ngrid].append(inputfiles[i])

    for i_flist in range(0,len(flistForEachJob)):
        flist = flistForEachJob[i_flist]
        input_list_name = 'inputs_%s.list'%(i_flist)
        input_list_path = MasterJobDir + '/' + input_list_name

        with open(input_list_path, 'w') as list_out:
            for f in flist:
                local_name = re.sub(r'^(?:root://[^/]+)?(?:/pnfs/fnal.gov/usr/)?', '', f)
                local_name = local_name.lstrip('/')
                local_name = local_name.replace('/', '__')
                list_out.write(local_name + '\n')

        out = open(MasterJobDir + '/run_%s.sh'%(i_flist),'w')
        out.write('#!/bin/bash\n')
        cmd = 'python run_df_maker.py -c ' + args.config + ' -o ' + args.output + '_%d'%i_flist + '.df -ncpu 7 -l ' + input_list_name
        for i_f in range(0,len(flist)):
            out.write('echo "[run_%s.sh] input %d : %s"\n'%(i_flist, i_f, flist[i_f]))
            local_name = re.sub(r'^(?:root://[^/]+)?(?:/pnfs/fnal.gov/usr/)?', '', flist[i_f])
            local_name = local_name.lstrip('/')
            local_name = local_name.replace('/', '__')
            out.write('xrdcp ' + flist[i_f] + ' ' + local_name + '\n') ## -- for checking auth
        out.write('echo "[run_%s.sh] input list file: %s"\n'%(i_flist, input_list_name))
        out.write('echo "[run_%s.sh] number of inputs:"\n'%(i_flist))
        out.write('wc -l ' + input_list_name + '\n')
        out.write('ls -alh\n')
        out.write(cmd)
        out.close()

    os.system('cp ./bin/grid_executable.sh %s' %MasterJobDir)

    # 5) prepare a package for xrootd
    CAFPYANA_WD = os.environ['CAFPYANA_WD']
    cp_XRootD = "cp -r " + CAFPYANA_WD + "/envs/xrootd-5.6.9/build/lib.linux-x86_64-cpython-310/XRootD " + MasterJobDir
    cp_pyxrootd = "cp -r " + CAFPYANA_WD + "/envs/xrootd-5.6.9/build/lib.linux-x86_64-cpython-310/pyxrootd " + MasterJobDir
    os.system(cp_XRootD)
    os.system(cp_pyxrootd)

    # 6) git archive of the current branch's last commit
    archive_repo = "git archive -o " + MasterJobDir + "/cafpyana.tar.gz HEAD"
    os.system(archive_repo)

    # 7) move to MasterJobDir to submit jobs
    os.chdir(MasterJobDir)
    tar_cmd = 'tar cf bin_dir.tar ./'
    os.system(tar_cmd)

    submitCMD = '''jobsub_submit \\
-G sbnd \\
--auth-methods="token" \\
-e LC_ALL=C \\
--role=Analysis \\
--resource-provides="usage_model=DEDICATED,OPPORTUNISTIC" \\
-l '+SingularityImage=\\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9\:9.7\\"' \\
--lines '+FERMIHTC_AutoRelease=True' --lines '+FERMIHTC_GraceMemory=2000' --lines '+FERMIHTC_GraceLifetime=3600' \\
--append_condor_requirements='(TARGET.HAS_SINGULARITY=?=true)' \\
--tar_file_name "dropbox://$(pwd)/bin_dir.tar" \\
-N %d \\
--disk 10GB \\
--cpu 7 \\
--memory 5GB \\
--expected-lifetime 1h \\
"file://$(pwd)/grid_executable.sh" \\
"%s" \\
"%s"'''%(ngrid,OutputDir,args.output)

    print(submitCMD)
    os.system(submitCMD)
    
    # go back to working dir
    os.chdir(CAFPYANA_WD)
    
if __name__ == "__main__":
    printhelp = ((args.inputfiles == "" and args.inputfilelist == "") or args.config == "" or args.output == "")
    if printhelp:
        parser.print_help()
        print(parser.epilog)
        sys.exit(1)
        
    else:
        ### Organize input list 
        InputSamples = []
        StringForHash = ""
        if args.inputfilelist != "":
            lines = open(args.inputfilelist)
            for line in lines:
                if "#" in line:
                    continue
                line = line.strip('\n')
                InputSamples.append(line)
                StringForHash += line
        else:
            split_inputfiles = args.inputfiles.split(",")
            for split_inputfile in split_inputfiles:
                InputSamples.append(split_inputfile)
                StringForHash += args.inputfiles

        if(args.NFiles > 0 and len(InputSamples) > args.NFiles):
            InputSamples = InputSamples[:args.NFiles]
                
        ### check if it is grid mode for pool mode
        if args.NGridJobs == 0:
            print("Running Pool mode");
            exec(open(args.config).read())
            run_pool(args.output, InputSamples, "auto" if args.NCPU < 0 else args.NCPU)

        elif args.NGridJobs > 0:
            print("Running Grid mode");
            run_grid(InputSamples)
            
        else:
            print("-ngrid must be greater than 0.");
