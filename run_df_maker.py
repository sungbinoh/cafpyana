#!/usr/bin/env python3 
import os,sys,time
import datetime
#from TimeTools import *
import argparse
from pyanalib.ntuple_glob import NTupleGlob
import pandas as pd
import warnings

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
# warnings.simplefilter(action='ignore', category=NaturalNameWarning)

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
parser.add_argument('-ngrid', dest='NGridJobs', default=0, type=int, help="Number of grid jobs. Default = 0, no grid submission.")
parser.add_argument('-nfile', dest='NFiles', default=0, type=int, help="Number of files to run. Default = 0, run all input files.")
args = parser.parse_args()

def run_pool(output, inputs):
    os.nice(10)
    ntuples = NTupleGlob(inputs, None)

    dfs = ntuples.dataframes(nproc="auto", fs=DFS)
    output = output + ".df"
    with pd.HDFStore(output) as hdf:
        for k,df in zip(reversed(NAMES), reversed(dfs)): # go in reverse order so we can delete along the way
            try:
                hdf.put(key=k, value=df, format="fixed")
            except Exception as e:
                print("Table %s failed to save, skipping. Exception: %s" % (k, str(e)))

            del df

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
        out = open(MasterJobDir + '/run_%s.sh'%(i_flist),'w')
        out.write('#!/bin/bash\n')
        cmd = 'python run_df_maker.py -c ' + args.config + ' -o ' + args.output + '_%d'%i_flist + '.df -i'
        for i_f in range(0,len(flist)):
            out.write('echo "[run_%s.sh] input %d : %s"\n'%(i_flist, i_f, flist[i_f]))
            if i_f == 0:
                cmd += ' ' + flist[i_f]
            else: 
                cmd += ',' + flist[i_f]
            #out.write('xrdcp ' + flist[i_f] + ' .\n') ## -- for checking auth
        out.write(cmd)
        out.close()

    os.system('cp ./bin/grid_executable.sh %s' %MasterJobDir)

    # 5) prepare a package for xrootd
    CAFPYANA_WD = os.environ['CAFPYANA_WD']
    cp_XRootD = "cp -r " + CAFPYANA_WD + "/envs/xrootd-5.6.1/build/lib.linux-x86_64-3.9/XRootD " + MasterJobDir
    cp_pyxrootd = "cp -r " + CAFPYANA_WD + "/envs/xrootd-5.6.1/build/lib.linux-x86_64-3.9/pyxrootd " + MasterJobDir
    os.system(cp_XRootD)
    os.system(cp_pyxrootd)

    os.chdir(MasterJobDir)
    tar_cmd = 'tar cf bin_dir.tar ./'
    os.system(tar_cmd)

    submitCMD = '''jobsub_submit \\
-G sbnd \\
--auth-methods="token,proxy" \\
-e LC_ALL=C \\
--role=Analysis \\
--resource-provides="usage_model=DEDICATED,OPPORTUNISTIC" \\
--lines '+FERMIHTC_AutoRelease=True' --lines '+FERMIHTC_GraceMemory=1000' --lines '+FERMIHTC_GraceLifetime=3600' \\
--append_condor_requirements='(TARGET.HAS_SINGULARITY=?=true)' \\
--tar_file_name "dropbox://$(pwd)/bin_dir.tar" \\
--email-to sungbin.oh555@gmail.com \\
-N %d \\
--disk 100GB \\
--expected-lifetime 10h \\
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
            print("Runing Pool mode");
            exec(open(args.config).read())
            run_pool(args.output, InputSamples)

        elif args.NGridJobs > 0:
            print("Runing Grid mode");
            run_grid(InputSamples)
            
        else:
            print("-ngrid must be greater than 0.");
