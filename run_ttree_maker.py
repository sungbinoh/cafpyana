#!/usr/bin/env python3 
import os,sys,time
import datetime
#from TimeTools import *
import argparse
from pyanalib.ntuple_glob import NTupleGlob
import pandas as pd
import warnings

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

## Arguments
parser = argparse.ArgumentParser(
    description="Data frame maker command: process input flatcaf files and generate output dataframes.",
    epilog="""\
Examples:

  -- Use Pool
  $ python run_ttree_maker.py -c ./configs/cohpi_ttree.py -o output.root -i input.df
""",
    formatter_class=argparse.RawTextHelpFormatter  # Ensures line breaks are preserved
)
parser.add_argument('-c', dest='config', default="", help="Path to the data frame maker configuration file in ./configs, i.e.) -c ./configs/mcnu.py.")
parser.add_argument('-o', dest='output', default="", help="output data frame name prefix")
parser.add_argument('-i', dest='inputfiles', default="", help="input root file path, you can submit multiple files using comma, i.e.) -i input_0.root,input_1.root")
parser.add_argument('-l', dest='inputfilelist', default="", help="a file of list for input root files")
parser.add_argument('-nfile', dest='NFiles', default=0, type=int, help="Number of files to run. Default = 0, run all input files.")
args = parser.parse_args()

def save(outf, recodf, truedf):
    with uproot.recreate(outf) as f:
        f["SelectedEvents"] = recodf

        if truedf is not None and not truedf.empty:
            f["TrueEvents"] = truedf.to_dict(orient="list")
        else:
            print("Warning: truedf is empty or None — skipping writing 'TrueEvents'")
        #f["TrueEvents"] = truedf

def run(output, inputs):

    ## for MC, expect a tuple of two data frames
    ## for data, expect only one data frame
    result = list(map(TTREEMKR, inputs))

    if all(isinstance(r, (list, tuple)) and len(r) == 2 and 
           isinstance(r[0], pd.DataFrame) and isinstance(r[1], pd.DataFrame) for r in result):
        recodfs, truedfs = zip(*result)
        recodf = pd.concat(recodfs)
        truedf = pd.concat(truedfs)
    elif all(isinstance(r, pd.DataFrame) for r in result):
        # Only recodf returned; truedf will be empty
        recodf = pd.concat(result)
        truedf = pd.DataFrame()
    else:
        raise ValueError("TTREEMKR must return either a DataFrame or a tuple of two DataFrames.")

    save(output, recodf, truedf)

    #recodfs, truedfs = zip(*list(map(TTREEMKR, inputs)))
    #recodf = pd.concat(recodfs)
    #truedf = pd.concat(truedfs)

    #save(output, recodf, truedf)

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
        
        exec(open(args.config).read())
        run(args.output, InputSamples)
