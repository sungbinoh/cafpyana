#!/usr/bin/env python3 
import os,sys,time
import datetime
#from TimeTools import *
import argparse
from pyanalib.ntuple_glob import NTupleGlob
from pyanalib.split_df_helpers import *
import pandas as pd
from tqdm.auto import tqdm
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

def run(output, inputs):
    first_fill = True

    # count total splits beforehand for nice tqdm bar
    total_splits = sum(get_n_split(inp) for inp in inputs)

    with uproot.recreate(output) as f, tqdm(total=total_splits, desc="Processing splits") as pbar:
        for input in inputs:
            this_n_split = get_n_split(input)
            for split in range(this_n_split):
                result = TTREEMKR(input, split)

                if isinstance(result, tuple) and len(result) == 2 \
                   and isinstance(result[0], pd.DataFrame) and isinstance(result[1], pd.DataFrame):
                    recodf, truedf = result
                elif isinstance(result, pd.DataFrame):
                    recodf, truedf = result, pd.DataFrame()
                else:
                    raise TypeError(
                        "TTREEMKR must return either a (recodf, truedf) tuple of DataFrames "
                        "or a single reco DataFrame."
                    )

                write_true = truedf is not None and not truedf.empty
                if first_fill:
                    f["SelectedEvents"] = recodf.to_dict(orient="list")
                    if write_true:
                        f["TrueEvents"] = truedf.to_dict(orient="list")
                    first_fill = False
                else:
                    f["SelectedEvents"].extend(recodf.to_dict(orient="list"))
                    if write_true:
                        f["TrueEvents"].extend(truedf.to_dict(orient="list"))

                # update progress bar
                pbar.update(1)

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
