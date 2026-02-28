import sys
import os
import pandas as pd
workspace_root = os.getcwd()
sys.path.insert(0, workspace_root + "/../../")
from pyanalib.split_df_helpers import *

def scale_pot(df_hdr):
    """Scale DataFrame by desired POT."""
    pot = sum(df_hdr.pot.tolist())
    return pot

def grab_pot(files, mc_bools, sep_bool=True):
    print(f"running: {files}")
    pot = []

    if isinstance(files, (str)):
        files = [files]
    elif not isinstance(files, (list)):
        print("Check input file data type!")
        sys.exit()

    if isinstance(mc_bools, (bool)):
        mc_bools = [mc_bools]*len(files)

    for file, mc_bool in zip(files, mc_bools):
        detector = pd.read_hdf(file, "evt_0").detector.iloc[0]
        if mc_bool:
            tot_pot = 0
            for n in range(get_n_split(file)):
                tot_pot += scale_pot(pd.read_hdf(file,"hdr_"+str(n)))
            pot.append(tot_pot)
            print(f"{file} sample pot: {tot_pot}")
        else:
            if(detector == "ICARUS"):
                N_GATES_ON_PER_5e12POT = 1.3886218026202426
            elif(detector == "SBND"):
                N_GATES_ON_PER_5e12POT = 1.05104

            ngates_OFF = 0
            for i in range(get_n_split(file)):
                if detector == "ICARUS": single_ngates_OFF = pd.read_hdf(file, "trig_"+str(i)).gate_delta.sum()*(1-1/20.)
                elif detector == "SBND": single_ngates_OFF = pd.read_hdf(file, "hdr_"+str(i)).noffbeambnb.sum()
                print(f"{single_ngates_OFF} for {i} num gates")
                ngates_OFF += single_ngates_OFF
            pot.append(5e12*ngates_OFF/N_GATES_ON_PER_5e12POT)
            print(f"{file} sample \"pot\": {pot[-1]}")

    if any(p < 0 for p in pot):
        print("Cannot have negative POT!!!")
        sys.exit()

        if sep_bool: 
            if len(pot) == 1:
                return pot[0]
            else:
                return pot
        else:
            return sum(pot)

def test():
    prefix = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/"

    cv_files = []
    for i in range(19):
        cv_files.append(prefix+f"SBND_SpringMC_rewgt_{i}.df")

    grab_pot(cv_files, True)
    grab_pot(prefix+"ICARUS_SpringMCOverlay_rewgt.df", True)
    grab_pot(prefix+"SBND_SpringLowEMC.df", True)
    grab_pot(prefix+"ICARUS_SpringMCDirt_slimwgt.df", True)
    grab_pot(prefix+"SBND_SpringBNBOffData_5000.df", False)
    grab_pot(prefix+"ICARUS_SpringRun2BNBOff_unblind_prescaled.df", False)

if __name__ == "__main__":
    test()
