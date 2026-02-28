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

def main():
    prefix = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/"

    cv_files = []
    for i in range(19):
        cv_files.append(prefix+f"SBND_SpringMC_rewgt_{i}.df")
    cv_files.append(prefix+"ICARUS_SpringMCOverlay_rewgt.df")

    tot_tot_pot = 0
    for input in cv_files:
        nsplits = get_n_split(input)
        tot_pot = 0
        for n in range(nsplits):
            tot_pot += scale_pot(pd.read_hdf(input,"hdr_"+str(n)))
        if "SBND" in input:
            tot_tot_pot += tot_pot
        print(f"{input} sample pot: {tot_pot}")
        print(f"Combined SBND pot: {tot_tot_pot}")

    prefix = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/"
    dirt_nd = prefix+"SBND_SpringLowEMC.df"
    dirt_fd = prefix+"ICARUS_SpringMCDirt_slimwgt.df"

    for input in [dirt_nd, dirt_fd]:
        nsplits = get_n_split(input)
        tot_pot = 0
        for n in range(nsplits):
            tot_pot += scale_pot(pd.read_hdf(input,"hdr_"+str(n)))
        print(f"{input} sample pot: {tot_pot}")

    offbeam_nd = prefix+"SBND_SpringBNBOffData_5000.df"
    offbeam_fd = prefix+"ICARUS_SpringRun2BNBOff_unblind_prescaled.df"
    N_GATES_ON_PER_5e12POT = 1.3886218026202426
    nsplits = get_n_split(offbeam_fd)

    ngates_OFF = 0
    for i in range(nsplits):
        single_ngates_OFF = pd.read_hdf(offbeam_fd, "trig_"+str(i)).gate_delta.sum()*(1-1/20.)
        #print(f"{single_ngates_OFF} for {i} trig \"pot\"")
        ngates_OFF += single_ngates_OFF
    print(f"{offbeam_fd} sample \"pot\": {5e12*ngates_OFF/N_GATES_ON_PER_5e12POT}")

    N_GATES_ON_PER_5e12POT = 1.05104
    nsplits = get_n_split(offbeam_nd)

    ngates_OFF = 0
    for i in range(nsplits):
        single_ngates_OFF = pd.read_hdf(offbeam_nd, "hdr_"+str(i)).noffbeambnb.sum()
        # print(f"{single_ngates_OFF} for {i} trig \"pot\"")
        ngates_OFF += single_ngates_OFF
    print(f"{offbeam_nd} sample pot: {5e12*ngates_OFF/N_GATES_ON_PER_5e12POT}")

if __name__ == "__main__":
    main()
