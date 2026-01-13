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
    dffile_nd_1 = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-3-fix-zexp-again-again/SBND_SpringMC_rewgt_1.df"
    dffile_nd_2 = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-3-fix-zexp-again-again/SBND_SpringMC_rewgt_2.df"
    dffile_nd_3 = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-3-fix-zexp-again-again/SBND_SpringMC_rewgt_3.df"
    dffile_nd_4 = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-3-fix-zexp-again-again/SBND_SpringMC_rewgt_4.df"
    dffile_nd_5 = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-3-fix-zexp-again-again/SBND_SpringMC_rewgt_5.df"
    dffile_nd_6 = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-3-fix-zexp-again-again/SBND_SpringMC_rewgt_6.df"
    dffile_nd_7 = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-3-fix-zexp-again-again/SBND_SpringMC_rewgt_7.df"
    dffile_fd = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-3-fix-zexp-again-again/ICARUS_SpringMC_Dev_rewgt.df"

    inputs = [dffile_nd_1, dffile_nd_2, dffile_nd_3, dffile_nd_4, dffile_nd_5, dffile_nd_6, dffile_nd_7, dffile_fd]

    tot_tot_pot = 0
    for input in inputs:
        nsplits = get_n_split(input)
        tot_pot = 0
        for n in range(nsplits):
            tot_pot += scale_pot(pd.read_hdf(input,"hdr_"+str(n)))
        if "SBND" in input:
            tot_tot_pot += tot_pot
        print(f"{input} sample pot: {tot_pot}")
        print(f"Combined SBND pot: {tot_tot_pot}")

    offbeam_nd = "/exp/sbnd/data/users/gputnam/GUMP/sbn-wgted/SBND_SPINE_SpringBNBOffData.df"
    offbeam_fd = "/exp/sbnd/data/users/gputnam/GUMP/sbn-wgted/ICARUS_Run2_BNBoff_uncalo_prescaled.df"

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
