import re, array, os
import argparse
import ROOT
from ROOT import TFile, TTree, gROOT, addressof, double
from array import array
import pandas as pd
import numpy as np
import awkward as ak


## Function to read in data from dataframe and fill the ROOT tree
def main(input_name, output_name):
    f = TFile(output_name, 'RECREATE' )

    df_var = pd.read_hdf(input_name, key='var')
    df_systs = pd.read_hdf(input_name, key='wgt')

    NEvents = len(df_var)

    selectedNutree = TTree('selectedNu', 'selectedNu')
    multisigmatree = TTree('multisigmaTree', 'multisigmaTree')

    systs = ['CoulombCCQE', 'FrAbs_N', 'FrCEx_N', 'MaCCRES', 'MaNCRES', 'MFP_N', 'MvCCRES',
             'NormCCMEC', 'RPA_CCQE', 'ZExpA1CCQE', 'ZExpA2CCQE', 'ZExpA3CCQE', 'ZExpA4CCQE']

    val_vecs = []
    sigma_vecs = []
    val_branches = []
    sigma_branches = []

    # selectedNu variables
    trueE = array('d', [0])
    recoE = array('d', [0])
    dPt = array('d', [0])
    truepdg = array('d', [0])
    iscc = array('d', [0])
    baseline = array('d', [0])

    selectedNutree.Branch("trueE", trueE, 'trueE/D')
    selectedNutree.Branch("recoE", recoE, 'recoE/D')
    selectedNutree.Branch("dPt", dPt, 'dPt/D')
    selectedNutree.Branch("truepdg", truepdg, 'truepdg/D')
    selectedNutree.Branch("iscc", iscc, 'iscc/D')
    selectedNutree.Branch("baseline", baseline, 'baseline/D')
    prefix = "GENIEReWeight_SBN_v1_multisigma_"

    str_to_sigma = {'ms3': -3,
                    'ms2': -2,
                    'ms1': -1,
                    'ps1':  1,
                    'ps2':  2,
                    'ps3':  3}

    for syst in systs:
        val_vecs.append(ROOT.std.vector('double')())
        sigma_vecs.append(ROOT.std.vector('double')())

        val_branches.append(multisigmatree.Branch(prefix+syst, val_vecs[-1]))
        sigma_branches.append(multisigmatree.Branch(prefix+syst+'_sigma', sigma_vecs[-1]))

    for eventi in range(NEvents):
        # Loop over systs
        for systi in range(len(systs)):
            val_vecs[systi].clear()
            sigma_vecs[systi].clear()
    
            for val in df_systs[prefix+systs[systi]].values[eventi]:
                val_vecs[systi].push_back(val)
            val_vecs[systi].push_back(1)

            for col in df_systs[prefix+systs[systi]].columns:
                sigma_vecs[systi].push_back(str_to_sigma[col])
            sigma_vecs[systi].push_back(0)

        multisigmatree.Fill()

        # Now do non-syst values
        trueE[0] = df_var.trueE.values[eventi]  
        recoE[0] = df_var.caloE.values[eventi]  
        dPt[0] = df_var.dPt.values[eventi]  
        truepdg[0] = df_var.truepdg.values[eventi]  
        iscc[0] = df_var.iscc.values[eventi]  
        baseline[0] = df_var.baseline.values[eventi]  
        selectedNutree.Fill()

    f.Write()
    f.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert dataframe object to sBruce tree for use in ProFit")
    parser.add_argument('-i', '--input_file', type=str, help="dataframe input.")
    parser.add_argument('-o', '--output_file', type=str, help="Output sBruce file.")
    args = parser.parse_args()

    assert args.input_file and args.output_file, "Must include dataframe input and sBruce output files."
    main(args.input_file, args.output_file)
