import pandas as pd
import os
import sys
from cycler import cycler

# Add the head direcoty to sys.path
workspace_root = os.getcwd()
sys.path.insert(0, workspace_root + "/../../")
from analysis_village.gump.plot_tools import *
import pyanalib.pandas_helpers as ph
from pyanalib.split_df_helpers import *
from makedf.util import *

from itertools import product 
from more_itertools import distinct_permutations

def get_gain_vars(recodf, plot=False):
    """Main analysis pipeline."""
    b_E = np.array([0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.25, 1.5, 3.0])
    b_p = np.array([0.0, 0.2, 0.4, 0.6])
    ret_dfs = []

    DETECTOR = recodf.detector.iloc[0]

    # make sure '' is first otherwise we lose it :(
    var_dict = {'':'CV', 'lo':'Low', 'hi':'High'}
    for s in var_dict.keys():
        if s != '': # 4 + 6 + 4 + 1 = 15 permutations for non-cv
            # generate permutatons
            for i in range(4):
                curr_iter = ['']*i + [s]*(4-i)
                for j in distinct_permutations(curr_iter, 4):
                    for (p1, p2), k in zip(list(product(['mu', 'prot'], repeat=2)), range(4)):
                        tempdf = recodf.copy()
                        tempdf[f'{p1}_chi2_of_{p2}_cand'] = tempdf[f'{p1}_chi2{j[k]}_of_{p2}_cand']
                    tempdf = all_cuts(tempdf, DETECTOR)
                    ret_dfs.append(tempdf)
                    if plot:
                        n, _, _ = plt.hist(tempdf.nu_E_calo, histtype='step', bins=b_E, alpha=0.5, color='orange')
                    #plt.hist(recodf_cut.nu_E_calo, label=' '.join([var_dict[chi] for chi in j]), histtype='step', bins=b_E)
    if plot:
        tempdf = all_cuts(recodf.copy(), DETECTOR)
        plt.hist(tempdf.nu_E_calo, label='CV', histtype='step', bins=b_E, color='black')
        plt.xlabel(r'Reconstructed Energy $E_{calo}$ [GeV]')
        plt.ylabel('Events')
        plt.title('Recombination Variations Demo', fontsize=20)
        plt.legend()
        plt.savefig('recomb_test.png')

    return ret_dfs

def get_smear_vars(recodf, plot=False):
    """Main analysis pipeline."""
    b_E = np.array([0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.25, 1.5, 3.0])
    b_p = np.array([0.0, 0.2, 0.4, 0.6])
    ret_dfs = []

    # make sure '' is first otherwise we lose it :(
    var_dict = {'':'CV', 'smear5':'smear5', 'smear13':'smear13'}
    for s in var_dict.keys():
        if s != '': # 4 + 6 + 4 + 1 = 15 permutations for non-cv
            # generate permutatons
            for i in range(4):
                curr_iter = ['']*i + [s]*(4-i)
                for j in distinct_permutations(curr_iter, 4):
                    for (p1, p2), k in zip(list(product(['mu', 'prot'], repeat=2)), range(4)):
                        tempdf = recodf.copy()
                        tempdf[f'{p1}_chi2_of_{p2}_cand'] = tempdf[f'{p1}_chi2{j[k]}_of_{p2}_cand']
                    tempdf = all_cuts(tempdf, tempdf.detector.iloc[0])
                    ret_dfs.append(tempdf)
                    if plot:
                        n, _, _ = plt.hist(tempdf.nu_E_calo, histtype='step', bins=b_E, alpha=0.5, color='orange')
                    #plt.hist(recodf_cut.nu_E_calo, label=' '.join([var_dict[chi] for chi in j]), histtype='step', bins=b_E)
    if plot:
        tempdf = all_cuts(recodf.copy(), tempdf.detector.iloc[0])
        plt.hist(tempdf.nu_E_calo, label='CV', histtype='step', bins=b_E, color='black')
        plt.xlabel(r'Reconstructed Energy $E_{calo}$ [GeV]')
        plt.ylabel('Events')
        plt.title('Recombination Variations Demo', fontsize=20)
        plt.legend()
        plt.savefig('recomb_test.png')

    return ret_dfs

def test():
    b_E = np.array([0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.25, 1.5, 3.0])
    b_p = np.array([0.0, 0.2, 0.4, 0.6])

    # dfs = []
    # for i in range(6):
    #     dfs.append(load_dfs(f"/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/SBND_SpringMC_rewgt_{i}.df", ['evt'])['evt'])
    # recodf = pd.concat(dfs)

    recodf = load_dfs(f"/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/ICARUS_SpringMCOverlay_rewgt.df", ['evt'])['evt']

    get_gain_vars(recodf, plot=True)

if __name__ == "__main__":
    test()
    #get_smear_vars()
    #get_gain_vars()
