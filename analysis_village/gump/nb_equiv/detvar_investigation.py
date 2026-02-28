import pandas as pd
import os
import sys
from cycler import cycler
from plot_tools import *
workspace_root = os.getcwd()
sys.path.insert(0, workspace_root + "/../../")
import pyanalib.pandas_helpers as ph
from pyanalib.split_df_helpers import *
from makedf.util import *
from rwt_map import *

def main():
    prefix = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/"
    b_E = np.array([0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.25, 1.5])#, 3.0])

    dataframes, pot = match_across_detvar([f"{prefix}SBND_SpringMC_Nom.df", f"{prefix}SBND_SpringMC_WMXThetaXW.df"])

    print(pot)
 
    recodf = dataframes[0]
    wm_df = dataframes[1]

    recodf = all_cuts(recodf, recodf.detector.iloc[0])
    wm_df = all_cuts(wm_df, wm_df.detector.iloc[0])

    nom_n, _ = np.histogram(recodf.nu_E_calo, bins=b_E, weights=np.ones(len(recodf.nu_E_calo))/pot[0])
    wm_n, _  = np.histogram(wm_df.nu_E_calo, bins=b_E, weights=np.ones(len(wm_df.nu_E_calo))/pot[1])

    recodf_rwgt = apply_map(recodf, 'XThetaXW.txt', 'WireMod_SBN_v1_multisigma_WMXThetaXW')
    rwgt_n, _ = np.histogram(recodf.nu_E_calo, weights=recodf_rwgt['WireMod_SBN_v1_multisigma_WMXThetaXW']['ps1'].values, bins=b_E)

    bin_centers = (b_E[:-1] + b_E[1:]) / 2


    diag = [2.45024255e+03, 1.74269777e-02, 3.49286335e+02, 1.73215343e+01,
            5.33963286e+01, 3.68239971e+03, 2.88429364e+02, 1.85497344e+00,
            2.08044201e+03, 2.99063334e+02, 3.96952341e+02, 1.11363614e+03,
            1.07783442e+03, 2.19530232e+03, 2.16877292e+02,]

    SCV = [1678.61650228, 2584.47310732, 3446.55703653, 3848.30506537, 3988.38820891,
           3844.80516827, 3480.171428,   3042.64577638, 2598.51184074, 2150.30845446,
           1738.17203275, 1455.08260177, 1150.79105022, 2796.5478703,   489.42385818]

    diff = [1629.11657755, 2584.34109604, 3427.86783284, 3844.14315064, 3981.08092976,
            3784.12238895, 3463.1882197,  3041.28380228, 2552.89997815, 2133.01500676,
            1718.24836991, 1421.71141445, 1117.96066157, 2749.69381694,  474.6971039 ]

    var = [1698.27820814, 2534.62328929, 3522.40610367, 3919.16242926, 4082.87470857,
           3766.04230412, 3538.41568019, 3155.5911302 , 2593.47607966, 2194.80371138,
           1842.68083038, 1427.29248814, 1163.40298779, 2799.79073279,  470.20055886]

    for arr, l in zip([wm_n, rwgt_n], ['Sample', 'Reweight']):

        rel_diff = np.divide(np.abs(arr-nom_n), nom_n, out=np.zeros_like(nom_n, dtype=float), where=nom_n!=0)

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True,
                                   gridspec_kw={'height_ratios': [3, 1]})
        fig.subplots_adjust(hspace=0.05)

        ax1.step(bin_centers, nom_n, where='mid', label='Nominal', lw=2)
        ax1.step(bin_centers, arr, where='mid', label='WireMod', lw=2)
        #ax1.step(bin_centers, SCV, where='mid', label='GrayNominal', lw=2)
        #ax1.step(bin_centers, var, where='mid', label='GrayWireMod', lw=2)
        ax1.set_ylabel('Counts')
        ax1.legend()
        ax1.set_title('Detector Variation Investigation')

        ax2.axhline(0, color='black', linestyle='--', alpha=0.5) # Reference line at 0
        ax2.step(bin_centers, np.abs(rel_diff), where='mid', color='red', lw=2)
        ax2.set_ylabel('Rel. Diff')
        ax2.set_xlabel('Value')
        ax2.set_ylim(0.0, 0.125) # Adjust based on your data spread
        plt.savefig(f'debug/{l}invest.png')

if __name__ == "__main__":
    main()
