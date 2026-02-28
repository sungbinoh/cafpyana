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
    b_E = np.array([0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.25, 1.5, 3.0])


    recodf = dataframes[0]
    wm_df = dataframes[1]

    DF_DIR = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/"
    SCV_FILES = [DF_DIR + "SBND_SpringMC_rewgt_%i.df" % i for i in range(20)]
    recodf, Smatch, Spot = loaddf.loadl(SCV_FILES, njob=min(len(SCV_FILES), 10), preselection=FVSBND, reweight_aFF=True)

    nom_n, _ = np.histogram(recodf.nu_E_calo, bins=b_E)
    wm_n, _  = np.histogram(recodf.nu_E_calo, bins=b_E)

    bin_centers = (b_E[:-1] + b_E[1:]) / 2

    rel_diff = np.divide(np.sqrt(wm_n), nom_n, out=np.zeros_like(nom_n, dtype=float), where=nom_n!=0)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True,
                               gridspec_kw={'height_ratios': [3, 1]})
    fig.subplots_adjust(hspace=0.05)

    ax1.step(bin_centers, nom_n, where='mid', label='Nominal', lw=2)
    ax1.step(bin_centers, wm_n, where='mid', label='WireMod', lw=2)
    ax1.set_ylabel('Counts')
    ax1.legend()
    ax1.set_title('Detector Variation Investigation')

    ax2.axhline(0, color='black', linestyle='--', alpha=0.5) # Reference line at 0
    ax2.step(bin_centers, np.abs(rel_diff), where='mid', color='red', lw=2)
    ax2.set_ylabel('Rel. Difference')
    ax2.set_xlabel('Value')
    ax2.set_ylim(0.0, 0.125) # Adjust based on your data spread
    plt.savefig('debug/invest.png')

if __name__ == "__main__":
    main()
