# Constants for this analysis
# Today: 2026-02-04

import uproot
import matplotlib.pyplot as plt
from makedf.constants import *
from os import path
import pandas as pd
from pyanalib.split_df_helpers import *
from pyanalib.pandas_helpers import *

plot = False

DETECTOR = "SBND_nohighyz"
EPSILON = 1e-6 # for clipping distributions at bin ranges

def get_xsec_unit():
    # ==== xsec unit calculation ====
    # TODO: z-dependence?
    # flux file, units: /m^2/10^6 POT 
    # 50 MeV bins
    fluxfile = "/exp/sbnd/data/users/munjung/flux/sbnd_original_flux.root"
    flux = uproot.open(fluxfile)
    numu_flux = flux["flux_sbnd_numu"].to_numpy()
    bin_edges = numu_flux[1]
    flux_vals = numu_flux[0]

    if plot:
        fig, ax = plt.subplots()
        plt.hist(bin_edges[:-1], bins=bin_edges, weights=flux_vals, histtype="step", linewidth=2, color="C0")
        plt.xlim(0, 3)
        plt.xlabel("Neutrino Energy [GeV]")
        plt.ylabel("Flux [/m$^{2}$/10$^{6}$ POT]")
        plt.title("SBND $\\nu_\\mu$ Flux")
        plt.savefig("sbnd-flux.pdf", bbox_inches='tight')

    # get integrated flux
    # TODO: refactor this to be always consistent with what's used in the ana notebooks
    _file_dir = "/exp/sbnd/data/users/munjung/xsec/2025Spring_v10_06_00_09"
    _mc_file = path.join(_file_dir, "MC", "BNB_cosmics", "aa.df")
    _mc_split_df = pd.read_hdf(_mc_file, key="split")
    _mc_n_split = get_n_split(_mc_file)
    print("mc_n_split: %d" %(_mc_n_split))
    _n_max_concat = 3
    _mc_hdr_df = load_dfs(_mc_file, ["hdr"], n_max_concat=_n_max_concat)["hdr"]
    _mc_tot_pot = _mc_hdr_df["pot"].sum()
    print("mc_tot_pot: %.3e" %(_mc_tot_pot))

    integrated_flux = _mc_tot_pot * flux_vals.sum() / (1e4  * 1e6) # to cm2 # to POT
    print("Integrated flux: %.3e" % integrated_flux)

    V_SBND = 380 * 380 * 440 # cm3, the active volume of the detector 
    NTARGETS = RHO * V_SBND * N_A / M_AR
    print("# of targets: ", NTARGETS)

    xsec_unit = 1 / (integrated_flux * NTARGETS)
    # TODO: fix scalar overflow error in python v3.10+
    if xsec_unit == 0:
        print("XSEC_UNIT is 0, setting to 1e-38")
        xsec_unit = 1e-38
    print("xsec unit: ", xsec_unit)
    return xsec_unit

XSEC_UNIT = get_xsec_unit()