#!/usr/bin/env python3
"""
GUMP Processing Pipeline: Dataframe Processing -> Uproot Staging -> MakesBruce Macro
"""
import argparse
import os
import subprocess
import sys

import awkward as ak
import numpy as np
import pandas as pd
import uproot

workspace_root = os.getcwd()
sys.path.insert(0, workspace_root + "/../gump/")

from sbruce import *

sys.path.insert(0, workspace_root + "/../maple/")
# Import local physics framework utilities
import loaddf

from rwt_map import resolve_function

def get_cols_to_drop():
    """Returns the explicit lists of columns to drop to save memory."""

    return ['is_clear_cosmic', 'crlongtrkdiry', 'p_len', 'mu_E', 'mu_T', 
            'p_E', 'p_T', 'del_Tp', 'del_phi', 'has_stub', 
            'true_pcand_pdg', 'true_p_dir_x', 'true_p_dir_y', 
            'true_p_dir_z', 'true_pcand_dir_x', 'true_pcand_dir_y', 
            'true_pcand_dir_z', 'true_pcand_end_x', 'true_pcand_end_y', 
            'true_pcand_end_z', 'true_mucand_pdg', 'true_mu_dir_x', 
            'true_mu_dir_y', 'true_mu_dir_z', 'true_mucand_dir_x', 
            'true_mucand_dir_y', 'true_mucand_dir_z', 'true_mucand_end_x', 
            'true_mucand_end_y', 'true_mucand_end_z', 'stub_l0_5cm_dedx',
            'stub_l0_5cm_charge','stub_l1cm_dedx','stub_l1cm_charge',
            'stub_l2cm_dedx','stub_l2cm_charge','stub_l3cm_dedx',
            'stub_l3cm_charge','stub_l4cm_dedx','stub_l4cm_charge',
            'prot_chi2smear5_of_prot_cand', 'prot_chi2smear5_of_mu_cand', 
            'mu_chi2smear5_of_mu_cand', 'mu_chi2smear5_of_prot_cand', 
            'tmatch_pur', 'tmatch_eff', 'true_baseline', 
            'true_nu_pdg_x', 'true_nu_pdg_y','true_nmu_27MeV', 
            'true_np_20MeV', 'true_np_50MeV', 'true_npi_30MeV', 
            'is_cosmic', 'flash_sumpe', 'true_mucand_p', 'true_pcand_p', 
            'mu_true_p', 'p_true_p', 'true_mu_end_x', 
            'true_p_end_x', 'true_mu_end_y', 'true_p_end_y', 
            'true_mu_end_z', 'true_p_end_z','crthit', 'true_nu_E', 
            'p_true_pdg', 'mu_true_pdg', 'mu_chi22lo_of_mu_cand', 
            'mu_chi22hi_of_mu_cand', 'prot_chi22lo_of_mu_cand', 
            'prot_chi22hi_of_mu_cand', 'mu_chi22lo_of_prot_cand', 
            'mu_chi22hi_of_prot_cand', 'prot_chi22lo_of_prot_cand', 
            'prot_chi22hi_of_prot_cand', 'true_mu_p', 'true_p_p', 
            'pot_univ', "other_shw_length","other_trk_length",
            "slc_vtx_x", "slc_vtx_y", "slc_vtx_z", "nu_score",
            "mu_chi2_of_mu_cand", "mu_chi2_of_prot_cand",
            "prot_chi2_of_mu_cand", "prot_chi2_of_prot_cand",
            "mu_chi2smear13_of_mu_cand", "mu_chi2smear13_of_prot_cand",
            "prot_chi2smear13_of_mu_cand", "prot_chi2smear13_of_prot_cand",
            "mu_chi2_lo_of_mu_cand", "mu_chi2_lo_of_prot_cand",
            "prot_chi2_lo_of_mu_cand", "prot_chi2_lo_of_prot_cand",
            "mu_chi2_hi_of_mu_cand", "mu_chi2_hi_of_prot_cand",
            "prot_chi2_hi_of_mu_cand", "prot_chi2_hi_of_prot_cand",
            "mu_len", "mu_end_x", "mu_end_y", "mu_end_z",
            "p_end_x", "p_end_y", "p_end_z",
            "true_nu_vtx_x", "true_nu_vtx_y", "true_nu_vtx_z",
            "flash_maxpe"
        ]

def main():
    parser = argparse.ArgumentParser(
        description="GUMP Command Line TTree Processing Pipeline"
    )

    VALID_CONFIGS = ["mc", "data"]

    parser.add_argument(
        "-c",
        "--config",
        required=True,
        choices=VALID_CONFIGS,
        help="specify config (e.g., mc or data)",
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Input .df dataframe file path"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Final output .root file path"
    )
    parser.add_argument(
        "-w", "--weights", action="store_true", help="include weights"
    )
    parser.add_argument(
        "-j", "--cores", type=int, default=1, help="Number of cores for loadl"
    )
    parser.add_argument(
        "-f", "--dvsplines", type=str, default="rwt_outputs", help="directory of detvar splines"
    )

    parser.add_argument(
        "-s", "--selection",
        type=resolve_function,
        default=None,
        help="Selection function to run (e.g., 'gc.all_cuts' or 'gc.coworker_cuts')"
    )

    args = parser.parse_args()

    if args.config == "mc":
        is_mc = True
    else:
        is_mc = False

    # Determine detector configuration context from string names
    if "ICARUS" in args.input:
        if "Run4" in args.input:
            detector_context = "ICARUS Run4"
        else:
            detector_context = "ICARUS Run2"
    else:
        detector_context = "SBND"

    print(f"=== Initializing Pipeline ===")
    print(f"Input File:    {args.input}")
    print(f"Output Target: {args.output}")
    print(f"Detector:      {detector_context} | Mode: {'MC' if is_mc else 'Data'}")

    # Generate transient intermediate root name
    temp_stage1_root = args.output.replace(".root", "_stage1_tmp.root")

    # Step 1: Execute loaddf framework loader
    print("\n--- Step 1: Loading Dataframe via loaddf ---")
    # Stub placeholder cuts if not importing global cut library modules directly
    # Assumes 'gc.all_cuts' equivalent logic is packaged or ignored for flat passes

    df, _, _ = loaddf.loadl(
        [args.input],
        njob=args.cores,
        xsec_univ=False,
        flux_univ=False,
        sep_flux_univ=args.weights,  # Enable multisim loops only for MC configurations
        sep_g4_univ=args.weights,  # Enable multisim loops only for MC configurations
        xsec_spline=args.weights,
        pot_spline=args.weights,
        match_Enu=is_mc,
        load_truth=is_mc,
        load_evtrec=is_mc,
        detvar_spline=args.weights,
        spline_dir=args.dvsplines,
        include_syst=args.weights,
        reweight_aFF=args.weights,
        preselection=args.selection,
        detector=detector_context,
        lightmem=True,
    )

    print(f"args.input: {args.input}")
    # Step 2: Export to intermediate uproot template structure
    print("\n--- Step 2: Exporting to Staging Uproot Structure ---")
    export_dataframe_to_uproot(df, temp_stage1_root)

    # Step 3: Run C++ ROOT Macros to balance structures
    print("\n--- Step 3: Processing ROOT Branch Rebalancing Macro ---")
    run_makesbruce_macro(temp_stage1_root, args.output)

    os.remove(temp_stage1_root)
    print(f"\n[+] Pipeline execution completed successfully for: {args.output}\n")

if __name__ == "__main__":
    main()
