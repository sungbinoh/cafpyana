import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from plot_tools import *
# Add the head direcoty to sys.path
workspace_root = os.getcwd()  
sys.path.insert(0, workspace_root + "/../../")

# import this repo's classes
from pyanalib.split_df_helpers import *
import pyanalib.pandas_helpers as ph
from makedf.util import *

import kinematics

def load_data(file):
    """Load event, header, and mcnu data from HDF file."""
    nfiles = get_n_split(file)
    for s in range(nfiles):
        print("df index:"+str(s))
        df_evt = pd.read_hdf(file, "evt_"+str(s))
        df_hdr = pd.read_hdf(file, "hdr_"+str(s))
        df_mcnu = pd.read_hdf(file, "mcnu_"+str(s))
        df_stub = pd.read_hdf(file, "stub_"+str(s))

        matchdf = df_evt.copy()
        matchdf.columns = pd.MultiIndex.from_tuples([(col, '') for col in matchdf.columns])
        df_evt = ph.multicol_merge(matchdf.reset_index(), df_mcnu.reset_index(),
                                    left_on=[("__ntuple", ""), ("entry", ""), ("tmatch_idx", "")],
                                    right_on=[("__ntuple", ""), ("entry", ""), ("rec.mc.nu..index", "")],
                                    how="left") ## -- save all sllices
        cols_to_drop = []
        for c in df_evt.columns:
            if 'GENIE' in c[0] or 'Flux' in c[0]:
                cols_to_drop.append(c)

        df_evt.drop(columns=cols_to_drop, inplace=True)
        del df_mcnu

        if s == 0:
            res_df_evt = df_evt
            res_df_hdr = df_hdr
            res_df_stub = df_stub
        else:
            res_df_evt = pd.concat([res_df_evt, df_evt])
            res_df_hdr = pd.concat([res_df_hdr, df_hdr])
            res_df_stub = pd.concat([res_df_stub, df_stub])
        del df_evt
        del df_hdr
        del df_stub

    return res_df_evt, res_df_hdr, res_df_stub

def scale_pot(df, df_hdr, desired_pot):
    """Scale DataFrame by desired POT."""
    pot = sum(df_hdr.pot.tolist())
    if pot == 0:
        print("Warning, no POT found, using scale of 1.")
        df['glob_scale'] = 1
        return 1, 1
    print(f"POT: {pot}\nScaling to: {desired_pot}")
    scale = desired_pot / pot
    df['glob_scale'] = scale
    return pot, scale

def apply_cuts(df_nd, df_nd_stub, df_fd, df_fd_stub, nd_POT, fd_POT, dffile_nd, dffile_fd):
    comp_sbnd = []
    comp_icarus = []
    cut_labels = ["Initial Sample", "FV Cut", "NuScore Cut", "Two Prong", "PID", "CRT Veto"]
    mode_labels = ['QE', 'MEC', 'RES', 'SIS/DIS', 'COH', "other"]
    det_labels = ("SBND", "ICARUS")
    
    df_nd['og_sig_ct'] = len(df_nd.is_sig[df_nd.is_sig == True])
    df_fd['og_sig_ct'] = len(df_fd.is_sig[df_fd.is_sig == True])

    # Grab some plots w/o cuts
    sbnd_comp, icarus_comp = make_all_plots(df_nd, df_fd, "No Cut", mode_labels, top_labels, det_labels)
    comp_sbnd.append(sbnd_comp)
    comp_icarus.append(icarus_comp)

    # FV Cut
    nd_vtx = pd.DataFrame({'x': df_nd.slc_vtx_x,
                           'y': df_nd.slc_vtx_y,
                           'z': df_nd.slc_vtx_z})
    fd_vtx = pd.DataFrame({'x': df_fd.slc_vtx_x,
                           'y': df_fd.slc_vtx_y,
                           'z': df_fd.slc_vtx_z})

    df_nd = df_nd[fv_cut(nd_vtx, "SBND")]
    df_fd = df_fd[fv_cut(fd_vtx, "ICARUS")]
    sbnd_comp, icarus_comp = make_all_plots(df_nd, df_fd, "FV Cut", mode_labels, top_labels, det_labels)

    comp_sbnd.append(sbnd_comp)
    comp_icarus.append(icarus_comp)

    # NuScore cut
    plot_nuscore_cut([
        df_nd.nu_score[cosmic_cut(df_nd)], df_nd.nu_score[np.invert(cosmic_cut(df_nd))]
    ], [0.5], [
        df_fd.nu_score[cosmic_cut(df_fd)], df_fd.nu_score[np.invert(cosmic_cut(df_fd))]
    ], [0.5], "NuScore", "NuScore.png", [0, 1],
        ['SBND Neutrino', 'SBND Cosmic', 'ICARUS Neutrino', 'ICARUS Cosmic'], 'right', 'select as neutrino')

    df_nd = df_nd[cosmic_cut(df_nd)]
    df_fd = df_fd[cosmic_cut(df_fd)]
    sbnd_comp, icarus_comp = make_all_plots(df_nd, df_fd, "NuScore Cut", mode_labels, top_labels, det_labels)
    comp_sbnd.append(sbnd_comp)
    comp_icarus.append(icarus_comp)

    # Two prong cut
    df_nd = df_nd[twoprong_cut(df_nd)]
    df_fd = df_fd[twoprong_cut(df_fd)]
    sbnd_comp, icarus_comp = make_all_plots(df_nd, df_fd, "Two Prong Cut", mode_labels, top_labels, det_labels)
    comp_sbnd.append(sbnd_comp)
    comp_icarus.append(icarus_comp)

    # PID cuts
    plot_PID_cut([[df_nd.mu_chi2_of_prot_cand, df_nd.mu_chi2_of_mu_cand],
                  [df_fd.mu_chi2_of_prot_cand, df_fd.mu_chi2_of_mu_cand]],
                  [25, 25], "MuScore", "MuScore.png", [0, 65],
                  ['SBND, Protons', 'SBND, Muons', 'ICARUS, Protons', 'ICARUS, Muons'],
                  'left', 'select as muon')

    plot_PID_cut([[df_nd.prot_chi2_of_prot_cand, df_nd.prot_chi2_of_mu_cand],
                  [df_fd.prot_chi2_of_prot_cand, df_fd.prot_chi2_of_mu_cand]],
                  [100, 100], "ProtonScore", "ProtonScore.png", [0, 200],
                  ['SBND, Protons', 'SBND, Muons', 'ICARUS, Protons', 'ICARUS, Muons'],
                  'right', 'select as muon')
    
    df_nd = df_nd[pid_cut(df_nd.mu_chi2_of_mu_cand, df_nd.mu_chi2_of_prot_cand, 
                          df_nd.prot_chi2_of_mu_cand, df_nd.prot_chi2_of_prot_cand, df_nd.mu_len)]
    df_fd = df_fd[pid_cut(df_fd.mu_chi2_of_mu_cand, df_fd.mu_chi2_of_prot_cand, 
                          df_fd.prot_chi2_of_mu_cand, df_fd.prot_chi2_of_prot_cand, df_fd.mu_len)]
    sbnd_comp, icarus_comp = make_all_plots(df_nd, df_fd, "1mu+1p Cut", mode_labels, top_labels, det_labels)
    comp_sbnd.append(sbnd_comp)
    comp_icarus.append(icarus_comp)

    ### crthitveto cut
    df_fd = df_fd[crthitveto_cut(df_fd)]
    sbnd_comp, icarus_comp = make_all_plots(df_nd, df_fd, "CRT Veto", mode_labels, top_labels, det_labels)
    comp_sbnd.append(sbnd_comp)
    comp_icarus.append(icarus_comp)

    # Composition plots
    plot_composition(np.array(comp_sbnd), cut_labels, top_labels, 'SBND', 'SBNDComp.png')
    plot_composition(np.array(comp_icarus), cut_labels, top_labels, 'ICARUS', 'ICARUSComp.png')

    return df_nd, df_fd

def main():
    """Main analysis pipeline."""
    dffile_nd = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-3-fix-zexp-again-again/SBND_SpringMC_rewgt_10.df"
    dffile_fd = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-3-fix-zexp-again-again/ICARUS_SpringMC_Dev_wgt.df"

    df_nd, df_nd_hdr, df_nd_stub = load_data(dffile_nd)
    df_fd, df_fd_hdr, df_fd_stub = load_data(dffile_fd)

    des_nd_POT = 1e20
    des_fd_POT = 5e20

    nd_POT = scale_pot(df_nd, df_nd_hdr, des_nd_POT)
    fd_POT = scale_pot(df_fd, df_fd_hdr, des_fd_POT)
    df_nd, df_fd = apply_cuts(df_nd, df_nd_stub, df_fd, df_fd_stub, nd_POT, fd_POT, dffile_nd, dffile_fd)

if __name__ == "__main__":
    main()
