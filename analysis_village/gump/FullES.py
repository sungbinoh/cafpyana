import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from plot_tools import *
# Add the head direcoty to sys.path
workspace_root = os.getcwd()  
sys.path.insert(0, workspace_root + "/../../")

# import this repo's classes
from makedf.util import *

import kinematics

def load_data(file, start=0, stop=10000):
    """Load event, header, and mcnu data from HDF file."""
    df_evt = pd.read_hdf(file, "evt_0")
    df_hdr = pd.read_hdf(file, "hdr_0")
    return df_evt, df_hdr

def scale_pot(df, df_hdr, desired_pot):
    """Scale DataFrame by desired POT."""
    pot = sum(df_hdr.pot.tolist())
    print(f"POT: {pot}\nScaling to: {desired_pot}")
    scale = desired_pot / pot
    df['glob_scale'] = scale
    return pot, scale


def fv_cut(df_nd, df_fd):
    plot_cut(df_nd.slc.vertex.x, [SBNDFVCuts['x']['min'], SBNDFVCuts['x']['max']],
             df_fd.slc.vertex.x, [ICARUSFVCuts['C0']['x']['min'], ICARUSFVCuts['C0']['x']['max'],
                                 ICARUSFVCuts['C1']['x']['min'], ICARUSFVCuts['C1']['x']['max']],
             "Vertex X", "FVCutVarX.png", [-300, 300])
    plot_cut(df_nd.slc.vertex.y, [SBNDFVCuts['y']['min'], SBNDFVCuts['y']['max']],
             df_fd.slc.vertex.y, [ICARUSFVCuts['C0']['y']['min'], ICARUSFVCuts['C0']['y']['max'],
                                 ICARUSFVCuts['C1']['y']['min'], ICARUSFVCuts['C1']['y']['max']],
             "Vertex Y", "FVCutVarY.png", [-250, 250])
    plot_cut(df_nd.slc.vertex.z, [SBNDFVCuts['z']['min'], SBNDFVCuts['z']['max']],
             df_fd.slc.vertex.z, [ICARUSFVCuts['C0']['z']['min'], ICARUSFVCuts['C0']['z']['max'],
                                 ICARUSFVCuts['C1']['z']['min'], ICARUSFVCuts['C1']['z']['max']],
             "Vertex Z", "FVCutVarZ.png", [-1000, 1000])
    df_nd = df_nd[SelFV(df_nd.slc.vertex, "SBND")]
    df_fd = df_fd[SelFV(df_fd.slc.vertex, "ICARUS")]
    return df_nd, df_fd

def cosmic_cut(df_nd, df_fd):
    plot_nuscore_cut([
        df_nd.slc.nu_score[is_cosmic(df_nd)], df_nd.slc.nu_score[np.invert(is_cosmic(df_nd))]
    ], [0.5], [
        df_fd.slc.nu_score[is_cosmic(df_fd)], df_fd.slc.nu_score[np.invert(is_cosmic(df_fd))]
    ], [0.5], "NuScore", "NuScore.png", [0, 1],
        ['SBND Cosmic', 'SBND Neutrino', 'ICARUS Cosmic', 'ICARUS Neutrino'], 'right', 'select as neutrino')
    df_nd = df_nd[df_nd.slc.nu_score > 0.5]
    df_fd = df_fd[df_fd.slc.nu_score > 0.5]
    return df_nd, df_fd

def twoprong_cut(df_nd, df_fd):
    df_nd = df_nd[np.isnan(df_nd.other_shw_length) & np.isnan(df_nd.other_trk_length)]
    df_fd = df_fd[np.isnan(df_fd.other_shw_length) & np.isnan(df_fd.other_trk_length)]
    return df_nd, df_fd

def pid_cut(df_nd, df_fd):

    plot_PID_cut([[df_nd.p.pfp.trk.chi2pid.I2.chi2_muon, df_nd.mu.pfp.trk.chi2pid.I2.chi2_muon], 
                  [df_fd.p.pfp.trk.chi2pid.I2.chi2_muon, df_fd.mu.pfp.trk.chi2pid.I2.chi2_muon]], 
                  [25, 25], "MuScore", "MuScore.png", [0, 65],
                  ['SBND, Protons', 'SBND, Muons', 'ICARUS, Protons', 'ICARUS, Muons'], 
                  'left', 'select as muon')

    plot_PID_cut([[df_nd.p.pfp.trk.chi2pid.I2.chi2_proton, df_nd.mu.pfp.trk.chi2pid.I2.chi2_proton], 
                  [df_fd.p.pfp.trk.chi2pid.I2.chi2_proton, df_fd.mu.pfp.trk.chi2pid.I2.chi2_proton]],
                  [100, 100], "ProtonScore", "ProtonScore.png", [0, 200],
                  ['SBND, Protons', 'SBND, Muons', 'ICARUS, Protons', 'ICARUS, Muons'], 
                  'right', 'select as muon')

    plot_cut(df_nd.mu.pfp.trk.len, [50], df_fd.mu.pfp.trk.len, [50], "MuLen", "MuonMuLen.png", [0, 500])
    MUSEL_MUSCORE_TH, MUSEL_PSCORE_TH, MUSEL_LEN_TH = 25, 100, 50
    mu_cut_nd = (df_nd.mu.pfp.trk.chi2pid.I2.chi2_muon < MUSEL_MUSCORE_TH) & \
                (df_nd.mu.pfp.trk.chi2pid.I2.chi2_proton > MUSEL_PSCORE_TH) & \
                (df_nd.mu.pfp.trk.len > MUSEL_LEN_TH)
    mu_cut_fd = (df_fd.mu.pfp.trk.chi2pid.I2.chi2_muon < MUSEL_MUSCORE_TH) & \
                (df_fd.mu.pfp.trk.chi2pid.I2.chi2_proton > MUSEL_PSCORE_TH) & \
                (df_fd.mu.pfp.trk.len > MUSEL_LEN_TH)
    PSEL_MUSCORE_TH, PSEL_PSCORE_TH = 0, 90
    p_cut_nd = (df_nd.p.pfp.trk.chi2pid.I2.chi2_muon > PSEL_MUSCORE_TH) & \
               (df_nd.p.pfp.trk.chi2pid.I2.chi2_proton < PSEL_PSCORE_TH)
    p_cut_fd = (df_fd.p.pfp.trk.chi2pid.I2.chi2_muon > PSEL_MUSCORE_TH) & \
               (df_fd.p.pfp.trk.chi2pid.I2.chi2_proton < PSEL_PSCORE_TH)
    df_nd = df_nd.loc[mu_cut_nd & p_cut_nd]
    df_fd = df_fd.loc[mu_cut_fd & p_cut_fd]
    return df_nd, df_fd

def stub_cut(df_nd, df_fd):
    l_cuts = [5.5e5, 4e5, 3e5, 0]
    plot_stub_2d(df_nd, l_cuts, "SBND.png", "SBND Stub Cut")
    plot_stub_2d(df_nd[np.invert(is_signal(df_nd, "SBND"))], l_cuts, "FPSBND.png", "False Positive SBND Stub Cut")
    plot_stub_2d(df_fd, l_cuts, "ICARUS.png", "ICARUS Stub Cut")
    plot_stub_2d(df_fd[np.invert(is_signal(df_fd, "ICARUS"))], l_cuts, "FPICARUS.png", "False Positive ICARUS Stub Cut")
    for i, l in enumerate(["0_5", "1", "2", "3"]):
        nd_cut = np.invert((df_nd.stub[f"l{l}cm"].Q / df_nd.stub[f"l{l}cm"].length) > l_cuts[i])
        fd_cut = np.invert((df_fd.stub[f"l{l}cm"].Q / df_fd.stub[f"l{l}cm"].length) > l_cuts[i])
        df_nd = df_nd.loc[nd_cut]
        df_fd = df_fd.loc[fd_cut]
    return df_nd, df_fd

def apply_cuts(df_nd, df_fd, nd_POT, fd_POT, dffile_nd, dffile_fd):
    comp_sbnd = []
    comp_icarus = []
    cut_labels = ["FV Cut", "Cosmic Cut", "Two Prong", "PID", "Stub"]
    mode_labels = ['QE', 'MEC', 'RES', 'SIS/DIS', 'COH', "other"]

    # FV Cut
    df_nd, df_fd = fv_cut(df_nd, df_fd)
    comp_sbnd.append(plot_int(df_nd, 'del_p', "FV Cut", "int0delp_fv_sbnd.png", r"$\delta p$ [GeV/c]", mode_labels, 'SBND', eff_bool=True))
    comp_icarus.append(plot_int(df_fd, 'del_p', "FV Cut", "int0delp_fv_icarus.png", r"$\delta p$ [GeV/c]", mode_labels, 'ICARUS', eff_bool=True))
    plot_int(df_nd, 'del_Tp', "FV Cut", "int0delpT_fv_sbnd.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'SBND')
    plot_int(df_fd, 'del_Tp', "FV Cut", "int0delpT_fv_icarus.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'ICARUS')
    plot_fs(df_nd, 'del_p', "FV Cut", "fs0_fv_sbnd.png", 'SBND')
    plot_fs(df_fd, 'del_p', "FV Cut", "fs0_fv_icarus.png", 'ICARUS')

    # Cosmic cut
    df_nd, df_fd = cosmic_cut(df_nd, df_fd)
    comp_sbnd.append(plot_int(df_nd, 'del_p', "Cosmic Rejection Cut", "int1delp_cosmic_rej_sbnd.png", r"$\delta p$ [GeV/c]", mode_labels, 'SBND', eff_bool=True))
    comp_icarus.append(plot_int(df_fd, 'del_p', "Cosmic Rejection Cut", "int1delp_cosmic_rej_icarus.png", r"$\delta p$ [GeV/c]", mode_labels, 'ICARUS', eff_bool=True))
    plot_int(df_nd, 'del_Tp', "Cosmic Rejection Cut", "int1delpT_cosmic_rej_sbnd.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'SBND')
    plot_int(df_fd, 'del_Tp', "Cosmic Rejection Cut", "int1delpT_cosmic_rej_icarus.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'ICARUS')
    plot_fs(df_nd, 'del_p', "Cosmic Rejection Cut", "fs1_cosmic_rej_sbnd.png", 'SBND')
    plot_fs(df_fd, 'del_p', "Cosmic Rejection Cut", "fs1_cosmic_rej_icarus.png", 'ICARUS')

    # Two prong cut
    df_nd, df_fd = twoprong_cut(df_nd, df_fd)
    comp_sbnd.append(plot_int(df_nd, 'del_p', "Two Prong Cut", "int2delp_two_prong_sbnd.png", r"$\delta p$ [GeV/c]", mode_labels, 'SBND', eff_bool=True))
    comp_icarus.append(plot_int(df_fd, 'del_p', "Two Prong Cut", "int2delp_two_prong_icarus.png", r"$\delta p$ [GeV/c]", mode_labels, 'ICARUS', eff_bool=True))
    plot_int(df_nd, 'del_Tp', "Two Prong Cut", "int2delpT_two_prong_sbnd.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'SBND')
    plot_int(df_fd, 'del_Tp', "Two Prong Cut", "int2delpT_two_prong_icarus.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'ICARUS')
    plot_fs(df_nd, 'del_p', "Two Prong Cut", "fs2_two_prong_sbnd.png", 'SBND')
    plot_fs(df_fd, 'del_p', "Two Prong Cut", "fs2_two_prong_icarus.png", 'ICARUS')

    # PID cuts
    df_nd, df_fd = pid_cut(df_nd, df_fd)
    comp_sbnd.append(plot_int(df_nd, 'del_p', "1mu+1p Cut", "int3delp_1mu1p_sbnd.png", r"$\delta p$ [GeV/c]", mode_labels, 'SBND', eff_bool=True))
    comp_icarus.append(plot_int(df_fd, 'del_p', "1mu+1p Cut", "int3delp_1mu1p_icarus.png", r"$\delta p$ [GeV/c]", mode_labels, 'ICARUS', eff_bool=True))
    plot_int(df_nd, 'del_Tp', "1mu+1p Cut", "int3delpT_1mu1p_sbnd.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'SBND')
    plot_int(df_fd, 'del_Tp', "1mu+1p Cut", "int3delpT_1mu1p_icarus.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'ICARUS')
    plot_fs(df_nd, 'del_p', "1mu+1p Cut", "fs3_1mu1p_sbnd.png", 'SBND')
    plot_fs(df_fd, 'del_p', "1mu+1p Cut", "fs3_1mu1p_icarus.png", 'ICARUS')

    # Stub cut
    df_nd, df_fd = stub_cut(df_nd, df_fd)
    comp_sbnd.append(plot_int(df_nd, 'del_p', "Stub Cut", "int4delp_stub_sbnd.png", r"$\delta p$ [GeV/c]", mode_labels, 'SBND', eff_bool=True))
    comp_icarus.append(plot_int(df_fd, 'del_p', "Stub Cut", "int4delp_stub_icarus.png", r"$\delta p$ [GeV/c]", mode_labels, 'ICARUS', eff_bool=True))
    plot_int(df_nd, 'del_Tp', "Stub Cut", "int4delpT_stub_sbnd.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'SBND')
    plot_int(df_fd, 'del_Tp', "Stub Cut", "int4delpT_stub_icarus.png", r"$\delta p_{T}$ [GeV/c]", mode_labels, 'ICARUS')
    plot_fs(df_nd, 'del_p', "Stub Cut", "fs4_stub_sbnd.png", 'SBND')
    plot_fs(df_fd, 'del_p', "Stub Cut", "fs4_stub_icarus.png", 'ICARUS')

    # Composition plots
    comp_sbnd = np.array(comp_sbnd)
    comp_icarus = np.array(comp_icarus)
    plot_composition(comp_sbnd, cut_labels, mode_labels, 'SBND', 'SBNDComp.png')
    plot_composition(comp_icarus, cut_labels, mode_labels, 'ICARUS', 'ICARUSComp.png')

    return df_nd, df_fd

def main():
    """Main analysis pipeline."""
    dffile_nd = "no_cuts.df"
    dffile_fd = "no_cuts.df"
    df_nd, df_nd_hdr = load_data(dffile_nd)
    df_fd, df_fd_hdr = load_data(dffile_fd)
    df_nd['og_sig_ct'] = len(df_nd[is_signal(df_nd, "SBND")])
    df_fd['og_sig_ct'] = len(df_fd[is_signal(df_fd, "ICARUS")])
    des_nd_POT = 1e20
    des_fd_POT = 5e20
    nd_POT, _ = scale_pot(df_nd, df_nd_hdr, des_nd_POT)
    fd_POT, _ = scale_pot(df_fd, df_fd_hdr, des_fd_POT)
    df_nd, df_fd = apply_cuts(df_nd, df_fd, nd_POT, fd_POT, dffile_nd, dffile_fd)

if __name__ == "__main__":
    main()
