# Standard library imports
import sys
import os

# Third-party imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Add the head direcoty to sys.path
workspace_root = os.getcwd()
sys.path.insert(0, workspace_root + "/../../")
# Local imports
import kinematics
from makedf.util import *
from gump_cuts import *

import matplotlib as mpl

class PlotObj:
  def __init__(self, title, label):
    self.title = title
    self.label = label

# Colors for plots
HAWKS_COLORS = ["#315031", "#d54c28", "#1e3f54", "#c89648", "#43140b", "#95af8b"]
FONTSIZE = 14
plt.style.use('dune.mplstyle')

def make_all_plots(df_nd, df_fd, cut_stage, mode_labels, top_labels, det_labels):
    sbnd_title = f"{cut_stage}"
    icarus_title = f"{cut_stage}"

    sbnd_file_p = f"intdelp_{cut_stage.lower().replace(' ', '_')}_sbnd.png"
    icarus_file_p = f"intdelp_{cut_stage.lower().replace(' ', '_')}_icarus.png"
    plot_int(df_nd, 'del_p', sbnd_title, sbnd_file_p, r"$\delta p$ [GeV/c]", mode_labels, det_labels[0])
    plot_int(df_fd, 'del_p', icarus_title, icarus_file_p, r"$\delta p$ [GeV/c]", mode_labels, det_labels[1])

    sbnd_file_pT = f"intdelpT_{cut_stage.lower().replace(' ', '_')}_sbnd.png"
    icarus_file_pT = f"intdelpT_{cut_stage.lower().replace(' ', '_')}_icarus.png"
    plot_int(df_nd, 'del_Tp', sbnd_title, sbnd_file_pT, r"$\delta p_{T}$ [GeV/c]", mode_labels, det_labels[0])
    plot_int(df_fd, 'del_Tp', icarus_title, icarus_file_pT, r"$\delta p_{T}$ [GeV/c]", mode_labels, det_labels[1])

    sbnd_file_E = f"intE_{cut_stage.lower().replace(' ', '_')}_sbnd.png"
    icarus_file_E = f"intE_{cut_stage.lower().replace(' ', '_')}_icarus.png"
    plot_int(df_nd, 'nu_E_calo', sbnd_title, sbnd_file_E, r"$E_{reco}$ [GeV]", mode_labels, det_labels[0])
    plot_int(df_fd, 'nu_E_calo', icarus_title, icarus_file_E, r"$E_{reco}$ [GeV]", mode_labels, det_labels[1])

    sbnd_fs_file = f"fs_{cut_stage.lower().replace(' ', '_')}_sbnd.png"
    icarus_fs_file = f"fs_{cut_stage.lower().replace(' ', '_')}_icarus.png"
    plot_fs(df_nd, sbnd_title, sbnd_fs_file, det_labels[0])
    plot_fs(df_fd, icarus_title, icarus_fs_file, det_labels[1])

    sbnd_file_top = f"topE_{cut_stage.lower().replace(' ', '_')}_sbnd.png"
    icarus_file_top = f"topE_{cut_stage.lower().replace(' ', '_')}_icarus.png"
    sbnd_comp = plot_top(df_nd, 'nu_E_calo', sbnd_title, sbnd_file_top, r"$E_{reco}$ [GeV]", top_labels, det_labels[0], eff_bool=True)
    icarus_comp = plot_top(df_fd, 'nu_E_calo', icarus_title, icarus_file_top, r"$E_{reco}$ [GeV]", top_labels, det_labels[1], eff_bool=True)
    return sbnd_comp, icarus_comp

def plot_top(df, var, title, outfile, label, mode_labels, det, eff_bool=False):
    var_data = df[var]
    og_sig_ct = df['og_sig_ct'].to_numpy()[0]
    glob_scale = df['glob_scale'].to_numpy()[0]
    tvar = breakdown_top(var_data, df)

    # Do purity/efficiency calculation
    ret = [len(t) for t in tvar]
    sig_ct = ret[0] 
    tot_ct = sum(ret)
    eff = sig_ct/og_sig_ct
    pur = sig_ct/tot_ct

    fig, ax = plt.subplots(1, 1, figsize=(6, 4.5))
    plt.xlabel(label)
    plt.ylabel('Events')

    n, bins, _ = plt.hist(tvar, bins=np.linspace(0,2.5,21), stacked=True, label=top_labels, 
                        color=HAWKS_COLORS, weights=[glob_scale*np.ones_like(t) for t in tvar])

    plt.title(f"$\\bf{{{det}}}$  {title}")
    plt.legend()

    if eff_bool:
        ax.text(0.6, 0.33, f"Purity {{:.2f}}%".format(100*pur), 
                transform=ax.transAxes, fontsize=FONTSIZE)
        ax.text(0.6, 0.4, f"Efficiency {{:.2f}}%".format(100*eff), 
                transform=ax.transAxes, fontsize=FONTSIZE)

    plt.savefig("EventSelectionPlots/"+outfile, bbox_inches='tight')
    plt.clf()
    plt.close()
    return ret/np.sum(ret)

def plot_int(df, var, title, outfile, label, mode_labels, det):
    var_data = df[var]
    og_sig_ct = df['og_sig_ct'].to_numpy()[0]
    glob_scale = df['glob_scale'].to_numpy()[0]
    pvar = breakdown_mode(var_data, df)

    fig, ax = plt.subplots(1, 1, figsize=(6, 4.5))
    plt.xlabel(label)
    plt.ylabel('Events')

    n, bins, _ = plt.hist(pvar, bins=np.linspace(0,2,21), stacked=True, label=mode_labels, 
                      color=HAWKS_COLORS, weights=[glob_scale*np.ones_like(p) for p in pvar])

    plt.title(f"$\\bf{{{det}}}$  {title}")
    plt.legend()

    plt.savefig("EventSelectionPlots/"+outfile, bbox_inches='tight')
    plt.clf()
    plt.close()

def plot_fs(df, title, outfile, det):
    og_sig_ct = df['og_sig_ct'].to_numpy()[0]
    glob_scale = df['glob_scale'].to_numpy()[0]

    fig, ax = plt.subplots(1, 1, figsize=(6, 4.5))

    plt.xlabel("Final State Particle Count")
    plt.ylabel("Final State Particle Count")

    prim_v = ['np', 'nn', 'nmu', 'npi', 'npi0']#, 'ng']
    prim_v_labels = ['Protons', 'Neutrons', 'Muons', 'Charged Pions', 'Neutral Pions']#, 'Gammas']

    b = np.arange(1, 11, 1)
    bin_centers = 0.5 * (b[:-1] + b[1:])
    bin_numbers = range(1, len(bin_centers)+1)

    df_fs = [df[p] for p in prim_v]
    ax.set_xticks(bin_centers, [str(i) for i in bin_numbers])
    plt.hist(df_fs, bins=b, stacked=True, label=prim_v_labels, color=HAWKS_COLORS[:-1], 
            weights=[glob_scale*np.ones_like(p) for p in df_fs])

    plt.title(f"$\\bf{{{det}}}$  {title}")
    plt.legend()

    plt.savefig("EventSelectionPlots/"+outfile, bbox_inches='tight')
    plt.clf()
    plt.close()

def plot_PID_cut(var, cut_vals, title, outfile, xlims=[0, 100], 
                 labels=['SBND', 'ICARUS'], arrow_dir='None', arrow_txt='None'):

    var_nd, var_fd = var
    cut_val_nd, cut_val_fd = cut_vals
    b = np.linspace(xlims[0], xlims[1], 40)


    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 4.5), sharex=True)
    ax1.set_title(title+" Cut")
    plt.subplots_adjust(hspace=0)

    n1, bins, _ = ax1.hist(var_nd[0], bins=b, histtype='step', color=HAWKS_COLORS[0], 
                           label=labels[0], stacked=False, weights=[1/len(var_nd[0])]*len(var_nd[0]), linestyle='-')
    n2, bins, _ = ax1.hist(var_nd[1], bins=b, histtype='step', color=HAWKS_COLORS[0], 
                           label=labels[1], stacked=False, weights=[1/len(var_nd[1])]*len(var_nd[1]), linestyle='--')
    n3, bins, _ = ax2.hist(var_fd[0], bins=b, histtype='step', color=HAWKS_COLORS[1], 
                           label=labels[2], stacked=False, weights=[1/len(var_fd[0])]*len(var_fd[0]), linestyle='-')
    n4, bins, _ = ax2.hist(var_fd[1], bins=b, histtype='step', color=HAWKS_COLORS[1], 
                           label=labels[3], stacked=False, weights=[1/len(var_fd[1])]*len(var_fd[1]), linestyle='--')

    for ax in [ax1, ax2]:
        ax.set_ylabel('Events')
    ax2.set_xlabel(title)
    ax1.legend()
    ax2.legend()

    max_1 = np.max([n1, n2])*1.2
    max_2 = np.max([n3, n4])*1.2
    ax1.vlines(cut_val_nd, 0, max_1, colors='black', linestyle='--')
    ax2.vlines(cut_val_fd, 0, max_2, colors='black', linestyle='--')
    ax1.set_ylim(0, max_1)
    ax2.set_ylim(0, max_2)

    if arrow_dir == 'right':
       ax1.arrow(cut_val_nd, max_1/2, (xlims[1]-xlims[0])/8, 0, head_width=0.1*max_1, 
                 head_length=(xlims[1]-xlims[0])/40, fc='red', ec='red')
       ax2.arrow(cut_val_fd, max_2/2, (xlims[1]-xlims[0])/8, 0, head_width=0.1*max_2, 
                 head_length=(xlims[1]-xlims[0])/40, fc='red', ec='red')
       ax1.text(cut_val_nd+10, 0.35*max_1, arrow_txt, color='red')
       ax2.text(cut_val_fd+10, 0.35*max_2, arrow_txt, color='red')
    elif arrow_dir == 'left':           
       ax1.arrow(cut_val_nd, max_1/2, -(xlims[1]-xlims[0])/8, 0, head_width=0.1*max_1, 
                 head_length=(xlims[1]-xlims[0])/40, fc='red', ec='red')
       ax2.arrow(cut_val_fd, max_2/2, -(xlims[1]-xlims[0])/8, 0, head_width=0.1*max_2, 
                 head_length=(xlims[1]-xlims[0])/40, fc='red', ec='red')
       ax1.text(cut_val_nd-20, 0.8*max_1, arrow_txt, color='red')
       ax2.text(cut_val_fd-20, 0.8*max_2, arrow_txt, color='red')

    plt.tight_layout()
    plt.savefig("EventSelectionPlots/"+outfile, bbox_inches='tight')
    plt.clf()
    plt.close()

def plot_nuscore_cut(var_nd, cut_vals_nd, var_fd, cut_vals_fd, title, outfile, xlims=[0, 100], 
                     labels=['SBND', 'ICARUS'], arrow_dir='None', arrow_txt='None'):
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 4.5), sharex=True)
    plt.subplots_adjust(hspace=0)
    
    b = np.linspace(xlims[0], xlims[1], 40)
    n1, bins, _ = ax1.hist(var_nd[0], bins=b, histtype='step', color=HAWKS_COLORS[0], 
                           label=labels[0], stacked=False, weights=[1/len(var_nd[0])]*len(var_nd[0]), linestyle='-')
    n2, bins, _ = ax1.hist(var_nd[1], bins=b, histtype='step', color=HAWKS_COLORS[0], 
                           label=labels[1], stacked=False, weights=[1/len(var_nd[1])]*len(var_nd[1]), linestyle='--')
    n3, bins, _ = ax2.hist(var_fd[0], bins=b, histtype='step', color=HAWKS_COLORS[1], 
                           label=labels[2], stacked=False, weights=[1/len(var_fd[0])]*len(var_fd[0]), linestyle='-')

    n4, bins, _ = ax2.hist(var_fd[1], bins=b, histtype='step', color=HAWKS_COLORS[1], 
                           label=labels[3], stacked=False, weights=[1/len(var_fd[1])]*len(var_fd[1]), linestyle='--')

    # Style
    for ax in [ax1, ax2]:
        ax.tick_params(axis='both', which='both', direction='in', length=6, width=1.5,  top=True, right=True)
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)
        ax.set_ylabel('Events')
    ax2.set_xlabel('Nu Score')
    ax1.legend()
    ax2.legend()

    max_1 = np.max([n1, n2])*3.0
    max_2 = np.max(n4)*12.0
    ax1.vlines(cut_vals_nd, 0, max_1, colors='black', linestyle='--')
    ax2.vlines(cut_vals_fd, 0, max_2, colors='black', linestyle='--')
    ax1.set_ylim(0, max_1)
    ax2.set_ylim(0, max_2)

    ax1.arrow(cut_vals_nd[0], max_1/2, (xlims[1]-xlims[0])/8, 0, head_width=0.1*max_1, 
              head_length=(xlims[1]-xlims[0])/40, fc='red', ec='red')
    ax2.arrow(cut_vals_fd[0], max_2/2, (xlims[1]-xlims[0])/8, 0, head_width=0.1*max_2, 
              head_length=(xlims[1]-xlims[0])/40, fc='red', ec='red')
    ax1.text(cut_vals_nd[0]+(xlims[1]-xlims[0])/20, 0.35*max_1, arrow_txt, color='red')
    ax2.text(cut_vals_fd[0]+(xlims[1]-xlims[0])/20, 0.35*max_2, arrow_txt, color='red')

    plt.tight_layout()
    plt.savefig("EventSelectionPlots/"+outfile, bbox_inches='tight')
    plt.clf()
    plt.close()

def plot_stub_2d(length, dqdx, outfile, title='test'):
    fig, ax = plt.subplots(figsize=(8, 6), dpi=80)
    dqdx_cuts = np.array([5.5e5, 3.5e5, 3e5, 2e5, 0])
    length_cuts = np.array([0, 0.5, 1, 2, 4])
    plt.hlines(dqdx_cuts[:-1], length_cuts[1:], length_cuts[:-1], color='red', lw=1.0, linestyle='--')
    plt.vlines(length_cuts[1:], dqdx_cuts[1:], dqdx_cuts[:-1], color='red', lw=1.0, linestyle='--')

    plt.text(0.5, 2e5, 'Select', color='red', fontsize=20)
    plt.text(1.8, 6e5, 'Remove', color='red', fontsize=20)

    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

    plt.title(title)
    plt.tick_params(axis='both', which='both', direction='in', length=6, width=1.5, top=True, right=True)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel('Stub Length [cm]')
    plt.ylabel('dQ/dx [electrons/cm]')
    h = plt.hist2d(length, dqdx, bins=[8, 16], range=[[0, 4],[0, 8e5]])
    plt.gca().yaxis.get_offset_text().set_fontsize(FONTSIZE)
    cbar = plt.colorbar(h[3])

    plt.tight_layout()
    plt.savefig("EventSelectionPlots/"+outfile)
    plt.clf()
    plt.close()

def plot_composition(percentages, time_labels=None, components=None, title='', outfile='comp.png'):
    percentages = np.flip(np.array(percentages), axis=0)
    num_time_points, num_components = percentages.shape

    y = np.arange(num_time_points)

    left = np.zeros(num_time_points)
    for i in range(num_components):
        plt.barh(y, percentages[:, i], left=left, label=components[i])
        left += percentages[:, i]

    plt.yticks(y, labels=reversed(time_labels))
    plt.tick_params(axis='both', which='both', direction='in', length=6, width=1.5, top=True, right=True)
    plt.xlabel('Percentage')
    plt.xlim(0, 1)
    plt.legend(frameon=True)
    plt.title('Purity, '+title)
    plt.tight_layout()
    plt.savefig("EventSelectionPlots/"+outfile)
    plt.clf()
