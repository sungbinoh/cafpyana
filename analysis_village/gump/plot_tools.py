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
from sel_tools import *

import matplotlib as mpl

class PlotObj:
  def __init__(self, title, label):
    self.title = title
    self.label = label

# Colors for plots
HAWKS_COLORS = ["#315031", "#d54c28", "#1e3f54", "#c89648", "#43140b", "#95af8b"]

top_labels = ["Signal",
              "Other numu CC",
              "NC",
              "Out of FV",
              "Cosmic",
              "Other"]

FONTSIZE = 14

def plot_int(df, var, title, outfile, label, mode_labels, det, eff_bool=False):
    var_data = df[var]
    og_sig_ct = df['og_sig_ct'].to_numpy()[0]
    glob_scale = df['glob_scale'].to_numpy()[0]

    # Use SBND style: thick axes, ticks in, bold label, no 'Work in Progress', detector name bold
    plt.style.use('default')
    fig, ax = plt.subplots(1, 1, figsize=(6, 4.5))
    pvar = breakdown_mode(var_data, df)
    n, bins, _ = ax.hist(pvar, bins=np.linspace(0,2,21), stacked=True, label=mode_labels, 
                        color=HAWKS_COLORS, weights=[glob_scale*np.ones_like(p) for p in pvar])

    # Style
    ax.tick_params(axis='both', which='both', direction='in', length=6, width=1.5, labelsize=FONTSIZE, top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    ax.set_xlabel(label, fontsize=FONTSIZE, fontweight='bold')
    ax.set_ylabel('Events', fontsize=FONTSIZE, fontweight='bold')
    ax.set_title(f"$\\bf{{{det}}}$  {title}", fontsize=FONTSIZE+2)
    plt.legend(fontsize=FONTSIZE)

    sig_ct = len(var_data[is_signal(df, det)])
    if eff_bool:
        ax.text(0.5, 0.33, f"Purity {{:.2f}}%".format(100*sig_ct/len(df)), 
                transform=ax.transAxes, fontsize=FONTSIZE)
        ax.text(0.5, 0.4, f"Efficiency {{:.2f}}%".format(100*sig_ct/og_sig_ct), 
                transform=ax.transAxes, fontsize=FONTSIZE)

    plt.tight_layout()
    plt.gca().yaxis.get_offset_text().set_fontsize(FONTSIZE)
    plt.savefig("EventSelectionPlots/"+outfile, bbox_inches='tight')
    plt.clf()
    plt.close()

    tot = sum(n[-1])
    ret = [sum(n[0])/tot]
    for i in range(1, len(n)):
        ret.append((sum(n[i])-sum(n[i-1]))/tot)
    return ret

def plot_fs(df, var, title, outfile, det):
    og_sig_ct = df['og_sig_ct'].to_numpy()[0]
    glob_scale = df['glob_scale'].to_numpy()[0]

    plt.style.use('default')
    fig, ax = plt.subplots(1, 1, figsize=(6, 4.5))
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    prim_v = ['np', 'nn', 'nmu', 'npi', 'npi0', 'ng']
    prim_v_labels = ['Protons', 'Neutrons', 'Muons', 'Charged Pions', 'Neutral Pions', 'Gammas']
    b = np.arange(1, 11, 1)
    bin_centers = 0.5 * (b[:-1] + b[1:])
    bin_numbers = range(1, len(bin_centers)+1)

    df_fs = [df[p] for p in prim_v]
    ax.set_xticks(bin_centers, [str(i) for i in bin_numbers])
    ax.hist(df_fs, bins=b, stacked=True, label=prim_v_labels, color=HAWKS_COLORS, 
            weights=[glob_scale*np.ones_like(p) for p in df_fs])

    # Style
    ax.tick_params(axis='both', which='both', direction='in', length=6, width=1.5, labelsize=FONTSIZE, top=True, right=True)
    ax.set_xlabel("Final State Particle Count", fontsize=FONTSIZE, fontweight='bold')
    ax.set_ylabel('Events', fontsize=FONTSIZE, fontweight='bold')
    ax.set_title(f"$\\bf{{{det}}}$  {title}", fontsize=FONTSIZE+2)
    plt.legend(fontsize=FONTSIZE)
    plt.gca().yaxis.get_offset_text().set_fontsize(FONTSIZE)
    plt.tight_layout()
    plt.savefig("EventSelectionPlots/"+outfile, bbox_inches='tight')
    plt.clf()
    plt.close()

def plot_PID_cut(var, cut_vals, title, outfile, xlims=[0, 100], 
                 labels=['SBND', 'ICARUS'], arrow_dir='None', arrow_txt='None'):
    
    var_nd, var_fd = var
    cut_val_nd, cut_val_fd = cut_vals
    b = np.linspace(xlims[0], xlims[1], 40)
    plt.style.use('default')
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 4.5), sharex=True)
    plt.subplots_adjust(hspace=0)

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
        ax.tick_params(axis='both', which='both', direction='in', length=6, width=1.5, labelsize=FONTSIZE, top=True, right=True)
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)
        ax.set_ylabel('Events', fontsize=FONTSIZE, fontweight='bold')
    ax2.set_xlabel('x label', fontsize=FONTSIZE, fontweight='bold')
    ax1.legend(fontsize=FONTSIZE)
    ax2.legend(fontsize=FONTSIZE)

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

    fig.suptitle(title, fontsize=FONTSIZE+2, fontweight='bold', x=0.01, ha='left')
    plt.tight_layout()
    plt.savefig("EventSelectionPlots/"+outfile, bbox_inches='tight')
    plt.clf()
    plt.close()

def plot_nuscore_cut(var_nd, cut_vals_nd, var_fd, cut_vals_fd, title, outfile, xlims=[0, 100], 
                     labels=['SBND', 'ICARUS'], arrow_dir='None', arrow_txt='None'):
    
    plt.style.use('default')
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 4.5), sharex=True)
    plt.subplots_adjust(hspace=0)
    
    b = np.linspace(xlims[0], xlims[1], 40)
    n1, bins, _ = ax1.hist(var_nd[0], bins=b, histtype='step', color=HAWKS_COLORS[0], 
                           label=labels[0], stacked=False, weights=[1/len(var_nd[0])]*len(var_nd[0]), linestyle='-')
    n2, bins, _ = ax1.hist(var_nd[1], bins=b, histtype='step', color=HAWKS_COLORS[0], 
                           label=labels[1], stacked=False, weights=[1/len(var_nd[1])]*len(var_nd[1]), linestyle='--')
    n4, bins, _ = ax2.hist(var_fd[1], bins=b, histtype='step', color=HAWKS_COLORS[1], 
                           label=labels[3], stacked=False, weights=[1/len(var_fd[1])]*len(var_fd[1]), linestyle='--')

    # Style
    for ax in [ax1, ax2]:
        ax.tick_params(axis='both', which='both', direction='in', length=6, width=1.5, labelsize=FONTSIZE, top=True, right=True)
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)
        ax.set_ylabel('Events', fontsize=FONTSIZE, fontweight='bold')
    ax2.set_xlabel('x label', fontsize=FONTSIZE, fontweight='bold')
    ax1.legend(fontsize=FONTSIZE)
    ax2.legend(fontsize=FONTSIZE)

    max_1 = np.max([n1, n2])*3.0
    max_2 = np.max(n4)*3.0
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

    fig.suptitle(title, fontsize=FONTSIZE+2, fontweight='bold', x=0.01, ha='left')
    plt.tight_layout()
    plt.savefig("EventSelectionPlots/"+outfile, bbox_inches='tight')
    plt.clf()
    plt.close()

def plot_cut(var_nd, cut_vals_nd, var_fd, cut_vals_fd, title, outfile, xlims=[0, 100], 
             labels=['SBND', 'ICARUS'], arrow_dir='None', arrow_txt='None'):
    plt.style.use('default')
    b = np.linspace(xlims[0], xlims[1], 40)
    fig, ax1 = plt.subplots(1, 1, figsize=(6, 4.5))
    n, bins, _ = ax1.hist(var_nd, bins=b, histtype='step', color=HAWKS_COLORS[0], label='SBND')
    n, bins, _ = ax1.hist(var_fd, bins=b, histtype='step', color=HAWKS_COLORS[1], label='ICARUS')
    ax1.vlines(cut_vals_fd, 0, np.max(n), colors='black', linestyle='--')
    # Style
    ax1.tick_params(axis='both', which='both', direction='in', length=6, width=1.5, labelsize=FONTSIZE, top=True, right=True)
    for spine in ax1.spines.values():
        spine.set_linewidth(1.5)
    ax1.set_xlabel('x label', fontsize=FONTSIZE, fontweight='bold')
    ax1.set_ylabel('Events', fontsize=FONTSIZE, fontweight='bold')
    fig.suptitle(title, fontsize=FONTSIZE+2, fontweight='bold', x=0.01, ha='left')
    plt.legend(fontsize=FONTSIZE)
    plt.tight_layout()
    plt.savefig("EventSelectionPlots/"+outfile, bbox_inches='tight')
    plt.clf()
    plt.close()

def plot_stub_2d(df, cuts, outfile, title='test'):
    plt.style.use('default')
    dQdx = []
    x = []
    for l in ["0_5", "1", "2", "3"]:
        x.extend(df.stub["l"+l+"cm"].length)
        dQdx.extend(df.stub["l"+l+"cm"].Q/df.stub["l"+l+"cm"].length)

    fig, ax = plt.subplots(figsize=(8, 6), dpi=80)
    plt.hlines(cuts, [0, 0.5, 1.0, 2.0], [0.5, 1.0, 2.0, 3.0], color='red', lw=1.0, linestyle='--')
    plt.vlines([0.5, 1.0, 2.0], cuts[1:], cuts[:-1], color='red', lw=1.0, linestyle='--')

    plt.text(0.5, 2e5, 'Select', color='red', fontsize=FONTSIZE*2)
    plt.text(1.8, 6e5, 'Remove', color='red', fontsize=FONTSIZE*2)

    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    plt.title(title, fontsize=FONTSIZE+5, fontweight='bold', loc='left')
    plt.tick_params(axis='both', which='both', direction='in', length=6, width=1.5, labelsize=FONTSIZE, top=True, right=True)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel('Stub Length [cm]', fontsize=FONTSIZE, fontweight='bold')
    plt.ylabel('dQ/dx [electrons/cm]', fontsize=FONTSIZE, fontweight='bold')
    h = plt.hist2d(x, dQdx, bins=[6, 16], range=[[0, 3],[0, 8e5]])
    plt.gca().yaxis.get_offset_text().set_fontsize(FONTSIZE)
    cbar = plt.colorbar(h[3])
    cbar.ax.tick_params(labelsize=FONTSIZE)

    plt.tight_layout()
    plt.savefig("EventSelectionPlots/"+outfile)
    plt.clf()
    plt.close()

def plot_composition(percentages, time_labels=None, components=None, title='', outfile='comp.png'):
    percentages = np.flip(np.array(percentages), axis=0)
    num_time_points, num_components = percentages.shape

    if components is None:
        components = [f"Component {i+1}" for i in range(num_components)]
    if time_labels is None:
        time_labels = [f"Time {i+1}" for i in range(num_time_points)]

    plt.style.use('default')
    y = np.arange(num_time_points)

    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

    left = np.zeros(num_time_points)
    for i in range(num_components):
        plt.barh(y, percentages[:, i], left=left, label=components[i])
        left += percentages[:, i]

    plt.yticks(y, labels=reversed(time_labels), fontsize=FONTSIZE)
    plt.tick_params(axis='both', which='both', direction='in', length=6, width=1.5, labelsize=FONTSIZE, top=True, right=True)
    plt.xlabel('Percentage', fontsize=FONTSIZE, fontweight='bold')
    plt.xlim(0, 1)
    plt.legend(fontsize=FONTSIZE)
    plt.title('Purity, '+title, fontsize=FONTSIZE+2, fontweight='bold', loc='left')
    plt.tight_layout()
    plt.savefig("EventSelectionPlots/"+outfile)
    plt.clf()
