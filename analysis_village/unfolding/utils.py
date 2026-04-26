import numpy as np
import pandas as pd
from scipy.stats import chi2
from tqdm import tqdm

import sys
sys.path.append('../../')
from pyanalib.split_df_helpers import *
from pyanalib.stat_helpers import *
from makedf.constants import *
from analysis_village.unfolding.unfolding_inputs import *
from analysis_village.unfolding.wienersvd import *

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
last_sys = sys.path[-1]
from analysis_village.plot_style.sbnd_style import *

##### List of functions and their purposes #####################################################
# This file is for functions preparing histograms in multiverses to be used for covariance matrices measurements
#
# get_clipped_evts: for a column and bins, clip values to include over/underflow events.
#                   This is for catching all possible smearing for the unfolding
#
# get_univ_rates: make histograms for CV + multiverses
#                 - var_config: a class defined in analysis_village.unfolding.variable_configs.py.
#                               contains binning info
#
# get_genie_univs: make histograms for CV + multiverses, but dedicated to GENIE with two modes, "rate" vs. "xsec"
#                  - "rate": true signal rate is also shifted to multiverse's rate
#                  - "xsec": true signel rate is kept as CV, but response matrix of each multiverse is multiplied to produce expected signal rate.
#                  For covariance matrices including expected signal rate, the latter should be used for x-sec extraction
#
# some plotting functions
# ...
#
# signal_hists: produce histograms after clipping under/overflow events within binnings 
#
################################################################################################ 

EPSILON = 1e-6 # clip overflow values greater than bins[-1] - EPSILON

#  ====== calculation functions ======
def get_clipped_evts(df, var_col, bins):
    var = df[var_col]
    var = np.clip(var, bins[0], bins[-1] - EPSILON)

    if 'pot_weight' in df.columns:
        # print(f"pot scale: {df.pot_weight.unique()[0]}")
        weights = df.loc[:, 'pot_weight']
    else:
        print("No pot_weight column found, return 1 as pot scale (expected for data)")
        weights = np.ones_like(var)
    return var, weights


def get_univ_rates(evtdf, 
                   var_config, 
                   syst_name,
                   n_univ_max=100,
                   ):
    var = evtdf[var_config.var_evt_reco_col]
    n_univ = len(evtdf[syst_name].columns)
    if n_univ > n_univ_max:
        n_univ = n_univ_max
    univ_events = np.zeros((n_univ, len(var_config.bin_centers)))
    print(n_univ)
    print(univ_events.shape)
    for uidx in range(n_univ):
        #weights = evtdf[syst_name]["univ_{}".format(uidx)]
        weights = evtdf[syst_name][evtdf[syst_name].columns[uidx]]
        weights = np.where(np.isnan(weights), 1, weights) # nan are non-neutrino events
        print(weights)
        n_univ, _ = np.histogram(var, bins=var_config.bins, weights=weights)
        univ_events[uidx, :] = n_univ
    cv_events, _ = np.histogram(var, bins=var_config.bins)
    print(univ_events.shape)
    return univ_events, cv_events


def get_genie_univs(cov_type="rate", 
                    evtdf=None,
                    nudf=None,
                    var_config=None, 
                    syst_name="",
                    pot=None,
                    vol=None,
                    flux=None,
                    topology_list=None,
                    n_univ=100, 
                    plot=False):
    """
    for the GENIE uncertainty on the xsec measurement
    
    Notes for variables:
      - vol: TPC volume in cm^3 to measure total number of argon targets
      - flux: Total flux of your interested neutrinos in (cm^2 POT)^-1, note that usuall SBND flux is given in (m^2 10^6 POT)^-1
    """

    if cov_type == "xsec":
        print("getting {} universes for {} uncertainty on the xsec".format(n_univ, syst_name))

        if pot is None:
            raise ValueError("pot is None, cannot make integrated flux")
        if flux is None:
            raise ValueError("flux is None, cannot make integrated flux")
        if vol is None:
            raise ValueError("vol is None, cannot make total targets")

        INTEGRATED_FLUX = pot * flux
        NTARGETS = vol * RHO * N_A / M_AR
        XSEC_UNIT = 1 / (INTEGRATED_FLUX * NTARGETS)
        scale_factor = XSEC_UNIT

        print("XSEC_UNIT: {%.3e}" %XSEC_UNIT)
        print("NTARGETS: {%.3e}" %NTARGETS)
    elif cov_type == "rate":
        print("getting {} universes for {} uncertainty on the event rate".format(n_univ, syst_name))
        scale_factor = 1.0
    else:
        raise ValueError("Invalid covariance type: {}, choose in [xsec, rate]".format(cov_type))

    if topology_list is None:
        raise ValueError("topology_list is None, can't categorize")

    bins = var_config.bins

    ret = signal_hists(evtdf, nudf, var_config, return_data=True, plot=plot)
    signal_allmc_cv = ret["nevts_allmc"]
    nevts_signal_sel_truth = ret["nevts_sel_truth"]
    signal_sel_reco_cv = ret["nevts_sel_reco"]
    signal_sel_reco_cv *= scale_factor # = Response @ true_signal

    evtdf_signal = evtdf[evtdf.nuint_categ == 1]
    nudf_signal = nudf[nudf.nuint_categ == 1]

    # reco variable histogram, topology breakdown
    evtdf_div_topo = [evtdf[evtdf.nuint_categ == mode]for mode in topology_list]

    signal_univ_events = []
    #signal_univ_effs   = []
    #signal_univ_smears = []
    bkg_univ_events    = []
    for uidx in range(n_univ):
        univ_col = syst_name + ("univ_{}".format(uidx),)

        # ---- uncertainty on the signal rate ----
        # only consider effect on the response matrix for the signal channel
        if cov_type == "xsec":
            signal_allmc_univ, _ = np.histogram(ret["var_allmc"],
                                               weights=ret["wgt_allmc"]*nudf_signal[univ_col],
                                               bins=bins)
            
            reco_vs_true = get_smear_matrix(ret["var_sel_truth"], 
                                            ret["var_sel_reco"], 
                                            weights=evtdf_signal[univ_col],
                                            bins_2d=[bins, bins],
                                            plot=plot)
            #signal_univ_smears.append(reco_vs_true)

            eff = get_eff(reco_vs_true, signal_allmc_univ) 
            #signal_univ_effs.append(eff)

            response_univ = get_response_matrix(reco_vs_true, eff, bins, plot=plot)
            signal_univ = response_univ @ signal_allmc_cv # note that we multiply the CV signal rate!
            # signal_univ = signal_cv

        elif cov_type == "rate":
            signal_univ, _ = np.histogram(ret["var_sel_reco"], 
                                          weights=evtdf_signal[univ_col],
                                          bins=bins)

        else:
            raise ValueError("Invalid covariance type: {}, choose xsec or rate".format(cov_type))

        # ---- uncertainty on the background rate ----
        # loop over background categories
        # + univ background - cv background
        # note: cv background subtraction cancels out with the cv background subtraction for the cv event rate. 
        #       doing it anyways for the plot of universes on background subtracted event rate.
        bkg_events_this_univ = []
        for this_evtdf in evtdf_div_topo[1:]:
            var, wgt = get_clipped_evts(this_evtdf, var_config.var_evt_reco_col, bins)
            univ_wgt = this_evtdf[univ_col].copy()
            univ_wgt[np.isnan(univ_wgt)] = 1 ## IMPORTANT: make nan univ_wgt to 1. to ignore them
            wgt *= univ_wgt
            #background_cv, _   = np.histogram(var, bins=bins)
            background_univ, _ = np.histogram(var, bins=bins, weights=wgt)
            #signal_univ += background_univ - background_cv
            background_univ *= scale_factor
            bkg_events_this_univ.append(background_univ)

        signal_univ *= scale_factor
        signal_univ_events.append(signal_univ)
        bkg_univ_events.append(bkg_events_this_univ) ## for each univ, save array for each bkg evt category

    ## collect bkg cv arrays
    bkg_sec_rec_cv = []
    for this_evtdf in evtdf_div_topo[1:]:
        var, wgt = get_clipped_evts(this_evtdf, var_config.var_evt_reco_col, bins)
        background_cv, _   = np.histogram(var, bins=bins)
        background_cv = background_cv.astype(np.float64)
        background_cv *= scale_factor
        bkg_sec_rec_cv.append(background_cv)

    signal_univ_events = np.array(signal_univ_events)
    bkg_univ_events = np.array(bkg_univ_events)
    bkg_sec_rec_cv = np.array(bkg_sec_rec_cv)
    
    return signal_univ_events, signal_sel_reco_cv, bkg_univ_events, bkg_sec_rec_cv

# ====== linear algebra
def collect_inv_cov(cov):
    ## covriance matrix should be positive-definite, so we can use Cholesky decomposition to get the inverse
    L = np.linalg.cholesky(cov)      # C = L L^T
    L_inv = np.linalg.inv(L)
    C_inv = L_inv.T @ L_inv
    return C_inv

# ====== plotting functions ======

def get_textloc_x(values, bins, textloc=[0.05, 0.55]):
    textloc_x, _ = textloc
    n_firsthalf = np.sum(values[:len(bins)//2])
    n_secondhalf = np.sum(values[len(bins)//2:])
    if n_firsthalf < n_secondhalf:
        textloc_x, textloc_ha = textloc_x, 'left'
    else:
        textloc_x, textloc_ha = 1-textloc_x, 'right'
    return textloc_x, textloc_ha

def add_approval_text(approval, textloc_x, textloc_y, textloc_ha):
    # SBND approval rank
    if approval == "internal":
        approval_text = r"$\mathbf{SBND}$ Internal"
        textcolor = 'rosybrown'
    elif approval == "preliminary":
        approval_text = r"$\mathbf{SBND}$ Preliminary"
        textcolor = 'gray'
    else:
        return  # no approval or unknown: do not add

    ax = plt.gcf().axes[0]  # get the first axes of the current figure
    ax.text(
        textloc_x, textloc_y, 
        approval_text, 
        transform=ax.transAxes, 
        ha=textloc_ha, va='top',
        fontsize=20, color=textcolor
    )

def overlay_hists(mc_df=None,
                  data_df=None,
                  categ_col="nuint_categ",
                  categ_list=None,
                  categ_labels=None,
                  categ_colors=None,
                  intime_df=None,
                  var_config="",
                  plot_labels=["", "", ""],
                  ax_ylim_ratio=1.5,
                  ratio = False,
                  syst = False,
                  vline = None,
                  textloc=[0.05, 0.55],
                  approval="internal",
                  save_fig=False, 
                  save_name=None): 

    ## input vars
    ## -- categs are for mc_df
    ##    -- categ_list should be an array like [-1, 1, 3, ... ]
    ##    -- categ_labels should be an array like [r'$\nu_{\mu}$ CC QE', r'$\nu_{\mu}$ CC MEC', r'$\nu_{\mu}$ CC RES', ... ]
    ##    -- categ_colors should be an array like ["#9b5580", "#390C1E", "#2c7c94", ... ]

    # ==== prepare dfs for plotting ====

    # MC
    if mc_df is not None:
        vardf, _        = get_clipped_evts(mc_df, var_config.var_evt_reco_col, var_config.bins)
        var_categ = [vardf[mc_df[categ_col] == categ] for categ in categ_list]
        weights_categ = [list(mc_df.loc[mc_df[categ_col] == categ, 'pot_weight']) for categ in categ_list]

        # MC stat err
        each_mc_hist_data = []
        each_mc_hist_err2 = []  # sum of squared weights for error
        for v, w in zip(var_categ, weights_categ):
            hist_vals, _ = np.histogram(v, weights=w, bins=var_config.bins)
            hist_err2, _ = np.histogram(v, weights=np.square(w), bins=var_config.bins)
            each_mc_hist_data.append(hist_vals)
            each_mc_hist_err2.append(hist_err2)
        total_mc = np.sum(each_mc_hist_data, axis=0)
        total_mc_err2 = np.sum(each_mc_hist_err2, axis=0)
        mc_stat_err = np.sqrt(total_mc_err2)

        # if topology breakdown, add background CV
        if breakdown_type == "topology":
            total_mc_bkgd = total_mc - each_mc_hist_data[-1]
        else:
            total_mc_bkgd = None

    else:
        vardf = None
        var_categ = None
        total_mc = None
        print("No MC data provided")

    # Data
    if data_df is not None:
        vardf_data, _   = get_clipped_evts(data_df, var_config.var_evt_reco_col, var_config.bins)
        total_data, _ = np.histogram(vardf_data, bins=var_config.bins, weights=data_df.pot_weight)
        data_eylow, data_eyhigh = return_data_stat_err(total_data)

        # data/MC
        if total_mc is not None:
            data_ratio = total_data / total_mc
            data_ratio_eylow = data_eylow / total_mc
            data_ratio_eyhigh = data_eyhigh / total_mc
            data_ratio = np.nan_to_num(data_ratio, nan=-999.)
            data_ratio_eylow = np.nan_to_num(data_ratio_eylow, nan=0.)
            data_ratio_eyhigh = np.nan_to_num(data_ratio_eyhigh, nan=0.)
        
    else:
        vardf_data = None
        total_data = None
        print("No data data provided")

    # Intime cosmics
    if intime_df is not None:
        vardf_intime, _ = get_clipped_evts(intime_df, var_config.var_evt_reco_col, var_config.bins)
        var_categ = [vardf_intime] + var_categ
        weights_categ = [list(intime_df.pot_weight)] + weights_categ
        colors = colors + ["silver"]
        labels = labels + ["In-time\nCosmic"]

    else:
        vardf_intime = None
        print("No intime cosmics provided")


    # the order of cuts from get_*_category is reversed from the order of labels and colors
    colors, labels = colors[::-1], labels[::-1]

    # ========================================================

    # ==== plot template ====
    if ratio:
        fig, axs = plt.subplots(2, 1, figsize=(8.5, 8.5), 
                               sharex=True, gridspec_kw={'height_ratios': [4, 1]})
        ax, ax_r = axs[0], axs[1]
        fig.subplots_adjust(hspace=0.1)
        ax_r.axhline(1.0, color='red', linestyle='--', linewidth=1)
        ax_r.set_ylim(0.4, 1.6)
        ax_r.set_xlabel(plot_labels[0])
        ax_r.set_ylabel("Data/MC")
        ax_r.grid(True)
        ax_r.grid(which='minor', linestyle=':', linewidth=0.5, color='gray', alpha=0.5)
        ax_r.minorticks_on()

    else:
        fig, ax = plt.subplots(figsize=(8.5, 7))
        ax.set_xlabel(plot_labels[0])

    # common formatting
    ax.set_xlim(var_config.bins[0], var_config.bins[-1])
    ax.set_ylabel(plot_labels[1])
    ax.set_title(plot_labels[2])

    # ==== Plot histograms ====

    # == rate panel ==
    # MC
    if var_categ is not None:
        mc_stack, _, _ = ax.hist(var_categ,
                                 weights=weights_categ,
                                 bins=var_config.bins,
                                 stacked=True,
                                 color=colors,
                                 label=labels,
                                 linewidth=0,
                                 edgecolor='none',
                                 histtype='stepfilled')

        breakdown_accum = [np.sum(this_mode) for this_mode in mc_stack]
        breakdown_fractions = [breakdown_accum[0]] + [(breakdown_accum[i+1] - breakdown_accum[i]) for i in range(len(breakdown_accum) - 1)]
        breakdown_fractions = [frac / breakdown_accum[-1] for frac in breakdown_fractions]

        if breakdown_type == "genie_sb":
            # hatch background portion
            bottom = np.zeros(len(var_config.bins) - 1)
            for i, (v, w, h) in enumerate(zip(var_categ, weights_categ, hatches)):
                hist_vals, _ = np.histogram(v, weights=w, bins=var_config.bins)
                ax.bar(
                    var_config.bin_centers,
                    hist_vals,
                    width=np.diff(var_config.bins),
                    bottom=bottom,
                    color='none',
                    hatch=h,
                    edgecolor='white',
                    linewidth=0.0,
                    align='center'
                )
                bottom += hist_vals

    if syst:
        syst_err = mc_stat_err
        # TODO: other syst sources
        # if 
        # mc_stat_err_frac = 
        # syst_err = np.sqrt(syst_err**2 + this_syst_err**2)

        ax.bar(
            var_config.bin_centers,
            2 * syst_err,
            width=np.diff(var_config.bins),
            bottom=total_mc - syst_err,
            facecolor='none',             # transparent fill
            hatch='xxx',                 # hatch pattern similar to ROOT's 3004
            linewidth=0.0,
            edgecolor='dimgray',            # outline color of the hatching
            label='MC Stat. Unc.'
        )
 
    else:
        print("no syst provided")

    # Data
    if vardf_data is not None:
        ax.errorbar(var_config.bin_centers, 
                    total_data, 
                    yerr=np.vstack((data_eylow, data_eyhigh)),
                    color='black', 
                    fmt='o', markersize=5, capsize=3, linewidth=1.5,
                    label='Data')

    # == ratio panel ==
    if ratio:
        # MC 
        if syst:
            mc_content_ratio = total_mc / total_mc # dummy
            mc_stat_err_ratio = syst_err / total_mc
            mc_stat_err_ratio = np.nan_to_num(mc_stat_err_ratio, nan=0.)
            ax_r.bar(
                var_config.bin_centers,
                2*mc_stat_err_ratio,
                width=np.diff(var_config.bins),
                bottom=mc_content_ratio - mc_stat_err_ratio,
                facecolor='none',
                edgecolor='dimgray',
                hatch='xxx',
                linewidth=0.0,
                label='MC Stat. Unc.'
            )
        else:
            print("No syst provided")

        # data/MC 
        ax_r.errorbar(var_config.bin_centers, data_ratio, 
                        yerr=np.vstack((data_ratio_eylow, data_ratio_eyhigh)),
                        fmt='o', color='black',
                        markersize=5, capsize=3, linewidth=1.5)

    # ===============================

    # ==== Legend ====
    # legend order: data, mc, syst
    handles, labels_orig = ax.get_legend_handles_labels()
    ordered_handles = []
    ordered_labels = []

    if data_df is not None:
        data_handle_index = labels_orig.index('Data')
        data_handle = handles[data_handle_index]
        ordered_handles.extend([data_handle])
        ordered_labels.extend(['Data'])

    if mc_df is not None:
        if breakdown_type == "genie_sb":
            # collapse S and B into a single combined legend
            for i in range(len(genie_mode_colors)):
                ordered_handles.append(Patch(facecolor=genie_mode_colors[i], edgecolor='none'))

            legend_labels = []
            i, n = 0, len(labels)
            while i < n:
                # assume paired S/B in the order of MC legend entries
                label_base = labels[i]
                if (i + 1 < n) and (labels[i + 1] == label_base): # paired S/B
                    frac1 = breakdown_fractions[i]
                    frac2 = breakdown_fractions[i+1]
                    legend_label = f"{label_base}\n({frac2*100:.1f}%/{frac1*100:.1f}%)"
                    i += 2
                else: # single
                    frac1 = breakdown_fractions[i]
                    legend_label = f"{label_base}\n({frac1*100:.1f}%)"
                    i += 1
                legend_labels.append(legend_label)
            ordered_labels.extend(legend_labels[::-1])

        else:
            mc_handles = [h for i, h in enumerate(handles) if i != data_handle_index and 'Unc.' not in labels_orig[i]]
            mc_labels = [f"{label} ({frac*100:.1f}%)"
                                for label, frac in zip(labels, breakdown_fractions)]
            ordered_handles.extend(mc_handles)
            ordered_labels.extend(mc_labels[::-1]) # note the reverse order of mc_labels

    if syst:
        unc_handle = [h for i, h in enumerate(handles) if 'Unc.' in labels_orig[i]]
        unc_label = [l for l in labels_orig if 'Unc.' in l]
        ordered_handles.extend(unc_handle)
        ordered_labels.extend(unc_label)

    ax.legend(
        ordered_handles,
        ordered_labels,
        loc='upper left',
        fontsize=12,
        frameon=False,
        ncol=3,
        bbox_to_anchor=(0.015, 0.985)
    )

    # ===============================

    # ==== additional fomatting etc. ====

    # == set y-axis limit ==
    ax.set_ylim(0., ax_ylim_ratio* np.max(total_mc))

    # == option to plot vertical lines ==
    if vline is not None:
        for v in vline:
            ymax = ax.get_ylim()[1]
            ax.vlines(x=v, ymin=0, ymax=ymax*0.75, color='red', linestyle='--')

    textloc_x, textloc_ha = get_textloc_x(total_mc, var_config.bins, textloc)
    textloc_y = textloc[1]
    add_approval_text(approval, textloc_x, textloc_y, textloc_ha)

    # GEINE version
    ax.text(textloc_x, textloc_y-0.08, r"GENIE v3.4.0 AR23_00i_00_000", transform=ax.transAxes, 
            ha=textloc_ha, va='top',
            fontsize=11.5, color='gray')


    # ===============================

    # == save figure ==
    if save_fig:
        plt.savefig(save_name+".pdf", bbox_inches='tight')
    plt.show()

    return {"cuts": cuts, 
            "total_mc": total_mc, 
            "total_mc_bkgd": total_mc_bkgd,
            "total_data": total_data}

def plot_univ_hists(univ_events, 
                    cv_events,
                    syst_name,
                    var_config,
                    categ_name = "",
                    approval="internal",
                    textloc=[0.05, 0.55],
                    save_fig=False, 
                    save_name=None,
                    ylabel="Events / Bin"): 

    assert univ_events.shape[1] == len(cv_events) 
    n_univ = univ_events.shape[0]

    if (n_univ > 10):
        colors = ["orange", "green", "red"]

        sorted_univs = np.sort(univ_events, axis=0)

        n_68 = int(0.68 * n_univ)
        start_68 = (n_univ - n_68) // 2
        end_68 = start_68 + n_68

        n_95 = int(0.95 * n_univ)
        start_95 = (n_univ - n_95) // 2
        end_95 = start_95 + n_95

        # Define bins & groupings for easy handling: [(range, color, label, skip_68)]
        segs = [
            (range(start_68, end_68), colors[0], "Universe (68%)", False),
            (range(start_95, end_95), colors[1], "Universe (95%)", True),
            ( (i for i in range(n_univ) if i not in range(start_95, end_95)), colors[2], "Universe (100%)", False)
        ]

        plotted = set()  # tracks which labels were plotted already
        for r, color, label, skip_68 in segs:
            for i in r:
                if skip_68 and i in range(start_68, end_68): 
                    continue
                show_label = label if label not in plotted else None
                plt.hist(var_config.bin_centers, bins=var_config.bins, weights=sorted_univs[i],
                        histtype="step", color=color, alpha=0.7, label=show_label)
                plotted.add(label)

    # if too few universes, just plot all of them in gray
    else:
        for i in range(n_univ):
            show_label = "Syst. #" + str(i) # if i == 0 else None
            plt.hist(var_config.bin_centers, bins=var_config.bins, weights=univ_events[i], histtype="step", label=show_label)

    # plot CV last so that it is on top of the universes
    plt.hist(var_config.bin_centers, bins=var_config.bins, weights=cv_events, histtype="step", color="black", label="Central Value", linestyle="--", linewidth=2)

    plt.xlim(var_config.bins[0], var_config.bins[-1])
    plt.xlabel(var_config.var_labels[1])
    #plt.ylim(-1, )
    plt.ylabel(ylabel)
    plt.title(syst_name + " syst.: " + categ_name)

    plt.legend(frameon=False, fontsize=12)

    # == textbox to add approval rank ==
    # decide text location based on the distribution
    textloc_x, textloc_ha = get_textloc_x(cv_events, var_config.bins, textloc)
    textloc_y = textloc[1]
    add_approval_text(approval, textloc_x, textloc_y, textloc_ha)


    if save_fig:
        plt.savefig(save_name+".pdf", bbox_inches='tight')

    plt.show()


def plot_frac_unc(frac_unc, 
                  var_config, 
                  plot_labels=["", "", ""],
                  textloc=[0.05, 0.55],
                  approval="internal",
                  save_fig=False, 
                  save_name=None):

    plt.hist(var_config.bin_centers, bins=var_config.bins, weights=frac_unc, histtype="step", color="black")

    plt.xlim(var_config.bins[0], var_config.bins[-1])
    plt.xlabel(var_config.var_labels[0])
    plt.ylabel("Fractional Uncertainty")
    plt.title(plot_labels[2])
    plt.grid(True)

    textloc_x, textloc_ha = get_textloc_x(frac_unc, var_config.bins, textloc)
    textloc_y = textloc[1]
    add_approval_text(approval, textloc_x, textloc_y, textloc_ha)

    if save_fig:
        plt.savefig(save_name+".pdf", bbox_inches='tight')
    plt.show()


def plot_heatmap(matrix, 
                 var_config, 
                 title="", 
                 save_fig=False, 
                 save_fig_name=None):

    nbins = len(var_config.bins)
    assert nbins-1 == matrix.shape[0] == matrix.shape[1]
    unif_bin = np.linspace(0., float(nbins - 1), nbins)
    extent = [unif_bin[0], unif_bin[-1], unif_bin[0], unif_bin[-1]]

    x_edges, y_edges = np.array(var_config.bins), np.array(var_config.bins)
    x_tick_positions, y_tick_positions = (unif_bin[:-1] + unif_bin[1:]) / 2, (unif_bin[:-1] + unif_bin[1:]) / 2
    x_labels, y_labels = bin_range_labels(x_edges), bin_range_labels(y_edges)

    fig, ax = plt.subplots(figsize=(10, 10))
    plt.imshow(matrix, extent=extent, origin="lower")

    for i in range(nbins-1):      # rows (y)
        for j in range(nbins-1):  # columns (x)
            value = matrix[i, j]
            if not np.isnan(value):  # skip NaNs
                plt.text(
                    j + 0.5, i + 0.5,
                    f"{value:.2e}",
                    ha="center", va="center",   
                    color=get_text_color(value),
                    fontsize=10
                )

    plt.colorbar(shrink=0.7)
    plt.xticks(x_tick_positions, x_labels, rotation=45, ha="right")
    plt.yticks(y_tick_positions, y_labels)
    plt.xlabel(var_config.var_labels[0], fontsize=20)
    plt.ylabel(var_config.var_labels[1], fontsize=20)
    plt.title(title, fontsize=20)

    if save_fig:
        plt.savefig("{}.png".format(save_fig_name), bbox_inches='tight', dpi=300)
    plt.show()

def plot_overlay_with_cov(
    data,
    cov_data,
    signal,
    cov_signal,
    background,
    cov_background,
    varcfg,
    title="",
    draw_ratio=True,
    logy=False,
    ylims=None,
    bkg_label="Background",
    sig_label="Signal",
    sig_p_bkg_label="Signal + Background",
    pred_unc_label="Pred. Unc.",
    data_label="Data",
    y_label_top="Events",
    ratio_y_min=0.,
    ratio_y_max=2.,
):
    """
    Overlay data vs (signal+background) with:
      - data error bars from cov_data
      - prediction uncertainty band from (cov_signal + cov_background)
      - stacked step-filled components for signal and background

    Args
    ----
    data, signal, background : array-like, shape (n_bins,)
    cov_data : array-like, shape (n_bins,) [variances] or (n_bins, n_bins)
    cov_signal, cov_background : array-like, shape (n_bins, n_bins)
    varcfg : VariableConfig with .bins and .var_labels
    """

    data = np.asarray(data)
    signal = np.asarray(signal)
    background = np.asarray(background)
    pred = signal + background

    edges = np.asarray(varcfg.bins)
    centers = 0.5 * (edges[:-1] + edges[1:])
    xerr = 0.5 * (edges[1:] - edges[:-1])

    n_bins = len(centers)
    if data.shape[0] != n_bins or pred.shape[0] != n_bins:
        raise ValueError(
            f"Bin mismatch: expected n_bins={n_bins} from varcfg.bins, "
            f"got data={data.shape}, pred={pred.shape}"
        )

    # --- data y-errors from cov_data ---
    cov_data = np.asarray(cov_data)
    if cov_data.ndim == 1:
        if cov_data.shape[0] != n_bins:
            raise ValueError(f"cov_data length {cov_data.shape[0]} != n_bins {n_bins}")
        yerr_data = np.sqrt(np.clip(cov_data, 0, None))
    elif cov_data.ndim == 2:
        if cov_data.shape != (n_bins, n_bins):
            raise ValueError(f"cov_data shape {cov_data.shape} != {(n_bins, n_bins)}")
        yerr_data = np.sqrt(np.clip(np.diag(cov_data), 0, None))
    else:
        raise ValueError("cov_data must be 1D (variances) or 2D (cov matrix).")

    # --- prediction uncertainty from cov_signal + cov_background ---
    cov_pred = np.asarray(cov_signal) + np.asarray(cov_background)
    if cov_pred.shape != (n_bins, n_bins):
        raise ValueError(f"cov_pred shape {cov_pred.shape} != {(n_bins, n_bins)}")

    pred_err = np.sqrt(np.clip(np.diag(cov_pred), 0, None))

    # Step arrays on edges (n_bins+1) to cover full first/last bin widths
    pred_step = np.r_[pred, pred[-1]]
    low_step  = np.r_[pred - pred_err, (pred - pred_err)[-1]]
    high_step = np.r_[pred + pred_err, (pred + pred_err)[-1]]

    # Stacked components as step arrays on edges
    bkg_step = np.r_[background, background[-1]]
    sig_step = np.r_[signal, signal[-1]]
    tot_step = np.r_[background + signal, (background + signal)[-1]]

    # --- figure layout ---
    if draw_ratio:
        fig, (ax, rax) = plt.subplots(
            2, 1, figsize=(7, 6),
            gridspec_kw={"height_ratios": [3, 1]},
            sharex=True
        )
    else:
        fig, ax = plt.subplots(figsize=(7, 4))
        rax = None

    # --- stacked histogram look (filled steps) ---
    # Background fill from 0 -> bkg
    ax.fill_between(
        edges, 0.0, bkg_step,
        step="post",
        alpha=0.35,
        label=bkg_label,
    )
    # Signal fill stacked on top: bkg -> bkg+sig
    ax.fill_between(
        edges, bkg_step, tot_step,
        step="post",
        alpha=0.35,
        label=sig_label,
    )

    # --- total prediction line + uncertainty band ---
    ax.step(
        edges, pred_step,
        where="post",
        linewidth=2.0,
        linestyle="--",
        color="blue",
        label=sig_p_bkg_label,
    )
    ax.fill_between(
        edges,
        low_step,
        high_step,
        step="post",
        facecolor="none",
        hatch="///",
        edgecolor="blue",
        linewidth=0.0,
        label=pred_unc_label,
    )

    # --- data on top ---
    ax.errorbar(
        centers, data,
        xerr=0,
        yerr=yerr_data,
        fmt="o",
        color="black",
        label=data_label,
        capsize=2,
        zorder=10,
    )

    ax.set_ylabel(y_label_top)
    ax.set_title(title)
    ax.legend(fontsize='x-small')

    # --- auto y-range from what's plotted (data±err, pred±err, plus stack bounds) ---
    y_all = np.concatenate([
        data - yerr_data, data + yerr_data,
        pred - pred_err, pred + pred_err,
        np.array([0.0]),
    ])
    y_all = y_all[np.isfinite(y_all)]

    if len(y_all) > 0:
        min_y = np.min(y_all)
        max_y = np.max(y_all)

        if logy:
            pos = y_all[y_all > 0]
            if len(pos) > 0:
                ax.set_yscale("log")
                ax.set_ylim(np.min(pos) * 0.8, max_y * 1.2)
        else:
            ymin = 0 if min_y > 0 else min_y * 1.1
            ax.set_ylim(ymin, max_y * 1.1)

    if ylims is not None:
        ax.set_ylim(ylims[0], ylims[1])

    # --- x-limits (full bin width) ---
    x_min, x_max = edges[0], edges[-1]
    ax.set_xlim(x_min, x_max)

    # --- ratio panel ---
    if draw_ratio:
        denom = np.where(pred == 0, np.nan, pred)
        ratio = data / denom

        # ratio y-errors from data covariance only
        ratio_yerr = yerr_data / denom
        ratio_step = np.r_[ratio, ratio[-1]]

        ratio_pred = pred / denom
        ratio_pred_yerr = pred_err / denom
        rlow_step = np.r_[ratio_pred - ratio_pred_yerr, (ratio_pred - ratio_pred_yerr)[-1]]
        rhigh_step = np.r_[ratio_pred + ratio_pred_yerr, (ratio_pred + ratio_pred_yerr)[-1]]

        rax.errorbar(
            centers, ratio,
            xerr=0,
            yerr=ratio_yerr,
            fmt="o",
            color="black",
            capsize=2,
        )
        rax.fill_between(edges, rlow_step, rhigh_step, step="post", alpha=0.5, facecolor="none", hatch="///", edgecolor="blue", linewidth=0.0)
        rax.axhline(1.0, linestyle="--", color="blue")

        rax.set_ylabel("Data/Pred")
        rax.set_xlabel(varcfg.var_labels[0])
        rax.set_xlim(x_min, x_max)
        rax.set_ylim(ratio_y_min, ratio_y_max)
    else:
        ax.set_xlabel(varcfg.var_labels[0])

    plt.tight_layout()
    plt.show()

def plot_unfolded_result(unfold, 
                         measured, 
                         models,
                         var_config, 
                         textloc=[0.05, 0.55],
                         approval="internal",
                         save_fig=False, 
                         save_name=None,
                         closure_test=False):

    bins = var_config.bins
    bin_centers = var_config.bin_centers
    bin_widths = np.diff(bins)

    # unfolded result
    Unfolded = unfold['unfold']
    UnfoldedCov = unfold["UnfoldCov"]
    Unfolded_perwidth = Unfolded / bin_widths

    # --- stat uncertainties
    UnfoldCov_stat = unfold['StatUnfoldCov']
    Unfold_uncert_stat = np.diag(UnfoldCov_stat)
    # --- syst uncertainties
    UnfoldCov_syst = unfold['SystUnfoldCov']
    Unfold_uncert_syst = np.diag(UnfoldCov_syst)

    # --- decompose into norm and shape components
    # the first item in models dict is the nominal input model
    norm_model = list(models.keys())[0]
    SystUnfoldCov_norm, SystUnfoldCov_shape = Matrix_Decomp(models[norm_model], UnfoldCov_syst)
    Unfold_uncert_norm = np.sqrt(np.abs(np.diag(SystUnfoldCov_norm)))
    Unfold_uncert_shape = np.sqrt(np.abs(np.diag(SystUnfoldCov_shape)))

    # --- plot
    fig, ax = plt.subplots(figsize=(8.5, 7))
    # set err to 0 for closure test
    if closure_test:
        dummy_err = np.zeros_like(Unfolded_perwidth)
        bar_handle = plt.errorbar(bin_centers, Unfolded_perwidth, yerr=dummy_err, fmt='o', color='black')

    else:
        # plot shape uncertainty as error bars
        Unfold_uncert_stat_perwidth = Unfold_uncert_stat / bin_widths
        Unfold_uncert_shape_perwidth = Unfold_uncert_shape / bin_widths
        # Unfold_uncert_stat_shape_perwidth = Unfold_uncert_stat_perwidth + Unfold_uncert_shape_perwidth
        Unfold_uncert_stat_shape_perwidth = Unfold_uncert_shape_perwidth
        bar_handle = plt.errorbar(bin_centers, Unfolded_perwidth, yerr=Unfold_uncert_stat_shape_perwidth, fmt='o', color='black')

        # plot syst norm component as histogram at the bottom
        Unfold_uncert_norm_perwidth = Unfold_uncert_norm / bin_widths
        norm_handle = plt.bar(bin_centers, Unfold_uncert_norm_perwidth, width=bin_widths, label='Syst. error (norm)', alpha=0.5, color='gray')

    # divide measured & model by bin width
    measured_perwidth = measured / bin_widths
    reco_handle, = plt.step(bins, np.append(measured_perwidth, measured_perwidth[-1]), where='post', label='Meausred Signal (Input)')

    # --- get chi2 values for each model to compare
    chi2_vals = []
    p_values = []
    model_handles = []
    model_labels = []
    for midx, mkey in enumerate(models.keys()):
        model_smeared = unfold['AddSmear'] @ models[mkey]

        chi2_val, p_val = get_chi2(Unfolded, model_smeared, UnfoldCov_syst)
        chi2_vals.append(chi2_val)
        p_values.append(p_val)

        model_smeared_perwidth = model_smeared / bin_widths
        model_handle, = plt.step(bins, np.append(model_smeared_perwidth, model_smeared_perwidth[-1]), where='post')
        model_handles.append(model_handle)
        model_labels.append(f'$A_c \\otimes$ {mkey} ($\chi^2$ = {chi2_vals[midx]:.2f}/{len(bins)-1}), p-value = {p_values[midx]:.3f}')

    # legend
    if closure_test:
        handles = [bar_handle, reco_handle] + model_handles
        labels = ['Unfolded Asimov Data', 'Measured Signal'] + model_labels
    else:
        handles = [bar_handle, norm_handle, reco_handle] + model_handles
        labels = ['Unfolded', 'Norm. Syst. Unc.', 'Measured Signal'] + model_labels
    plt.legend(handles, labels, 
               loc='upper left', fontsize=12, frameon=False, ncol=1, bbox_to_anchor=(0.02, 0.98))

    textloc_x, textloc_ha = get_textloc_x(Unfolded_perwidth, var_config.bins, textloc)
    textloc_y = textloc[1]
    add_approval_text(approval, textloc_x, textloc_y, textloc_ha)

    plt.xlabel(var_config.var_labels[0])
    plt.ylabel(var_config.xsec_label)
    plt.xlim(bins[0], bins[-1])
    plt.ylim(0., np.max(Unfolded_perwidth)*1.7)

    if save_fig:
        plt.savefig(save_name, bbox_inches='tight')
    plt.show()


# Plotters for detector variation analysis
def variation_hists(evtdfs, 
                    var_name, 
                    bins,
                    var_colors, 
                    var_labels,
                    plot_labels=["", "", ""],
                    vline = None,
                    textloc=[0.05, 0.55],
                    approval="internal",
                    save_fig=False, save_name=None): 

    vardfs, wgtdfs = [], []
    for df in evtdfs:
        vardf, wgtf = get_clipped_evts(df, var_name, bins)
        vardfs.append(vardf)
        wgtdfs.append(wgtf)

    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    fig, axs = plt.subplots(2, 1, figsize=(7.5, 7), 
                            sharex=True, gridspec_kw={'height_ratios': [4, 1]})
    fig.subplots_adjust(hspace=0.05)
    ax = axs[0]
    ax_r = axs[1]

    total_mc_list = []
    mc_stat_err_list = []
    for sidx in range(len(vardfs)):
        total_mc, _ = np.histogram(vardfs[sidx], bins=bins, weights=wgtdfs[sidx])
        total_mc_err2, _ = np.histogram(vardfs[sidx], bins=bins, weights=wgtdfs[sidx]**2)
        mc_stat_err = np.sqrt(total_mc_err2)
        total_mc_list.append(total_mc)
        mc_stat_err_list.append(mc_stat_err)

        ax.hist(vardfs[sidx], 
                weights=wgtdfs[sidx],
                bins=bins, 
                histtype="step" , 
                color=var_colors[sidx], 
                label=var_labels[sidx])

    ax.set_ylabel(plot_labels[1])
    ax.set_title(plot_labels[2])
    ax.set_xlim(bins[0], bins[-1])
    ax.legend()

    # == option to plot vertical lines ==
    if vline is not None:
        for v in vline:
            ymax = ax.get_ylim()[1]
            ax.vlines(x=v, ymin=0, ymax=ymax*0.75, color='red', linestyle='--')

    textloc_x, textloc_ha = get_textloc_x(total_mc, bins, textloc)
    textloc_y = textloc[1]
    add_approval_text(approval, textloc_x, textloc_y, textloc_ha)


    # ==== var/CV ratio panel

    for i, (vardf, wgtf) in enumerate(zip(vardfs, wgtdfs)):
        if i == 0:
            continue
        # Avoid division by zero: ignore bins where denominator is 0
        ratio = np.full_like(total_mc_list[0], np.nan, dtype=float)
        nonzero_mask = total_mc_list[0] != 0
        ratio[nonzero_mask] = total_mc_list[i][nonzero_mask] / total_mc_list[0][nonzero_mask]
        # this_err = np.sqrt(
        #     (mc_stat_err_list[0] / total_mc_list[0])**2 + 
        #     (mc_stat_err_list[i] / total_mc_list[i])**2
        # )
        # Only plot nonzero, non-nan elements in the ratio
        valid_mask = (~np.isnan(ratio)) & (ratio != 0)
        ax_r.hist(bin_centers[valid_mask], bins=bins, weights=ratio[valid_mask], linewidth=1, histtype="step")

    ax_r.axhline(1.0, color='red', linestyle='--', linewidth=1)
    ax_r.set_xlim(bins[0], bins[-1])
    ax_r.set_ylim(0.9, 1.1)
    ax_r.set_xlabel(plot_labels[0])
    ax_r.set_ylabel("Variation / CV")
    ax_r.grid(True)
    ax_r.minorticks_on()
    ax_r.grid(which='minor', linestyle=':', linewidth=0.5, color='gray', alpha=0.5)


    # == save figure ==
    if save_fig:
        plt.savefig(save_name+".pdf", bbox_inches='tight')
    plt.show()

    return total_mc_list

def signal_hists(evtdf=None,  # df with selected & reco'ed events
                 nudf=None,   # df with all MC truth
                 var_config=None,
                 return_data=False,
                 plot=True,
                 save_fig=False, 
                 save_name=None):
    """
    plot generate / selected / reco'ed signal events

    nuint_categ == 1 : signal
    topology_list has "1" as the first item, corresponding to the signal topology
    """

    bins = var_config.bins
    reco_col = var_config.var_evt_reco_col
    truth_col = var_config.var_evt_truth_col
    nu_col = var_config.var_nu_col

    # ===== all selected events =====
    # reco'ed
    var_allsel_reco, wgt_allsel_reco = get_clipped_evts(evtdf, reco_col, bins)
    var_allsel_truth, wgt_allsel_truth = get_clipped_evts(evtdf, truth_col, bins)

    # ===== true signal events =====
    evtdf_signal = evtdf[evtdf.nuint_categ == 1]
    # selected, reco'ed 
    var_sel_reco, wgt_sel_reco = get_clipped_evts(evtdf_signal, reco_col, bins)
    # selected, truth (for response matrix)
    var_sel_truth, wgt_sel_truth = get_clipped_evts(evtdf_signal, truth_col, bins)
    # all MC, truth (for efficiency vector)
    nudf_signal = nudf[nudf.nuint_categ == 1]
    var_allmc, wgt_allmc = get_clipped_evts(nudf_signal, nu_col, bins)

    nevts_allmc, _, _ = plt.hist(var_allmc, weights=wgt_allmc, bins=bins, histtype="step", label="All Signal in MC")
    nevts_allsel_truth, _, _ = plt.hist(var_allsel_truth, weights=wgt_allsel_truth, bins=bins, histtype="step", label="All Selected Events, True", color="orange")
    nevts_allsel_reco, _, _ = plt.hist(var_allsel_reco,  weights=wgt_allsel_reco,  bins=bins, histtype="step", label="All Selected Events, Reco", color="orange", linestyle='--')
    nevts_sel_truth, _, _ = plt.hist(var_sel_truth, weights=wgt_sel_truth, bins=bins, histtype="step", label="Selected Signal, True", color="blue")
    nevts_sel_reco, _, _ = plt.hist(var_sel_reco,  weights=wgt_sel_reco,  bins=bins, histtype="step", label="Selected Signal, Reco", color="blue", linestyle='--')

    plt.legend()
    plt.xlabel(var_config.var_labels[0])
    plt.ylabel("Events / Bin")
    plt.xlim(bins[0], bins[-1])

    if save_fig:
        plt.savefig(save_name, bbox_inches='tight')

    if plot == False:
        plt.close()
    else:
        plt.show()

    if return_data:
        return {
            "var_allmc": var_allmc,
            "var_sel_truth": var_sel_truth,
            "var_sel_reco": var_sel_reco,
            "var_allsel_truth": var_allsel_truth,
            "var_allsel_reco": var_allsel_reco,
            "wgt_allmc": wgt_allmc,
            "wgt_sel_truth": wgt_sel_truth,
            "wgt_sel_reco": wgt_sel_reco,
            "wgt_allsel_truth": wgt_allsel_truth,
            "wgt_allsel_reco": wgt_allsel_reco,
            "nevts_allmc": nevts_allmc,
            "nevts_sel_truth": nevts_sel_truth,
            "nevts_sel_reco": nevts_sel_reco,
            "nevts_allsel_truth": nevts_allsel_truth,
            "nevts_allsel_reco": nevts_allsel_reco,
        }
