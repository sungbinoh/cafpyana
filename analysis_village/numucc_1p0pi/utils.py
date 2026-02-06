import numpy as np
import pandas as pd
from tqdm import tqdm
import string

import sys
sys.path.append('../../')
from pyanalib.split_df_helpers import *
from pyanalib.stat_helpers import *
from makedf.constants import *
from analysis_village.unfolding.wienersvd import *
from analysis_village.numucc_1p0pi.categories import *
from analysis_village.numucc_1p0pi.constants import *

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib as mpl
plt.style.use("presentation.mplstyle")
cmap = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)


#  ====== calculation functions ======
def generate_tags(end_tag=""):
    tags = []
    for first in string.ascii_lowercase:
        for second in string.ascii_lowercase:
            tag = first + second
            if tag == end_tag:
                break
            tags.append(tag)
        if tag == end_tag:
            break
    return tags


def get_clipped_evts(df, var_col, bins, verbose=False):
    var = df[var_col]
    var = np.clip(var, bins[0], bins[-1] - EPSILON)

    if 'pot_weight' in df.columns:
        weights = df.loc[:, 'pot_weight']
    else:
        if verbose:
            print("No pot_weight column found, return 1 as pot scale (expected for data)")
        weights = np.ones_like(var)
    return var, weights


def get_univ_rates(cov_type="rate", 
                    evtdf=None, 
                    nudf=None, 
                    var_config=None, 
                    syst_name="", 
                    n_univ=100, 
                    bkgd_subtract=True,
                    plot=False):
    """
    for the GENIE uncertainty on the xsec measurement
    """

    if cov_type == "xsec":
        print("getting {} universes for {} uncertainty on the xsec".format(n_univ, syst_name))
        scale_factor = XSEC_UNIT
    elif cov_type == "rate":
        print("getting {} universes for {} uncertainty on the event rate".format(n_univ, syst_name))
        scale_factor = 1.0
    else:
        raise ValueError("Invalid covariance type: {}, choose in [xsec, rate]".format(cov_type))

    bins = var_config.bins

    evtdf_signal = evtdf[evtdf.nuint_categ == 1]
    # reco variable histogram, topology breakdown
    evtdf_div_topo = [evtdf[evtdf.nuint_categ == mode]for mode in topology_list]

    if nudf is not None:
        nudf_signal = nudf[nudf.nuint_categ == 1]

    ret = signal_hists(evtdf, nudf, var_config, return_data=True, plot=plot)

    univ_events = []
    univ_effs   = []
    univ_smears = []

    for uidx in range(n_univ):
        univ_col = "univ_{}".format(uidx)

        # ---- uncertainty on the signal rate ----
        # only consider effect on the response matrix for the signal channel
        if cov_type == "xsec":
            # smearing matrix
            reco_vs_true = get_smear_matrix(ret["var_sel_truth"], 
                                            ret["var_sel_reco"], 
                                            weights=ret["wgt_sel_truth"]*evtdf_signal[syst_name][univ_col],
                                            bins_2d=[bins, bins],
                                            plot=plot)
            univ_smears.append(reco_vs_true)

            # efficiency
            signal_allmc_univ, _ = np.histogram(ret["var_allmc"],
                                               weights=ret["wgt_allmc"]*nudf_signal[syst_name][univ_col],
                                               bins=bins)
            signal_sel_univ, _ = np.histogram(ret["var_sel_truth"],
                                               weights=ret["wgt_sel_truth"]*evtdf_signal[syst_name][univ_col],
                                               bins=bins)
            eff = signal_sel_univ / signal_allmc_univ
            univ_effs.append(eff)

            response_univ = get_response_matrix(reco_vs_true, eff, bins, plot=plot)
            signal_univ = response_univ @ ret["nevts_allmc"] # note that we multiply the CV signal rate!
            # signal_univ = signal_cv

        elif cov_type == "rate":
            signal_univ, _ = np.histogram(ret["var_sel_reco"], 
                                          weights=ret["wgt_sel_reco"]*evtdf_signal[syst_name][univ_col],
                                          bins=bins)

        else:
            raise ValueError("Invalid covariance type: {}, choose xsec or rate".format(cov_type))

        # TODO: this isn't computationally efficient, but it's useful for debugging
        # ---- uncertainty on the background rate ----
        # loop over background categories
        # + univ background - cv background
        # note: cv background subtraction cancels out with the cv background subtraction for the cv event rate. 
        #       doing it anyways for the plot of universes on background subtracted event rate.
        for this_evtdf in evtdf_div_topo[1:]:
            var, wgt = get_clipped_evts(this_evtdf, var_config.var_evt_reco_col, bins)
            univ_wgt = this_evtdf[syst_name][univ_col].copy()
            univ_wgt[np.isnan(univ_wgt)] = 1 ## IMPORTANT: make nan univ_wgt to 1. to ignore them
            background_cv, _   = np.histogram(var, bins=bins, weights=wgt)
            background_univ, _ = np.histogram(var, bins=bins, weights=wgt*univ_wgt)

            if bkgd_subtract:
                signal_univ += background_univ - background_cv
            else:
                signal_univ += background_univ

        signal_univ *= scale_factor
        univ_events.append(signal_univ)

    univ_events = np.array(univ_events)

    if bkgd_subtract:
        cv_events = ret["nevts_sel_reco"]
        cv_events *= scale_factor # = Response @ true_signal
    else:
        cv_events = ret["nevts_allsel_reco"]
        cv_events *= scale_factor # = Response @ true_signal

    return univ_events, cv_events


def get_response_matrix(reco_vs_true, 
                        eff):
    denom = reco_vs_true.T.sum(axis=0)
    num = reco_vs_true.T
    response = np.divide(
        num * eff, denom,
        out=np.zeros_like(num, dtype=float),  # fill with 0 where invalid
        where=denom != 0
    )
    return response



# ====== plotting functions ======

# ==== plot additions ====
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
    if approval == "internal":
        approval_text = r"$\mathbf{SBND}$ Internal"
        textcolor = 'rosybrown'

    elif approval == "preliminary":
        approval_text = r"$\mathbf{SBND}$ Preliminary"
        textcolor = 'gray'

    else:
        return # don't add anything

    ax = plt.gcf().axes[0]  # get the first axes of the current figure
    ax.text(
        textloc_x, textloc_y, 
        approval_text, 
        transform=ax.transAxes, 
        ha=textloc_ha, va='top',
        fontsize=20, color=textcolor
    )


def add_genie_version_text(textloc_x, textloc_y, textloc_ha):
    ax = plt.gcf().axes[0]  # get the first axes of the current figure
    ax.text(textloc_x, textloc_y, 
            r"GENIE v3.4.0 AR23_00i_00_000", 
            transform=ax.transAxes, 
            ha=textloc_ha, va='top',
            fontsize=11.5, color='gray')


# ==== bar plot ====
def bar_plot(breakdown_type="topology", 
             generator=None,
             evtdf=None,
             show_plot=True, 
             plot_labels=["", "", ""],
             save_fig=False, 
             save_name=None): 

    if breakdown_type == "nu_cosmics":
        ncateg = 3
        labels = nu_cosmics_labels
        colors = nu_cosmics_colors

        cut_cosmic = IsCosmic(evtdf)
        cut_nu_outfv = IsNuOutFV(evtdf)
        cut_nu_intfv = IsNu(evtdf) & InFV(evtdf.slc.truth.position, det=DETECTOR)
        cuts = [cut_cosmic, cut_nu_outfv, cut_nu_intfv]

        if breakdown_type == "topology":
            labels = topology_labels
            colors = topology_colors
            cuts = get_int_category(evtdf, ret_cuts=True)

        elif breakdown_type == "genie":
            labels = genie_mode_labels
            colors = genie_mode_colors
            cuts = get_genie_category(evtdf, ret_cuts=True)

        else:
            raise ValueError("Invalid breakdown_type: %s, please choose between [topology, genie, or genie_sb]" % breakdown_type)

    fig, ax = plt.subplots(figsize = (6, ncateg*0.6))

    scale = evtdf.pot_weight.unique()[0]
    # size = [scale*len(evtdf[i]) for i in cuts]

    if generator == "GiBUU":
        wgt_means = []
        for i in cuts:
            wgt_means.append(np.nan_to_num(evtdf[i].genweight.mean(), nan=1.))
        size = [scale*len(evtdf[i])*wgt_mean for i, wgt_mean in zip(cuts, wgt_means)]
    else:
        size = [scale*len(evtdf[i]) for i in cuts]

    bars = plt.barh(labels, size, align='center', color = colors)
    tot_count = np.array(size).sum()
    
    perc_list = []
    for bar in bars:
        width = bar.get_width()
        label_y_pos = bar.get_y() + bar.get_height() / 2
        perc = 100*(width+0.)/(tot_count+0.)
        ax.text(width+1, label_y_pos, s= ("%0.1f"%(100*(width+0.)/(tot_count+0.)) + "%"), va='center')
        perc_list.append(perc)

    plt.xlabel(plot_labels[0])
    plt.xlim(0, 1.12 * np.max(size))

    if not show_plot:
        plt.close()
        # instead print breakdowns as a table
        print("Breakdown ", breakdown_type, " ", generator)
        print("-"*20)
        for label, perc in zip(labels, perc_list):
            print(f"{label}: {perc:.1f}%")
        print("-"*20)
        print(f"Total: {np.sum(perc_list):.1f}%")

    if save_fig:
        plt.savefig(save_name, bbox_inches="tight")

    ret_dict = {"cuts": cuts,
                "perc_list": perc_list}
    return ret_dict


# ==== histograms ====
def overlay_hists(breakdown_type="topology",
                  mc_df=None,
                  data_df=None,
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

    # ==== prepare dfs for plotting ====

    # MC
    if mc_df is not None:
        vardf, _        = get_clipped_evts(mc_df, var_config.var_evt_reco_col, var_config.bins)

        # breakdown MC events into truth categories
        if breakdown_type == "topology":
            labels = topology_labels
            colors = topology_colors
            cuts = get_int_category(mc_df, ret_cuts=True)

        elif breakdown_type == "genie":
            labels = genie_mode_labels
            colors = genie_mode_colors
            cuts = get_genie_category(mc_df, ret_cuts=True)

        elif breakdown_type == "genie_sb":
            labels = genie_sb_mode_labels
            colors = genie_sb_mode_colors
            cuts = get_genie_sb_category(mc_df, ret_cuts=True)
            # hatches for marking S/B on plot
            hatches = [None] * len(labels)
            for i in range(5, len(labels), 2):
                hatches[i] = '////'
        else:
            raise ValueError("Invalid breakdown_type: %s, please choose between [topology, genie, or genie_sb]" % breakdown_type)
        var_categ = [vardf[i] for i in cuts]
        weights_categ = [list(mc_df.loc[cuts[i], 'pot_weight']) for i in range(len(cuts))] 

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

    # adjust fontsize so that legend fits in the figure
    fontsize = 11.3
    ncol = 3
    if breakdown_type == "genie_sb":
        fontsize = 10.5
        ncol = 4

    ax.legend(
        ordered_handles,
        ordered_labels,
        loc='upper left',
        fontsize=fontsize,
        frameon=False,
        ncol=ncol,
        bbox_to_anchor=(0.01, 0.99)
    )

    # ===============================

    # ==== plot additions ====

    # y-axis limit
    ax.set_ylim(0., ax_ylim_ratio* np.max(total_mc))

    # vertical lines
    if vline is not None:
        for v in vline:
            ymax = ax.get_ylim()[1]
            ax.vlines(x=v, ymin=0, ymax=ymax*0.75, color='red', linestyle='--')

    # textboxes
    textloc_x, textloc_ha = get_textloc_x(total_mc, var_config.bins, textloc)
    textloc_y = textloc[1]
    add_approval_text(approval, textloc_x, textloc_y, textloc_ha)
    add_genie_version_text(textloc_x, textloc_y-0.08, textloc_ha)

    # ===============================

    # == save figure ==
    if save_fig:
        plt.savefig(save_name+".pdf", bbox_inches='tight')
    plt.show()

    return {"cuts": cuts, 
            "total_mc": total_mc, 
            "total_mc_bkgd": total_mc_bkgd,
            "total_data": total_data}


def plot_univ_hists(
                univ_events, 
                cv_events,
                syst_name, 
                var_config, 
                approval="internal",
                textloc=[0.05, 0.55],
                save_fig=False, 
                save_name=None): 

    assert univ_events.shape[1] == len(cv_events) 
    n_univ = univ_events.shape[0]

    if (n_univ > 10):
        colors = ["C0", "C1", "C2"]

        sorted_univs = np.sort(univ_events, axis=0)

        n_68 = int(0.68 * n_univ)
        start_68 = (n_univ - n_68) // 2
        end_68 = start_68 + n_68

        n_95 = int(0.95 * n_univ)
        start_95 = (n_univ - n_95) // 2
        end_95 = start_95 + n_95

        # Define bins & groupings of universes: [(range, color, label, skip_68)]
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
            show_label = "Universe" if i == 0 else None
            plt.hist(var_config.bin_centers, bins=var_config.bins, weights=univ_events[i], histtype="step", color="gray", label=show_label)

    # plot CV last so that it is on top of the universes
    plt.hist(var_config.bin_centers, bins=var_config.bins, weights=cv_events, histtype="step", color="k", label="Central Value")

    plt.xlim(var_config.bins[0], var_config.bins[-1])
    plt.xlabel(var_config.var_labels[1])
    plt.ylabel("Events / Bin")
    plt.title(syst_name)

    plt.legend(frameon=False)

    # ==== plot additions ====
    textloc_x, textloc_ha = get_textloc_x(cv_events, var_config.bins, textloc)
    textloc_y = textloc[1]
    add_approval_text(approval, textloc_x, textloc_y, textloc_ha)


    if save_fig:
        plt.savefig(save_name+".pdf", bbox_inches='tight')

    plt.show();


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

    plt.xlabel(var_config.var_labels[0])
    plt.ylabel(var_config.xsec_label)
    plt.xlim(bins[0], bins[-1])
    plt.ylim(0., np.max(Unfolded_perwidth)*1.7)

    # ==== plot additions
    textloc_x, textloc_ha = get_textloc_x(Unfolded_perwidth, var_config.bins, textloc)
    textloc_y = textloc[1]
    add_approval_text(approval, textloc_x, textloc_y, textloc_ha)


    if save_fig:
        plt.savefig(save_name, bbox_inches='tight')
    plt.show()


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

    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    vardfs, wgtdfs = [], []
    for df in evtdfs:
        vardf, wgtf = get_clipped_evts(df, var_name, bins)
        vardfs.append(vardf)
        wgtdfs.append(wgtf)

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

    # ==== plot additions
    # vertical lines on main panel
    if vline is not None:
        for v in vline:
            ymax = ax.get_ylim()[1]
            ax.vlines(x=v, ymin=0, ymax=ymax*0.75, color='red', linestyle='--')

    # approval rank
    textloc_x, textloc_ha = get_textloc_x(total_mc, bins, textloc)
    textloc_y = textloc[1]
    add_approval_text(approval, textloc_x, textloc_y, textloc_ha)

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
                 textloc=[0.05, 0.55],
                 approval="internal",
                 save_fig=False, 
                 save_name=None):
    """
    generate / selected / reco'ed signal events

    nuint_categ == 1 : signal
    topology_list has "1" as the first item, corresponding to the signal topology

    items with "sel" tags are selected events that are signal in truth 
    item with "allsel" tags are all selected events, signal + background
    """

    bins = var_config.bins
    bin_centers = var_config.bin_centers
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

    nevts_allsel_truth, _ = np.histogram(var_allsel_truth, weights=wgt_allsel_truth, bins=bins)
    nevts_allsel_reco, _  = np.histogram(var_allsel_reco,  weights=wgt_allsel_reco,  bins=bins)
    nevts_sel_truth, _    = np.histogram(var_sel_truth,    weights=wgt_sel_truth,    bins=bins)
    nevts_sel_reco, _     = np.histogram(var_sel_reco,     weights=wgt_sel_reco,     bins=bins)

    if nudf is not None:
        # all MC, truth (for efficiency vector)
        nudf_signal = nudf[nudf.nuint_categ == 1]
        var_allmc, wgt_allmc = get_clipped_evts(nudf_signal, nu_col, bins)
        nevts_allmc, _ = np.histogram(var_allmc, weights=wgt_allmc, bins=bins)
    else:
        var_allmc = None
        wgt_allmc = None
        nevts_allmc = None

    if plot:
        plt.hist(bin_centers, weights=nevts_allsel_truth, bins=bins, histtype="step", label="Selected Events, True", color="C0")
        plt.hist(bin_centers, weights=nevts_allsel_reco,  bins=bins, histtype="step", label="Selected Events, Reco", color="C0", linestyle="--")
        plt.hist(bin_centers, weights=nevts_sel_truth,    bins=bins, histtype="step", label="Selected Signal, True", color="C1")
        plt.hist(bin_centers, weights=nevts_sel_reco,     bins=bins, histtype="step", label="Selected Signal, Reco", color="C1", linestyle="--")

        if nevts_allmc is not None:
            nevts_allmc, _, _ = plt.hist(var_allmc,     weights=wgt_allmc,     bins=bins, histtype="step", label="All Signal in MC", color="black")
        else:
            nevts_allmc = None

        plt.xlabel(var_config.var_labels[0])
        plt.ylabel("Events / Bin")
        plt.xlim(bins[0], bins[-1])
        plt.legend()

        # ==== plot additions ====
        textloc_x, textloc_ha = get_textloc_x(nevts_allsel_truth, bins, textloc)
        textloc_y = textloc[1]
        add_approval_text(approval, textloc_x, textloc_y, textloc_ha)

        if save_fig:
            plt.savefig(save_name, bbox_inches='tight')
        plt.show()

    if return_data:
        return {
            "var_allmc": var_allmc,
            "wgt_allmc": wgt_allmc,
            "nevts_allmc": nevts_allmc,

            "var_sel_truth": var_sel_truth,
            "wgt_sel_truth": wgt_sel_truth,
            "nevts_sel_truth": nevts_sel_truth,

            "var_sel_reco": var_sel_reco,
            "wgt_sel_reco": wgt_sel_reco,
            "nevts_sel_reco": nevts_sel_reco,

            "var_allsel_truth": var_allsel_truth,
            "wgt_allsel_truth": wgt_allsel_truth,
            "nevts_allsel_truth": nevts_allsel_truth,

            "var_allsel_reco": var_allsel_reco,
            "wgt_allsel_reco": wgt_allsel_reco,
            "nevts_allsel_reco": nevts_allsel_reco,
        }


# ==== fractional uncertainty plot ====
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
    plt.show();


# ==== 2D plots ====

def get_text_color(value):
    rgba = cmap(norm(value))
    # Compute luminance (perceived brightness)
    luminance = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
    return "black" if luminance > 0.5 else "white"


def bin_range_labels(edges):
    return [f"{edges[i]:.2f}–{edges[i+1]:.2f}" for i in range(len(edges)-1)]


def plot_heatmap(matrix, 
                 bins,
                 plot_labels=["", "", ""],
                 approval="internal",
                 save_fig=False, 
                 save_name=None):

    nbins = len(bins)
    assert nbins-1 == matrix.shape[0] == matrix.shape[1]
    unif_bin = np.linspace(0., float(nbins - 1), nbins)
    extent = [unif_bin[0], unif_bin[-1], unif_bin[0], unif_bin[-1]]

    x_edges, y_edges = np.array(bins), np.array(bins)
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
                    f"{value:.2f}",
                    ha="center", va="center",   
                    color=get_text_color(value),
                    fontsize=10
                )

    plt.colorbar(shrink=0.7, label=plot_labels[2])
    plt.xticks(x_tick_positions, x_labels, rotation=45, ha="right")
    plt.yticks(y_tick_positions, y_labels)
    plt.xlabel(plot_labels[0], fontsize=20)
    plt.ylabel(plot_labels[1], fontsize=20)
    # plt.title(plot_labels[2], fontsize=20)

    # ===== plot additions =====
    add_approval_text(approval, 0.95, 1.05, "right")

    if save_fig:
        plt.savefig("{}.png".format(save_name), bbox_inches='tight', dpi=300)
    plt.show();

