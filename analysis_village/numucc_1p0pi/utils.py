import numpy as np
from makedf.constants import *
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnchoredText
from matplotlib.offsetbox import AnchoredOffsetbox, DrawingArea, HPacker, VPacker, TextArea
from matplotlib.legend import Legend

from scipy.stats import chi2

import sys
sys.path.append('../../')
from analysis_village.numucc_1p0pi.selection_definitions import *
from analysis_village.numucc_1p0pi.variable_configs import VariableConfig
from analysis_village.numucc_1p0pi.utils import *
from pyanalib.split_df_helpers import *
from pyanalib.stat_helpers import *
from pyanalib.variable_calculator import get_cc1p0pi_tki
from pyanalib.pandas_helpers import pad_column_name
from analysis_village.unfolding.wienersvd import *
from analysis_village.unfolding.unfolding_inputs import *
from analysis_village.unfolding.wienersvd import *

DETECTOR = "SBND_nohighyz"
# DETECTOR = "SBND"

eps = 1e-6


# import pickle
# with open("/exp/sbnd/app/users/munjung/xsec/cafpyana/analysis_village/numucc1p0pi/NuINT_uncs/tot_uncert_dict.pkl", "rb") as f:
#     syst_uncert_dict = pickle.load(f)
# syst_uncert_dict.keys()


# KE <-> p conversion for providing the sig threshold in terms of KE in technote
def p_to_KE(p, mass):
    return np.sqrt(p**2 + mass**2) - mass

def KE_to_p(KE, mass):
    E = KE + mass
    return np.sqrt(E**2 - mass**2)


### Plotters

# def hist_plot(plot_type,
#               evtdf, vardf, 
#               vardf_data, var_intime,
#               bins,
#               generator="GENIE",
#               plot_labels=["", "", ""],
#               ratio = True,
#               vline = None,
#               save_fig=False, save_name=None): 

#     assert len(evtdf) == len(vardf)

#     if plot_type == "nu_cosmics":
#         labels = nu_cosmics_labels
#         colors = nu_cosmics_colors

#         cut_cosmic = IsCosmic(evtdf)
#         cut_nu_outfv = IsNuOutFV(evtdf)
#         cut_nu_infv = IsNuInFV(evtdf)
#         cuts = [cut_cosmic, cut_nu_outfv, cut_nu_infv]

#     elif plot_type == "topology":
#         labels = topology_labels[::-1]
#         colors = topology_colors[::-1]

#         cut_cosmic = IsCosmic(evtdf)
#         cut_nu_outfv = IsNuOutFV(evtdf)
#         cut_nu_infv_nu_other = IsNuInFV_NuOther(evtdf)
#         cut_nu_infv_numu_nc = IsNuInFV_NumuNC(evtdf)
#         cut_nu_infv_numu_cc_other = IsNuInFV_NumuCC_Other(evtdf)
#         cut_nu_infv_numu_cc_np0pi = IsNuInFV_NumuCC_Np0pi(evtdf)
#         cut_nu_infv_numu_cc_1p0pi = IsNuInFV_NumuCC_1p0pi(evtdf)
#         cuts = [cut_cosmic, cut_nu_outfv, cut_nu_infv_nu_other, cut_nu_infv_numu_nc, 
#                 cut_nu_infv_numu_cc_other, cut_nu_infv_numu_cc_np0pi, cut_nu_infv_numu_cc_1p0pi]

#     elif plot_type == "genie" and generator == "GENIE":
#         labels = genie_mode_labels[::-1]
#         colors = genie_mode_colors[::-1]

#         cut_cosmic = IsCosmic(evtdf)
#         cut_nu_outfv = IsNuOutFV(evtdf)
#         cut_nu_infv_nu_other = IsNuInFV_NuOther(evtdf)
#         cut_nu_infv_numu_nc = IsNuInFV_NumuNC(evtdf)
#         # print("numu NC", cut_nu_infv_numu_nc.sum())
#         # cut_nu_infv_numu_coh = IsNuInFV_NumuCC_COH(evtdf)
#         cut_nu_infv_numu_othermode = IsNuInFV_NumuCC_OtherMode(evtdf)
#         cut_nu_infv_numu_cc_dis = IsNuInFV_NumuCC_DIS(evtdf)
#         cut_nu_infv_numu_cc_res = IsNuInFV_NumuCC_RES(evtdf)
#         cut_nu_infv_numu_cc_me = IsNuInFV_NumuCC_MEC(evtdf)
#         cut_nu_infv_numu_cc_qe = IsNuInFV_NumuCC_QE(evtdf)
#         cuts = [cut_cosmic, cut_nu_outfv, cut_nu_infv_nu_other, cut_nu_infv_numu_nc, 
#                 cut_nu_infv_numu_othermode, cut_nu_infv_numu_cc_dis, cut_nu_infv_numu_cc_res, 
#                 cut_nu_infv_numu_cc_me, cut_nu_infv_numu_cc_qe]

#     elif plot_type == "genie" and generator == "GiBUU":
#         labels = gibuu_mode_labels[::-1]
#         colors = gibuu_mode_colors[::-1]

#         cut_cosmic = IsCosmic(evtdf)
#         cut_nu_outfv = IsNuOutFV(evtdf)
#         cut_nu_infv_nu_other = IsNuInFV_NuOther(evtdf)
#         cut_nu_infv_numu_nc = IsNuInFV_NumuNC(evtdf)
#         print("numu NC", cut_nu_infv_numu_nc.sum())
#         # cut_nu_infv_numu_coh = IsNuInFV_NumuCC_COH(evtdf)
#         cut_nu_infv_numu_othermode = IsNuInFV_NumuCC_OtherMode_GiBUU(evtdf)
#         cut_nu_infv_numu_cc_dis = IsNuInFV_NumuCC_DIS_GiBUU(evtdf)
#         cut_nu_infv_numu_cc_res = IsNuInFV_NumuCC_RES_GiBUU(evtdf)
#         cut_nu_infv_numu_cc_me = IsNuInFV_NumuCC_MEC_GiBUU(evtdf)
#         cut_nu_infv_numu_cc_qe = IsNuInFV_NumuCC_QE_GiBUU(evtdf)
#         cuts = [cut_cosmic, cut_nu_outfv, cut_nu_infv_nu_other, cut_nu_infv_numu_nc, 
#                 cut_nu_infv_numu_othermode, cut_nu_infv_numu_cc_dis, cut_nu_infv_numu_cc_res, 
#                 cut_nu_infv_numu_cc_me, cut_nu_infv_numu_cc_qe]
        

#     else:
#         raise ValueError("Invalid plot_type: %s, please choose between [nu_cosmics, topolgy, or genie]" % plot_type)

#     # --- Plot template
#     if ratio:
#         fig, axs = plt.subplots(2, 1, figsize=(8.5, 8), 
#                                sharex=True, gridspec_kw={'height_ratios': [3, 1]})
#         fig.subplots_adjust(hspace=0.05)
#         ax = axs[0]
#         ax_r = axs[1]
#     else:
#         fig, ax = plt.subplots()

#     bin_centers = 0.5 * (bins[:-1] + bins[1:])

#     # --- Data
#     total_data, bins = np.histogram(vardf_data, bins=bins)
#     data_eylow, data_eyhigh = return_data_stat_err(total_data)
#     ax.errorbar(bin_centers, total_data, yerr=np.vstack((data_eylow, data_eyhigh)), 
#                 fmt='o', color='black')

#     # --- MC
#     # collect all MC + intime
#     var_categ = [var_intime] + [vardf[i] for i in cuts]
#     scale_mc = evtdf.pot_weight.unique()[0]
#     # TODO:
#     scale_intime_to_lightdata = 0.073
#     scale_intime = scale_intime_to_lightdata
#     weights_categ = [scale_intime*np.ones_like(var_intime)] + [scale_mc*np.ones_like(vardf[i]) for i in cuts] 
#     if generator == "GiBUU":
#         # TODO: GiBUU
#         weights_categ = [scale_intime*np.ones_like(var_intime)] + [np.nan_to_num(scale_mc*evtdf[i].genweight, nan=1.) for i in cuts] 
#         weights_categ = [scale_intime*np.ones_like(var_intime)] + [scale_mc*np.ones_like(vardf[i]) for i in [cuts[0]]] + [np.nan_to_num(scale_mc*evtdf[i].genweight, nan=1.) for i in cuts[1:]] 
#     colors = ["black"] + colors
#     labels = ["Cosmic\n(In-time)"] + labels

#     mc_stack, _, _ = ax.hist(var_categ,
#                                 bins=bins,
#                                 weights=weights_categ,
#                                 stacked=True,
#                                 color=colors,
#                                 label=labels,
#                                 edgecolor='none',
#                                 linewidth=0,
#                                 density=False,
#                                 histtype='stepfilled')


#     # ---- MC stat err
#     each_mc_hist_data = []
#     each_mc_hist_err2 = []  # sum of squared weights for error

#     for data, w in zip(var_categ, weights_categ):
#         hist_vals, _ = np.histogram(data, bins=bins, weights=w)
#         hist_err2, _ = np.histogram(data, bins=bins, weights=np.square(w))
#         each_mc_hist_data.append(hist_vals)
#         each_mc_hist_err2.append(hist_err2)

#     total_mc = np.sum(each_mc_hist_data, axis=0)
#     total_mc_err2 = np.sum(each_mc_hist_err2, axis=0)
#     mc_stat_err = np.sqrt(total_mc_err2)

#     ax.bar(
#        bin_centers,
#         2 * mc_stat_err,
#         width=np.diff(bins),
#         bottom=total_mc - mc_stat_err,
#         facecolor='none',             # transparent fill
#         edgecolor='dimgray',            # outline color of the hatching
#         hatch='xxxx',                 # hatch pattern similar to ROOT's 3004
#         linewidth=0.0,
#         label='MC Stat. Unc.'
#     )

#     # ax.errorbar(bin_centers, total_mc, yerr=np.sqrt(total_mc), fmt='o', color='red')

#     ax.set_xlim(bins[0], bins[-1])
#     if ratio == False: # only plot xlabel if we're not plotting the ratio panel
#         ax_r.set_xlabel(plot_labels[0])
#     ax.set_ylabel(plot_labels[1])

#     if vline is not None:
#         ax.axvline(x=vline, color='red', linestyle='--')


#     if ratio:
#         # MC stat err
#         mc_stat_err_ratio = mc_stat_err / total_mc
#         mc_content_ratio = total_mc / total_mc
#         mc_stat_err_ratio = np.nan_to_num(mc_stat_err_ratio, nan=0.)
#         mc_content_ratio = np.nan_to_num(mc_content_ratio, nan=-999.)
#         ax_r.bar(
#             bin_centers,
#             2*mc_stat_err_ratio,
#             width=np.diff(bins),
#             bottom=mc_content_ratio - mc_stat_err_ratio,
#             facecolor='none',             # transparent fill
#             edgecolor='dimgray',          # outline color of the hatching
#             hatch='xxxx',                 # hatch pattern similar to ROOT's 3004
#             linewidth=0.0,
#             label='MC Stat. Unc.'
#         )

#         # data/MC ratio err
#         data_ratio = total_data / total_mc
#         data_ratio_eylow = data_eylow / total_mc
#         data_ratio_eyhigh = data_eyhigh / total_mc
#         data_ratio = np.nan_to_num(data_ratio, nan=-999.)
#         data_ratio_eylow = np.nan_to_num(data_ratio_eylow, nan=0.)
#         data_ratio_eyhigh = np.nan_to_num(data_ratio_eyhigh, nan=0.)
        
#         #data_ratio_errors = data_ratio_eylow + data_ratio_eyhigh
#         #ax_ratio.errorbar(bin_centers, data_ratio, yerr=data_ratio_errors,
#         #                 fmt='o', color='black', label='Data',
#         #                 markersize=5, capsize=3, linewidth=1.5)

#         ax_r.errorbar(bin_centers, data_ratio,
#                   yerr=np.vstack((data_ratio_eylow, data_ratio_eyhigh)),
#                   fmt='o', color='black')
#                 #   , label='Data')
#                 #   markersize=5, capsize=3, linewidth=1.5)
        
#         # if highest value is greater than 2.0, set ylim to 2.0
#         # if np.max(data_ratio) > 2.0:
#         ax_r.set_ylim(0.5, 1.5)

#         ax_r.axhline(1.0, color='red', linestyle='--', linewidth=1)
        
#         ax_r.grid(True)
#         ax_r.minorticks_on()
#         ax_r.grid(which='minor', linestyle=':', linewidth=0.5, color='gray', alpha=0.5)

#         ax_r.set_xlabel(plot_labels[0])
#         ax_r.set_ylabel("Data/MC")

#     # --- Legend
#     accum_sum = [np.sum(data) for data in mc_stack]
#     accum_sum = [0.] + accum_sum
#     total_sum = accum_sum[-1]
#     individual_sums = [accum_sum[i + 1] - accum_sum[i] for i in range(len(accum_sum) - 1)]
#     fractions = [(count / total_sum) * 100 for count in individual_sums]
#     legend_labels = [f"{label} ({frac:.1f}%)" for label, frac in zip(labels[::-1], fractions[::-1])]
#     legend_labels += ["Data", "MC Stat. Unc."]
#     leg = ax.legend(legend_labels, 
#                     loc='upper left', 
#                     fontsize=10, 
#                     frameon=False, 
#                     ncol=3, 
#                     bbox_to_anchor=(0.02, 0.98))
#     leg_height = leg.get_bbox_to_anchor().height
#     max_data_with_err = np.max(total_data + data_eyhigh)
#     # ax.set_ylim(0., 1.05 * max_data_with_err + leg_height)
#     ax.set_ylim(0., 1.4 * max_data_with_err)

#     # textbox
#     # textbox = ax.text(0.95, 0.7, "GiBUU", transform=ax.transAxes, fontsize=20,
#     #                   verticalalignment='top', horizontalalignment='right',
#     #                   bbox=dict(facecolor='white', alpha=0.5))

#     if save_fig:
#         plt.savefig(save_name, bbox_inches='tight') #, dpi=300)
#     plt.show()

#     # bolder figure lines?
#     # ax.tick_params(width=2, length=10)
#     # for spine in ax.spines.values():
#     #     spine.set_linewidth(2)
    
#     ret_dict = {"cuts": cuts}
#     return ret_dict



def bar_plot(breakdown_type, generator,
             evtdf,
             show_plot=True, plot_labels=["", "", ""],
             save_fig=False, save_name=None): #, scale, stage):

    if breakdown_type == "nu_cosmics":
        ncateg = 3
        labels = nu_cosmics_labels
        colors = nu_cosmics_colors

        cut_cosmic = IsCosmic(evtdf)
        cut_nu_outfv = IsNuOutFV(evtdf)
        cut_nu_intfv = IsNu(evtdf) & InFV(evtdf.slc.truth.position, det=DETECTOR)
        cuts = [cut_cosmic, cut_nu_outfv, cut_nu_intfv]

    elif breakdown_type == "topology":
        ncateg = 7
        labels = topology_labels[::-1]
        colors = topology_colors[::-1]

        cut_cosmic = IsCosmic(evtdf)
        cut_nu_outfv = IsNuOutFV(evtdf)
        cut_nu_infv_nu_other = IsNuInFV_NuOther(evtdf)
        cut_nu_infv_numu_nc = IsNuInFV_NumuNC(evtdf)
        cut_nu_infv_numu_cc_other = IsNuInFV_NumuCC_Other(evtdf)
        cut_nu_infv_numu_cc_np0pi = IsNuInFV_NumuCC_Np0pi(evtdf)
        cut_nu_infv_numu_cc_1p0pi = IsNuInFV_NumuCC_1p0pi(evtdf)

        cuts = [cut_cosmic, cut_nu_outfv, cut_nu_infv_nu_other, cut_nu_infv_numu_nc, 
                cut_nu_infv_numu_cc_other, cut_nu_infv_numu_cc_np0pi, cut_nu_infv_numu_cc_1p0pi]

    elif breakdown_type == "genie":
        ncateg = 9

        if generator == "GENIE":
            labels = genie_mode_labels[::-1]
            colors = genie_mode_colors[::-1]

            cut_cosmic = IsCosmic(evtdf)
            cut_nu_outfv = IsNuOutFV(evtdf)
            cut_nu_infv_nu_other = IsNuInFV_NuOther(evtdf)
            cut_nu_infv_numu_nc = IsNuInFV_NumuNC(evtdf)
            # cut_nu_infv_numu_coh = IsNuInFV_NumuCC_COH(evtdf)
            cut_nu_infv_numu_othermode = IsNuInFV_NumuCC_OtherMode(evtdf)
            cut_nu_infv_numu_cc_dis = IsNuInFV_NumuCC_DIS(evtdf)
            cut_nu_infv_numu_cc_res = IsNuInFV_NumuCC_RES(evtdf)
            cut_nu_infv_numu_cc_me = IsNuInFV_NumuCC_MEC(evtdf)
            cut_nu_infv_numu_cc_qe = IsNuInFV_NumuCC_QE(evtdf)
            cuts = [cut_cosmic, cut_nu_outfv, cut_nu_infv_nu_other, cut_nu_infv_numu_nc, 
                    cut_nu_infv_numu_othermode, cut_nu_infv_numu_cc_dis, cut_nu_infv_numu_cc_res, 
                    cut_nu_infv_numu_cc_me, cut_nu_infv_numu_cc_qe]

        elif generator == "GiBUU":
            labels = gibuu_mode_labels[::-1]
            colors = gibuu_mode_colors[::-1]

            cut_cosmic = IsCosmic(evtdf)
            cut_nu_outfv = IsNuOutFV(evtdf)
            cut_nu_infv_nu_other = IsNuInFV_NuOther(evtdf)
            cut_nu_infv_numu_nc = IsNuInFV_NumuNC(evtdf)
            print("numu NC", cut_nu_infv_numu_nc.sum())
            # cut_nu_infv_numu_coh = IsNuInFV_NumuCC_COH(evtdf)
            cut_nu_infv_numu_othermode = IsNuInFV_NumuCC_OtherMode_GiBUU(evtdf)
            cut_nu_infv_numu_cc_dis = IsNuInFV_NumuCC_DIS_GiBUU(evtdf)
            cut_nu_infv_numu_cc_res = IsNuInFV_NumuCC_RES_GiBUU(evtdf)
            cut_nu_infv_numu_cc_me = IsNuInFV_NumuCC_MEC_GiBUU(evtdf)
            cut_nu_infv_numu_cc_qe = IsNuInFV_NumuCC_QE_GiBUU(evtdf)
            cuts = [cut_cosmic, cut_nu_outfv, cut_nu_infv_nu_other, cut_nu_infv_numu_nc, 
                    cut_nu_infv_numu_othermode, cut_nu_infv_numu_cc_dis, cut_nu_infv_numu_cc_res, 
                    cut_nu_infv_numu_cc_me, cut_nu_infv_numu_cc_qe]
        else:
            raise ValueError("Invalid generator: %s, please choose between [GENIE, GiBUU]" % generator)

    else:
        raise ValueError("Invalid breakdown_type: %s, please choose between [nu_cosmics, topolgy, or genie]" % breakdown_type)


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

    # # make sure that the categories don't overlap
    # for i in range(len(cuts)):
    #     for j in range(i+1, len(cuts)):
    #         if (cuts[i] & cuts[j]).sum() > 0:
    #             print("Categories overlap:", labels[i], labels[j])
    # # and check if the categories cover all events
    # if not np.array(size).sum() == len(evtdf):
    #     print("Categories do not cover all events")
    #     print("Total events:", len(evtdf))
    #     print("Sum of categories:", np.array(size).sum())

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


#TODO: move to selection_definitions.py
nu_cosmics_labels = ["Cosmic", r"Out-FV $\nu$", r"FV $\nu$"]
# nu_cosmics_colors = ["#ED5564", "#FFCE54", "#A0D568"]
nu_cosmics_colors = ["gray", "C0", "C1"]

topology_labels = ["Cosmic", r"Out-FV $\nu$", r"Other $\nu$", r"$\nu$ NC",  
                  r"$\nu_{\mu}$ CC Other", r"$\nu_{\mu}$ CC Np0$\pi$", r"$\nu_{\mu}$ CC 1p0$\pi$"]
topology_colors = ["gray", "sienna", "crimson", "darkgreen", 
                  "coral", "darkslateblue", "mediumslateblue"]

genie_labels = ["Cosmic", r"Out-FV $\nu$", r"Other $\nu$", r"$\nu$ NC",
                r"$\nu_{\mu}$ CC Other", r"$\nu_{\mu}$ CCDIS",
                r"$\nu_{\mu}$ CCRES", r"$\nu_{\mu}$ CCMEC",
                r"$\nu_{\mu}$ CCQE"]
genie_colors = ["gray", "sienna", "crimson", "darkgreen",
                "#BFB17C", "#D88A3B", "#2c7c94",
                "#390C1E", "#9b5580"]


# TODO: combine with "plot_topology_breakdown" and "plot_genie_breakdown" functions below
def hist_plot(breakdown_type="topology",
              mc_df=None,
              data_df=None,
              intime_df=None,
              var_config="",
              plot_labels=["", "", ""],
              ax_ylim_ratio=1.5,
              ratio = False,
              syst = False,
              vline = None,
              approval="internal",
              save_fig=False, 
              save_name=None): 

    # ==== prepare dfs for plotting ====

    # MC
    if mc_df is not None:
        vardf, _        = get_clipped_evts(mc_df, var_config.var_evt_reco_col, var_config.bins)
        # TODO: get from column?
        scale_mc = mc_df.pot_weight.unique()[0]

        # breakdown MC events into truth categories
        if breakdown_type == "topology":
            labels = topology_labels
            colors = topology_colors
            cuts = get_int_category(mc_df, ret_cuts=True)

        elif breakdown_type == "genie":
            labels = genie_labels
            colors = genie_colors
            cuts = get_genie_category(mc_df, ret_cuts=True)

        else:
            raise ValueError("Invalid breakdown_type: %s, please choose between [topology, or genie]" % breakdown_type)
        var_categ = [vardf[i] for i in cuts]
        weights_categ = [scale_mc*np.ones_like(vardf[i]) for i in cuts] 

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

    else:
        vardf = None
        var_categ = None
        total_mc = None
        print("No MC data provided")

    # Data
    if data_df is not None:
        vardf_data, _   = get_clipped_evts(data_df, var_config.var_evt_reco_col, var_config.bins)
        total_data, _ = np.histogram(vardf_data, bins=var_config.bins)
        data_err = np.sqrt(total_data)
    else:
        vardf_data = None
        total_data = None
        print("No data data provided")

    # Intime cosmics
    if intime_df is not None:
        vardf_intime, _ = get_clipped_evts(intime_df, var_config.var_evt_reco_col, var_config.bins)
        var_categ = [vardf_intime] + var_categ
        # TODO: get the correct cosmic scale factor
        scale_intime_to_lightdata = 0.008
        scale_intime = scale_intime_to_lightdata
        weights_categ = [scale_intime*np.ones_like(vardf_intime)] + weights_categ
        colors = ["silver"] + colors
        labels = ["In-time\nCosmic"] + labels

    else:
        vardf_intime = None
        print("No intime cosmics provided")

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
                    yerr=data_err,
                    color='black', 
                    fmt='o', markersize=5, capsize=3, linewidth=1.5,
                    label='Data')


    if ratio:
        # MC 
        if syst:
            mc_content_ratio = total_mc / total_mc # dummy
            # TODO: other syst sources
            syst_err = mc_stat_err
            mc_stat_err_ratio = mc_stat_err / total_mc
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
        data_eylow, data_eyhigh = return_data_stat_err(total_data)

        data_ratio = total_data / total_mc
        data_ratio_eylow = data_eylow / total_mc
        data_ratio_eyhigh = data_eyhigh / total_mc
        data_ratio = np.nan_to_num(data_ratio, nan=-999.)
        data_ratio_eylow = np.nan_to_num(data_ratio_eylow, nan=0.)
        data_ratio_eyhigh = np.nan_to_num(data_ratio_eyhigh, nan=0.)
        
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
        mc_handles = [h for i, h in enumerate(handles) if i != data_handle_index and 'Unc.' not in labels_orig[i]]
        mc_labels = [f"{label} ({frac*100:.1f}%)"
                            for label, frac in zip(labels[::-1], breakdown_fractions[::-1])]
        ordered_handles.extend(mc_handles)
        ordered_labels.extend(mc_labels)

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

    # ==== set y-axis limit ====
    max_data_with_err = np.max(total_mc + mc_stat_err)
    ax.set_ylim(0., ax_ylim_ratio* max_data_with_err)

    # ==== option to plot vertical lines ====
    if vline is not None:
        for v in vline:
            ymax = ax.get_ylim()[1]
            ax.vlines(x=v, ymin=0, ymax=ymax*0.75, color='red', linestyle='--')

    # ==== textbox to add approval rank and GENIE version ====
    # decide text location based on the distribution
    n_firsthalf = np.sum(total_mc[:len(var_config.bins)//2])
    n_secondhalf = np.sum(total_mc[len(var_config.bins)//2:])
    if n_firsthalf > n_secondhalf:
        textloc_x, textloc_ha = 0.95, 'right'
    else:
        textloc_x, textloc_ha = 0.05, 'left'

    # GEINE version
    ax.text(textloc_x, 0.57, r"GENIE v3.4.0 AR23_00i_00_000", transform=ax.transAxes, 
            fontsize=11.5, color='gray',
            ha=textloc_ha, va='top')

    # SBND approval rank
    if approval == "internal":
        approval_text = r"$\mathbf{SBND}$ Internal"
    elif approval == "preliminary":
        approval_text = r"$\mathbf{SBND}$ Preliminary"

    ax.text(textloc_x, 0.65, 
            approval_text, 
            transform=ax.transAxes, 
            ha=textloc_ha, va='top',
            fontsize=20, color='gray')

    if save_fig:
        plt.savefig(save_name+".pdf", bbox_inches='tight')
    plt.show()

    return {"cuts": cuts, 
            "total_mc": total_mc, 
            "total_data": total_data}



def hist_plot_genie_sb(type,
              evtdf, vardf, 
              vardf_data, var_intime,
              bins,
              var_config="",
              plot_labels=["", "", ""], legendloc="right",
              ax_ylim_ratio=1.5,
              ratio = False,
              syst = False,
              vline = None,
              approval="internal",
              save_fig=False, save_name=None): 

    assert len(evtdf) == len(vardf)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    tot_mc_count = len(evtdf)
    tot_data_count = len(vardf_data)

    # repeat later numu labels twice
    labels = genie_labels[:4] 
    for i in range(4, len(genie_labels)):
        labels += [genie_labels[i]] * 2
    colors = genie_colors[:4]
    for i in range(4, len(genie_colors)):
        colors += [genie_colors[i]] * 2
 
    cut_cosmic = IsCosmic(evtdf)
    cut_nu_outfv = IsNuOutFV(evtdf)
    cut_nu_infv_nu_other = IsNuInFV_NuOther(evtdf)
    cut_nu_infv_numu_nc = IsNuInFV_NumuNC(evtdf)
    print("numu NC", cut_nu_infv_numu_nc.sum())
    # cut_nu_infv_numu_coh = IsNuInFV_NumuCC_COH(evtdf)
    # break down numu modes into signal and background
    cut_nu_infv_numu_othermode_s = IsNuInFV_NumuCC_OtherMode(evtdf) & IsNuInFV_NumuCC_1p0pi(evtdf)
    cut_nu_infv_numu_othermode_b = IsNuInFV_NumuCC_OtherMode(evtdf) & ~IsNuInFV_NumuCC_1p0pi(evtdf)
    cut_nu_infv_numu_cc_dis_s = IsNuInFV_NumuCC_DIS(evtdf) & IsNuInFV_NumuCC_1p0pi(evtdf)
    cut_nu_infv_numu_cc_dis_b = IsNuInFV_NumuCC_DIS(evtdf) & ~IsNuInFV_NumuCC_1p0pi(evtdf)
    cut_nu_infv_numu_cc_res_s = IsNuInFV_NumuCC_RES(evtdf) & IsNuInFV_NumuCC_1p0pi(evtdf)
    cut_nu_infv_numu_cc_res_b = IsNuInFV_NumuCC_RES(evtdf) & ~IsNuInFV_NumuCC_1p0pi(evtdf)
    cut_nu_infv_numu_cc_mec_s = IsNuInFV_NumuCC_MEC(evtdf) & IsNuInFV_NumuCC_1p0pi(evtdf)
    cut_nu_infv_numu_cc_mec_b = IsNuInFV_NumuCC_MEC(evtdf) & ~IsNuInFV_NumuCC_1p0pi(evtdf)
    cut_nu_infv_numu_cc_qe_s = IsNuInFV_NumuCC_QE(evtdf) & IsNuInFV_NumuCC_1p0pi(evtdf)
    cut_nu_infv_numu_cc_qe_b = IsNuInFV_NumuCC_QE(evtdf) & ~IsNuInFV_NumuCC_1p0pi(evtdf)
    cuts = [cut_cosmic, cut_nu_outfv, cut_nu_infv_nu_other, cut_nu_infv_numu_nc, 
            cut_nu_infv_numu_othermode_b, cut_nu_infv_numu_othermode_s,
            cut_nu_infv_numu_cc_dis_b, cut_nu_infv_numu_cc_dis_s,
            cut_nu_infv_numu_cc_res_b, cut_nu_infv_numu_cc_res_s,
            cut_nu_infv_numu_cc_mec_b, cut_nu_infv_numu_cc_mec_s,
            cut_nu_infv_numu_cc_qe_b, cut_nu_infv_numu_cc_qe_s]


    # ratio = False

    # --- Plot template
    if ratio:
        fig, axs = plt.subplots(2, 1, figsize=(8.5, 8.5), 
                               sharex=True, gridspec_kw={'height_ratios': [4, 1]})
        fig.subplots_adjust(hspace=0.1)
        ax = axs[0]
        ax_r = axs[1]

        # # --- Data
        total_data, bins = np.histogram(vardf_data, bins=bins)
        data_err = np.sqrt(total_data)
        ax.errorbar(bin_centers, total_data, yerr=data_err, 
                    fmt='o', color='black')  # error bars

    else:
        fig, ax = plt.subplots(figsize=(8.5, 7))

    # --- MC
    # collect all MC + intime
    var_categ = [var_intime] + [vardf[i] for i in cuts]
    scale_mc = evtdf.pot_weight.unique()[0]
    # TODO:
    scale_intime_to_lightdata = 0.008
    scale_intime = scale_intime_to_lightdata
    weights_categ = [scale_intime*np.ones_like(var_intime)] + [scale_mc*np.ones_like(vardf[i]) for i in cuts] 
    colors = ["silver"] + colors
    labels = ["In-time\nCosmic"] + labels
    hatches = [None] * len(labels)
    for i in range(5, len(labels), 2):
        hatches[i] = '////'
    # hatches = ['*'] * len(colors)

    mc_stack, _, _ = ax.hist(var_categ,
                                bins=bins,
                                weights=weights_categ,
                                stacked=True,
                                color=colors,
                                label=labels,
                                edgecolor='none',
                                linewidth=0,
                                density=False,
                                histtype='stepfilled')

    # Switch to stacked bar plots using ax.bar for hatches
    bottom = np.zeros(len(bins) - 1)
    for i, (data, w, color, label, hatch) in enumerate(zip(var_categ, weights_categ, colors, labels, hatches)):
        hist_vals, _ = np.histogram(data, bins=bins, weights=w)
        ax.bar(
            bin_centers,
            hist_vals,
            width=np.diff(bins),
            bottom=bottom,
            color='none',
            hatch=hatch,
            edgecolor='white',
            linewidth=0.0,
            align='center'
        )
        bottom += hist_vals


    # ---- MC stat err
    each_mc_hist_data = []
    each_mc_hist_err2 = []  # sum of squared weights for error

    for data, w in zip(var_categ, weights_categ):
        hist_vals, _ = np.histogram(data, bins=bins, weights=w)
        hist_err2, _ = np.histogram(data, bins=bins, weights=np.square(w))
        each_mc_hist_data.append(hist_vals)
        each_mc_hist_err2.append(hist_err2)

    total_mc = np.sum(each_mc_hist_data, axis=0)
    total_mc_err2 = np.sum(each_mc_hist_err2, axis=0)
    mc_stat_err = np.sqrt(total_mc_err2)

    if syst:
        print("No syst provided")
        # syst_err = syst_uncert_dict[var_config.var_save_name]

        # vardf = np.clip(vardf, bins[0], bins[-1] - eps)
        # n_univ = []
        # for i in range(100):
        #     weights = evtdf.mc.GENIE["univ_{}".format(i)]
        #     weights = np.nan_to_num(weights, nan=1)
        #     n, bins = np.histogram(vardf, bins=bins, weights=weights)
        #     n_univ.append(n)
        # n_cv, _ = np.histogram(vardf, bins=bins)
        # genie_unc = np.std(n_univ, axis=0)
        # genie_unc_frac = genie_unc / n_cv
        # bin_centers = 0.5 * (bins[:-1] + bins[1:])
        # # plt.errorbar(bin_centers, n_cv, yerr=unc, fmt='o', color='black')

        # print(genie_unc_frac)
        # print(syst_err)
        # syst_err = np.sqrt(genie_unc_frac**2 + syst_err**2)

        # ax.bar(
        #     bin_centers,
        #     2 * syst_err * total_mc,
        #     width=np.diff(bins),
        #     bottom=total_mc - syst_err * total_mc,
        #     facecolor='none',             # transparent fill
        #     edgecolor='dimgray',            # outline color of the hatching
        #     hatch='xxx',                 # hatch pattern similar to ROOT's 3004
        #     linewidth=0.0,
        #     label='Syst. Unc.'
        # )

    else:
        ax.bar(
        bin_centers,
            2 * mc_stat_err,
            width=np.diff(bins),
            bottom=total_mc - mc_stat_err,
            facecolor='none',             # transparent fill
            edgecolor='dimgray',            # outline color of the hatching
            hatch='xxx',                 # hatch pattern similar to ROOT's 3004
            linewidth=0.0,
            label='MC Stat. Unc.'
        )


    ax.set_xlim(bins[0], bins[-1])
    if ratio == False: # only plot xlabel if we're not plotting the ratio panel
        ax.set_xlabel(plot_labels[0])
    ax.set_ylabel(plot_labels[1])
    ax.set_title(plot_labels[2])

    if ratio:
        # MC stat err
        mc_stat_err_ratio = mc_stat_err / total_mc
        mc_content_ratio = total_mc / total_mc
        mc_stat_err_ratio = np.nan_to_num(mc_stat_err_ratio, nan=0.)
        mc_content_ratio = np.nan_to_num(mc_content_ratio, nan=-999.)
        if syst:
            ax_r.bar(
                bin_centers,
                2 * syst_err,
                width=np.diff(bins),
                bottom=mc_content_ratio - syst_err,
                facecolor='none',             # transparent fill
                edgecolor='dimgray',            # outline color of the hatching
                hatch='xxxx',                 # hatch pattern similar to ROOT's 3004
                linewidth=0.0,
                label='Syst. Unc.'
            )

        else:
            ax_r.bar(
                bin_centers,
                2*mc_stat_err_ratio,
                width=np.diff(bins),
                bottom=mc_content_ratio - mc_stat_err_ratio,
                facecolor='none',             # transparent fill
                edgecolor='dimgray',          # outline color of the hatching
                hatch='xxxx',                 # hatch pattern similar to ROOT's 3004
                linewidth=0.0,
                label='MC Stat. Unc.'
            )

        # data/MC ratio err
        data_eylow, data_eyhigh = return_data_stat_err(total_data)

        data_ratio = total_data / total_mc
        data_ratio_eylow = data_eylow / total_mc
        data_ratio_eyhigh = data_eyhigh / total_mc
        data_ratio = np.nan_to_num(data_ratio, nan=-999.)
        data_ratio_eylow = np.nan_to_num(data_ratio_eylow, nan=0.)
        data_ratio_eyhigh = np.nan_to_num(data_ratio_eyhigh, nan=0.)
        
        data_ratio_errors = data_ratio_eylow + data_ratio_eyhigh
        ax_ratio.errorbar(bin_centers, data_ratio, yerr=data_ratio_errors,
                        fmt='o', color='black', label='Data',
                        markersize=5, capsize=3, linewidth=1.5)

        ax_r.errorbar(bin_centers, data_ratio,
                  yerr=np.vstack((data_ratio_eylow, data_ratio_eyhigh)),
                  fmt='o', color='black')
                #   , label='Data')
                #   markersize=5, capsize=3, linewidth=1.5)
        
        # if highest value is greater than 2.0, set ylim to 2.0
        ax_r.axhline(1.0, color='red', linestyle='--', linewidth=1)
        
        ax_r.grid(True)
        ax_r.minorticks_on()
        ax_r.grid(which='minor', linestyle=':', linewidth=0.5, color='gray', alpha=0.5)

        ax_r.set_xlabel(plot_labels[0])
        ax_r.set_ylabel("Data/MC")

        # if np.max(data_ratio + data_ratio_eyhigh) > 2.0:
        ax_r.set_ylim(0.4, 1.6)


    # # --- Legend
    # # dummy plot for legend
    accum_sum = [np.sum(data) for data in mc_stack]
    accum_sum = [0.] + accum_sum
    total_sum = accum_sum[-1]
    individual_sums = [accum_sum[i + 1] - accum_sum[i] for i in range(len(accum_sum) - 1)]
    fractions = [(count / total_sum) * 100 for count in individual_sums]
    legend_labels = []
    for lidx in range(5):
        legend_labels.append("{}\n({:0.1f}%/{:0.1f}%)".format(genie_labels[::-1][lidx], fractions[::-1][lidx*2], fractions[::-1][lidx*2+1]))
    for label, frac in zip(genie_labels[::-1][5:], fractions[::-1][10:]):
        legend_labels.append("{}\n({:0.1f}%)".format(label, frac))
    legend_labels += ["In-time\nCosmic\n({:0.1f}%)".format(fractions[0])]

    # --- Legend
    # 1. Data handle (black marker)
    # Data handle with errorbar for legend
    data_handle = plt.errorbar([], [], yerr=[], marker='o', color='black', linestyle='None')
    # data_err_handle = Line2D([0], [0], marker='o', color='black', linestyle='None', markersize=0,
    #                          markerfacecolor='none', markeredgecolor='none', 
    #                          linewidth=1.5, capsize=3, label='Data error', 
    #                          alpha=0.0)  # invisible marker, just for errorbar

    # 2. MC handles
    mc_handles = []
    mc_labels = []

    # First 5 genie categories (with two fractions each)
    genie_colors_reversed = genie_colors[::-1]
    genie_labels_reversed = genie_labels[::-1]
    for i in range(5):
        mc_handles.append(Patch(facecolor=genie_colors_reversed[i], edgecolor='none'))
        mc_labels.append("{}\n({:0.1f}%/{:0.1f}%)".format(genie_labels_reversed[i],
                                                        fractions[::-1][i*2],
                                                        fractions[::-1][i*2+1]))

    # Remaining genie categories (single fraction)
    for i, label in enumerate(genie_labels_reversed[5:]):
        mc_handles.append(Patch(facecolor=genie_colors_reversed[i+5], edgecolor='none'))
        mc_labels.append("{}\n({:0.1f}%)".format(label, fractions[::-1][10+i]))

    # In-time cosmic
    mc_handles.append(Patch(facecolor="silver", edgecolor='none'))
    mc_labels.append("In-time\nCosmic\n({:0.1f}%)".format(fractions[0]))

    # 3. Uncertainty handle
    if syst:
        unc_handle = Patch(facecolor='none', edgecolor='dimgray', hatch='xxxx', linewidth=0, label='Syst. Unc.')
    else:
        unc_handle = Patch(facecolor='none', edgecolor='dimgray', hatch='xxxx', linewidth=0, label='MC Stat. Unc.')

    # Combine handles and labels in order: Data, MC, Unc
    if ratio:
        handles = [data_handle] + mc_handles + [unc_handle]
        # labels_for_legend = [h.get_label() for h in handles]
        labels_for_legend = ["Data"] + mc_labels + ["Syst. Unc."]
    else:
        handles = mc_handles + [unc_handle]
        # labels_for_legend = [h.get_label() for h in handles]
        labels_for_legend = mc_labels + ["Syst. Unc."]

    # Draw legend
    leg = ax.legend(handles, labels_for_legend,
                    loc='upper left',
                    fontsize=11.5,
                    frameon=False,
                    ncol=4,
                    bbox_to_anchor=(0.007, 0.995))


    leg_height = leg.get_bbox_to_anchor().height
    max_data_with_err = np.max(total_mc + mc_stat_err)
    if ratio:
        max_data_with_err = np.max(total_data + data_eyhigh)
    ax.set_ylim(0., ax_ylim_ratio * max_data_with_err)

    if vline is not None:
        for v in vline:
            ymax = ax.get_ylim()[1]
            ax.vlines(x=v, ymin=0, ymax=ymax*0.7, color='red', linestyle='--')

    # --- "SBND Preliminary" textbox
    # decide if the distribution is tilted to the right or left
    n_firsthalf = np.sum(total_mc[:len(bins)//2])
    n_secondhalf = np.sum(total_mc[len(bins)//2:])
    if n_firsthalf > n_secondhalf:
        textloc_x = 0.97
        ha = 'right'
    else:
        textloc_x = 0.04
        ha = 'left'

    if approval == "internal":
        ax.text(textloc_x, 0.57, r"$\mathbf{SBND}$ Internal", transform=ax.transAxes, 
                fontsize=20, color='rosybrown',
                ha=ha, va='top')
        ax.text(textloc_x, 0.49, r"GENIE v3.4.0 AR23_00i_00_000", transform=ax.transAxes, 
                fontsize=11.5, color='gray',
                ha=ha, va='top')

    elif approval == "preliminary":
        ax.text(textloc_x, 0.57, r"$\mathbf{SBND}$ Preliminary", transform=ax.transAxes, 
                fontsize=20, color='gray',
                ha=ha, va='top')
        ax.text(textloc_x, 0.49, r"GENIE v3.4.0 AR23_00i_00_000", transform=ax.transAxes, 
                fontsize=11.5, color='gray',
                ha=ha, va='top')

    # separate legend box with S / B hatches
    example_signal = Patch(facecolor="black", edgecolor='white', label='Signal')
    example_background = Patch(facecolor="black", edgecolor='white', hatch='////', linewidth=0, label='Background')
    if legendloc == "left":
        box_ax = ax.inset_axes([0.17, 0.56, 0.13, 0.13], transform=ax.transAxes)
    elif legendloc == "right":
        box_ax = ax.inset_axes([0.7, 0.56, 0.13, 0.13], transform=ax.transAxes)
    box_ax.axis('off')
    mini_legend = Legend(
        box_ax,
        handles=[example_signal, example_background],
        labels=['Signal', 'Background'],
        loc='center',
        fontsize=11.3,
        frameon=False,
        borderpad=0.7,
        handlelength=2,
        handleheight=1.5,
        ncol=2,
        fancybox=True,
        framealpha=1.0
    )
    box_ax.add_artist(mini_legend)

    if save_fig:
        plt.savefig(save_name+".pdf", bbox_inches='tight')
    plt.show()

    # bolder figure lines?
    # ax.tick_params(width=2, length=10)
    # for spine in ax.spines.values():
    #     spine.set_linewidth(2)
    
    ret_dict = {"cuts": cuts}
    return ret_dict



def hist_plot_pdg_breakdown(
              trkdf, vardf, 
              vardf_data, 
              bins,
              plot_labels=["", "", ""],
              ratio = True,
              vline = None,
              save_fig=False, save_name=None): 

    assert len(trkdf) == len(vardf)

    tot_mc_count = len(trkdf)
    tot_data_count = len(vardf_data)

    labels = ["Muon", "Pion", "Proton", "Electron", "Photon", "Other"]
    colors = ["blue", "green", "red", "purple", "orange", "gray"]

    cut_muon = (np.abs(trkdf.pfp.trk.truth.p.pdg) == 13)
    cut_pion = (np.abs(trkdf.pfp.trk.truth.p.pdg) == 211)
    cut_proton = (np.abs(trkdf.pfp.trk.truth.p.pdg) == 2212)
    cut_electron = (np.abs(trkdf.pfp.trk.truth.p.pdg) == 11)
    cut_photon = (np.abs(trkdf.pfp.trk.truth.p.pdg) == 22)
    cut_other = ~cut_muon & ~cut_pion & ~cut_proton & ~cut_electron & ~cut_photon

    cuts = [cut_muon, cut_pion, cut_proton, cut_electron, cut_photon, cut_other]

    # --- Plot template
    if ratio:
        fig, axs = plt.subplots(2, 1, figsize=(8.5, 8), 
                               sharex=True, gridspec_kw={'height_ratios': [4, 1]})
        fig.subplots_adjust(hspace=0.05)
        ax = axs[0]
        ax_r = axs[1]
    else:
        fig, ax = plt.subplots()

    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    # --- Data
    total_data, bins = np.histogram(vardf_data, bins=bins)
    data_err = np.sqrt(total_data)
    ax.errorbar(bin_centers, total_data, yerr=data_err, 
                fmt='o', color='black')  # error bars

    # --- MC
    # collect all MC + intime
    var_categ = [vardf[i] for i in cuts]
    scale_mc = trkdf.pot_weight.unique()[0]
    # TODO:
    weights_categ = [scale_mc*np.ones_like(vardf[i]) for i in cuts] 
    colors = colors
    labels = labels

    mc_stack, _, _ = ax.hist(var_categ,
                                bins=bins,
                                weights=weights_categ,
                                stacked=True,
                                color=colors,
                                label=labels,
                                edgecolor='none',
                                linewidth=0,
                                density=False,
                                histtype='stepfilled')


    # ---- MC stat err
    each_mc_hist_data = []
    each_mc_hist_err2 = []  # sum of squared weights for error

    for data, w in zip(var_categ, weights_categ):
        hist_vals, _ = np.histogram(data, bins=bins, weights=w)
        hist_err2, _ = np.histogram(data, bins=bins, weights=np.square(w))
        each_mc_hist_data.append(hist_vals)
        each_mc_hist_err2.append(hist_err2)

    total_mc = np.sum(each_mc_hist_data, axis=0)
    total_mc_err2 = np.sum(each_mc_hist_err2, axis=0)
    mc_stat_err = np.sqrt(total_mc_err2)

    ax.bar(
       bin_centers,
        2 * mc_stat_err,
        width=np.diff(bins),
        bottom=total_mc - mc_stat_err,
        facecolor='none',             # transparent fill
        edgecolor='dimgray',            # outline color of the hatching
        hatch='xxxx',                 # hatch pattern similar to ROOT's 3004
        linewidth=0.0,
        label='MC Stat. Unc.'
    )


    ax.set_xlim(bins[0], bins[-1])
    if ratio == False: # only plot xlabel if we're not plotting the ratio panel
        ax_r.set_xlabel(plot_labels[0])
    ax.set_ylabel(plot_labels[1])
    ax.set_title(plot_labels[2])

    if ratio:
        # MC stat err
        mc_stat_err_ratio = mc_stat_err / total_mc
        mc_content_ratio = total_mc / total_mc
        mc_stat_err_ratio = np.nan_to_num(mc_stat_err_ratio, nan=0.)
        mc_content_ratio = np.nan_to_num(mc_content_ratio, nan=-999.)
        ax_r.bar(
            bin_centers,
            2*mc_stat_err_ratio,
            width=np.diff(bins),
            bottom=mc_content_ratio - mc_stat_err_ratio,
            facecolor='none',             # transparent fill
            edgecolor='dimgray',          # outline color of the hatching
            hatch='xxxx',                 # hatch pattern similar to ROOT's 3004
            linewidth=0.0,
            label='MC Stat. Unc.'
        )

        # data/MC ratio err
        data_eylow, data_eyhigh = return_data_stat_err(total_data)

        data_ratio = total_data / total_mc
        data_ratio_eylow = data_eylow / total_mc
        data_ratio_eyhigh = data_eyhigh / total_mc
        data_ratio = np.nan_to_num(data_ratio, nan=-999.)
        data_ratio_eylow = np.nan_to_num(data_ratio_eylow, nan=0.)
        data_ratio_eyhigh = np.nan_to_num(data_ratio_eyhigh, nan=0.)
        
        #data_ratio_errors = data_ratio_eylow + data_ratio_eyhigh
        #ax_ratio.errorbar(bin_centers, data_ratio, yerr=data_ratio_errors,
        #                 fmt='o', color='black', label='Data',
        #                 markersize=5, capsize=3, linewidth=1.5)

        ax_r.errorbar(bin_centers, data_ratio,
                  yerr=np.vstack((data_ratio_eylow, data_ratio_eyhigh)),
                  fmt='o', color='black')
                #   , label='Data')
                #   markersize=5, capsize=3, linewidth=1.5)
        
        # if highest value is greater than 2.0, set ylim to 2.0
        ax_r.axhline(1.0, color='red', linestyle='--', linewidth=1)
        
        ax_r.grid(True)
        ax_r.minorticks_on()
        ax_r.grid(which='minor', linestyle=':', linewidth=0.5, color='gray', alpha=0.5)

        ax_r.set_xlabel(plot_labels[0])
        ax_r.set_ylabel("Data/MC")

        if np.max(data_ratio + data_ratio_eyhigh) > 2.0:
            ax_r.set_ylim(0., 2.0)

    # --- Legend
    accum_sum = [np.sum(data) for data in mc_stack]
    accum_sum = [0.] + accum_sum
    total_sum = accum_sum[-1]
    individual_sums = [accum_sum[i + 1] - accum_sum[i] for i in range(len(accum_sum) - 1)]
    fractions = [(count / total_sum) * 100 for count in individual_sums]
    legend_labels = [f"{label} ({frac:.1f}%)" for label, frac in zip(labels[::-1], fractions[::-1])]
    legend_labels += ["Data ({})".format(tot_data_count), "MC Stat. Unc."]
    leg = ax.legend(legend_labels, 
                    loc='upper left', 
                    fontsize=10, 
                    frameon=False, 
                    ncol=3, 
                    bbox_to_anchor=(0.015, 0.985))
    leg_height = leg.get_bbox_to_anchor().height
    max_data_with_err = np.max(total_data + data_eyhigh)
    # ax.set_ylim(0., 1.05 * max_data_with_err + leg_height)
    ax.set_ylim(0., 1.4 * max_data_with_err)

    if vline is not None:
        # Draw a vertical line from y=0 up to the max value of the histogram
        for v in vline:
            ymax = ax.get_ylim()[1]
            ax.vlines(x=v, ymin=0, ymax=ymax*0.75, color='red', linestyle='--')

    # --- "SBND Preliminary" textbox
    # decide if the distribution is tilted to the right or left
    n_firsthalf = np.sum(total_data[:len(bins)//2])
    n_secondhalf = np.sum(total_data[len(bins)//2:])
    if n_firsthalf > n_secondhalf:
        textloc_x = 0.95
        textloc_ha = 'right'
    else:
        textloc_x = 0.05
        textloc_ha = 'left'

    ax.text(textloc_x, 0.65, "SBND Preliminary", 
            transform=ax.transAxes, 
            fontsize=12, 
            # fontweight='bold', 
            color='gray',
            ha=textloc_ha, 
            va='top')


    if save_fig:
        plt.savefig(save_name, bbox_inches='tight', dpi=300)
    plt.show()

    # bolder figure lines?
    # ax.tick_params(width=2, length=10)
    # for spine in ax.spines.values():
    #     spine.set_linewidth(2)
    
    ret_dict = {"cuts": cuts}
    return ret_dict

def match_trkdf_to_slcdf(trkdf, slcdf):
    # trkdf: df to match
    # slcdf: df to match to
    matched_trkdf = trkdf.reset_index(level=[3]).loc[slcdf.index].reset_index().set_index(trkdf.index.names)
    return matched_trkdf


def plot_univ_hist(evtdf, var_config, syst_name, univ_events, cv_events):
    # n_univ = len(evtdf[syst_name].columns)
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
            if i == 0:
                plt.hist(var_config.bin_centers, bins=var_config.bins, weights=univ_events[i], histtype="step", color="gray", label="Universe")
            else:
                plt.hist(var_config.bin_centers, bins=var_config.bins, weights=univ_events[i], histtype="step", color="gray")

    # CV universe
    plt.hist(var_config.bin_centers, bins=var_config.bins, weights=cv_events, histtype="step", color="k", label="Central Value")

    plt.xlim(var_config.bins[0], var_config.bins[-1])
    plt.xlabel(var_config.var_labels[1])
    plt.ylabel("Events / Bin")
    plt.title(syst_name)

    # plt.text(0.05, 0.05, "SBND Internal", fontsize=12, color="rosybrown", ha="left", va="bottom", transform=ax.transAxes)
    # plt.legend(reverse=True, frameon=False)
    plt.legend(frameon=False)
    plt.show();

def get_covariance(univ_events, cv_events):
    n_univ = univ_events.shape[0]
    n_bins = univ_events.shape[1]

    Covariance_Frac = np.zeros((n_bins, n_bins))
    Covariance = np.zeros((n_bins, n_bins))

    # I'm looping & calculating with the CV value for clarity, 
    # but techincally np.cov should also be fine under the assumption of gaussian universes that we're using
    for uidx in range(n_univ):
        for i in range(univ_events.shape[1]):
            for j in range(univ_events.shape[1]):
                nom_i = cv_events[i] 
                nom_j = cv_events[j] 

                univ_i = univ_events[uidx, i] 
                univ_j = univ_events[uidx, j] 

                cov_entry = (univ_i - nom_i) * (univ_j - nom_j)
                frac_cov_entry = ((univ_i - nom_i) / nom_i) * ( (univ_j - nom_j) / nom_j)

                # TODO: this clipping exists in the uboone code, but I'm not sure why..?
                # if cov_entry > 0:
                #     this_cov = max( cov_entry, eps * scale_factor)
                # else:
                #     this_cov = min( cov_entry, eps * scale_factor)

                # if frac_cov_entry > 0:
                #     this_frac_cov = max( frac_cov_entry, eps * scale_factor)
                # else:
                #     this_frac_cov = min( frac_cov_entry, eps * scale_factor)

                Covariance[i, j] += cov_entry
                Covariance_Frac[i, j] += frac_cov_entry

    Covariance = Covariance / n_univ
    Covariance_Frac = Covariance_Frac / n_univ
    Correlation = np.zeros_like(Covariance)
    for i in range(len(cv_events)):
        for j in range(len(cv_events)):
            Correlation[i, j] = Covariance[i, j] / (np.sqrt(Covariance[i, i]) * np.sqrt(Covariance[j, j]))

    return {"Covariance_Frac": Covariance_Frac, 
            "Covariance": Covariance,
            "Correlation": Correlation,
            }


def plot_heatmap(matrix, var_config, title="", save_fig=False, save_fig_name=None):

    unif_bin = np.linspace(0., float(len(var_config.bins) - 1), len(var_config.bins))
    extent = [unif_bin[0], unif_bin[-1], unif_bin[0], unif_bin[-1]]

    x_edges = np.array(var_config.bins)
    y_edges = np.array(var_config.bins)
    x_tick_positions = (unif_bin[:-1] + unif_bin[1:]) / 2
    y_tick_positions = (unif_bin[:-1] + unif_bin[1:]) / 2

    x_labels = bin_range_labels(x_edges)
    y_labels = bin_range_labels(y_edges)

    fig, ax = plt.subplots(figsize=(10, 10))
    plt.imshow(matrix, extent=extent, origin="lower")
    plt.colorbar(shrink=0.7)
    plt.xticks(x_tick_positions, x_labels, rotation=45, ha="right")
    plt.yticks(y_tick_positions, y_labels)

    plot_labels = var_config.var_labels
    plt.xlabel(plot_labels[0], fontsize=20)
    plt.ylabel(plot_labels[1], fontsize=20)
    for i in range(matrix.shape[0]):      # rows (y)
        for j in range(matrix.shape[1]):  # columns (x)
            value = matrix[i, j]
            if not np.isnan(value):  # skip NaNs
                plt.text(
                    j + 0.5, i + 0.5,
                    f"{value:.2f}",
                    ha="center", va="center",   
                    color=get_text_color(value),
                    fontsize=10
                )
    plt.title(title, fontsize=20)
    if save_fig:
        plt.savefig("{}.png".format(save_fig_name), bbox_inches='tight', dpi=300)
    plt.show();


def get_univ_rates(evtdf, var_config, syst_name):
    var = evtdf[var_config.var_evt_reco_col]
    cv_events, _ = np.histogram(var, bins=var_config.bins)

    univ_events = []
    for uidx in range(len(evtdf[syst_name].columns)):
        weights = evtdf[syst_name]["univ_{}".format(uidx)]
        # fill nan (these are non-neutrino events) with 1
        weights = np.where(np.isnan(weights), 1, weights)
        n, bins = np.histogram(var, bins=var_config.bins, weights=weights)
        univ_events.append(n)
    univ_events = np.array(univ_events)
    return univ_events, cv_events

def get_clipped_evts(df, var_col, bins):
    eps = 1e-8
    var = df[var_col]
    var = np.clip(var, bins[0], bins[-1] - eps)

    if 'pot_weight' in df.columns:
        weights = df.loc[:, 'pot_weight']
    else:
        print("No pot_weight column found, return 1 as pot scale")
        weights = np.ones_like(var)
    return var, weights

def get_distributions(mc_evt_df, mc_nu_df, var_config):
    # Total MC reco muon momentum: for fake data
    var_total_mc, weights_total_mc = get_clipped_evts(mc_evt_df, var_config.var_evt_reco_col, var_config.bins)

    # --- all events, selected ---
    # mc_evt_df divided into topology modes for subtraction from data in future
    # first item in list is the signal topology
    mc_evt_df_divided = [mc_evt_df[mc_evt_df.nuint_categ == mode]for mode in topology_list]

    # Reco variable distribution for each 'nuint_categ' for stack plot and subtraction from the fake data
    var_per_nuint_categ_mc, weights_per_categ = [], []
    for mode in topology_list:
        var, weights = get_clipped_evts(mc_evt_df[mc_evt_df.nuint_categ == mode], var_config.var_evt_reco_col, var_config.bins)
        var_per_nuint_categ_mc.append(var)
        weights_per_categ.append(weights)

    # Reco variable distribution for each genie mode
    var_per_genie_mode_mc, weights_per_genie_mode = [], []
    for mode in genie_mode_list:
        var, weights = get_clipped_evts(mc_evt_df[mc_evt_df.genie_categ == mode], var_config.var_evt_reco_col, var_config.bins)
        var_per_genie_mode_mc.append(var)
        weights_per_genie_mode.append(weights)

    # --- signal events ---
    # selected, for response matrix
    # Signal event's reco muon momentum after the event selection
    var_signal_sel_reco, weight_signal = get_clipped_evts(mc_evt_df[mc_evt_df.nuint_categ == 1], var_config.var_evt_reco_col, var_config.bins)

    # Signal event's true muon momentum after the event selection
    var_signal_sel_truth, weight_true_signal = get_clipped_evts(mc_evt_df[mc_evt_df.nuint_categ == 1], var_config.var_evt_truth_col, var_config.bins)

    # total generated, for efficiency vector
    # Signal event's true muon momentum without event selection
    var_truth_signal, weight_truth_signal = get_clipped_evts(mc_nu_df[mc_nu_df.nuint_categ == 1], var_config.var_nu_col, var_config.bins)

    return {
        "var_total_mc": var_total_mc,
        "weights_total_mc": weights_total_mc,
        "mc_evt_df_divided": mc_evt_df_divided,
        "var_per_nuint_categ_mc": var_per_nuint_categ_mc,
        "weights_per_categ": weights_per_categ,
        "var_per_genie_mode_mc": var_per_genie_mode_mc,
        "weights_per_genie_mode": weights_per_genie_mode,
        "var_signal_sel_reco": var_signal_sel_reco,
        "weight_signal": weight_signal,
        "var_signal_sel_truth": var_signal_sel_truth,
        "weight_true_signal": weight_true_signal,
        "var_truth_signal": var_truth_signal,
        "weight_truth_signal": weight_truth_signal,
    }


# TODO: load xsec unit factor as a global variable, instead of calculating it every time
def get_genie_univs(cov_type, mc_evt_df, mc_nu_df, var_config, syst_name, n_univ=100, plot=False):
    """
    for the GENIE uncertainty on the xsec measurement
    """
    XSEC_UNIT = 1e-38

    if cov_type == "xsec":
        print("generating covariance for xsec, using scale factor: {}".format(XSEC_UNIT))
        scale_factor = XSEC_UNIT
    elif cov_type == "rate":
        scale_factor = 1.0
    else:
        raise ValueError("Invalid covariance type: {}, choose xsec or rate".format(cov_type))

    ret_dists = get_distributions(mc_evt_df, mc_nu_df, var_config)
    nevts_signal_truth, _     = np.histogram(ret_dists["var_truth_signal"],     bins=var_config.bins, weights=ret_dists["weight_truth_signal"])
    nevts_signal_sel_truth, _ = np.histogram(ret_dists["var_signal_sel_truth"], bins=var_config.bins, weights=ret_dists["weight_signal"])
    nevts_signal_sel_reco, _  = np.histogram(ret_dists["var_signal_sel_reco"],  bins=var_config.bins, weights=ret_dists["weight_signal"])

    signal_cv = nevts_signal_sel_reco * scale_factor # = Response @ true_signal

    univ_events = []
    univ_effs   = []
    univ_smears = []
    for uidx in range(n_univ):
        # univ_col_evt = ("mc", syst_name, "univ_{}".format(uidx)) #, "", "", "", "") #, "", "")
        # univ_col_mc = ("mc", syst_name, "univ_{}".format(uidx), "")
        univ_col = ("mc", syst_name, "univ_{}".format(uidx))

        # ---- uncertainty on the signal rate ----
        # only consider effect on the response matrix for the signal channel
        if cov_type == "xsec":
            true_signal_univ, _ = np.histogram(ret_dists["var_truth_signal"], bins=var_config.bins, 
                                               weights=ret_dists["weight_truth_signal"]*mc_nu_df[mc_nu_df.nuint_categ == 1][univ_col])
            
            reco_vs_true = get_smear_matrix(ret_dists["var_signal_sel_truth"], ret_dists["var_signal_sel_reco"], [var_config.bins, var_config.bins],
                                            weights=mc_evt_df[mc_evt_df.nuint_categ == 1][univ_col], plot=plot)
            univ_smears.append(reco_vs_true)

            eff = get_eff(reco_vs_true, true_signal_univ) 
            univ_effs.append(eff)

            Response_univ = get_response_matrix(reco_vs_true, eff, var_config.bins, plot=plot)
            signal_univ = Response_univ @ nevts_signal_truth # note that we multiply the CV signal rate!
            # signal_univ = signal_cv

        elif cov_type == "rate":
            signal_univ, _ = np.histogram(ret_dists["var_signal_sel_reco"], bins=var_config.bins, 
                                          weights=mc_evt_df[mc_evt_df.nuint_categ == 1][univ_col])

        else:
            raise ValueError("Invalid covariance type: {}, choose xsec or rate".format(cov_type))

        # ---- uncertainty on the background rate ----
        # loop over background categories
        # + univ background - cv background
        # note: cv background subtraction cancels out with the cv background subtraction for the cv event rate. 
        #       doing it anyways for the plot of universes on background subtracted event rate.
        for this_mc_evt_df in ret_dists["mc_evt_df_divided"][1:]:
            weights = this_mc_evt_df[univ_col].copy()
            weights[np.isnan(weights)] = 1 ## IMPORTANT: make nan weights to 1. to ignore them
            this_var = this_mc_evt_df[var_config.var_evt_reco_col]
            this_var = np.clip(this_var, var_config.bins[0], var_config.bins[-1] - eps)
            background_univ, _ = np.histogram(this_var, bins=var_config.bins, weights=weights)
            background_cv, _   = np.histogram(this_var, bins=var_config.bins)
            signal_univ += background_univ - background_cv

        signal_univ *= scale_factor
        univ_events.append(signal_univ)

    univ_events = np.array(univ_events)
    return univ_events, signal_cv

def cov_from_fraccov(cov_frac, cv_vals):
    cov = np.zeros_like(cov_frac)
    for i in range(cov_frac.shape[0]):
        for j in range(cov_frac.shape[1]):
            cov[i, j] = cov_frac[i, j] * (cv_vals[i] * cv_vals[j])
    return cov

def corr_from_fraccov(cov_frac):
    corr = np.zeros_like(cov_frac)
    for i in range(cov_frac.shape[0]):
        for j in range(cov_frac.shape[1]):
            corr[i, j] = cov_frac[i, j] / np.sqrt(cov_frac[i, i] * cov_frac[j, j])
    return corr


def plot_unfolded_result(unfold, bins, measured, models,
                         plot_labels, model_names, 
                         save_fig=False, save_name=None,
                         closure_test=False):

    # need to divide by bin width for differential xsec
    bin_widths = np.diff(bins)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    fig, ax = plt.subplots(figsize=(8.5, 7))

    # --- stat uncertainties
    UnfoldCov_stat = unfold['StatUnfoldCov']
    Unfold_uncert_stat = np.diag(UnfoldCov_stat)

    # --- syst uncertainties
    UnfoldCov_syst = unfold['SystUnfoldCov']
    Unfold_uncert_syst = np.diag(UnfoldCov_syst)

    # decomose into norm and shape components
    # TODO: the first item in models is the true model
    SystUnfoldCov_norm, SystUnfoldCov_shape = Matrix_Decomp(models[0], UnfoldCov_syst)
    print(np.diag(SystUnfoldCov_norm))
    Unfold_uncert_norm = np.sqrt(np.abs(np.diag(SystUnfoldCov_norm)))
    Unfold_uncert_shape = np.sqrt(np.abs(np.diag(SystUnfoldCov_shape)))

    # --- plot
    # unfolded result
    Unfolded = unfold['unfold']
    UnfoldedCov = unfold["UnfoldCov"]
    Unfolded_perwidth = Unfolded / bin_widths

    # set err to 0 for closure test
    if closure_test:
        dummy_err = np.zeros_like(Unfolded_perwidth)
        bar_handle = plt.errorbar(bin_centers, Unfolded_perwidth, yerr=dummy_err, fmt='o', color='black')

    else:
        # plot shape syst and stat as error bars
        Unfold_uncert_stat_perwidth = Unfold_uncert_stat / bin_widths
        Unfold_uncert_shape_perwidth = Unfold_uncert_shape / bin_widths
        # Unfold_uncert_stat_shape_perwidth = Unfold_uncert_stat_perwidth + Unfold_uncert_shape_perwidth
        Unfold_uncert_stat_shape_perwidth = Unfold_uncert_shape_perwidth
        bar_handle = plt.errorbar(bin_centers, Unfolded_perwidth, yerr=Unfold_uncert_stat_shape_perwidth, fmt='o', color='black')

        # plot syst norm component as histogram at the bottom
        Unfold_uncert_norm_perwidth = Unfold_uncert_norm / bin_widths
        norm_handle = plt.bar(bin_centers, Unfold_uncert_norm_perwidth, width=bin_widths, label='Syst. error (norm)', alpha=0.5, color='gray')

    # Also divide measured & model by bin width
    measured_perwidth = measured / bin_widths
    reco_handle, = plt.step(bins, np.append(measured_perwidth, measured_perwidth[-1]), where='post', label='Meausred Signal (Input)')

    # --- get chi2 values for each model to compare
    chi2_vals = []
    p_values = []
    model_handles = []
    model_labels = []
    for midx, model in enumerate(models):

        model_smeared = unfold['AddSmear'] @ model
        # chi2_val = chi2(Unfolded, model_smeared, UnfoldCov_syst)
        # chi2_vals.append(chi2_val)
        # p_val = 0

        chi2_val, p_val = get_chi2(Unfolded, model_smeared, UnfoldCov_syst)
        chi2_vals.append(chi2_val)
        p_values.append(p_val)

        # # calculated the p-value for the chi2 test
        # import scipy.stats as stats
        # dof = len(bins) - 1
        # p_value = 1 - stats.chi2.cdf(chi2_val, dof)
        # p_values.append(p_value)

        # Example histograms
        # statistic, p_value = chisquare(Unfolded, f_exp=model_smeared)
        # chi2_vals.append(statistic)
        # p_values.append(p_value)

        model_smeared_perwidth = model_smeared / bin_widths
        model_handle, = plt.step(bins, np.append(model_smeared_perwidth, model_smeared_perwidth[-1]), where='post')
        model_handles.append(model_handle)
        model_labels.append(f'$A_c \\otimes$ {model_names[midx]} ($\chi^2$ = {chi2_vals[midx]:.2f}/{len(bins)-1}), p-value = {p_values[midx]:.3f}')

    # legend
    if closure_test:
        handles = [bar_handle, reco_handle] + model_handles
        labels = ['Unfolded Asimov Data', 'Measured Signal'] + model_labels
    else:
        handles = [bar_handle, norm_handle, reco_handle] + model_handles
        labels = ['Unfolded', 'Norm. Syst. Unc.', 'Measured Signal'] + model_labels
    plt.legend(handles, labels, 
               loc='upper left', fontsize=12, frameon=False, ncol=1, bbox_to_anchor=(0.02, 0.98))

    n_firsthalf = np.sum(Unfolded_perwidth[:len(bins)//2])
    n_secondhalf = np.sum(Unfolded_perwidth[len(bins)//2:])
    if n_firsthalf > n_secondhalf:
        textloc_x = 0.95
        ha = 'right'
    else:
        textloc_x = 0.05
        ha = 'left'

    ax.text(textloc_x, 0.65, r"$\mathbf{SBND}$ Preliminary", 
            transform=ax.transAxes, 
            fontsize=16, 
            color='gray',
            ha=ha, 
            va='top')
    # ax.text(textloc_x, 0.58, r"Work in Progress", 
    #         transform=ax.transAxes, 
    #         fontsize=16, 
    #         color='gray',
    #         ha=ha, 
    #         va='top')

    # leave space for legend
    plt.xlabel(plot_labels[0])
    plt.ylabel(plot_labels[1])
    plt.xlim(bins[0], bins[-1])
    plt.ylim(0., np.max(Unfolded_perwidth)*1.7)

    if save_fig:
        plt.savefig(save_name, bbox_inches='tight')
    plt.show()

def plot_unfolded_data(unfold, bins, measured, models,
                         plot_labels, model_names,
                         save_fig=False, save_name=None):

    # need to divide by bin width for differential xsec
    bin_widths = np.diff(bins)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    fig, ax = plt.subplots(figsize=(6, 5))

    Unfolded = unfold['unfold']

    # --- stat uncertainties
    UnfoldCov_stat = unfold['StatUnfoldCov'] 
    # print(UnfoldCov_stat)
    UnfoldCov_stat = np.diag(measured/XSEC_UNIT)
    # print(UnfoldCov_stat)
    Unfold_uncert_stat = np.diag(UnfoldCov_stat)  * XSEC_UNIT
    # print(Unfold_uncert_stat)
    # print(UnfoldCov_stat)

    # mask = Unfolded != 0
    UnfoldCov_stat_frac = np.zeros_like(UnfoldCov_stat)
    for i in range(len(Unfolded)):
        for j in range(len(Unfolded)):
            if Unfolded[i] != 0 and Unfolded[j] != 0:
                UnfoldCov_stat_frac[i, j] = UnfoldCov_stat[i, j] / (Unfolded[i] * Unfolded[j])
            else:
                UnfoldCov_stat_frac[i, j] = 0.0
    UnfoldCov_stat_frac = UnfoldCov_stat_frac * XSEC_UNIT**2
    # print(UnfoldCov_stat_frac)


    # --- syst uncertainties
    UnfoldCov_syst = unfold['SystUnfoldCov']
    Unfold_uncert_syst = np.diag(UnfoldCov_syst)

    UnfoldCov_syst_frac = np.zeros_like(UnfoldCov_syst)
    for i in range(len(Unfolded)):
        for j in range(len(Unfolded)):
            if Unfolded[i] != 0 and Unfolded[j] != 0:
                UnfoldCov_syst_frac[i, j] = UnfoldCov_syst[i, j] / (Unfolded[i] * Unfolded[j])
            else:
                UnfoldCov_syst_frac[i, j] = 0.0
    # print(UnfoldCov_syst_frac)

    UnfoldedCov_frac = UnfoldCov_stat_frac + UnfoldCov_syst_frac
    UnfoldedCov = np.zeros_like(UnfoldedCov_frac)
    for i in range(len(Unfolded)):
        for j in range(len(Unfolded)):
            UnfoldedCov[i, j] = Unfolded[i] * Unfolded[j] * UnfoldedCov_frac[i, j]


    # decomose into norm and shape components
    # TODO: the first item in models is the true model
    SystUnfoldCov_norm, SystUnfoldCov_shape = Matrix_Decomp(models[0], UnfoldCov_syst)
    Unfold_uncert_norm = np.sqrt(np.abs(np.diag(SystUnfoldCov_norm)))
    Unfold_uncert_shape = np.sqrt(np.abs(np.diag(SystUnfoldCov_shape)))
    # print(Unfold_uncert_stat)
    # print(Unfold_uncert_shape)

    # --- plot
    # unfolded result
    Unfolded_perwidth = Unfolded / bin_widths

    # set err to 0 for closure test
    # plot shape syst and stat as error bars
    Unfold_uncert_stat_perwidth = Unfold_uncert_stat / bin_widths
    Unfold_uncert_shape_perwidth = Unfold_uncert_shape / bin_widths
    # Plot errorbars so that inner cap is stat, outer cap is syst
    # First, plot stat-only errorbars with smaller capsize (inner)
    bar_handle_stat = plt.errorbar(
        bin_centers, Unfolded_perwidth, 
        yerr=Unfold_uncert_stat_perwidth, 
        fmt='o', color='black', linewidth=1.5, 
        capsize=3, elinewidth=1.5, markeredgewidth=1.5
    )
    # Then, plot total (stat+syst) errorbars with larger capsize (outer), but no marker
    Unfold_uncert_total = np.sqrt(Unfold_uncert_stat_perwidth**2 + Unfold_uncert_shape_perwidth**2)
    bar_handle = plt.errorbar(
        bin_centers, Unfolded_perwidth, 
        yerr=Unfold_uncert_total, 
        fmt='none', ecolor='black', elinewidth=1, capsize=7
    )

    # plot syst norm component as histogram at the bottom
    Unfold_uncert_norm_perwidth = Unfold_uncert_norm / bin_widths
    norm_handle = plt.bar(bin_centers, Unfold_uncert_norm_perwidth, width=bin_widths, label='Syst. error (norm)', alpha=0.5, color='gray')

    # Also divide measured & model by bin width
    measured_perwidth = measured / bin_widths

    # --- get chi2 values for each model to compare
    chi2_vals = []
    model_handles = []
    model_labels = []
    for midx, model in enumerate(models):

        model_smeared = unfold['AddSmear'] @ model
        chi2_val = chi2(Unfolded, model_smeared, UnfoldedCov)
        # chi2_val = chi2(Unfolded, model_smeared, UnfoldCov)
        chi2_vals.append(chi2_val)

        model_smeared_perwidth = model_smeared / bin_widths
        model_handle, = plt.step(bins, np.append(model_smeared_perwidth, model_smeared_perwidth[-1]), where='post', linewidth=2)
        model_handles.append(model_handle)
        model_labels.append(f'$A_c \\otimes$ {model_names[midx]} ($\chi^2$ = {chi2_vals[midx]:.0f}/{len(bins)-1})')

    # legend
    # Show both inner and outer cap in legend, and write legend accordingly
    # Only unfolded result and models
    # For errorbars, show both stat and syst errorbar handles in legend
    # Plot main legend for unfolded result
    handles = [bar_handle_stat, bar_handle]
    labels = [
        'Stat. Unc.',
        'Stat. $\oplus$ Shape Syst. Unc.'
    ]
    valign = 0.6
    main_legend = plt.legend(handles, labels, 
                                loc='upper left', fontsize=10, frameon=False, ncol=1, bbox_to_anchor=(valign, 0.88),
                                title="Unfolded Data", title_fontproperties={'weight': 'bold'})

    plt.gca().add_artist(main_legend)
    model_legend = plt.legend(model_handles, model_labels, 
                                loc='upper left', fontsize=10, frameon=False, ncol=1, bbox_to_anchor=(valign, 0.98))

    # leave space for legend
    plt.xlabel(plot_labels[0])
    plt.ylabel(plot_labels[1])
    plt.xlim(bins[0], bins[-1])
    plt.ylim(0., np.max(Unfolded_perwidth)*1.3)

    if save_fig:
        plt.savefig(save_name, bbox_inches='tight')
    plt.show()

def get_chi2(data, model, cov):
    chi2_value = (data - model) @ np.linalg.inv(cov) @ (data - model)
    ndof = len(data)
    p_value = 1 - chi2.cdf(chi2_value, df=ndof)
    return chi2_value, p_value



def plot_topology_breakdown(var_categ, weights_categ, var_total, weights_total, bins,
                            plot_labels, 
                            colors, labels, 
                            save_fig=False, save_name=None):

    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    plt.figure(figsize=(8.5, 6))

    # TODO: make this more general?
    # our categroy breakdown list have signal at the front
    # stack in reverse order so that signal is on top
    var_categ = var_categ[::-1]
    weights_categ = weights_categ[::-1]
    colors = colors[::-1]
    labels = labels[::-1]
    mc_stack, _, _ = plt.hist(var_categ,
                              bins=bins,
                              weights=weights_categ,
                              stacked=True,
                              color=colors,
                              edgecolor='none',
                              linewidth=0,
                              density=False,
                              histtype='stepfilled')

    # background_cv = mc_stack[-1] - mc_stack[0]
    background_cv = mc_stack[-2]

    # use MC as fake data for closure test
    totmc, _ = np.histogram(var_total, bins=bins, weights=weights_total)
    fake_data     = totmc
    fake_data_err = np.sqrt(totmc)
    plt.errorbar(bin_centers, fake_data, yerr=fake_data_err, 
                 fmt='o', color='black')

    # note the % breakdown in the legend
    accum_sum = [0.] + [np.sum(data) for data in mc_stack]
    total_sum = accum_sum[-1]
    individual_sums = [accum_sum[i + 1] - accum_sum[i] for i in range(len(accum_sum) - 1)]
    fractions = [(count / total_sum) * 100 for count in individual_sums]
    legend_labels = [f"{label} ({frac:.1f}%)" for label, frac in zip(labels[::-1], fractions[::-1])]
    legend_labels.append("Fake Data")
    plt.legend(legend_labels, 
                loc='upper left', 
                fontsize=10, 
                frameon=False, 
                ncol=3, 
                bbox_to_anchor=(0.02, 0.98))

    n_firsthalf = np.sum(totmc[:len(bins)//2])
    n_secondhalf = np.sum(totmc[len(bins)//2:])
    if n_firsthalf > n_secondhalf:
        textloc_x = 0.95
        ha = 'right'
    else:
        textloc_x = 0.05
        ha = 'left'

    plt.text(textloc_x, 0.5, "SBND Simulation\nSBND Preliminary", 
            transform=plt.gca().transAxes, 
            fontsize=12, 
            color='gray',
            ha=ha, 
            va='top')

    # leave whitespace at the top for the legend
    plt.xlabel(plot_labels[0])
    plt.ylabel(plot_labels[1])
    plt.xlim(bins[0], bins[-1])
    plt.ylim(0., 1.45 * fake_data.max())

    if save_fig:
        plt.savefig(save_name, bbox_inches='tight', dpi=300)
    plt.show()
    
    return fake_data, background_cv


def plot_genie_breakdown(var_categ, weights_categ, var_total, weights_total, bins,
                         plot_labels, 
                         colors, labels, 
                         save_fig=False, save_name=None):

    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    plt.figure(figsize=(8.5, 6))

    # stack in reverse order so that signal is on top 
    var_categ = var_categ[::-1]
    weights_categ = weights_categ[::-1]
    colors = colors[::-1]
    labels = labels[::-1]

    mc_stack, _, _ = plt.hist(var_categ,
                                bins=bins,
                                weights=weights_categ,
                                stacked=True,
                                color=colors,
                                edgecolor='none',
                                linewidth=0,
                                density=False,
                                histtype='stepfilled')

    # use MC as fake data for closure test
    totmc, bins = np.histogram(var_total, bins=bins, weights=weights_total)
    fake_data     = totmc
    fake_data_err = np.sqrt(totmc)
    plt.errorbar(bin_centers, fake_data, yerr=fake_data_err, 
                 fmt='o', color='black')

    # note the % breakdown in the legend
    accum_sum = [np.sum(data) for data in mc_stack]
    accum_sum = [0.] + accum_sum
    total_sum = accum_sum[-1]
    individual_sums = [accum_sum[i + 1] - accum_sum[i] for i in range(len(accum_sum) - 1)]
    fractions = [(count / total_sum) * 100 for count in individual_sums]
    legend_labels = [f"{label} ({frac:.1f}%)" for label, frac in zip(labels[::-1], fractions[::-1])]
    legend_labels.append("Fake Data")
    plt.legend(legend_labels, 
                loc='upper left', 
                fontsize=10, 
                frameon=False, 
                ncol=3, 
                bbox_to_anchor=(0.02, 0.98))

    n_firsthalf = np.sum(totmc[:len(bins)//2])
    n_secondhalf = np.sum(totmc[len(bins)//2:])
    if n_firsthalf > n_secondhalf:
        textloc_x = 0.95
        ha = 'right'
    else:
        textloc_x = 0.05
        ha = 'left'

    plt.text(textloc_x, 0.5, "SBND Simulation\nSBND Preliminary", 
            transform=plt.gca().transAxes, 
            fontsize=12, 
            color='gray',
            ha=ha, 
            va='top')

    # leave whitespace at the top for the legend
    plt.xlabel(plot_labels[0])
    plt.ylabel(plot_labels[1])
    plt.xlim(bins[0], bins[-1])
    plt.ylim(0., 1.45 * fake_data.max())

    if save_fig:
        plt.savefig(save_name, bbox_inches='tight', dpi=300)
    plt.show()


# Plotters for detector variation analysis
def overlay_hists(evtdfs, vardfs, bins,
              colors, labels,
              plot_labels=["", "", ""],
              vline = None,
              approval="internal",
              save_fig=False, save_name=None): 

    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    fig, axs = plt.subplots(2, 1, figsize=(7.5, 7), 
                            sharex=True, gridspec_kw={'height_ratios': [4, 1]})
    fig.subplots_adjust(hspace=0.05)
    ax = axs[0]
    ax_r = axs[1]

    total_mc_list = []
    mc_stat_err_list = []
    for i, (evtdf, vardf) in enumerate(zip(evtdfs, vardfs)):
        # scale_mc = evtdf.pot_weight.unique()[0]
        total_mc, _ = np.histogram(vardf, bins=bins)
        total_mc_err2, _ = np.histogram(vardf, bins=bins)
        mc_stat_err = np.sqrt(total_mc_err2)
        total_mc_list.append(total_mc)
        mc_stat_err_list.append(mc_stat_err)

        ax.hist(vardf, bins=bins, histtype="step", 
                 color=colors[i], label=labels[i])

    # ax.set_xlabel(plot_labels[0])
    ax.set_ylabel(plot_labels[1])
    ax.set_title(plot_labels[2])
    ax.set_xlim(bins[0], bins[-1])
    ax.legend()

    # ratio plot -- divide by the first df (CV)

    for i, (evtdf, vardf) in enumerate(zip(evtdfs, vardfs)):
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
    
    ax_r.grid(True)
    ax_r.minorticks_on()
    ax_r.grid(which='minor', linestyle=':', linewidth=0.5, color='gray', alpha=0.5)

    ax_r.set_xlabel(plot_labels[0])
    ax_r.set_ylabel("Variation / CV")
    ax_r.set_xlim(bins[0], bins[-1])
    ax_r.set_ylim(0.9, 1.1)

    if vline is not None:
        for v in vline:
            ax.axvline(x=v, color='red', linestyle='--')
            ax_r.axvline(x=v, color='red', linestyle='--')

    # --- approval textbox
    # decide if the distribution is tilted to the right or left
    n_firsthalf = np.sum(total_mc[:len(bins)//2])
    n_secondhalf = np.sum(total_mc[len(bins)//2:])
    if n_firsthalf > n_secondhalf:
        textloc_x = 0.95
        textloc_ha = 'right'
    else:
        textloc_x = 0.05
        textloc_ha = 'left'

    if approval == "internal":
        ax.text(textloc_x, 0.65, "SBND Internal", transform=ax.transAxes, 
                fontsize=14, color='rosybrown',
                ha=textloc_ha, va='top')

    elif approval == "preliminary":
        ax.text(textloc_x, 0.65, "SBND Preliminary", transform=ax.transAxes, 
                fontsize=14, color='gray',
                ha=textloc_ha, va='top')

    if save_fig:
        plt.savefig(save_name, bbox_inches='tight', dpi=300)
    plt.show()

    return total_mc_list