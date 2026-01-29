import numpy as np
from makedf.constants import *
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnchoredText
from matplotlib.offsetbox import AnchoredOffsetbox, DrawingArea, HPacker, VPacker, TextArea
from matplotlib.legend import Legend

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

DETECTOR = "SBND_nohighyz"
# DETECTOR = "SBND"

eps = 1e-6


import pickle
with open("/exp/sbnd/app/users/munjung/xsec/cafpyana/analysis_village/numucc1p0pi/NuINT_uncs/tot_uncert_dict.pkl", "rb") as f:
    syst_uncert_dict = pickle.load(f)
syst_uncert_dict.keys()


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


def hist_plot(type,
              evtdf, vardf, 
              vardf_data, var_intime,
              bins,
              var_config="",
              plot_labels=["", "", ""],
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

    if type == "nu_cosmics":
        labels = nu_cosmics_labels
        colors = nu_cosmics_colors

        cut_cosmic = IsCosmic(evtdf)
        cut_nu_outfv = IsNuOutFV(evtdf)
        cut_nu_infv = IsNuInFV(evtdf)
        cuts = [cut_cosmic, cut_nu_outfv, cut_nu_infv]

    elif type == "topology":
        labels = topology_labels
        colors = topology_colors

        cut_cosmic = IsCosmic(evtdf)
        cut_nu_outfv = IsNuOutFV(evtdf)
        cut_nu_infv_nu_other = IsNuInFV_NuOther(evtdf)
        cut_nu_infv_numu_nc = IsNuInFV_NumuNC(evtdf)
        cut_nu_infv_numu_cc_other = IsNuInFV_NumuCC_Other(evtdf)
        cut_nu_infv_numu_cc_np0pi = IsNuInFV_NumuCC_Np0pi(evtdf)
        cut_nu_infv_numu_cc_1p0pi = IsNuInFV_NumuCC_1p0pi(evtdf)
        cuts = [cut_cosmic, cut_nu_outfv, cut_nu_infv_nu_other, cut_nu_infv_numu_nc, 
                cut_nu_infv_numu_cc_other, cut_nu_infv_numu_cc_np0pi, cut_nu_infv_numu_cc_1p0pi]

    elif type == "genie":
        labels = genie_labels
        colors = genie_colors

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
        

    else:
        raise ValueError("Invalid type: %s, please choose between [nu_cosmics, topolgy, or genie]" % type)

    # ratio = False

    # --- Plot template
    if ratio:
        fig, axs = plt.subplots(2, 1, figsize=(8.5, 8.5), 
                               sharex=True, gridspec_kw={'height_ratios': [4, 1]})
        fig.subplots_adjust(hspace=0.1)
        ax = axs[0]
        ax_r = axs[1]

        # --- Data
        total_data, bins = np.histogram(vardf_data, bins=bins)
        data_err = np.sqrt(total_data)
        ax.errorbar(bin_centers, total_data, yerr=data_err, 
                    fmt='o', color='black', label='Data')  # error bars

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

    if syst:
        syst_err = syst_uncert_dict[var_config.var_save_name]

        n_univ = []
        for i in range(100):
            weights = evtdf.mc.GENIE["univ_{}".format(i)]
            weights = np.nan_to_num(weights, nan=1)
            n, bins = np.histogram(vardf, bins=bins, weights=weights)
            n_univ.append(n)
        n_cv, _ = np.histogram(vardf, bins=bins)
        genie_unc = np.std(n_univ, axis=0)
        print("genie_unc",genie_unc)
        print("n_cv",n_cv)
        genie_unc_frac = genie_unc / n_cv
        # bin_centers = 0.5 * (bins[:-1] + bins[1:])
        # # plt.errorbar(bin_centers, n_cv, yerr=unc, fmt='o', color='black')

        # print(genie_unc_frac)
        # print(syst_err)
        syst_err = np.sqrt(genie_unc_frac**2 + syst_err**2)

        ax.bar(
        bin_centers,
            2 * syst_err * total_mc,
            width=np.diff(bins),
            bottom=total_mc - syst_err * total_mc,
            facecolor='none',             # transparent fill
            edgecolor='dimgray',            # outline color of the hatching
            hatch='xxx',                 # hatch pattern similar to ROOT's 3004
            linewidth=0.0,
            label='Syst. Unc.'
        )
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
                hatch='xxx',                 # hatch pattern similar to ROOT's 3004
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
                hatch='xxx',                 # hatch pattern similar to ROOT's 3004
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

        # if np.max(data_ratio + data_ratio_eyhigh) > 2.0:
        ax_r.set_ylim(0.4, 1.6)

    # --- Legend
    accum_sum = [np.sum(data) for data in mc_stack]
    # need to add the intime contribution
    accum_sum = [0.] + accum_sum 
    total_sum = accum_sum[-1]
    individual_sums = [accum_sum[i + 1] - accum_sum[i] for i in range(len(accum_sum) - 1)]
    fractions = [(count / total_sum) * 100 for count in individual_sums]
    fractions[-1] = 92.1
    mc_legend_labels = [f"{label} ({frac:.1f}%)"
                        for label, frac in zip(labels[::-1], fractions[::-1])]

    handles, labels_orig = ax.get_legend_handles_labels()

    if ratio:
        # Find the Data handle (usually from errorbar) and put it first
        data_handle_index = labels_orig.index('Data')
        data_handle = handles[data_handle_index]

        # Collect MC handles (exclude Data and uncertainties)
        mc_handles = [h for i, h in enumerate(handles) if i != data_handle_index and 'Unc.' not in labels_orig[i]]
        mc_labels = [l for l in mc_legend_labels if l != 'Data' and 'Unc.' not in l]

        # Collect uncertainty handle
        unc_handle = [h for i, h in enumerate(handles) if 'Unc.' in labels_orig[i]]
        unc_label = [l for l in labels_orig if 'Unc.' in l]

        # Reorder handles and labels: Data first, MC next, uncertainty last
        ordered_handles = [data_handle] + mc_handles + unc_handle
        ordered_labels = ['Data'] + mc_labels + unc_label

    else:
        # Collect MC handles (exclude Data and uncertainties)
        mc_handles = [h for i, h in enumerate(handles) if 'Unc.' not in labels_orig[i]]
        mc_labels = [l for l in mc_legend_labels if 'Unc.' not in l]

        # Collect uncertainty handle
        unc_handle = [h for i, h in enumerate(handles) if 'Unc.' in labels_orig[i]]
        unc_label = [l for l in labels_orig if 'Unc.' in l]
        ordered_handles = mc_handles + unc_handle
        ordered_labels = mc_labels + unc_label

    leg = ax.legend(
        ordered_handles,
        ordered_labels,
        loc='upper left',
        fontsize=11.5,
        frameon=False,
        ncol=3,
        bbox_to_anchor=(0.015, 0.985)
    )

    leg_height = leg.get_bbox_to_anchor().height

    max_data_with_err = np.max(total_mc + mc_stat_err)
    if ratio:
        max_data_with_err = np.max(total_data + data_eyhigh)
        # ax.set_ylim(0., 1.05 * max_data_with_err + leg_height)
    ax.set_ylim(0., ax_ylim_ratio* max_data_with_err)

    if vline is not None:
        # Draw a vertical line from y=0 up to the max value of the histogram
        for v in vline:
            ymax = ax.get_ylim()[1]
            ax.vlines(x=v, ymin=0, ymax=ymax*0.75, color='red', linestyle='--')

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
        ax.text(textloc_x, 0.65, r"$\mathbf{SBND}$ Internal", transform=ax.transAxes, 
                fontsize=20, color='rosybrown',
                ha=textloc_ha, va='top')
        ax.text(textloc_x, 0.57, r"GENIE v3.4.0 AR23_00i_00_000", transform=ax.transAxes, 
                fontsize=11.5, color='gray',
                ha=textloc_ha, va='top')

    elif approval == "preliminary":
        ax.text(textloc_x, 0.65, r"$\mathbf{SBND}$ Preliminary", transform=ax.transAxes, 
                fontsize=20, color='gray',
                ha=textloc_ha, va='top')
        ax.text(textloc_x, 0.57, r"GENIE v3.4.0 AR23_00i_00_000", transform=ax.transAxes, 
                fontsize=11.5, color='gray',
                ha=textloc_ha, va='top')


    if save_fig:
        plt.savefig(save_name+".pdf", bbox_inches='tight')
    plt.show()

    # bolder figure lines?
    # ax.tick_params(width=2, length=10)
    # for spine in ax.spines.values():
    #     spine.set_linewidth(2)
    
    if ratio:
        ret_dict = {"cuts": cuts, "total_data": total_data}
    else:
        ret_dict = {"cuts": cuts}
    return ret_dict



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
        syst_err = syst_uncert_dict[var_config.var_save_name]

        vardf = np.clip(vardf, bins[0], bins[-1] - eps)
        n_univ = []
        for i in range(100):
            weights = evtdf.mc.GENIE["univ_{}".format(i)]
            weights = np.nan_to_num(weights, nan=1)
            n, bins = np.histogram(vardf, bins=bins, weights=weights)
            n_univ.append(n)
        n_cv, _ = np.histogram(vardf, bins=bins)
        genie_unc = np.std(n_univ, axis=0)
        genie_unc_frac = genie_unc / n_cv
        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        # plt.errorbar(bin_centers, n_cv, yerr=unc, fmt='o', color='black')

        print(genie_unc_frac)
        print(syst_err)
        syst_err = np.sqrt(genie_unc_frac**2 + syst_err**2)

        ax.bar(
            bin_centers,
            2 * syst_err * total_mc,
            width=np.diff(bins),
            bottom=total_mc - syst_err * total_mc,
            facecolor='none',             # transparent fill
            edgecolor='dimgray',            # outline color of the hatching
            hatch='xxx',                 # hatch pattern similar to ROOT's 3004
            linewidth=0.0,
            label='Syst. Unc.'
        )

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
    n_univ = len(evtdf[syst_name].columns)
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
    plt.legend(reverse=True, frameon=False)
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