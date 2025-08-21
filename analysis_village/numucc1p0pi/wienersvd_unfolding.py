# Wiener SBND Unfolding for numu CC 1p0pi 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib as mpl
from os import path
import sys
import uproot
from tqdm import tqdm

# local imports
from variable_configs import *

sys.path.append('/exp/sbnd/app/users/munjung/xsec/wienersvd/cafpyana')
from analysis_village.unfolding.wienersvd import *
from analysis_village.unfolding.unfolding_inputs import *
from analysis_village.numucc1p0pi.selection_definitions import *
from pyanalib.split_df_helpers import *
from makedf.mcstat import *
from makedf.geniesyst import regen_systematics_sbnd_multisigma, regen_systematics_sbnd_morph
from makedf.constants import *
from pyanalib.variable_calculator import get_cc1p0pi_tki
from selection_definitions import *

plt.style.use("presentation.mplstyle")
cmap = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)


# configure 
save_fig = False
save_fig_dir = "/exp/sbnd/data/users/munjung/plots/wiener_svd/1mu1p"


# ===== samples =====

# -- MC 
# mc_file = "/exp/sbnd/data/users/munjung/xsec/2025B/MC_test.df"
# mc_file = "/exp/sbnd/app/users/munjung/xsec/wienersvd/cafpyana/test.df"
mc_file = "/exp/sbnd/data/users/munjung/xsec/2025B/MCP2025B_bnb_sel_half.df"
mc_split_df = pd.read_hdf(mc_file, key="split")
mc_n_split = get_n_split(mc_file)
print("mc_n_split: %d" %(mc_n_split))
print_keys(mc_file)

N_MAX_CONCAT = 5
mc_keys2load = ['evt', 'hdr', 'mcnuwgtslim']
mc_dfs = load_dfs(mc_file, mc_keys2load, n_max_concat=N_MAX_CONCAT)

mc_evt_df = mc_dfs['evt']
mc_hdr_df = mc_dfs['hdr']
mc_nu_df = mc_dfs['mcnuwgtslim']


# --- light-triggered data
# data_file = "/exp/sbnd/data/users/munjung/xsec/2025B/devdata_raw_daught.df"
data_file = "/exp/sbnd/data/users/munjung/xsec/2025B/DevData_sel.df"
data_split_df = pd.read_hdf(data_file, key="split")
data_n_split = get_n_split(data_file)
print("data_n_split: %d" %(data_n_split))
print_keys(data_file)

data_hdr_file = "/exp/sbnd/data/users/munjung/xsec/2025B/DevData_sel_hdr.df"
data_hdr_df = pd.read_hdf(data_hdr_file, key="split")
data_n_split = get_n_split(data_hdr_file)
print("data_n_split: %d" %(data_n_split))
print_keys(data_hdr_file)

data_keys2load = ['evt', 'trk']
data_dfs = load_dfs(data_file, data_keys2load, n_max_concat=N_max_concat)
data_evt_df = data_dfs['evt']
data_trk_df = data_dfs['trk']

data_keys2load = ['hdr']
data_dfs = load_dfs(data_hdr_file, data_keys2load, n_max_concat=n_max_concat)
data_hdr_df = data_dfs['hdr']


# -- intime data
basepath = "/exp/sbnd/data/users/munjung/xsec/2025B"
fname = "trash/Data_intime_1000.df"
intime_hdr_df = pd.read_hdf(path.join(basepath, fname), "hdr")
intime_evt_df = pd.read_hdf(path.join(basepath, fname), "slc")
intime_trk_df = pd.read_hdf(path.join(basepath, fname), "trk")



# ===== absolute normalization =====

# --- light-triggered data
# nu
data_tot_pot = data_hdr_df['pot'].sum()
pot_str = "{:.2e} $\\times 10^{:.0f}$".format(data_tot_pot/np.power(10,np.floor(np.log10(data_tot_pot))), np.floor(np.log10(data_tot_pot)))
print("data_tot_pot: %.3e" %(data_tot_pot))
print("data_tot_pot: {}".format(pot_str))
data_evt_df["pot_weight"] = np.ones(len(data_evt_df))
data_trk_df["pot_weight"] = np.ones(len(data_trk_df))
# cosmics
data_gates = data_hdr_df.nbnbinfo.sum()
print("data tot gates : %.3e" %(data_gates))

# --- MC
# nu
mc_tot_pot = mc_hdr_df['pot'].sum()
mc_pot_scale = data_tot_pot / mc_tot_pot
print("mc_tot_pot: %.3e" %(mc_tot_pot))
print("mc_pot_scale: %.3e" %(mc_pot_scale))
mc_evt_df["pot_weight"] = mc_pot_scale * np.ones(len(mc_evt_df))
mc_trk_df["pot_weight"] = mc_pot_scale * np.ones(len(mc_trk_df))

# --- intime data
intime_gates = intime_hdr_df[intime_hdr_df['first_in_subrun'] == 1]['noffbeambnb'].sum()
print("intime cosmics data gates: {:.2e}".format(intime_gates))
# TODO: check this scaling factor
f = 0.08
scale_intime_to_lightdata = (1-f)*data_gates/intime_gates
print("goal scale: {:.2f}".format(scale_intime_to_lightdata))
intime_evt_df["gates_weight"] = scale_intime_to_lightdata * np.ones(len(intime_evt_df))
intime_trk_df["gates_weight"] = scale_intime_to_lightdata * np.ones(len(intime_trk_df))


# ===== calculate TKI =====
# TODO: do this in df creation
mudf = mc_evt_df.mu
pdf = mc_evt_df.p
P_mu_col = ("pfp", "trk", "rangeP", "p_muon" )
P_p_col = ("pfp", "trk", "rangeP", "p_proton" )
ret_tki = get_cc1p0pi_tki(mudf, pdf, P_mu_col, P_p_col)
mc_evt_df["del_Tp"] = ret_tki["del_Tp"]
mc_evt_df["del_p"] = ret_tki["del_p"]
mc_evt_df["del_alpha"] = ret_tki["del_alpha"]
mc_evt_df["del_phi"] = ret_tki["del_phi"]
mc_evt_df["del_alpha"] *= 180/np.pi
mc_evt_df["del_phi"] *= 180/np.pi

# truth
P_mu_col = ("pfp", "trk", "truth", "p" , "totp")
P_p_col = ("pfp", "trk", "truth", "p" , "totp")
ret_tki = get_cc1p0pi_tki(mudf, pdf, P_mu_col, P_p_col)
mc_evt_df["del_Tp_truth"] = ret_tki["del_Tp"]
mc_evt_df["del_p_truth"] = ret_tki["del_p"]
mc_evt_df["del_alpha_truth"] = ret_tki["del_alpha"]
mc_evt_df["del_phi_truth"] = ret_tki["del_phi"]
mc_evt_df["del_alpha_truth"] *= 180/np.pi
mc_evt_df["del_phi_truth"] *= 180/np.pi

# nu dfs
mudf = mc_nu_df.mu
pdf = mc_nu_df.p
P_mu_col = ("totp",)
P_p_col = ("totp",)
ret_tki = get_cc1p0pi_tki(mudf, pdf, P_mu_col, P_p_col)
mc_nu_df["del_Tp"] = ret_tki["del_Tp"]
mc_nu_df["del_p"] = ret_tki["del_p"]
mc_nu_df["del_alpha"] = ret_tki["del_alpha"]
mc_nu_df["del_phi"] = ret_tki["del_phi"]
mc_nu_df["del_alpha"] *= 180/np.pi
mc_nu_df["del_phi"] *= 180/np.pi


# Generate MCstat uncertainty universes
n_univ_mcstat = 500
mc_evt_df, MCstat_univ_events = mcstat(mc_evt_df, mc_hdr_df, n_universes=n_univ_mcstat)

# %% [markdown]
# # POT

# %%
## total pot
mc_tot_pot = mc_hdr_df['pot'].sum()
print("mc_tot_pot: %.3e" %(mc_tot_pot))

target_pot = 1e20
mc_pot_scale = target_pot / mc_tot_pot
print("mc_pot_scale: %.3e" %(mc_pot_scale))
mc_pot_scale = 1.

mc_evt_df["pot_weight"] = mc_pot_scale * np.ones(len(mc_evt_df))

# %% [markdown]
# # Constants

# %%
# TODO: z-dependence?
# flux file, units: /m^2/10^6 POT 
# 50 MeV bins
fluxfile = "/exp/sbnd/data/users/munjung/flux/sbnd_original_flux.root"
flux = uproot.open(fluxfile)
print(flux.keys())

# numu flux
numu_flux = flux["flux_sbnd_numu"].to_numpy()
bin_edges = numu_flux[1]
flux_vals = numu_flux[0]

plt.hist(bin_edges[:-1], bins=bin_edges, weights=flux_vals, histtype="step", linewidth=2)
plt.xlabel("Neutrino Energy [GeV]")
plt.ylabel("Flux [/m$^{2}$/10$^{6}$ POT]")
plt.title("SBND $\\nu_\\mu$ Flux")

# if save_fig:
#     plt.savefig("{}/sbnd-flux.pdf".format(save_fig_dir))
plt.savefig("sbnd-flux.pdf", bbox_inches='tight')
plt.show()

# get integrated flux
integrated_flux = flux_vals.sum()
integrated_flux /= 1e4 # to cm2
INTEGRATED_FLUX = integrated_flux * mc_tot_pot / 1e6 # POT
print("Integrated flux: %.3e" % INTEGRATED_FLUX)

# %%
V_SBND = 380 * 380 * 440 # cm3, the active volume of the detector 
NTARGETS = RHO * V_SBND * N_A / M_AR
print("# of targets: ", NTARGETS)

# %%
# set to 1 for event rates
XSEC_UNIT = 1 / (INTEGRATED_FLUX * NTARGETS)

# XSEC_UNIT = 1
print("xsec unit: ", XSEC_UNIT)

# %% [markdown]
# # Set up utils and selections according to target channel

# %%
# to make it work with the get_covariance function...

# Turn GENIE multisigma into multisims with 2 universes
for knob in regen_systematics_sbnd_multisigma:
    mc_evt_df[(knob, "univ_0", "", "", "", "", "", "")] = mc_evt_df[knob].ps1
    mc_evt_df[(knob, "univ_1", "", "", "", "", "", "")] = mc_evt_df[knob].ms1

# Turn GENIE unisim into multisims with 1 universe
for knob in regen_systematics_sbnd_morph:
    mc_evt_df[(knob, "univ_0", "", "", "", "", "", "")] = mc_evt_df[knob].morph

# %%
# classify events into categories
mc_evt_df.loc[:,'nuint_categ'] = get_int_category(mc_evt_df)
mc_nu_df.loc[:,'nuint_categ'] = get_int_category(mc_nu_df)

mc_evt_df.loc[:,'genie_categ'] = get_genie_category(mc_evt_df)
mc_nu_df.loc[:,'genie_categ'] = get_genie_category(mc_nu_df)

print(mc_evt_df.nuint_categ.value_counts())
print(mc_nu_df.nuint_categ.value_counts()) # won't have -1 because nudf is all nu events

# %% [markdown]
# # Choose Variable to Unfold

# %%
save_fig = False
# choose a variable to unfold, defined in variable_configs.py
# var_config = VariableConfig.muon_momentum()
# var_config = VariableConfig.muon_direction()
# var_config = VariableConfig.proton_momentum()
# var_config = VariableConfig.proton_direction()
# var_config = VariableConfig.tki_del_Tp()
# var_config = VariableConfig.tki_del_alpha()
var_config = VariableConfig.tki_del_phi()

# %% [markdown]
# # Make dfs for analysis

# %% [markdown]
# np.clip is for including underflow events into the first bin and overflow events into the last bin

# %%
# Total MC reco muon momentum: for fake data
eps = 1e-8
var_total_mc = mc_evt_df[var_config.var_evt_reco_col]
var_total_mc = np.clip(var_total_mc, var_config.bins[0], var_config.bins[-1] - eps)
weights_total_mc = mc_evt_df.loc[:, 'pot_weight']

# --- all events, selected ---
# mc_evt_df divided into topology modes for subtraction from data in future
# first item in list is the signal topology
mc_evt_df_divided = [mc_evt_df[mc_evt_df.nuint_categ == mode]for mode in mode_list]

# Reco variable distribution for each 'nuint_categ' for stack plot and subtraction from the fake data
var_per_nuint_categ_mc = [mc_evt_df[mc_evt_df.nuint_categ == mode][var_config.var_evt_reco_col]for mode in mode_list]
var_per_nuint_categ_mc = [s.clip(var_config.bins[0], var_config.bins[-1] - eps) for s in var_per_nuint_categ_mc]
weights_per_categ = [mc_evt_df.loc[mc_evt_df.nuint_categ == mode, 'pot_weight'] for mode in mode_list]

# Reco variable distribution for each genie mode
var_per_genie_mode_mc = [mc_evt_df[mc_evt_df.genie_categ == mode][var_config.var_evt_reco_col]for mode in genie_mode_list]
var_per_genie_mode_mc = [s.clip(var_config.bins[0], var_config.bins[-1] - eps) for s in var_per_genie_mode_mc]
weights_per_genie_mode = [mc_evt_df.loc[mc_evt_df.genie_categ == mode, 'pot_weight'] for mode in genie_mode_list]


# --- signal events ---
# selected, for response matrix
# Signal event's reco muon momentum after the event selection
var_signal_sel_reco = mc_evt_df[mc_evt_df.nuint_categ == 1][var_config.var_evt_reco_col]
var_signal_sel_reco = np.clip(var_signal_sel_reco, var_config.bins[0], var_config.bins[-1] - eps)
weight_signal = mc_evt_df.loc[mc_evt_df.nuint_categ == 1, 'pot_weight']

# Signal event's true muon momentum after the event selection
var_signal_sel_truth = mc_evt_df[mc_evt_df.nuint_categ == 1][var_config.var_evt_truth_col]
var_signal_sel_truth = np.clip(var_signal_sel_truth, var_config.bins[0], var_config.bins[-1] - eps)
weight_true_signal = mc_evt_df.loc[mc_evt_df.nuint_categ == 1, 'pot_weight']

# total generated, for efficiency vector
# Signal event's true muon momentum without event selection
var_truth_signal = mc_nu_df[mc_nu_df.nuint_categ == 1][var_config.var_nu_col]
var_truth_signal = np.clip(var_truth_signal, var_config.bins[0], var_config.bins[-1] - eps)
weight_truth_signal = np.full_like(var_truth_signal, mc_pot_scale, dtype=float)

# %% [markdown]
# # Response Matrix

# %% [markdown]
# Draw true (before event selection) and reco (after event selection) muon momentum distributions of signal events.
# Print entries for double check.

# %%
nevts_signal_truth, _, _ = plt.hist(var_truth_signal, bins=var_config.bins, weights=weight_truth_signal, histtype="step", label="True Signal")
nevts_signal_sel_reco, _, _ = plt.hist(var_signal_sel_reco, bins=var_config.bins, weights=weight_signal, histtype="step", label="Reco Selected Signal", color="k")
nevts_signal_sel_truth, _, _ = plt.hist(var_signal_sel_truth, bins=var_config.bins, weights=weight_signal, histtype="step", label="True Selected Signal")
print(nevts_signal_truth)
print(nevts_signal_sel_reco)
print(nevts_signal_sel_truth)
plt.legend()
plt.ylabel("Events / Bin")
plt.xlim(var_config.bins[0], var_config.bins[-1])
plt.xlabel(var_config.var_labels[0])
if save_fig:
    plt.savefig("{}/{}-sel_event_rates.pdf".format(save_fig_dir, var_config.var_save_name), bbox_inches='tight')
plt.show()

# %%
bins_2d = var_config.bins# = [np.array([0.2, 2]), np.array([0.2, 2])] # commented out lines for 1 bin MC closure test

save_fig_name = "{}/{}-reco_vs_true".format(save_fig_dir, var_config.var_save_name)
reco_vs_true = get_smear_matrix(var_signal_sel_truth, var_signal_sel_reco, bins_2d, var_labels=var_config.var_labels,
                                save_fig=save_fig, save_fig_name=save_fig_name)
eff = get_eff(reco_vs_true, nevts_signal_truth)
print("eff")
print(eff)

save_fig_name = "{}/{}-response_matrix".format(save_fig_dir, var_config.var_save_name)
Response = get_response_matrix(reco_vs_true, eff, var_config.bins, var_labels=var_config.var_labels,
                               save_fig=save_fig, save_fig_name=save_fig_name)

# %% [markdown]
# # Covariance

# %%
def get_covariance(cov_type, syst_name, n_univ, 
                   nevts_signal_sel_reco, var_signal_sel_truth, var_signal_sel_reco, bins, 
                   plot_labels, save_fig=False, save_fig_name=None):

    if cov_type == "xsec":
        scale_factor = XSEC_UNIT
        print("generating covariance for xsec, using scale factor: {}".format(scale_factor))

    elif cov_type == "event":
        print("generating covariance for event rate")
        scale_factor = 1

    else:
        raise ValueError("Invalid cov_type: {}".format(cov_type))

    
    signal_cv = nevts_signal_sel_reco * scale_factor # = Response @ true_signal

    Covariance_Frac = np.zeros((len(signal_cv), len(signal_cv)))
    Covariance = np.zeros((len(signal_cv), len(signal_cv)))

    # init figure to plot event rates
    fig, ax = plt.subplots()

    univ_events = []
    for uidx in range(n_univ):
        univ_col_evt = (syst_name, "univ_{}".format(uidx), "", "", "", "", "", "")
        univ_col_mc = (syst_name, "univ_{}".format(uidx), "")

        # ---- uncertainty on the signal rate ----
        # GENIE syst need special treatment
        # we don't want uncertainty on the xsec
        # only consider its effect on the response matrix
        if syst_name == "GENIE" and cov_type == "xsec":
            true_signal_univ, _ = np.histogram(var_truth_signal, bins=var_config.bins, 
                                            weights=weight_truth_signal*mc_nu_df[mc_nu_df.nuint_categ == 1][univ_col_mc])
            
            # new response matrix for univ
            reco_vs_true = get_smear_matrix(var_signal_sel_truth, var_signal_sel_reco, var_config.bins, 
                                            weights=mc_evt_df[mc_evt_df.nuint_categ == 1][univ_col_evt], plot=False)

            eff = get_eff(reco_vs_true, true_signal_univ) 

            Response_univ = get_response_matrix(reco_vs_true, eff, bins, plot=False)
            signal_univ = Response_univ @ nevts_signal_truth # note that we multiply the CV signal rate!

        # for other systs, we just take the univ signal event rate
        else:
            signal_univ, _ = np.histogram(var_signal_sel_reco, bins=var_config.bins, 
                                             weights=mc_evt_df[mc_evt_df.nuint_categ == 1][univ_col_evt])

        signal_univ = np.array(signal_univ) * scale_factor

        # ---- uncertainty on the background rate ----
        # loop over background categories
        # + univ background - cv background
        # note: cv background subtraction cancels out with the cv background subtraction for the cv event rate. 
        #       doing it anyways for the plot of universes on background subtracted event rate.
        for this_mc_evt_df in mc_evt_df_divided[1:]:
            weights = this_mc_evt_df[univ_col_evt].copy()
            weights[np.isnan(weights)] = 1 ## IMPORTANT: make nan weights to 1. to ignore them
            this_var = this_mc_evt_df[var_config.var_evt_reco_col]
            this_var = np.clip(this_var, var_config.bins[0], var_config.bins[-1] - eps)
            background_univ, _ = np.histogram(this_var, bins=var_config.bins, weights=weights)
            background_cv, _ = np.histogram(this_var, bins=var_config.bins)
            background_univ = np.array(background_univ) * scale_factor
            background_cv = np.array(background_cv) * scale_factor
            signal_univ += background_univ - background_cv


        univ_events.append(signal_univ)
        plt.hist(var_config.bin_centers, bins=var_config.bins, weights=signal_univ, histtype="step", color="steelblue", linewidth=1)

        # ---- covariance calculation for this universe ----
        # I'm looping & calculating with the CV value for clarity, 
        # but techincally np.cov should also be fine under the assumption of gaussian universes that we're using
        for i in range(len(signal_univ)):
            for j in range(len(signal_univ)):
                nom_i = signal_cv[i] 
                nom_j = signal_cv[j] 

                univ_i = signal_univ[i] 
                univ_j = signal_univ[j] 

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

    plt.hist(var_config.bin_centers, bins=var_config.bins, weights=signal_cv, histtype="step", color="black")

    plt.xlim(var_config.bins[0], var_config.bins[-1])
    plt.xlabel(plot_labels[0])
    plt.ylabel(plot_labels[1])
    plt.title(syst_name)

    # legend from dummies
    plt.hist([], bins=var_config.bins, histtype="step", color="black", label="Central Value")
    plt.hist([], bins=var_config.bins, histtype="step", color="steelblue", label="Universe Value")
    plt.legend()

    Covariance = Covariance / n_univ
    Covariance_Frac = Covariance_Frac / n_univ
    Correlation = np.zeros_like(Covariance)
    for i in range(len(signal_cv)):
        for j in range(len(signal_cv)):
            Correlation[i, j] = Covariance[i, j] / (np.sqrt(Covariance[i, i]) * np.sqrt(Covariance[j, j]))

    if save_fig:
        plt.savefig("{}.png".format(save_fig_name), bbox_inches='tight', dpi=300)
    plt.show()

    return {"Covariance_Frac": Covariance_Frac, 
            "Covariance": Covariance,
            "Correlation": Correlation,
            "cv_events": signal_cv,
            "univ_events": univ_events,
            }

# %%
# pretty heatmap plotter

unif_bin = np.linspace(0., float(len(var_config.bins) - 1), len(var_config.bins))
extent = [unif_bin[0], unif_bin[-1], unif_bin[0], unif_bin[-1]]

x_edges = np.array(var_config.bins)
y_edges = np.array(var_config.bins)
x_tick_positions = (unif_bin[:-1] + unif_bin[1:]) / 2
y_tick_positions = (unif_bin[:-1] + unif_bin[1:]) / 2

x_labels = bin_range_labels(x_edges)
y_labels = bin_range_labels(y_edges)

def plot_heatmap(matrix, title, plot_labels=var_config.var_labels, save_fig=False, save_fig_name=None):
    fig, ax = plt.subplots(figsize=(12, 12))
    plt.imshow(matrix, extent=extent, origin="lower")
    plt.colorbar(shrink=0.7)
    plt.xticks(x_tick_positions, x_labels, rotation=45, ha="right")
    plt.yticks(y_tick_positions, y_labels)
    plt.xlabel(plot_labels[0])
    plt.ylabel(plot_labels[1])
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
    plt.title(title)
    if save_fig:
        plt.savefig("{}.png".format(save_fig_name), bbox_inches='tight', dpi=300)
    plt.show();

# %% [markdown]
# ## Detector
# 
# - use difference between det var samples as a unisim syst

# %%
# currently unavailable... assigining a flat uncertainty

# shape should be (n_bins-1, n_bins-1)
Detector_Covariance_Frac = np.zeros((len(var_config.bins)-1, len(var_config.bins)-1))

det_unc = 0.1
# fully uncorrelated version
for i in range(len(var_config.bins)-1):
    for j in range(len(var_config.bins)-1):
        if i == j:
            Detector_Covariance_Frac[i, j] = det_unc**2
        else:
            Detector_Covariance_Frac[i, j] = 0

# %% [markdown]
# # Corsika
# 
# - use difference between off-beam MC and data as a unisim syst
# - using dfs with off-beam data / MC with the **same** nu event selection applied

# %%
data_keys2load = ['evt', 'hdr']
# data_offbeam_file = '/exp/sbnd/data/users/munjung/xsec/2025B/data_offbeam.df'
data_offbeam_file = "/exp/sbnd/data/users/munjung/xsec/2025B/data_offbeam_full.df"
data_offbeam_dfs = load_dfs(data_offbeam_file, data_keys2load, n_max_concat=n_max_concat)
data_offbeam_evt_df = data_offbeam_dfs['evt']
data_offbeam_hdr_df = data_offbeam_dfs['hdr']
data_offbeam_gates = data_offbeam_hdr_df[data_offbeam_hdr_df['first_in_subrun'] == 1]['noffbeambnb'].sum()
print("intime cosmics data gates: {:.2e}".format(data_offbeam_gates))

mc_keys2load = ['evt', 'hdr']
mc_offbeam_file = '/exp/sbnd/data/users/munjung/xsec/2025B/MCP2025B_offbeam.df'
mc_offbeam_dfs = load_dfs(mc_offbeam_file, mc_keys2load, n_max_concat=n_max_concat)
mc_offbeam_evt_df = mc_offbeam_dfs['evt']
mc_offbeam_hdr_df = mc_offbeam_dfs['hdr']
mc_offbeam_gates = mc_offbeam_hdr_df[mc_offbeam_hdr_df['first_in_subrun'] == 1]['ngenevt'].sum()
print("intime cosmics MC gates: {:.2e}".format(mc_offbeam_gates))

offbeam_scale = data_offbeam_gates/mc_offbeam_gates
print("offbeam scale: {:.2f}".format(offbeam_scale))

# convert del_alpha and del_phi unit to degrees
data_offbeam_evt_df["del_alpha"] = data_offbeam_evt_df["del_alpha"] * 180 / np.pi
data_offbeam_evt_df["del_phi"] = data_offbeam_evt_df["del_phi"] * 180 / np.pi
mc_offbeam_evt_df["del_alpha"] = mc_offbeam_evt_df["del_alpha"] * 180 / np.pi
mc_offbeam_evt_df["del_phi"] = mc_offbeam_evt_df["del_phi"] * 180 / np.pi


# %%
# check that the average of universes' weights is ~ CV value
fig, ax = plt.subplots(2,1, figsize=(6, 6), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
data_avg_events, _, _ = ax[0].hist(data_offbeam_evt_df[var_config.var_evt_reco_col], bins=var_config.bins, histtype="step", color="red", label="Off-Beam Data")
mc_weights = np.ones(len(mc_offbeam_evt_df))*offbeam_scale
mc_events, _, _ = ax[0].hist(mc_offbeam_evt_df[var_config.var_evt_reco_col], bins=var_config.bins, weights=mc_weights, histtype="step", color="black", label="Off-Beam MC")
ax[0].set_xlim(var_config.bins[0], var_config.bins[-1])
ax[0].set_ylabel("Events / Bin")
ax[0].legend()

bin_centers = (var_config.bins[1:] + var_config.bins[:-1]) / 2
ratio = data_avg_events/mc_events
ax[1].hist(bin_centers, bins=var_config.bins, weights=ratio, histtype="step", color="black")
ax[1].set_xlim(var_config.bins[0], var_config.bins[-1])
ax[1].set_xlabel(var_config.var_labels[0])
ax[1].set_ylabel("Data/MC")
ax[1].set_xlabel(var_config.var_labels[0])
ax[1].axhline(1.0, color="gray", linestyle="--")
ax[1].set_ylim(0., 2.)
ax[1].grid(True)

plt.show();

# %%
# use data/MC ratio as bin scale factors to make alt univ MC
cosmics_univ_weights = ratio
# fill inf and nan with 2
cosmics_univ_weights[np.isinf(cosmics_univ_weights)] = 2
cosmics_univ_weights[np.isnan(cosmics_univ_weights)] = 2
print(cosmics_univ_weights)


# %%
# directly calculate cov (can't convert to universes & plug into get_covariance because we're using bin scale factors)
# TODO: probably good to define the cov element calculation as a function called by both places (here and in get_covariance)

scale_factor = 1  # event rate
signal_cv = nevts_signal_sel_reco * scale_factor 

Cosmics_Covariance_Frac = np.zeros((len(signal_cv), len(signal_cv)))
Cosmics_Covariance = np.zeros((len(signal_cv), len(signal_cv)))
univ_events = []

# signal doesn't change
signal_univ, _ = np.histogram(var_signal_sel_reco, bins=var_config.bins)
signal_univ = np.array(signal_univ) * scale_factor

# only background changing is the cosmics
# cosmics are saved as the last item in mc_evt_df_divided
this_mc_evt_df = mc_evt_df_divided[-1]
this_var = this_mc_evt_df[var_config.var_evt_reco_col]
this_var = np.clip(this_var, var_config.bins[0], var_config.bins[-1] - eps)
background_cv, _ = np.histogram(this_var, bins=var_config.bins)
background_univ = np.array(background_cv) * np.array(cosmics_univ_weights) * scale_factor
background_cv = np.array(background_cv) * scale_factor
signal_univ = signal_univ.astype(float)
signal_univ += background_univ - background_cv

univ_events.append(signal_univ)
print(signal_univ)
print(signal_cv)
plt.hist(var_config.bin_centers, bins=var_config.bins, weights=signal_cv, histtype="step", color="black")
plt.hist(var_config.bin_centers, bins=var_config.bins, weights=signal_univ, histtype="step", color="steelblue")

# ---- covariance calculation for this universe ----
# I'm looping & calculating with the CV value for clarity, 
# but techincally np.cov should also be fine under the assumption of gaussian universes that we're using
for i in range(len(signal_univ)):
    for j in range(len(signal_univ)):
        nom_i = signal_cv[i] 
        nom_j = signal_cv[j] 

        univ_i = signal_univ[i] 
        univ_j = signal_univ[j] 

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

        Cosmics_Covariance[i, j] += cov_entry
        Cosmics_Covariance_Frac[i, j] += frac_cov_entry

plt.xlim(var_config.bins[0], var_config.bins[-1])
plt.xlabel(var_config.var_labels[0])
plt.ylabel(var_config.var_labels[1])
plt.title("Cosmics")
plt.show()

# %%
# sign of fractional covariance is different from the other two?
syst_name = "cosmics"
save_fig_name = "{}/{}-{}-cosmics_covariance".format(save_fig_dir, var_config.var_save_name, syst_name)
plot_heatmap(Cosmics_Covariance, "Covariance - cosmics",
             save_fig=save_fig, save_fig_name=save_fig_name)
save_fig_name = "{}/{}-{}-cosmics_covariance_frac".format(save_fig_dir, var_config.var_save_name, syst_name)
plot_heatmap(Cosmics_Covariance_Frac, "Fractional Covariance - cosmics",
             save_fig=save_fig, save_fig_name=save_fig_name)

# %% [markdown]
# ## MC Stat

# %%
# check that the average of universes' weights is ~ CV value
n_univ_mcstat = 500

univ_avg = mc_evt_df.MCstat.mean(axis=1)
fig, ax = plt.subplots(2,1, figsize=(6, 6), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
n_cv, bins = np.histogram(mc_evt_df[var_config.var_evt_reco_col], bins=var_config.bins)
cv_events, _, _ = ax[0].hist(var_config.bin_centers, bins=bins, weights=n_cv, histtype="step", color="black", label="CV")
n_univ_avg, bins = np.histogram(mc_evt_df[var_config.var_evt_reco_col], bins=var_config.bins, weights=univ_avg)
univ_avg_events, _, _ = ax[0].hist(var_config.bin_centers, bins=bins, weights=n_univ_avg, histtype="step", color="red", label="UV Averaged")
ax[0].set_xlim(var_config.bins[0], var_config.bins[-1])
ax[0].set_ylabel("Events / Bin")
ax[0].legend()

bin_centers = (var_config.bins[1:] + var_config.bins[:-1]) / 2
ratio = univ_avg_events/cv_events
ax[1].hist(bin_centers, bins=var_config.bins, weights=ratio, histtype="step", color="black")
ax[1].set_xlim(var_config.bins[0], var_config.bins[-1])
ax[1].set_xlabel(var_config.var_labels[0])
ax[1].set_ylabel("UV avg. / CV")
ax[1].set_xlabel(var_config.var_labels[0])
ax[1].set_ylim(0.9, 1.1)

save_fig_name = "{}/{}-{}-tolerance_test".format(save_fig_dir, var_config.var_save_name, syst_name)
if save_fig:
    plt.savefig("{}.png".format(save_fig_name), bbox_inches='tight', dpi=300)
plt.show()
# tolerance = 0.05
# if np.all(np.abs(ratio - 1) < tolerance):
#     print("The average of universes' weights is within the tolerance of the CV value")
# else:
#     print("The average of universes' weights is not within the tolerance of the CV value")

fig, ax = plt.subplots()
for uidx in range(n_univ_mcstat):
    plt.hist(mc_evt_df[var_config.var_evt_reco_col], bins=var_config.bins, weights=mc_evt_df["MCstat"][f"univ_{uidx}"], histtype="step", color="steelblue")
plt.hist(mc_evt_df[var_config.var_evt_reco_col], bins=var_config.bins, histtype="step", color="black", label="nominal_event")
plt.xlim(var_config.bins[0], var_config.bins[-1])
plt.show()

# %%
syst_name = "MCstat"
n_univ = n_univ_mcstat
cov_type = "event"
plot_labels = [var_config.var_labels[1], "Flux Integrated Event Rate"]
save_fig_name = "{}/{}-{}-univ_{}_rates".format(save_fig_dir, var_config.var_save_name, syst_name, cov_type)
ret_mcstat = get_covariance(cov_type, syst_name, n_univ, 
                            nevts_signal_sel_reco, var_signal_sel_truth, var_signal_sel_reco, var_config.bins, 
                            plot_labels, save_fig=save_fig, save_fig_name=save_fig_name)

cov_type = "xsec"
plot_labels = [var_config.var_labels[1], "xsec"]
save_fig_name = "{}/{}-{}-univ_{}_rates".format(save_fig_dir, var_config.var_save_name, syst_name, cov_type)
ret_mcstat = get_covariance(cov_type, syst_name, n_univ, 
                            nevts_signal_sel_reco, var_signal_sel_truth, var_signal_sel_reco, var_config.bins, 
                            plot_labels, save_fig=save_fig, save_fig_name=save_fig_name)

# %%
# sign of fractional covariance is different from the other two?
syst_name = "MCstat"
save_fig_name = "{}/{}-{}-covariance".format(save_fig_dir, var_config.var_save_name, syst_name)
plot_heatmap(ret_mcstat["Covariance"], "{} Covariance".format(syst_name),
             save_fig=save_fig, save_fig_name=save_fig_name)
save_fig_name = "{}/{}-{}-covariance_frac".format(save_fig_dir, var_config.var_save_name, syst_name)
plot_heatmap(ret_mcstat["Covariance_Frac"], "{} Fractional Covariance".format(syst_name),
             save_fig=save_fig, save_fig_name=save_fig_name)
save_fig_name = "{}/{}-{}-correlation".format(save_fig_dir, var_config.var_save_name, syst_name)
plot_heatmap(ret_mcstat["Correlation"], "{} Correlation".format(syst_name),
             save_fig=save_fig, save_fig_name=save_fig_name)

# %% [markdown]
# ## Flux

# %%
syst_name = "Flux"
n_univ_flux = 500

cov_type = "event"
plot_labels = [var_config.var_labels[1], "Flux Integrated Event Rate"]
save_fig_name = "{}/{}-{}-univ_{}_rates".format(save_fig_dir, var_config.var_save_name, syst_name, cov_type)
ret_flux = get_covariance(cov_type, syst_name, n_univ_flux, 
                          nevts_signal_sel_reco, var_signal_sel_truth, var_signal_sel_reco, var_config.bins, 
                          plot_labels, save_fig=save_fig, save_fig_name=save_fig_name)

cov_type = "xsec"
plot_labels = [var_config.var_labels[1], "xsec"]
save_fig_name = "{}/{}-{}-univ_{}_rates".format(save_fig_dir, var_config.var_save_name, syst_name, cov_type)
ret_flux = get_covariance(cov_type, syst_name, n_univ_flux, 
                          nevts_signal_sel_reco, var_signal_sel_truth, var_signal_sel_reco, var_config.bins, 
                          plot_labels, save_fig=save_fig, save_fig_name=save_fig_name)

# %%
syst_name = "Flux"
save_fig_name = "{}/{}-{}-covariance".format(save_fig_dir, var_config.var_save_name, syst_name)
plot_heatmap(ret_flux["Covariance"], "{} Covariance".format(syst_name),
             save_fig=save_fig, save_fig_name=save_fig_name)
save_fig_name = "{}/{}-{}-covariance_frac".format(save_fig_dir, var_config.var_save_name, syst_name)
plot_heatmap(ret_flux["Covariance_Frac"], "{} Fractional Covariance".format(syst_name),
             save_fig=save_fig, save_fig_name=save_fig_name)
save_fig_name = "{}/{}-{}-correlation".format(save_fig_dir, var_config.var_save_name, syst_name)
plot_heatmap(ret_flux["Correlation"], "{} Correlation".format(syst_name),
             save_fig=save_fig, save_fig_name=save_fig_name)

# %% [markdown]
# ## GENIE

# %% [markdown]
# ### Multisims

# %%
syst_name = "GENIE"
n_univ_genie_multisim = 100

cov_type = "event"
plot_labels = [var_config.var_labels[1], "Flux Integrated Event Rate"]
save_fig_name = "{}/{}-{}-genie_univ_{}_rates".format(save_fig_dir, var_config.var_save_name, syst_name, cov_type)
ret_genie = get_covariance(cov_type, syst_name, n_univ_genie_multisim, 
                          nevts_signal_sel_reco, var_signal_sel_truth, var_signal_sel_reco, var_config.bins, 
                          plot_labels, save_fig=save_fig, save_fig_name=save_fig_name)

cov_type = "xsec"
plot_labels = [var_config.var_labels[1], "xsec"]
save_fig_name = "{}/{}-{}-genie_univ_{}_rates".format(save_fig_dir, var_config.var_save_name, syst_name, cov_type)
ret_genie = get_covariance(cov_type, syst_name, n_univ_genie_multisim, 
                          nevts_signal_sel_reco, var_signal_sel_truth, var_signal_sel_reco, var_config.bins, 
                          plot_labels, save_fig=save_fig, save_fig_name=save_fig_name)

# %%
syst_name = "GENIE"
save_fig_name = "{}/{}-geniemultisim-covariance".format(save_fig_dir, var_config.var_save_name)
plot_heatmap(ret_genie["Covariance"], "{} Covariance".format(syst_name),
             save_fig=save_fig, save_fig_name=save_fig_name)
save_fig_name = "{}/{}-geniemultisim-covariance_frac".format(save_fig_dir, var_config.var_save_name)
plot_heatmap(ret_genie["Covariance_Frac"], "{} Fractional Covariance".format(syst_name),
             save_fig=save_fig, save_fig_name=save_fig_name)
save_fig_name = "{}/{}-geniemultisim-correlation".format(save_fig_dir, var_config.var_save_name)
plot_heatmap(ret_genie["Correlation"], "{} Correlation".format(syst_name),
             save_fig=save_fig, save_fig_name=save_fig_name)

GENIEmultisim_Frac_Covariance = ret_genie["Covariance_Frac"]

# %% [markdown]
# ### Unisims (morph)

# %%
n_univ_genie_morph = 1

cov_type = "event"
GENIEmorph_Frac_Covariance_event = np.zeros((len(var_config.bins)-1, len(var_config.bins)-1))
plot_labels = [var_config.var_labels[1], "Flux Integrated Event Rate"]
for kidx, knob in enumerate(regen_systematics_sbnd_morph):
    syst_name = knob
    save_fig_name = "{}/{}-{}-geniemorph_univ_{}_rates".format(save_fig_dir, var_config.var_save_name, syst_name, cov_type)
    ret_geniemorph = get_covariance(cov_type, syst_name, n_univ_genie_morph, 
                            nevts_signal_sel_reco, var_signal_sel_truth, var_signal_sel_reco, var_config.bins, 
                            plot_labels, save_fig=save_fig, save_fig_name=save_fig_name)
    GENIEmorph_Frac_Covariance_event = GENIEmorph_Frac_Covariance_event + ret_geniemorph["Covariance_Frac"]

cov_type = "xsec"
GENIEmorph_Frac_Covariance = np.zeros((len(var_config.bins)-1, len(var_config.bins)-1))
plot_labels = [var_config.var_labels[1], "xsec"]
for kidx, knob in enumerate(regen_systematics_sbnd_morph):
    syst_name = knob
    save_fig_name = "{}/{}-{}-geniemorph_univ_{}_rates".format(save_fig_dir, var_config.var_save_name, syst_name, cov_type)
    ret_geniemorph = get_covariance(cov_type, syst_name, n_univ_genie_morph, 
                            nevts_signal_sel_reco, var_signal_sel_truth, var_signal_sel_reco, var_config.bins, 
                            plot_labels, save_fig=save_fig, save_fig_name=save_fig_name)
    GENIEmorph_Frac_Covariance = GENIEmorph_Frac_Covariance + ret_geniemorph["Covariance_Frac"]

# %%
# plot the covariance for xsec 
save_fig_name = "{}/{}-geniemorph_covariance_frac".format(save_fig_dir, var_config.var_save_name)
plot_heatmap(GENIEmorph_Frac_Covariance, "GENIE Morph Fractional Covariance",
             save_fig=save_fig, save_fig_name=save_fig_name)

GENIEmorph_correlation = np.zeros_like(GENIEmorph_Frac_Covariance)
for i in range(GENIEmorph_Frac_Covariance.shape[0]):
    for j in range(GENIEmorph_Frac_Covariance.shape[1]):
        GENIEmorph_correlation[i, j] = GENIEmorph_Frac_Covariance[i, j] / np.sqrt(GENIEmorph_Frac_Covariance[i, i] * GENIEmorph_Frac_Covariance[j, j])
save_fig_name = "{}/{}-geniemorph_correlation".format(save_fig_dir, var_config.var_save_name)
plot_heatmap(GENIEmorph_correlation, "GENIE Morph Correlation",
             save_fig=save_fig, save_fig_name=save_fig_name)

# %% [markdown]
# ### Multisigma
# 
# - we treat this as a unisim using 1 sigma 

# %% [markdown]
# #### Inspection

# %%
# plot multisigma event rates for a few knobs
for knob in regen_systematics_sbnd_multisigma[:3]:
    plt.hist(var_signal_sel_reco, bins=var_config.bins, weights=mc_evt_df[mc_evt_df.nuint_categ == 1][knob].ps1, histtype="step", 
            color="C0", label="ps1")
    plt.hist(var_signal_sel_reco, bins=var_config.bins, weights=mc_evt_df[mc_evt_df.nuint_categ == 1][knob].ps2, histtype="step", 
            color="C0", label="ps2", alpha=0.75)
    plt.hist(var_signal_sel_reco, bins=var_config.bins, weights=mc_evt_df[mc_evt_df.nuint_categ == 1][knob].ps3, histtype="step", 
            color="C0", label="ps3", alpha=0.5)
    plt.hist(var_signal_sel_reco, bins=var_config.bins, weights=mc_evt_df[mc_evt_df.nuint_categ == 1][knob].ms1, histtype="step", 
            color="C0", label="ms1", linestyle="--")
    plt.hist(var_signal_sel_reco, bins=var_config.bins, weights=mc_evt_df[mc_evt_df.nuint_categ == 1][knob].ms2, histtype="step", 
            color="C0", label="ms2", linestyle="--", alpha=0.75)
    plt.hist(var_signal_sel_reco, bins=var_config.bins, weights=mc_evt_df[mc_evt_df.nuint_categ == 1][knob].ms3, histtype="step", 
            color="C0", label="ms3", linestyle="--", alpha=0.5)

    # TODO: is there any purpose of using cv? seems redundant...
    # plt.hist(var_signal_sel_reco, bins=var_config.bins, histtype="step", color="black", label="nominal_event")
    plt.hist(var_signal_sel_reco, bins=var_config.bins, weights=mc_evt_df[mc_evt_df.nuint_categ == 1][knob].cv, histtype="step", color="black", label="nominal_event")

    plt.xlim(var_config.bins[0], var_config.bins[-1])
    plt.xlabel(var_config.var_labels[0])
    plt.ylabel(var_config.var_labels[1])
    plt.title(knob)
    plt.legend()
    plt.show()

# %%
# cov_type = "event"
# syst_name = knob
# knob = regen_systematics_sbnd_multisigma[1]
# n_univ = 2

# labels = [var_config.var_labels[1], "Flux Integrated Event Rate"]

# save_fig_name = "{}/{}-{}-genie_univ_events".format(save_fig_dir, var_config.var_save_name, syst_name)
# ret_geniemultisigma = get_covariance(cov_type, syst_name, n_univ, 
#                           nevts_signal_sel_reco, var_signal_sel_truth, var_signal_sel_reco, var_config.bins, 
#                           labels, save_fig=save_fig, save_fig_name=save_fig_name)

# %%
# save_fig_name = "{}/{}-{}-genie_covariance".format(save_fig_dir, var_config.var_save_name, syst_name)
# plot_heatmap(ret_geniemultisigma["Covariance"], "Covariance - genie",
#              save_fig=save_fig, save_fig_name=save_fig_name)
# save_fig_name = "{}/{}-{}-genie_covariance_frac".format(save_fig_dir, var_config.var_save_name, syst_name)
# plot_heatmap(ret_geniemultisigma["Covariance_Frac"], "Fractional Covariance - genie",
#              save_fig=save_fig, save_fig_name=save_fig_name)
# save_fig_name = "{}/{}-{}-genie_correlation".format(save_fig_dir, var_config.var_save_name, syst_name)
# plot_heatmap(ret_geniemultisigma["Correlation"], "Correlation - genie",
#              save_fig=save_fig, save_fig_name=save_fig_name)

# %%
# ---- justification for the unisim treatment
# throw multivariate normal distribution to the bin center using 1 sigma,
# show that the 2 and 3 sigma are close to the 2 and 3 sigmas of the distribution, per bin

# test one knob 
knob = regen_systematics_sbnd_multisigma[1]
labels = [var_config.var_labels[1], "Flux Integrated Event Rate"]
save_fig_name = None
ret_geniemultisigma = get_covariance(cov_type, knob, 2, 
                          nevts_signal_sel_reco, var_signal_sel_truth, var_signal_sel_reco, var_config.bins, 
                          labels, save_fig=save_fig, save_fig_name=save_fig_name)

# plot event rate variations from thrown universes vs. multisigma
n_cv, _ = np.histogram(var_signal_sel_reco, bins=var_config.bins, weights=mc_evt_df[mc_evt_df.nuint_categ == 1][knob].cv)
universes_from_1sigma = np.random.multivariate_normal(n_cv, ret_geniemultisigma["Covariance"], size=5000)
for i in range(len(universes_from_1sigma)):
    plt.hist(var_config.bin_centers, bins=var_config.bins, weights=universes_from_1sigma[i], histtype="step", color="gray")

n_cv , _, _ = plt.hist(var_config.bin_centers, bins=var_config.bins, weights=n_cv, histtype="step", color="black", label="CV")
n_ps1, _, _ = plt.hist(var_signal_sel_reco, bins=var_config.bins, weights=mc_evt_df[mc_evt_df.nuint_categ == 1][knob].ps1, histtype="step", 
        color="C0", label="ps1")
n_ps2, _, _ = plt.hist(var_signal_sel_reco, bins=var_config.bins, weights=mc_evt_df[mc_evt_df.nuint_categ == 1][knob].ps2, histtype="step", 
        color="C0", label="ps2", alpha=0.75)
n_ps3, _, _ = plt.hist(var_signal_sel_reco, bins=var_config.bins, weights=mc_evt_df[mc_evt_df.nuint_categ == 1][knob].ps3, histtype="step", 
        color="C0", label="ps3", alpha=0.5)
n_ms1, _, _ = plt.hist(var_signal_sel_reco, bins=var_config.bins, weights=mc_evt_df[mc_evt_df.nuint_categ == 1][knob].ms1, histtype="step", 
        color="C0", label="ms1", linestyle="--")
n_ms2, _, _ = plt.hist(var_signal_sel_reco, bins=var_config.bins, weights=mc_evt_df[mc_evt_df.nuint_categ == 1][knob].ms2, histtype="step", 
        color="C0", label="ms2", linestyle="--", alpha=0.75)
n_ms3, _, _ = plt.hist(var_signal_sel_reco, bins=var_config.bins, weights=mc_evt_df[mc_evt_df.nuint_categ == 1][knob].ms3, histtype="step", 
        color="C0", label="ms3", linestyle="--", alpha=0.5)

plt.xlim(var_config.bins[0], var_config.bins[-1])
plt.xlabel(var_config.var_labels[0])
plt.ylabel(var_config.var_labels[1])
plt.title(knob)
plt.legend()
plt.show()

# look at the distributinon of weights in one bin
bin_idx = 2
this_wgts = universes_from_1sigma[:, bin_idx]
plt.hist(this_wgts, bins=21, histtype="step", color="black", label="weights from CAF ps1/ms1")

# sig levels from thrown universes
labels = ["3$\\sigma$", "2$\\sigma$", "1$\\sigma$"]
wgts_sig = np.percentile(this_wgts, np.array([(1-0.997)/2., (1-0.95)/2., (1-0.68)/2.,
                                      0.5+(0.68)/2., 0.5+(0.95)/2., 0.5+(0.997)/2.])*100)
wgts_cv = np.percentile(this_wgts, 50)
plt.axvline(wgts_cv, color="green", label="CV from weights")
plt.text(wgts_cv*1.01, plt.ylim()[1]*0.9, "$\\mu$", color="green", fontsize=12)
for sidx, s in enumerate(wgts_sig):
    if sidx <= 2:
        linestyle = "-"
        plt.axvline(s, color="C2", linestyle=linestyle, alpha=0.5+0.25*(sidx%3))
        plt.text(s, plt.ylim()[1]*0.9, labels[::-1][sidx%3], color="C2", fontsize=12)
    else:
        linestyle = "--"
        plt.axvline(s, color="C2", linestyle=linestyle, alpha=1-0.25*(sidx%3))
        plt.text(s, plt.ylim()[1]*0.9, labels[sidx%3], color="C2", fontsize=12)

# sig levels from cafs
plt.axvline(n_cv[bin_idx], color="blue")
plt.text(n_cv[bin_idx], plt.ylim()[1]*0.8, "cv", color="blue", fontsize=12)
plt.axvline(n_ps1[bin_idx], color="C0")
plt.text(n_ps1[bin_idx], plt.ylim()[1]*0.8, "ps1", color="C0", fontsize=12)
plt.axvline(n_ps2[bin_idx], color="C0", alpha=0.75)
plt.text(n_ps2[bin_idx], plt.ylim()[1]*0.8, "ps2", color="C0", fontsize=12)
plt.axvline(n_ps3[bin_idx], color="C0", alpha=0.5)
plt.text(n_ps3[bin_idx], plt.ylim()[1]*0.8, "ps3", color="C0", fontsize=12)
plt.axvline(n_ms1[bin_idx], color="C0", linestyle="--")
plt.text(n_ms1[bin_idx], plt.ylim()[1]*0.8, "ms1", color="C0", fontsize=12)
plt.axvline(n_ms2[bin_idx], color="C0", linestyle="--", alpha=0.75)
plt.text(n_ms2[bin_idx], plt.ylim()[1]*0.8, "ms2", color="C0", fontsize=12)
plt.axvline(n_ms3[bin_idx], color="C0", linestyle="--", alpha=0.5)
plt.text(n_ms3[bin_idx], plt.ylim()[1]*0.8, "ms3", color="C0", fontsize=12)
plt.xlabel(var_config.var_labels[0] + ", Bin " + str(bin_idx))
plt.ylabel("Events / Bin")

plt.title(knob)
save_fig_name = "{}/{}-{}-geniemultisigma_univ_events.pdf".format(save_fig_dir, var_config.var_save_name, syst_name)
if save_fig:
    plt.savefig(save_fig_name, bbox_inches='tight')
plt.show()

# %% [markdown]
# #### Add cov from all multisigma knobs

# %%
n_univ_genie_multisigma = 2

cov_type = "event"
GENIEmultisigma_Frac_Covariance_event = np.zeros((len(var_config.bins)-1, len(var_config.bins)-1))
plot_labels = [var_config.var_labels[1], "Flux Integrated Event Rate"]
for kidx, knob in enumerate(regen_systematics_sbnd_multisigma):
    syst_name = knob
    save_fig_name = "{}/{}-{}-geniemultisigma_univ_{}_rates".format(save_fig_dir, var_config.var_save_name, syst_name, cov_type)
    ret_geniemultisigma = get_covariance(cov_type, syst_name, n_univ_genie_multisigma, 
                            nevts_signal_sel_reco, var_signal_sel_truth, var_signal_sel_reco, var_config.bins, 
                            plot_labels, save_fig=save_fig, save_fig_name=save_fig_name)
    GENIEmultisigma_Frac_Covariance_event = GENIEmultisigma_Frac_Covariance_event + ret_geniemultisigma["Covariance_Frac"]

cov_type = "xsec"
GENIEmultisigma_Frac_Covariance = np.zeros((len(var_config.bins)-1, len(var_config.bins)-1))
plot_labels = [var_config.var_labels[1], "xsec"]
for kidx, knob in enumerate(regen_systematics_sbnd_multisigma):
    syst_name = knob
    save_fig_name = "{}/{}-{}-geniemultisigma_univ_{}_rates".format(save_fig_dir, var_config.var_save_name, syst_name, cov_type)
    ret_geniemultisigma = get_covariance(cov_type, syst_name, n_univ_genie_multisigma, 
                            nevts_signal_sel_reco, var_signal_sel_truth, var_signal_sel_reco, var_config.bins, 
                            plot_labels, save_fig=save_fig, save_fig_name=save_fig_name)
    GENIEmultisigma_Frac_Covariance = GENIEmultisigma_Frac_Covariance + ret_geniemultisigma["Covariance_Frac"]

# %%
# plot the covariance for xsec 
save_fig_name = "{}/{}-geniemultisigma_covariance_frac".format(save_fig_dir, var_config.var_save_name)
plot_heatmap(GENIEmultisigma_Frac_Covariance, "GENIE multisigma Fractional Covariance",
             save_fig=save_fig, save_fig_name=save_fig_name)

GENIEmultisigma_correlation = np.zeros_like(GENIEmultisigma_Frac_Covariance)
for i in range(GENIEmultisigma_Frac_Covariance.shape[0]):
    for j in range(GENIEmultisigma_Frac_Covariance.shape[1]):
        GENIEmultisigma_correlation[i, j] = GENIEmultisigma_Frac_Covariance[i, j] / np.sqrt(GENIEmultisigma_Frac_Covariance[i, i] * GENIEmultisigma_Frac_Covariance[j, j])
save_fig_name = "{}/{}-geniemultisigma_correlation".format(save_fig_dir, var_config.var_save_name)
plot_heatmap(GENIEmultisigma_correlation, "GENIE multisigma Correlation",
             save_fig=save_fig, save_fig_name=save_fig_name)

# %% [markdown]
# ### Add all GENIE cov

# %%
GENIEtotal_Frac_Covariance = GENIEmultisim_Frac_Covariance + GENIEmultisigma_Frac_Covariance + GENIEmorph_Frac_Covariance

save_fig_name = "{}/{}-genie_total_covariance_frac".format(save_fig_dir, var_config.var_save_name)
plot_heatmap(GENIEtotal_Frac_Covariance, "GENIE Fractional Covariance",
             save_fig=save_fig, save_fig_name=save_fig_name)

GENIEtotal_correlation = np.zeros_like(GENIEtotal_Frac_Covariance)
for i in range(GENIEtotal_Frac_Covariance.shape[0]):
    for j in range(GENIEtotal_Frac_Covariance.shape[1]):
        GENIEtotal_correlation[i, j] = GENIEtotal_Frac_Covariance[i, j] / np.sqrt(GENIEtotal_Frac_Covariance[i, i] * GENIEtotal_Frac_Covariance[j, j])
save_fig_name = "{}/{}-genie_total_correlation".format(save_fig_dir, var_config.var_save_name)
plot_heatmap(GENIEtotal_correlation, "GENIE Total Correlation",
             save_fig=save_fig, save_fig_name=save_fig_name)

# %% [markdown]
# ## Total

# %%
# add fractional covariance of all systs, then multiply by the CV value to get the covariance
Total_Covariance_Frac = ret_mcstat["Covariance_Frac"] + ret_flux["Covariance_Frac"] + GENIEtotal_Frac_Covariance + Cosmics_Covariance_Frac + Detector_Covariance_Frac

save_fig_name = "{}/{}-total_covariance_frac".format(save_fig_dir, var_config.var_save_name)
plot_heatmap(Total_Covariance_Frac, "Total Fractional Covariance",
             save_fig=save_fig, save_fig_name=save_fig_name)

# covariance in xsec scale
Total_Covariance = np.zeros_like(Total_Covariance_Frac)
for i in range(len(var_config.bins)-1):
    for j in range(len(var_config.bins)-1):
        Total_Covariance[i, j] = Total_Covariance_Frac[i, j] * (nevts_signal_sel_reco[i] * nevts_signal_sel_reco[j]) * XSEC_UNIT**2

save_fig_name = "{}/{}-total_covariance".format(save_fig_dir, var_config.var_save_name)
plot_heatmap(Total_Covariance, "Total Covariance",
             save_fig=save_fig, save_fig_name=save_fig_name)

# correlation
Total_correlation = np.zeros_like(Total_Covariance_Frac)
for i in range(Total_Covariance_Frac.shape[0]):
    for j in range(Total_Covariance_Frac.shape[1]):
        Total_correlation[i, j] = Total_Covariance_Frac[i, j] / np.sqrt(Total_Covariance_Frac[i, i] * Total_Covariance_Frac[j, j])
save_fig_name = "{}/{}-total_correlation".format(save_fig_dir, var_config.var_save_name)
plot_heatmap(Total_correlation, "Total Correlation",
             save_fig=save_fig, save_fig_name=save_fig_name)

# %%
save_fig = True

# %%
# fractional uncertainty
frac_uncert_mcstat   = np.sqrt(np.diag(ret_mcstat["Covariance_Frac"]))
frac_uncert_flux     = np.sqrt(np.diag(ret_flux["Covariance_Frac"]))
frac_uncert_genie    = np.sqrt(np.diag(GENIEtotal_Frac_Covariance))
frac_uncert_cosmics  = np.sqrt(np.diag(Cosmics_Covariance_Frac))
frac_uncert_detector = np.sqrt(np.diag(Detector_Covariance_Frac))
frac_uncert_total    = np.sqrt(np.diag(Total_Covariance_Frac))
# 2% flat uncertainty for POT
frac_uncert_pot = np.zeros_like(frac_uncert_mcstat)
frac_uncert_pot[:] = 0.02
# 1% for ntargets
frac_uncert_ntargets = np.zeros_like(frac_uncert_mcstat)
frac_uncert_ntargets[:] = 0.01
# add to total
frac_uncert_total = np.sqrt(frac_uncert_total**2 + frac_uncert_pot**2 + frac_uncert_ntargets**2)

plt.hist(var_config.bin_centers, bins=var_config.bins, weights=frac_uncert_mcstat * 1e2,   histtype="step", color="C0", label="MCstat")
plt.hist(var_config.bin_centers, bins=var_config.bins, weights=frac_uncert_flux * 1e2,     histtype="step", color="C1", label="Flux")
plt.hist(var_config.bin_centers, bins=var_config.bins, weights=frac_uncert_genie * 1e2,    histtype="step", color="C2", label="GENIE")
plt.hist(var_config.bin_centers, bins=var_config.bins, weights=frac_uncert_cosmics * 1e2,  histtype="step", color="C3", label="Cosmics")
plt.hist(var_config.bin_centers, bins=var_config.bins, weights=frac_uncert_detector * 1e2, histtype="step", color="C4", label="Detector")
plt.hist(var_config.bin_centers, bins=var_config.bins, weights=frac_uncert_pot * 1e2,    histtype="step", color="C5", label="POT")
plt.hist(var_config.bin_centers, bins=var_config.bins, weights=frac_uncert_ntargets * 1e2, histtype="step", color="C6", label="NTargets")
plt.hist(var_config.bin_centers, bins=var_config.bins, weights=frac_uncert_total * 1e2,    histtype="step", color="k",  label="Total")

plt.xlim(var_config.bins[0], var_config.bins[-1])
plt.xlabel(var_config.var_labels[1])
plt.ylabel("Uncertainty [%]")
plt.legend(fontsize=12, ncol=3)
plt.ylim(0, max(frac_uncert_total*1e2) * 1.35)

if save_fig:
    plt.savefig("{}/uncertainty_breakdown-{}.pdf".format(save_fig_dir, var_config.var_save_name), bbox_inches='tight')
plt.show()


# %% [markdown]
# Singal distribution with error bars from diagonal components of covariance matrix

# %%
# Compute bin centers for error bars
plt.errorbar(bin_centers, nevts_signal_sel_reco*XSEC_UNIT, yerr=frac_uncert_total*nevts_signal_sel_reco*XSEC_UNIT, 
             fmt='o', color='black', label='Subtracted (syst. error)', capsize=3)
plt.xlim(var_config.bins[0], var_config.bins[-1])
plt.xlabel(var_config.var_labels[0])
plt.ylabel("Events / Bin")
plt.legend()

if save_fig:
    plt.savefig("{}/{}-bkg_subtracted_event_rates.pdf".format(save_fig_dir, var_config.var_save_name), bbox_inches='tight')
plt.show()

# %% [markdown]
# # Unfolding

# %%
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

    # leave whitespace at the top for the legend
    plt.xlabel(plot_labels[0])
    plt.ylabel(plot_labels[1])
    plt.xlim(bins[0], bins[-1])
    plt.ylim(0., 1.35 * fake_data.max())

    if save_fig:
        plt.savefig(save_name, bbox_inches='tight', dpi=300)
    plt.show()
    
    return fake_data, background_cv

# %%
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

    # leave whitespace at the top for the legend
    plt.xlabel(plot_labels[0])
    plt.ylabel(plot_labels[1])
    plt.xlim(bins[0], bins[-1])
    plt.ylim(0., 1.35 * fake_data.max())

    if save_fig:
        plt.savefig(save_name, bbox_inches='tight', dpi=300)
    plt.show()

# %%
# --- config for Wiener-SVD unfolding ---
C_type = 2
Norm_type = 0.5

# %% [markdown]
# ## Closure test 
# - use MC signal as fake data

# %%
def chi2(data, model, cov):
    return (data - model) @ np.linalg.inv(cov) @ (data - model)

# %%
def plot_unfolded_result(unfold, bins, measured, models,
                         plot_labels, model_names, 
                         save_fig=False, save_name=None,
                         closure_test=False):

    # need to divide by bin width for differential xsec
    bin_widths = np.diff(bins)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    fig, ax = plt.subplots(figsize=(8.5, 6))

    # --- stat uncertainties
    UnfoldCov_stat = unfold['StatUnfoldCov']
    Unfold_uncert_stat = np.diag(UnfoldCov_stat)

    # --- syst uncertainties
    UnfoldCov_syst = unfold['SystUnfoldCov']
    Unfold_uncert_syst = np.diag(UnfoldCov_syst)

    # decomose into norm and shape components
    # TODO: the first item in models is the true model
    SystUnfoldCov_norm, SystUnfoldCov_shape = Matrix_Decomp(models[0], UnfoldCov_syst)
    Unfold_uncert_norm = np.sqrt(np.diag(SystUnfoldCov_norm))
    Unfold_uncert_shape = np.sqrt(np.diag(SystUnfoldCov_shape))

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
    model_handles = []
    model_labels = []
    for midx, model in enumerate(models):

        model_smeared = unfold['AddSmear'] @ model
        chi2_val = chi2(Unfolded, model_smeared, UnfoldCov_syst)
        chi2_vals.append(chi2_val)

        model_smeared_perwidth = model_smeared / bin_widths
        model_handle, = plt.step(bins, np.append(model_smeared_perwidth, model_smeared_perwidth[-1]), where='post')
        model_handles.append(model_handle)
        model_labels.append(f'$A_c \\otimes$ {model_names[midx]} ($\chi^2$ = {chi2_vals[midx]:.2f}/{len(bins)-1})')

    # legend
    if closure_test:
        handles = [bar_handle, reco_handle] + model_handles
        labels = ['Unfolded', 'Measured Signal (Input)'] + model_labels
    else:
        handles = [bar_handle, norm_handle, reco_handle] + model_handles
        labels = ['Unfolded', 'Norm. Syst. Unc.', 'Measured Signal (Input)'] + model_labels
    plt.legend(handles, labels, 
               loc='upper left', fontsize=10, frameon=False, ncol=2, bbox_to_anchor=(0.02, 0.98))

    # leave space for legend
    plt.xlabel(plot_labels[0])
    plt.ylabel(plot_labels[1])
    plt.xlim(bins[0], bins[-1])
    plt.ylim(0., np.max(Unfolded_perwidth)*1.4)

    if save_fig:
        plt.savefig(save_name, bbox_inches='tight')
    plt.show()

# %%
var_total = var_total_mc
weights_total = weights_total_mc
bins = var_config.bins
plot_labels = [var_config.var_labels[1], "Events / Bin"]

var_categ = var_per_nuint_categ_mc
weights_categ = weights_per_categ
colors = mode_colors
labels = mode_labels
save_name = "{}/{}-closure_test_topology_breakdown.pdf".format(save_fig_dir, var_config.var_save_name)
fake_data, background_cv = plot_topology_breakdown(var_categ, weights_categ, var_total, weights_total, bins,
                            plot_labels, 
                            colors=colors, labels=labels, 
                            save_fig=save_fig, save_name=save_name)

var_categ = var_per_genie_mode_mc
weights_categ = weights_per_genie_mode
colors = genie_mode_colors
labels = genie_mode_labels
save_name = "{}/{}-closure_test_genie_mode_breakdown.pdf".format(save_fig_dir, var_config.var_save_name)
plot_genie_breakdown(var_categ, weights_categ, var_total, weights_total, bins,
                     plot_labels, 
                     colors=colors, labels=labels,
                     save_fig=save_fig, save_name=save_name)

# %%
Measured = nevts_signal_sel_reco * XSEC_UNIT 
Model = nevts_signal_truth * XSEC_UNIT
Covariance = ret_mcstat["Covariance"] # closure test -- just use MC stat uncertainty
unfold = WienerSVD(Response, Model, Measured, Covariance, C_type, Norm_type)
print(unfold.keys())
# decomp_cov = Matrix_Decomp(Model, unfold['SystUnfoldCov'])

# %%
measured = Measured
bins = var_config.bins
models = [Model]
model_names = ["AR23"]
plot_labels = [var_config.var_labels[0], var_config.xsec_label]
save_name = "{}/{}-closure_test.pdf".format(save_fig_dir, var_config.var_save_name)
plot_unfolded_result(unfold, bins, measured, models, 
                      plot_labels, model_names, 
                      save_fig=save_fig, save_name=save_name,
                      closure_test=True)

# %%
save_fig_name = "{}/{}-{}-add_smear".format(save_fig_dir, var_config.var_save_name, "closure_test")
plot_labels = [var_config.var_labels[2], var_config.var_labels[1]]
plot_heatmap(unfold["AddSmear"], "$A_c$", 
             plot_labels=plot_labels,
             save_fig=save_fig, save_fig_name=save_fig_name)

# %% [markdown]
# ## Fake Data Tests
# 
# - use alternate MC as fake data

# %%
# use only stat unc and xsec unc for fake data tests
Covariance_Frac = ret_mcstat["Covariance_Frac"] + GENIEtotal_Frac_Covariance
Covariance = np.zeros_like(Covariance_Frac)
for i in range(len(var_config.bins)-1):
    for j in range(len(var_config.bins)-1):
        Covariance[i, j] = Covariance_Frac[i, j] * (nevts_signal_sel_reco[i] * nevts_signal_sel_reco[j]) * XSEC_UNIT**2


# %% [markdown]
# ### MEC scale

# %%
test_name = "mec_test"

mec_scale = 2
weights_fake_data = np.ones(len(var_total_mc))
weights_fake_data[mc_evt_df.genie_mode == 10] *= mec_scale

# %%
var_total = var_total_mc
weights_total = weights_total_mc * weights_fake_data
bins = var_config.bins
plot_labels = [var_config.var_labels[1], "Events / Bin"]

var_categ = var_per_nuint_categ_mc
weights_categ = weights_per_categ
colors = mode_colors
labels = mode_labels
save_name = "{}/{}-{}-topology_breakdown.pdf".format(save_fig_dir, test_name, var_config.var_save_name)
fake_data, background_cv = plot_topology_breakdown(var_categ, weights_categ, var_total, weights_total, bins,
                            plot_labels, 
                            colors=colors, labels=labels, 
                            save_fig=save_fig, save_name=save_name)

var_categ = var_per_genie_mode_mc
weights_categ = weights_per_genie_mode
colors = genie_mode_colors
labels = genie_mode_labels
save_name = "{}/{}-{}-genie_mode_breakdown.pdf".format(save_fig_dir, test_name, var_config.var_save_name)
plot_genie_breakdown(var_categ, weights_categ, var_total, weights_total, bins,
                     plot_labels, 
                     colors=colors, labels=labels,
                     save_fig=save_fig, save_name=save_name)

# %%
Measured = fake_data * XSEC_UNIT - background_cv * XSEC_UNIT
Model = nevts_signal_truth * XSEC_UNIT
unfold = WienerSVD(Response, Model, Measured, Covariance, C_type, Norm_type)

# %%
weight_fakedata_signal_truth = weight_truth_signal.copy()

# nom MC
nevts_fakedata_signal_truth, _, _ = plt.hist(var_truth_signal, bins=var_config.bins, weights=weight_fakedata_signal_truth, 
                                             histtype="step", label="nominal MC")

# alt MC --> fake data
weight_fakedata_signal_truth[mc_nu_df[mc_nu_df.nuint_categ == 1].genie_mode == 10] = mec_scale
nevts_fakedata_signal_truth, _, _ = plt.hist(var_truth_signal, bins=var_config.bins, weights=weight_fakedata_signal_truth, 
                                             histtype="step", label="alternative MC")
plt.legend()
plt.show()

# %%
unfold = unfold
measured = Measured
bins = var_config.bins
models = [Model,  nevts_fakedata_signal_truth*XSEC_UNIT]
model_names = ["AR23", "Fake Data ($\\times${} MEC)".format(mec_scale)]
plot_labels = [var_config.var_labels[0], var_config.xsec_label]
save_name = "{}/{}-{}-unfolded_event_rates.pdf".format(save_fig_dir, test_name, var_config.var_save_name)
plot_unfolded_result(unfold, bins, measured, models, 
                     plot_labels, model_names, save_fig=save_fig, save_name=save_name)

# %%
save_fig_name = "{}/{}-{}-add_smear".format(save_fig_dir, var_config.var_save_name, syst_name)
plot_labels = [var_config.var_labels[2], var_config.var_labels[1]]
plot_heatmap(unfold["AddSmear"], "$A_c$", plot_labels=plot_labels,
             save_fig=save_fig, save_fig_name=save_fig_name)                                                    

# %% [markdown]
# ### QE scale

# %%
test_name = "qe_test"

qe_scale = 1.2
weights_fake_data = np.ones(len(var_total_mc))
weights_fake_data[mc_evt_df.genie_mode == 0] *= qe_scale

# %%
var_total = var_total_mc
weights_total = weights_total_mc * weights_fake_data
bins = var_config.bins
plot_labels = [var_config.var_labels[1], "Events / Bin"]

var_categ = var_per_nuint_categ_mc
weights_categ = weights_per_categ
colors = mode_colors
labels = mode_labels
save_name = "{}/{}-{}-topology_breakdown.pdf".format(save_fig_dir, test_name, var_config.var_save_name)
fake_data, background_cv = plot_topology_breakdown(var_categ, weights_categ, var_total, weights_total, bins,
                            plot_labels, 
                            colors=colors, labels=labels, 
                            save_fig=save_fig, save_name=save_name)

var_categ = var_per_genie_mode_mc
weights_categ = weights_per_genie_mode
colors = genie_mode_colors
labels = genie_mode_labels
save_name = "{}/{}-{}-genie_mode_breakdown.pdf".format(save_fig_dir, test_name, var_config.var_save_name)
plot_genie_breakdown(var_categ, weights_categ, var_total, weights_total, bins,
                     plot_labels, 
                     colors=colors, labels=labels,
                     save_fig=save_fig, save_name=save_name)

# %%
Measured = fake_data * XSEC_UNIT - background_cv * XSEC_UNIT
Model = nevts_signal_truth * XSEC_UNIT
unfold = WienerSVD(Response, Model, Measured, Covariance, C_type, Norm_type)

# %%
weight_fakedata_signal_truth = weight_truth_signal.copy()

# nom MC
nevts_fakedata_signal_truth, _, _ = plt.hist(var_truth_signal, bins=var_config.bins, weights=weight_fakedata_signal_truth, 
                                             histtype="step", label="nominal MC")

# alt MC --> fake data
weight_fakedata_signal_truth[mc_nu_df[mc_nu_df.nuint_categ == 1].genie_mode == 0] = qe_scale
nevts_fakedata_signal_truth, _, _ = plt.hist(var_truth_signal, bins=var_config.bins, weights=weight_fakedata_signal_truth, 
                                             histtype="step", label="alternative MC")
plt.legend()
plt.show()

# %%
unfold = unfold
measured = Measured
bins = var_config.bins
models = [Model,  nevts_fakedata_signal_truth*XSEC_UNIT]
model_names = ["AR23", "Fake Data ($\\times${} QE)".format(qe_scale)]
plot_labels = [var_config.var_labels[0], var_config.xsec_label]
save_name = "{}/{}-{}-unfolded_event_rates.pdf".format(save_fig_dir, test_name, var_config.var_save_name)
plot_unfolded_result(unfold, bins, measured, models, 
                     plot_labels, model_names, 
                     save_fig=save_fig, save_name=save_name)

# %%
save_fig_name = "{}/{}-{}-add_smear".format(save_fig_dir, var_config.var_save_name, test_name)
plot_labels = [var_config.var_labels[2], var_config.var_labels[1]]
plot_heatmap(unfold["AddSmear"], "$A_c$", plot_labels=plot_labels,
             save_fig=save_fig, save_fig_name=save_fig_name)

# %% [markdown]
# ### Np background scale

# %%
Np_scale = 2
weights_fake_data = np.ones(len(var_total_mc))
weights_fake_data[mc_evt_df.nuint_categ == 2] = Np_scale
test_name = "np_test"

# %%
var_categ = var_per_nuint_categ_mc
weights_categ = weights_per_categ
colors = mode_colors
labels = mode_labels
var_total = var_total_mc
weights_total = weights_total_mc * weights_fake_data
bins = var_config.bins

plot_labels = [var_config.var_labels[1], "Events / Bin"]
save_name = "{}/{}-{}-topology_breakdown.pdf".format(save_fig_dir, test_name, var_config.var_save_name)
fake_data, background_cv = plot_topology_breakdown(var_categ, weights_categ, var_total, weights_total, bins,
                            plot_labels, 
                            colors=colors, labels=labels, 
                            save_fig=save_fig, save_name=save_name)

var_categ = var_per_genie_mode_mc[::-1]
weights_categ = weights_per_genie_mode[::-1]
colors = genie_mode_colors[::-1]
labels = genie_mode_labels[::-1]
save_name = "{}/{}-{}-genie_mode_breakdown.pdf".format(save_fig_dir, test_name, var_config.var_save_name)
plot_genie_breakdown(var_categ, weights_categ, var_total, weights_total, bins,
                     plot_labels, 
                     colors=colors, labels=labels,
                     save_fig=save_fig, save_name=save_name)

# %%
Measured = fake_data * XSEC_UNIT - background_cv * XSEC_UNIT
Model = nevts_signal_truth * XSEC_UNIT
unfold = WienerSVD(Response, Model, Measured, Covariance, C_type, Norm_type)

# %%
weight_fakedata_signal_truth = weight_truth_signal.copy()

# nom MC
nevts_fakedata_signal_truth, _, _ = plt.hist(var_truth_signal, bins=var_config.bins, weights=weight_fakedata_signal_truth, 
                                             histtype="step", label="nominal MC")

# alt MC --> fake data
weight_fakedata_signal_truth[mc_nu_df[mc_nu_df.nuint_categ == 1].nuint_categ == 2] = Np_scale
nevts_fakedata_signal_truth, _, _ = plt.hist(var_truth_signal, bins=var_config.bins, weights=weight_fakedata_signal_truth, 
                                             histtype="step", label="alternative MC")
plt.legend()
plt.show()

# %%
unfold = unfold
measured = Measured
bins = var_config.bins
models = [Model,  nevts_fakedata_signal_truth*XSEC_UNIT]
model_names = ["AR23", "Fake Data ($\\times${} Np)".format(Np_scale)]
plot_labels = [var_config.var_labels[0], var_config.xsec_label]
save_name = "{}/{}-{}-unfolded_event_rates.pdf".format(save_fig_dir, test_name, var_config.var_save_name)
plot_unfolded_result(unfold, bins, measured, models, plot_labels, model_names, save_fig=save_fig, save_name=save_name)

# %%
save_fig_name = "{}/{}-{}-add_smear".format(save_fig_dir, var_config.var_save_name, test_name)
plot_labels = [var_config.var_labels[2], var_config.var_labels[1]]
plot_heatmap(unfold["AddSmear"], "$A_c$", plot_labels=plot_labels,
             save_fig=save_fig, save_fig_name=save_fig_name)

# %% [markdown]
# ### Signal scale

# %%
sig_scale = 1.2
weights_fake_data = np.ones(len(var_total_mc))
weights_fake_data[mc_evt_df.nuint_categ == 1] = sig_scale
test_name = "sig_test"

# %%
var_categ = var_per_nuint_categ_mc
weights_categ = weights_per_categ
colors = mode_colors
labels = mode_labels
var_total = var_total_mc
weights_total = weights_total_mc * weights_fake_data
bins = var_config.bins
plot_labels = [var_config.var_labels[1], "Events / Bin"]
save_name = "{}/{}-{}-topology_breakdown.pdf".format(save_fig_dir, test_name, var_config.var_save_name)
fake_data, background_cv = plot_topology_breakdown(var_categ, weights_categ, var_total, weights_total, bins,
                            plot_labels, 
                            colors=colors, labels=labels, 
                            save_fig=save_fig, save_name=save_name)

var_categ = var_per_genie_mode_mc
weights_categ = weights_per_genie_mode
colors = genie_mode_colors
labels = genie_mode_labels
save_name = "{}/{}-{}-genie_mode_breakdown.pdf".format(save_fig_dir, test_name, var_config.var_save_name)
plot_genie_breakdown(var_categ, weights_categ, var_total, weights_total, bins,
                     plot_labels, 
                     colors=colors, labels=labels,
                     save_fig=save_fig, save_name=save_name)

# %%
Measured = fake_data * XSEC_UNIT - background_cv * XSEC_UNIT
Model = nevts_signal_truth * XSEC_UNIT
unfold = WienerSVD(Response, Model, Measured, Covariance, C_type, Norm_type)

# %%
weight_fakedata_signal_truth = weight_truth_signal.copy()

# nom MC
nevts_fakedata_signal_truth, _, _ = plt.hist(var_truth_signal, bins=var_config.bins, weights=weight_fakedata_signal_truth, 
                                             histtype="step", label="nominal MC")

# alt MC --> fake data
weight_fakedata_signal_truth[mc_nu_df[mc_nu_df.nuint_categ == 1].nuint_categ == 1] = sig_scale
nevts_fakedata_signal_truth, _, _ = plt.hist(var_truth_signal, bins=var_config.bins, weights=weight_fakedata_signal_truth, 
                                             histtype="step", label="alternative MC")
plt.legend()
plt.show()

# %%
unfold = unfold
measured = Measured
bins = var_config.bins
models = [Model,  nevts_fakedata_signal_truth*XSEC_UNIT]
model_names = ["AR23", "Fake Data ($\\times${} sig)".format(sig_scale)]
plot_labels = [var_config.var_labels[0], var_config.xsec_label]
save_name = "{}/{}-{}-unfolded_event_rates.pdf".format(save_fig_dir, test_name, var_config.var_save_name)
plot_unfolded_result(unfold, bins, measured, models, plot_labels, model_names, save_fig=save_fig, save_name=save_name)

# %%
save_fig_name = "{}/{}-{}-add_smear".format(save_fig_dir, var_config.var_save_name, test_name)
plot_labels = [var_config.var_labels[2], var_config.var_labels[1]]
plot_heatmap(unfold["AddSmear"], "$A_c$", plot_labels=plot_labels,
             save_fig=save_fig, save_fig_name=save_fig_name)

# %%



