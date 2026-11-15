import numpy as np
from analysis_village.unfolding.variable_configs import VariableConfig

varcfg_p_mu = VariableConfig(
    var_save_name="muon-p",
    var_plot_name="P_\mu",
    var_labels=[r"$\mathrm{P_\mu}$ (GeV/c)", 
    r"$\mathrm{P_\mu^{reco.}}$ (GeV/c)", 
    r"$\mathrm{P_\mu^{true}}$ (GeV/c)"],
    bins=np.array([0., 0.399, 0.505, 0.657, 0.764, 0.906, 6.]),
    var_evt_reco_col='range_p_mu',
    var_evt_truth_col='true_p_mu',
    var_nu_col='true_p_mu',
    xsec_label=r"$\frac{d\sigma}{dP_\mu}$ ($\mathrm{cm}^2$ / GeV/c)"
)

varcfg_cos_mu = VariableConfig(
    var_save_name="muon-cos",
    var_plot_name="\cos\theta_{\mu^\mp}",
    var_labels=[r"$\cos\theta_{\mu^\mp}$",
    r"$\cos\theta_{\mu^\mp}^{reco.}$", 
    r"$\cos\theta_{\mu^\mp}^{true}$"],
    bins=np.array([-1., 0.929, 0.951, 0.974, 0.986, 0.994, 0.997, 1.]),
    var_evt_reco_col='long_dirz',
    var_evt_truth_col='true_cos_theta_mu',
    var_nu_col='true_cos_theta_mu',
    xsec_label=r"$\frac{d\sigma}{d\cos\theta_{\mu^\mp}}$ ($\mathrm{cm}^2$)"
)

varcfg_cos_pi = VariableConfig(
    var_save_name="pion-cos",
    var_plot_name="\cos\theta_{\pi^\pm}",
    var_labels=[r"$\cos\theta_{\pi^\pm}$", 
    r"$\cos\theta_{\pi^\pm}^{reco.}$", 
    r"$\cos\theta_{\pi^\pm}^{true}$"],
    bins=np.array([-1., 0.815, 0.863, 0.898, 0.937, 0.964, 0.977, 0.987, 1.]),
    var_evt_reco_col='short_dirz',
    var_evt_truth_col='true_cos_theta_pi',
    var_nu_col='true_cos_theta_pi',
    xsec_label=r"$\frac{d\sigma}{d\cos\theta_{\pi^\pm}}$ ($\mathrm{cm}^2$)"
)