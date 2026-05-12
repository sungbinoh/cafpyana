import numpy as np
from analysis_village.unfolding.variable_configs import VariableConfig

varcfg_p_mu = VariableConfig(
    var_save_name="muon-p",
    var_plot_name="P_\mu",
    var_labels=[r"$\mathrm{P_\mu}$ (GeV/c)", 
    r"$\mathrm{P_\mu^{reco.}}$ (GeV/c)", 
    r"$\mathrm{P_\mu^{true}}$ (GeV/c)"],
    bins=np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]),
    var_evt_reco_col='range_p_mu',
    var_evt_truth_col='true_p_mu',
    var_nu_col='true_p_mu',
    xsec_label=r"$\frac{d\sigma}{dP_\mu}$ ($\mathrm{cm}^2$ / GeV/c)"
)
