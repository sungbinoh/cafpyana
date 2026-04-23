import numpy as np
import inspect

class VariableConfig:
    """
    A configurable class for setting up unfolding variable configurations.
    Choose a configuration using one of the provided class methods,
    or instantiate directly with custom parameters.
    """
    def __init__(self, var_save_name, var_plot_name, var_labels, bins, var_evt_reco_col, var_evt_truth_col, var_nu_col, xsec_label):
        self.var_save_name = var_save_name
        self.var_plot_name = var_plot_name
        self.var_labels = var_labels
        self.bins = bins
        self.bin_centers = (bins[:-1] + bins[1:]) / 2.
        self.var_evt_reco_col = var_evt_reco_col
        self.var_evt_truth_col = var_evt_truth_col
        self.var_nu_col = var_nu_col
        self.xsec_label = xsec_label


    # for a integrated single-bin measurement of all events
    @classmethod
    def all_events(cls):
        return cls(
            var_save_name="all-events",
            var_plot_name="All Events",
            var_labels=[r"All Events", 
            r"All Events", ""],
            # use slc.producer as the dummy variable
            bins=np.array([-1, 1]),
            var_evt_reco_col=('slc', 'producer', '', '', '', '', ''),
            var_evt_truth_col=('slc', 'producer', '', '', '', '', ''),
            var_nu_col=('slc', 'producer', ''),
            xsec_label=r"$\frac{d\sigma}{dAll Events}$ ($\mathrm{cm}^2$)"
        )

    ## -- An example for P_mu
    #@classmethod
    #def muon_momentum(cls):
    #    return cls(
    #        var_save_name="muon-p",
    #        var_plot_name="P_\mu",
    #        var_labels=[r"$\mathrm{P_\mu}$ (GeV/c)", 
    #        r"$\mathrm{P_\mu^{reco.}}$ (GeV/c)", 
    #        r"$\mathrm{P_\mu^{true}}$ (GeV/c)"],
    #        bins=np.array([0.22, 0.27, 0.32, 0.37, 0.42, 0.47, 0.52, 0.57, 0.62, 0.67, 0.72, 0.77, 0.82, 0.9, 1.0]),
    #        var_evt_reco_col=('mu', 'pfp', 'trk', 'P', 'p_muon', '', ''),
    #        var_evt_truth_col=('mu', 'pfp', 'trk', 'truth', 'p', 'totp', ''),
    #        var_nu_col=('mc', 'mu', 'totp'),
    #        xsec_label=r"$\frac{d\sigma}{dP_\mu}$ ($\mathrm{cm}^2$ / GeV/c)"
    #    )
