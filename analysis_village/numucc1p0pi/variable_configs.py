import numpy as np

class VariableConfig:
    """
    A configurable class for setting up unfolding variable configurations.
    Choose a configuration using one of the provided class methods,
    or instantiate directly with custom parameters.
    """
    def __init__(self, var_save_name, var_plot_name, var_unit, bins, var_evt_reco_col, var_evt_truth_col, var_nu_col):
        self.var_save_name = var_save_name
        self.var_plot_name = var_plot_name
        self.var_unit = var_unit
        unit_suffix = f"~[{var_unit}]" if len(var_unit) > 0 else ""
        self.var_labels = [r"$\mathrm{" + var_plot_name + unit_suffix + "}$", 
                           r"$\mathrm{" + var_plot_name + "^{reco.}" + unit_suffix + "}$", 
                           r"$\mathrm{" + var_plot_name + "^{true}" + unit_suffix + "}$"]
        self.bins = bins
        self.bin_centers = (bins[:-1] + bins[1:]) / 2.
        self.var_evt_reco_col = var_evt_reco_col
        self.var_evt_truth_col = var_evt_truth_col
        self.var_nu_col = var_nu_col


    @classmethod
    def muon_momentum(cls):
        return cls(
            var_save_name="muon-p",
            var_plot_name="P_\mu",
            var_unit="GeV/c",
            bins=np.linspace(0.15, 1.2, 11),
            var_evt_reco_col=('mu', 'pfp', 'trk', 'P', 'p_muon', '', '', ''),
            var_evt_truth_col=('mu', 'pfp', 'trk', 'truth', 'p', 'totp', '', ''),
            var_nu_col=('mu', 'totp', '')
        )

    @classmethod
    def muon_direction(cls):
        return cls(
            var_save_name="muon-dir_z",
            var_plot_name="cos(\theta_\mu)",
            var_unit="",
            bins=np.linspace(-1, 1, 6),
            var_evt_reco_col=('mu', 'pfp', 'trk', 'dir', 'z', '', '', ''),
            var_evt_truth_col=('mu', 'pfp', 'trk', 'truth', 'p', 'dir', 'z', ''),
            var_nu_col=('mu', 'dir', 'z')
        )

    @classmethod
    def proton_momentum(cls):
        return cls(
            var_save_name="proton-p",
            var_plot_name="P_p",
            var_unit="GeV/c",
            bins=np.linspace(0.2, 2, 6),
            var_evt_reco_col=('p', 'pfp', 'trk', 'P', 'p_proton', '', '', ''),
            var_evt_truth_col=('p', 'pfp', 'trk', 'truth', 'p', 'totp', '', ''),
            var_nu_col=('p', 'totp', '')
        )

    @classmethod
    def proton_direction(cls):
        return cls(
            var_save_name="proton-dir_z",
            var_plot_name="cos(\theta_p)",
            var_unit="",
            bins=np.linspace(-1, 1, 6),
            var_evt_reco_col=('p', 'pfp', 'trk', 'dir', 'z', '', '', ''),
            var_evt_truth_col=('p', 'pfp', 'trk', 'truth', 'p', 'dir', 'z', ''),
            var_nu_col=('p', 'dir', 'z')
        )