import numpy as np
import inspect

# ===== References
# MicroBooNE tki bins: https://arxiv.org/abs/2301.03700

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


    # ==== variables for xsec measurement ====
    # for a integrated single-bin measurement of all events
    @classmethod
    def all_events(cls):
        return cls(
            var_save_name="integrated",
            var_plot_name="All Events",
            var_labels=[r"All Events", 
            r"All Events", ""],
            # use slc.producer as the dummy variable
            bins=np.linspace(0, 1, 2),
            var_evt_reco_col=('mc', 'iscc', '', '', '', '', ''),
            var_evt_truth_col=('mc', 'iscc', '', '', '', '', ''),
            var_nu_col=('mc', 'iscc', ''),
            xsec_label=r"$\frac{d\sigma}{dAll Events}$ ($\mathrm{cm}^2$)"
        )

    @classmethod
    def muon_momentum(cls):
        return cls(
            var_save_name="muon-p",
            var_plot_name="P_\mu",
            var_labels=[r"$\mathrm{P_\mu}$ (GeV/c)", 
            r"$\mathrm{P_\mu^{reco.}}$ (GeV/c)", 
            r"$\mathrm{P_\mu^{true}}$ (GeV/c)"],
            bins=np.array([0.22, 0.27, 0.32, 0.37, 0.42, 0.47, 0.52, 0.57, 0.62, 0.7, 0.8, 0.9, 1.0]),
            var_evt_reco_col=('mu', 'pfp', 'trk', 'P', 'p_muon', '', ''),
            var_evt_truth_col=('mu', 'pfp', 'trk', 'truth', 'p', 'totp', ''),
            var_nu_col=('mc', 'mu', 'totp'),
            xsec_label=r"$\frac{d\sigma}{dP_\mu}$ ($\mathrm{cm}^2$ / GeV/c)"
        )

    @classmethod
    def muon_direction(cls):
        return cls(
            var_save_name="muon-dir_z",
            var_plot_name="cos(\\theta_\mu)",
            var_labels=[r"$\mathrm{cos(\theta_\mu)}$", 
            r"$\mathrm{cos(\theta_\mu^{reco.})}$", 
            r"$\mathrm{cos(\theta_\mu^{true})}$"],
            bins=np.array([-1, -0.75, -0.6, -0.45, -0.3, -0.15, 0, 0.15, 0.3, 0.45, 0.6, 0.7, 0.8, 0.9, 1]),
            var_evt_reco_col=('mu', 'pfp', 'trk', 'dir', 'z', '', ''),
            var_evt_truth_col=('mu', 'pfp', 'trk', 'truth', 'p', 'dir', 'z'),
            var_nu_col=('mc', 'mu', 'dir', 'z'),
            xsec_label=r"$\frac{d\sigma}{dcos(\theta_\mu)}$ ($\mathrm{cm}^2$)"
        )

    @classmethod
    def proton_momentum(cls):
        return cls(
            var_save_name="proton-p",
            var_plot_name="P_p",
            var_labels=[r"$\mathrm{P_p}$ (GeV/c)", 
            r"$\mathrm{P_p^{reco.}}$ (GeV/c)", 
            r"$\mathrm{P_p^{true}}$ (GeV/c)"],
            bins=np.array([0.3, 0.37, 0.44, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.]),
            var_evt_reco_col=('p', 'pfp', 'trk', 'P', 'p_proton', '', ''),
            var_evt_truth_col=('p', 'pfp', 'trk', 'truth', 'p', 'totp', ''),
            var_nu_col=('mc', 'p', 'totp'),
            xsec_label=r"$\frac{d\sigma}{dP_p}$ ($\mathrm{cm}^2$ / GeV/c)"
        )

    @classmethod
    def proton_direction(cls):
        return cls(
            var_save_name="proton-dir_z",
            var_plot_name="cos(\\theta_p)",
            var_labels=[r"$\mathrm{cos(\theta_p)}$", 
            r"$\mathrm{cos(\theta_p^{reco.})}$", 
            r"$\mathrm{cos(\theta_p^{true})}$"],
            bins=np.array([-1,-0.7, -0.4, -0.1, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 1]),
            var_evt_reco_col=('p', 'pfp', 'trk', 'dir', 'z', '', ''),
            var_evt_truth_col=('p', 'pfp', 'trk', 'truth', 'p', 'dir', 'z'),
            var_nu_col=('mc', 'p', 'dir', 'z'),
            xsec_label=r"$\frac{d\sigma}{dcos(\theta_p)}$ ($\mathrm{cm}^2$)"
        )

    @classmethod
    def tki_del_p(cls):
        return cls(
            var_save_name="tki-del_p",
            var_plot_name="$\\delta p$",
            var_labels=[r"$\mathrm{\delta p}$ (GeV/c)", 
            r"$\mathrm{\delta p^{reco.}}$ (GeV/c)", 
            r"$\mathrm{\delta p^{true}}$ (GeV/c)"],
            bins=np.array([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.35, 0.45, 0.55, 0.65, 0.8, 1.0]),
            var_evt_reco_col=('del_p', '', '', '', '', '', ''),
            var_evt_truth_col=('mc_del_p', '', '', '', '', '', ''),
            var_nu_col=('mc', 'del_p', ''),
            xsec_label=r"$\frac{d\sigma}{d\delta p}$ ($\mathrm{cm}^2$ / GeV/c)"
        )
    
    @classmethod
    def tki_del_Tp(cls):
        return cls(
            var_save_name="tki-del_Tp",
            var_plot_name="$\\delta p_T$",
            var_labels=[r"$\mathrm{\delta p_T}$ (GeV/c)", 
            r"$\mathrm{\delta p_T^{reco.}}$ (GeV/c)", 
            r"$\mathrm{\delta p_T^{true}}$ (GeV/c)"],
            bins=np.array([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.35, 0.45, 0.55, 0.65, 0.8, 1.0]),
            var_evt_reco_col=('del_Tp', '', '', '', '', '', ''),
            var_evt_truth_col=('mc_del_Tp', '', '', '', '', '', ''),
            var_nu_col=('del_Tp', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\delta p_T}$ ($\mathrm{cm}^2$ / GeV/c)"
        )

    @classmethod
    def tki_del_Tp_x(cls):
        return cls(
            var_save_name="tki-del_Tp_x",
            var_plot_name="$\\delta p_T^x$",
            var_labels=[r"$\mathrm{\delta p_{T, x}}$ (GeV/c)", 
            r"$\mathrm{\delta p_{T, x}^{reco.}}$ (GeV/c)", 
            r"$\mathrm{\delta p_{T, x}^{true}}$ (GeV/c)"],
            bins=np.array([-0.7, -0.55, -0.45, -0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.7]),
            var_evt_reco_col=('del_Tp_x', '', '', '', '', '', ''),
            var_evt_truth_col=('mc_del_Tp_x', '', '', '', '', '', ''),
            var_nu_col=('del_Tp_x', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\delta p_T^x}$ ($\mathrm{cm}^2$ / GeV/c)"
        )

    @classmethod
    def tki_del_Tp_y(cls):
        return cls(
            var_save_name="tki-del_Tp_y",
            var_plot_name="$\\delta p_T^y$",
            var_labels=[r"$\mathrm{\delta p_{T, y}}$ (GeV/c)", 
            r"$\mathrm{\delta p_{T, y}^{reco.}}$ (GeV/c)", 
            r"$\mathrm{\delta p_{T, y}^{true}}$ (GeV/c)"],
            bins=np.array([-0.7, -0.55, -0.45, -0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.7]),
            var_evt_reco_col=('del_Tp_y', '', '', '', '', '', ''),
            var_evt_truth_col=('mc_del_Tp_y', '', '', '', '', '', ''),
            var_nu_col=('del_Tp_y', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\delta p_T^y}$ ($\mathrm{cm}^2$ / GeV/c)"
        )
    
    @classmethod
    def tki_del_alpha(cls):
        return cls(
            var_save_name="tki-del_alpha",
            var_plot_name="$\\delta \\alpha_T$",
            var_labels=[r"$\mathrm{\delta \alpha_T}$ (deg)", 
            r"$\mathrm{\delta \alpha_T^{reco.}}$ (deg)", 
            r"$\mathrm{\delta \alpha_T^{true}}$ (deg)"],
            bins=np.array([0,25,50,75,100,120,140,160,180]),
            var_evt_reco_col=('del_alpha', '', '', '', '', '', ''),
            var_evt_truth_col=('mc_del_alpha', '', '', '', '', '', ''),
            var_nu_col=('del_alpha', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\delta \alpha_T}$ ($\mathrm{cm}^2$ / deg)"
        )
    
    @classmethod
    def tki_del_phi(cls):
        return cls(
            var_save_name="tki-del_phi",
            var_plot_name="$\\delta \\phi_T$",
            var_labels=[r"$\mathrm{\delta \phi_T}$ (deg)", 
            r"$\mathrm{\delta \phi_T^{reco.}}$ (deg)", 
            r"$\mathrm{\delta \phi_T^{true}}$ (deg)"],
            bins=np.array([0,10,20,30,40,55,70,90,110,130,150,180]),
            var_evt_reco_col=('del_phi', '', '', '', '', '', ''),
            var_evt_truth_col=('mc_del_phi', '', '', '', '', '', ''),
            var_nu_col=('del_phi', '', ''),
            xsec_label=r"$\frac{d\sigma}{d\delta \phi_T}$ ($\mathrm{cm}^2$ / deg)"
        )

    @classmethod
    def opening_angle(cls):
        return cls(
            var_save_name="opening_angle",
            var_plot_name="$\\theta_{\\mu, p}$",
            var_labels=[r"$\mathrm{\theta_{\mu, p}}$ (deg)", 
            r"$\mathrm{\theta_{\mu, p}^{reco.}}$ (deg)", 
            r"$\mathrm{\theta_{\mu, p}^{true}}$ (deg)"],
            bins=np.linspace(0., np.pi, 21),
            var_evt_reco_col=('theta_mu_p', '', '', '', '', '', ''),
            var_evt_truth_col=('mc_theta_mu_p', '', '', '', '', '', ''),
            var_nu_col=('mc', 'theta_mu_p', ''),
            xsec_label=r"$\frac{d\sigma}{d\theta_{\\mu p}}$ ($\mathrm{cm}^2$ / deg)"
        )

    # ==== additional variables for efficiency inspection ====
    # only use var_nu_col
    @classmethod
    def neutrino_energy(cls):
        return cls(
            var_save_name="E_nu",
            var_plot_name="$E_{\\nu}$",
            var_labels=["$\mathrm{E_{\\nu}}$ (GeV)", 
            "", 
            ""],
            bins=np.linspace(0.2, 1.7, 16),
            var_evt_reco_col=('', '', '', '', '', '', ''),
            var_evt_truth_col=('', '', '', '', '', '', ''),
            var_nu_col=('mc', 'E', ''),
            xsec_label=r"$\frac{d\sigma}{dE_{\nu}}$ ($\mathrm{cm}^2$ / GeV)"
        )

    @classmethod
    def muon_direction_x(cls):
        return cls(
            var_save_name="muon-dir_x",
            var_plot_name="cos(\\theta_\mu^x)",
            var_labels=[r"$\mathrm{cos(\theta_\mu^x)}$", 
            r"$\mathrm{cos(\theta_\mu^x^{reco.})}$", 
            r"$\mathrm{cos(\theta_\mu^x^{true})}$"],
            bins=np.array([-1, -0.75, -0.6, -0.45, -0.3, -0.15, 0, 0.15, 0.3, 0.45, 0.6, 0.7, 0.8, 0.9, 1]),
            var_evt_reco_col=('mu', 'pfp', 'trk', 'dir', 'x', '', ''),
            var_evt_truth_col=('mu', 'pfp', 'trk', 'truth', 'p', 'dir', 'x'),
            var_nu_col=('mc', 'mu', 'dir', 'x'),
            xsec_label=r"$\frac{d\sigma}{dcos(\theta_{\\mu}^x)}$ ($\mathrm{cm}^2$)"
        )

    @classmethod
    def muon_direction_y(cls):
        return cls(
            var_save_name="muon-dir_y",
            var_plot_name="cos(\\theta_\mu^y)",
            var_labels=[r"$\mathrm{cos(\theta_\mu^y)}$", 
            r"$\mathrm{cos(\theta_\mu^y^{reco.})}$", 
            r"$\mathrm{cos(\theta_\mu^y^{true})}$"],
            bins=np.array([-1, -0.75, -0.6, -0.45, -0.3, -0.15, 0, 0.15, 0.3, 0.45, 0.6, 0.7, 0.8, 0.9, 1]),
            var_evt_reco_col=('mu', 'pfp', 'trk', 'dir', 'y', '', ''),
            var_evt_truth_col=('mu', 'pfp', 'trk', 'truth', 'p', 'dir', 'y'),
            var_nu_col=('mc', 'mu', 'dir', 'y'),
            xsec_label=r"$\frac{d\sigma}{dcos(\theta_{\\mu}^y)}$ ($\mathrm{cm}^2$)"
        )

    @classmethod
    def proton_direction_x(cls):
        return cls(
            var_save_name="proton-dir_x",
            var_plot_name="cos(\\theta_p^x)",
            var_labels=[r"$\mathrm{cos(\theta_p^x)}$", 
            r"$\mathrm{cos(\theta_p^x^{reco.})}$", 
            r"$\mathrm{cos(\theta_p^x^{true})}$"],
            bins=np.array([-1, -0.75, -0.6, -0.45, -0.3, -0.15, 0, 0.15, 0.3, 0.45, 0.6, 0.7, 0.8, 0.9, 1]),
            var_evt_reco_col=('p', 'pfp', 'trk', 'dir', 'x', '', ''),
            var_evt_truth_col=('p', 'pfp', 'trk', 'truth', 'p', 'dir', 'x'),
            var_nu_col=('mc', 'p', 'dir', 'x'),
            xsec_label=r"$\frac{d\sigma}{dcos(\theta_{\\p}^x)}$ ($\mathrm{cm}^2$)"
        )

    @classmethod
    def proton_direction_y(cls):
        return cls(
            var_save_name="proton-dir_y",
            var_plot_name="cos(\\theta_p^y)",
            var_labels=[r"$\mathrm{cos(\theta_p^y)}$", 
            r"$\mathrm{cos(\theta_p^y^{reco.})}$", 
            r"$\mathrm{cos(\theta_p^y^{true})}$"],
            bins=np.array([-1, -0.75, -0.6, -0.45, -0.3, -0.15, 0, 0.15, 0.3, 0.45, 0.6, 0.7, 0.8, 0.9, 1]),
            var_evt_reco_col=('p', 'pfp', 'trk', 'dir', 'y', '', ''),
            var_evt_truth_col=('p', 'pfp', 'trk', 'truth', 'p', 'dir', 'y'),
            var_nu_col=('mc', 'p', 'dir', 'y'),
            xsec_label=r"$\frac{d\sigma}{dcos(\theta_{\\p}^y)}$ ($\mathrm{cm}^2$)"
        )

    @classmethod
    def list_all_configs(cls, print_summary=True):
        members = inspect.getmembers(cls, predicate=inspect.isroutine)
        config_methods = [name for name, func in members
                          if getattr(func, "__self__", None) == cls
                          and not name.startswith("_")
                          and name not in ("list_all_configs")]
        print("Available VariableConfig options:")
        for name in config_methods:
            config = getattr(cls, name)()
            if print_summary:
                print(f"  {name:20}")
                print(f"    var_save_name   : {config.var_save_name}")
                print(f"    var_plot_name   : {config.var_plot_name}")
                print(f"    var_labels      : {config.var_labels}")
                print(f"    bins            : {np.array2string(config.bins, separator=', ')}")
                print(f"    bin_centers     : {np.array2string(config.bin_centers, separator=', ')}")
                print(f"    var_evt_reco_col: {config.var_evt_reco_col}")
                print(f"    var_evt_truth_col: {config.var_evt_truth_col}")
                print(f"    var_nu_col      : {config.var_nu_col}")
                print(f"    xsec_label      : {config.xsec_label}")
                print()
            else:
                print(f"  {name:20}")
        return config_methods