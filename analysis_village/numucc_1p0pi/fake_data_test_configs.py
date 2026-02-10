import numpy as np

class FakeDataWeights:
    def __init__(self, mc_evt_df, mc_nu_df, var_config):
        self.mc_evt_df = mc_evt_df
        self.mc_nu_df = mc_nu_df
        self.var_config = var_config

    def get_weights(self, test_name, **kwargs):
        # Always start from ones (unless otherwise needed)
        # TODO: place normalization here?
        weights_fake_data = np.ones(len(self.mc_evt_df))
        weight_fakedata_signal_truth = np.ones(len(self.mc_nu_df[self.mc_nu_df.nuint_categ == 1]))
        
        # MEC normalization
        if test_name == "mec_test":
            scale_factor = kwargs.get("scale_factor", 0.5)
            weights_fake_data[self.mc_evt_df.mc.genie_mode == 10] *= scale_factor
            weight_fakedata_signal_truth[self.mc_nu_df[self.mc_nu_df.nuint_categ == 1].mc.genie_mode == 10] *= scale_factor

        # QE normalization
        elif test_name == "qe_test":
            scale_factor = kwargs.get("scale_factor", 1.2)
            weights_fake_data[self.mc_evt_df.mc.genie_mode == 0] *= scale_factor
            weight_fakedata_signal_truth[self.mc_nu_df[self.mc_nu_df.nuint_categ == 1].mc.genie_mode == 0] *= scale_factor

        # np normalization
        elif test_name == "np_test":
            Np_scale = kwargs.get("Np_scale", 2)
            weights_fake_data[self.mc_evt_df.nuint_categ == 2] *= Np_scale
            weight_fakedata_signal_truth[self.mc_nu_df[self.mc_nu_df.nuint_categ == 1].nuint_categ == 2] *= Np_scale

        # sig normalization
        elif test_name == "sig_test":
            sig_scale = kwargs.get("sig_scale", 1.2)
            weights_fake_data[self.mc_evt_df.nuint_categ == 1] *= sig_scale
            weight_fakedata_signal_truth[self.mc_nu_df[self.mc_nu_df.nuint_categ == 1].nuint_categ == 1] *= sig_scale

        # Q2 tilt
        elif test_name.startswith("q2_test_alpha_"):
            alpha = kwargs.get("alpha", 0.3)
            Q2 = 2 * self.mc_evt_df.mc.E * self.mc_evt_df.mu.pfp.trk.truth.p.startE * (1 - self.mc_evt_df.mu.pfp.trk.truth.p.dir.z)
            weights_fake_data *= np.ones(len(self.mc_evt_df)) + alpha * (Q2 - Q2.mean())/Q2.mean()
            weights_fake_data[np.isnan(weights_fake_data)] = 1
            assert np.isnan(weights_fake_data).sum() == 0
            Q2_nu = self.mc_nu_df[self.mc_nu_df.nuint_categ == 1].mc.Q2
            weight_fakedata_signal_truth *= np.ones(len(Q2_nu)) + alpha * (Q2_nu - Q2_nu.mean())/Q2_nu.mean()
        
        # cos(theta) scale
        elif test_name.startswith("costh_weight_scale_"):
            scale_factor = kwargs.get("scale_factor", 0.7)
            weights_fake_data[self.mc_evt_df.mu.pfp.trk.truth.p.dir.z > 0.9] *= scale_factor
            weight_fakedata_signal_truth[self.mc_nu_df[self.mc_nu_df.nuint_categ == 1].mc.mu.dir.z > 0.9] *= scale_factor

        # Proton P tilt
        elif test_name.startswith("proton_P_tilt_alpha_"):
            alpha = kwargs.get("alpha", 0.3)
            P_p_evt = self.mc_evt_df.p.pfp.trk.truth.p.totp
            weights_fake_data = np.ones(len(self.mc_evt_df)) + alpha * (P_p_evt - P_p_evt.mean())/P_p_evt.mean()
            P_p_nu = self.mc_nu_df[self.mc_nu_df.nuint_categ == 1].mc.p.totp
            weight_fakedata_signal_truth = np.ones(len(P_p_nu)) + alpha * (P_p_nu - P_p_nu.mean())/P_p_nu.mean()

        # Bump
        elif test_name.startswith("bump_"):
            bump_pos = kwargs.get("bump_pos", 0.6)
            bump_width = kwargs.get("bump_width", 0.0015)
            bump_height = kwargs.get("bump_height", 0.001)
            bump_height = bump_height*len(self.mc_evt_df)/len(self.var_config.bin_centers)
            bump_var_evt = self.mc_evt_df[self.var_config.var_evt_truth_col]
            weights_fake_data = np.ones(len(self.mc_evt_df)) + bump_height * np.exp(-0.5 * (bump_var_evt - bump_pos)**2 / bump_width**2)
            weights_fake_data[np.isnan(weights_fake_data)] = 1.
            bump_var_nu = self.mc_nu_df[self.mc_nu_df.nuint_categ == 1][self.var_config.var_nu_col]
            weight_fakedata_signal_truth = np.ones(len(bump_var_nu)) + bump_height * np.exp(-0.5 * (bump_var_nu - bump_pos)**2 / bump_width**2)
        
        # TODO
        # # Enhance tail
        # elif test_name.startswith("enhance_tail_alpha_"):
        #     alpha = kwargs.get("alpha", 0.7)
        #     weights_fake_data = np.ones(len(self.mc_evt_df))
        #     weights_fake_data[self.mc_evt_df.del_p > 0.25] = alpha
        #     weights_fake_data[np.isnan(weights_fake_data)] = 1.
        #     weight_fakedata_signal_truth = np.ones(len(self.mc_nu_df[self.mc_nu_df.nuint_categ == 1]))
        #     weight_fakedata_signal_truth[self.mc_nu_df[self.mc_nu_df.nuint_categ == 1].mc.del_p > 0.25] = alpha
        #     weight_fakedata_signal_truth[np.isnan(weight_fakedata_signal_truth)] = 1.

        else:
            raise ValueError(f"Unknown test_name '{test_name}' provided.")

        return weights_fake_data, weight_fakedata_signal_truth