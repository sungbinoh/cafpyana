import numpy as np
from scipy.stats import gamma
from scipy.stats import chi2

alpha = 0.3173  # for 68.27% CL, alpha = 1 - 0.6827

def return_data_stat_err(data_array):
    L = np.where(data_array == 0, 0, gamma.ppf(alpha / 2, data_array))
    U = np.where(data_array == 0,
             gamma.isf(alpha, data_array + 1),
             gamma.isf(alpha / 2, data_array + 1))

    yerr_low = np.where(data_array == 0, 0.1, data_array - L)
    yerr_high = np.where(data_array == 0, 1.8, U - data_array)

    return yerr_low, yerr_high


def get_chi2(data, model, cov):
    chi2_value = (data - model) @ np.linalg.inv(cov) @ (data - model)
    ndof = len(data)
    p_value = 1 - chi2.cdf(chi2_value, df=ndof)
    return chi2_value, p_value