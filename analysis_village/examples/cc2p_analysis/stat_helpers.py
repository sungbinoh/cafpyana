import numpy as np
from scipy.stats import gamma

alpha = 0.3173  # for 68.27% CL, alpha = 1 - 0.6827

def return_data_stat_err(data_array):
    L = np.where(data_array == 0, 0, gamma.ppf(alpha / 2, data_array))
    U = np.where(data_array == 0,
             gamma.isf(alpha, data_array + 1),
             gamma.isf(alpha / 2, data_array + 1))

    yerr_low = np.where(data_array == 0, 0.1, data_array - L)
    yerr_high = np.where(data_array == 0, 1.8, U - data_array)

    return yerr_low, yerr_high