import numpy as np

def cov_from_fraccov(cov_frac, cv_vals):
    cov = np.zeros_like(cov_frac)
    for i in range(cov_frac.shape[0]):
        for j in range(cov_frac.shape[1]):
            cov[i, j] = cov_frac[i, j] * (cv_vals[i] * cv_vals[j])
    return cov

def corr_from_fraccov(cov_frac):
    corr = np.zeros_like(cov_frac)
    for i in range(cov_frac.shape[0]):
        for j in range(cov_frac.shape[1]):
            corr[i, j] = cov_frac[i, j] / np.sqrt(cov_frac[i, i] * cov_frac[j, j])
    return corr

def get_covariance_matrix(univ_events_a, cv_events_a,
                          univ_events_b, cv_events_b):

    # The two input array sets should have exactly the same number of bins and n_univ
    n_univ, n_bins = univ_events_a.shape

    cov_frac = np.zeros((n_bins, n_bins))
    cov = np.zeros((n_bins, n_bins))

    # looping & calculating with the CV value for clarity, 
    # but techincally np.cov should also be fine under the assumption of gaussian universes that we're using
    for uidx in range(n_univ):
        for i in range(univ_events_a.shape[1]):
            for j in range(univ_events_b.shape[1]):
                nom_i = cv_events_a[i]
                nom_j = cv_events_b[j]

                univ_i = univ_events_a[uidx, i] 
                univ_j = univ_events_b[uidx, j] 

                cov_entry = (univ_i - nom_i) * (univ_j - nom_j)
                frac_cov_entry = ((univ_i - nom_i) / nom_i) * ( (univ_j - nom_j) / nom_j)

                # TODO: uboone code has clipping that I'm not sure why.. investigate later
                # if cov_entry > 0:
                #     this_cov = max( cov_entry, eps * scale_factor)
                # else:
                #     this_cov = min( cov_entry, eps * scale_factor)

                # if frac_cov_entry > 0:
                #     this_frac_cov = max( frac_cov_entry, eps * scale_factor)
                # else:
                #     this_frac_cov = min( frac_cov_entry, eps * scale_factor)

                cov[i, j] += cov_entry
                cov_frac[i, j] += frac_cov_entry

    cov = cov / n_univ
    cov_frac = cov_frac / n_univ
    corr = np.zeros_like(cov)
    for i in range(len(cv_events_a)):
        for j in range(len(cv_events_b)):
            corr[i, j] = cov[i, j] / (np.sqrt(cov[i, i]) * np.sqrt(cov[j, j]))

    return {"cov_frac": cov_frac, 
            "cov": cov,
            "corr": corr,
            }

def get_covariance_matrix_self(univ_events,
                               cv_event):
    return get_covariance_matrix(univ_events, cv_event, univ_events, cv_event)
