import numpy as np

def Matrix_C(n, matrix_type):
    """
    Create the matrix C based on the specified type.
      - Type 0: Unit matrix
      - Type 1: First derivative matrix
      - Type 2: Second derivative matrix
      - Otherwise: Third derivative matrix
    """
    epsilon = 1e-6   # needed for 2nd derivative matrix inversion
    epsilon2 = 1e-2  # needed for 3rd derivative matrix inversion
    C = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if matrix_type == 0:
                if i == j:
                    C[i, j] = 1
            elif matrix_type == 1:
                if (j - i) == 1:
                    C[i, j] = 1
                if i == j:
                    C[i, j] = -1
            elif matrix_type == 2:
                if abs(i - j) == 1:
                    C[i, j] = 1
                if i == j:
                    C[i, j] = -2 + epsilon
                    if i == 0 or i == n - 1:
                        C[i, j] = -1 + epsilon
            else:  # third derivative matrix
                if (i - j) == 2:
                    C[i, j] = -1
                if (i - j) == 1:
                    C[i, j] = 2
                    if i == 1:
                        C[i, j] = 0
                if (i - j) == -1:
                    C[i, j] = -2
                    if i == n - 2:
                        C[i, j] = 0
                if (i - j) == -2:
                    C[i, j] = 1
                if i == j:
                    C[i, j] = epsilon2
                    if i == 0 or i == 1:
                        C[i, j] = 1 + epsilon2
                    if i == n - 1 or i == n - 2:
                        C[i, j] = -1 + epsilon2
    return C


def WienerSVD(Response, Signal, Measure, Covariance, C_type, Norm_type):
    """
    Perform Wiener-SVD unfolding.

    Parameters:
      Response   : 2D numpy array (m x n) - the response matrix.
      Signal     : 1D numpy array (n,)   - the signal prediction.
      Measure    : 1D numpy array (m,)   - the measured data.
      Covariance : 2D numpy array (m x m) - covariance matrix.
      C_type     : int - type specifier for the smoothness matrix.
      Norm_type  : float - normalization exponent for Signal.

    Returns:
      A dictionary containing:
         'unfold'      : unfolded spectrum (1D numpy array),
         'AddSmear'    : additive smearing matrix (2D numpy array),
         'WF'          : Wiener filter factors (1D numpy array),
         'UnfoldCov'   : covariance matrix of the unfolded spectrum (2D numpy array),
         'CovRotation' : covariance rotation matrix (2D numpy array).
    """
    m, n = Response.shape  # m measure, n signal bins

    # Decomposition of the Covariance matrix to obtain Q.
    U_cov, s_cov, Vh_cov = np.linalg.svd(Covariance)
    # Q0 is the transpose of V from the SVD (numpy's Vh is already V^T).
    Q0 = Vh_cov
    # Build a diagonal matrix of 1/sqrt(s) (with protection against division by zero)
    err_diag = np.array([1/np.sqrt(val) if val != 0 else 0 for val in s_cov])
    err = np.diag(err_diag)
    Q = err @ Q0

    # Transform Measure and Response
    M_trans = Q @ Measure
    R = Q @ Response

    # Build the smoothness matrix
    C0 = Matrix_C(n, C_type)
    normsig = np.zeros((n, n))
    for i in range(n):
        normsig[i, i] = 1.0 / (Signal[i] ** Norm_type)
    C0 = C0 @ normsig

    # Copy and invert the smoothness matrix
    C = C0.copy()
    C_inv = np.linalg.inv(C0)
    Signal_mod = C @ Signal
    R = R @ C_inv

    # SVD decomposition of R
    U, D, Vh = np.linalg.svd(R, full_matrices=False)
    U_t = U.T
    V = Vh.T
    # Construct D_t (an n x m matrix with diagonal elements set to D)
    D_t = np.zeros((n, m), dtype=float)
    for i in range(min(n, m)):
        D_t[i, i] = D[i]

    # Compute S = V_t * Signal_mod, where V_t is V^T (Vh in numpy)
    S_vec = Vh @ Signal_mod

    # Wiener Filter
    W = np.zeros((n, n))
    W0 = np.zeros((n, n))
    WF = np.zeros(n)
    for i in range(n):
        W[i, i] = (S_vec[i]**2) / ((D[i]**2 * S_vec[i]**2) + 1)
        WF[i] = D[i]**2 * W[i, i]
        W0[i, i] = WF[i]

    # Compute unfolded spectrum
    unfold = C_inv @ V @ W @ D_t @ U_t @ M_trans
    AddSmear = C_inv @ V @ W0 @ Vh @ C

    # Covariance matrix of the unfolded spectrum.
    covRotation = C_inv @ V @ W @ D_t @ U_t @ Q
    CovRotation = covRotation
    UnfoldCov = covRotation @ Covariance @ covRotation.T

    return {
        'unfold': unfold,
        'AddSmear': AddSmear,
        'WF': WF,
        'UnfoldCov': UnfoldCov,
        'CovRotation': CovRotation
    }