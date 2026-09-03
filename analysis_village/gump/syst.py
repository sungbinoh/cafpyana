import numpy as np
import kinematics
import pandas as pd

# Muon range->momentum, replicating the calculator that filled rangeP.p_muon in
# the CAFs: sbncode TrackMomentumCalculator.cxx @35baeab (LArReco). CSDA range
# table [g/cm^2] / 1.396 -> cm vs KE [MeV], float32 as in the C++ std::array;
# ROOT's TSpline3 default boundary conditions correspond to scipy's
# "not-a-knot" (verified to <5e-7 relative against the dataframes).
_MUON_RANGE_GRAMPERCM = np.array([
    9.833E-1, 1.786E0, 3.321E0, 6.598E0, 1.058E1, 3.084E1, 4.250E1, 6.732E1,
    1.063E2,  1.725E2, 2.385E2, 4.934E2, 6.163E2, 8.552E2, 1.202E3, 1.758E3,
    2.297E3,  4.359E3, 5.354E3, 7.298E3, 1.013E4, 1.469E4, 1.910E4, 3.558E4,
    4.326E4,  5.768E4, 7.734E4, 1.060E5, 1.307E5], dtype=np.float32)
_MUON_KE_MEV = np.array([
    10, 14, 20, 30, 40, 80, 100, 140, 200, 300, 400, 800, 1000, 1400, 2000,
    3000, 4000, 8000, 10000, 14000, 20000, 30000, 40000, 80000, 100000,
    140000, 200000, 300000, 400000], dtype=np.float32)
_MUON_M_MEV = 105.7 # sbncode constant; NOT kinematics.MUON_MASS

def muon_range_momentum(trkrange):
    """Range [cm] -> momentum [GeV] for muons, as in the CAFs.

    Invalid ranges and KE < 0 give NaN (sbncode returns -1/-999; NaN is the
    dataframe-native equivalent and falls out of histograms).
    """
    from scipy.interpolate import CubicSpline
    spline = CubicSpline(_MUON_RANGE_GRAMPERCM/np.float32(1.396),
                         _MUON_KE_MEV, bc_type="not-a-knot")
    trkrange = np.asarray(trkrange, dtype=float)
    with np.errstate(invalid="ignore"):
        KE = spline(trkrange)
        p = np.sqrt(KE**2 + 2*_MUON_M_MEV*KE)/1000.
    return np.where(np.isfinite(trkrange) & (trkrange >= 0) & (KE >= 0), p, np.nan)

def recompute_kinematics(s, mu_p=None, p_p=None, BE=None):
    """Recompute the downstream reco kinematics (nu_E_calo, nu_E_ccqe, del_p,
    del_Tp, del_phi, mu_E, mu_T) on the flat df `s`, in place.

    mu_p defaults to sqrt(mu_E^2 - m_mu^2) -- exact, since the reco momentum
    is range-based; BE defaults to kinematics.BE. Pass either to build a
    shifted-kinematics universe from the CV.

    The proton system is the SUMMED candidate-proton system (psum): requires
    the stored `psum_E` (summed energy) and `psum_dir_*` (unit direction of
    the vector-sum momentum) columns, so the kinematics generalize to the Np
    case. The momentum magnitude comes from the stored `psum_p` column when
    present; otherwise the frame is treated as a single on-shell proton and
    p_p = sqrt(psum_E^2 - m_p^2). Pass `p_p` to override the stored value.
    """

    if BE is None:
        BE = kinematics.BE

    if mu_p is None:
        if 'mu_p' in s.columns:
            mu_p = s.mu_p
        elif 'mu_E' in s.columns:
            mu_p = pd.Series(np.sqrt(np.maximum(s.mu_E.to_numpy()**2 - kinematics.MUON_MASS**2, 0)),
                             index=s.index)
        else:
           raise ValueError("You don't have the necessary columns to recompute kinematics.")

    if "psum_E" not in s.columns:
        raise KeyError("recompute_kinematics requires stored 'psum_E'/'psum_dir_*' "
                       "(summed proton system) columns.")

    mu_dir = pd.DataFrame({"x": s.mu_dir_x, "y": s.mu_dir_y, "z": s.mu_dir_z})
    p_dir = pd.DataFrame({"x": s.psum_dir_x, "y": s.psum_dir_y, "z": s.psum_dir_z})
    p_E = s.psum_E
    if p_p is None:
        if "psum_p" in s.columns:
            p_p = s.psum_p
        else:
            p_p = np.sqrt(np.maximum(p_E**2 - kinematics.PROTON_MASS**2, 0))

    # number of candidate protons = candidate pfps minus the muon; the TKI
    # residual mass and the calo BE term both scale with it
    n_proton = s.n_pfp - 1
    tki = kinematics.transverse_kinematics(mu_p, mu_dir, p_p, p_dir, p_E,
                                           n_proton=n_proton, BE=BE)
    # Calorimetric energy in the PRODUCTION convention (maple recoE):
    #     nu_E_calo = E_mu + sum_i ke_i + n*BE = E_mu + (psum_E - n*m_p) + n*BE
    # computed at the NOMINAL BE (the mode-scaled shift is applied below).
    # NB: deliberately NOT kinematics.neutrino_energy -- its KE term
    # (psum_E - M_inv, with M_inv the summed-system invariant mass) equals the
    # per-proton KE sum only at n = 1 and understates it by O(100 MeV) for
    # non-collinear multi-proton systems, which dwarfed the 25 MeV BE shift.
    # psum_E - n*m_p reproduces the stored psum_ke to float precision.
    s["nu_E_calo"] = tki["mu_E"] + (p_E - n_proton*kinematics.PROTON_MASS) \
        + n_proton*kinematics.BE

    # Binding-energy shift on the calorimetric estimator: subtract the actual
    # per-event shift (BE - nominal), which already carries the caller's
    # per-mode scaling (shift_binding_energy: x1 QE/RES/DIS, x sqrt(2) MEC,
    # 0 COH / non-neutrino). This keeps nu_E_calo consistent with nu_E_ccqe
    # and the TKI variables, which receive the same BE directly. At nominal BE
    # (the CV / track-splitting path) the term is identically zero.
    s["nu_E_calo"] = s["nu_E_calo"] - (BE - kinematics.BE)

    # muon-only CCQE energy estimator. Recomputed here (not just carried over as
    # a derived column) so both universes built on this function are consistent:
    # the binding-energy shift moves it via BE, and the track-splitting universe
    # moves it via the truncated mu_p.
    s["nu_E_ccqe"] = kinematics.neutrino_energy_ccqe(mu_p, s.mu_dir_z, BE=BE)
    s["del_p"] = tki["del_p"]
    s["del_Tp"] = tki["del_Tp"]
    s["del_phi"] = tki["del_phi"]
    s["mu_E"] = tki["mu_E"]
    s["mu_T"] = tki["mu_E"] - kinematics.MUON_MASS

    return s

# Columns recompute_kinematics writes -- the ones the BE universe overwrites
# on the shifted rows.
_BE_RECOMPUTED_COLS = ["nu_E_calo", "nu_E_ccqe", "del_p", "del_Tp", "del_phi",
                       "mu_E", "mu_T"]

def shift_binding_energy(df, dBE, fraction=0.5, scale="glob_scale"):
    """Universe df for a binding-energy shift: copy of the CV with the reco
    kinematics recomputed under BE -> BE + dBE, scaled per interaction mode:
    the nominal dBE for QE/RES/DIS (genie_mode 0/1/2), dBE*sqrt(2) for MEC
    (genie_mode 10, two-nucleon initial state), and no shift for COH and
    non-neutrino rows (genie_mode 3/NaN -- no struck bound nucleon).

    `fraction` applies the shift to only that fraction of events, IN PLACE on
    the copied frame: a deterministic, evenly-interleaved subset of row
    positions gets the shifted kinematics; every other row keeps its stored CV
    values (no recompute) and its full weight. The returned frame therefore
    has the same rows, index and weights as the input. Reproducible for a
    given row order; `scale` is unused and kept for API compatibility.

    The selection cuts use no recomputed column, so 'selected' etc. carry over
    from the CV unchanged.
    """
    s = df.copy()
    if fraction <= 0:
        return s

    n = len(s)
    if fraction >= 1.0:
        shifted = np.ones(n, dtype=bool)
    else:
        # evenly-interleaved deterministic selection of a fraction of the row
        # positions (fraction=0.5 -> every other row)
        shifted = (np.arange(n)*fraction) % 1.0 < fraction

    mode = s.genie_mode.to_numpy()[shifted]
    mode_scl = np.where(np.isin(mode, [0, 1, 2]), 1.,
               np.where(mode == 10, np.sqrt(2.), 0.))

    sub = s.iloc[shifted].copy()
    recompute_kinematics(sub, BE=kinematics.BE + dBE*mode_scl)

    # mu_T is not stored in the flat dfs; derive the unshifted rows' value
    # from the stored mu_E (exact -- the reco momentum is range-based).
    if "mu_T" not in s.columns:
        s["mu_T"] = s["mu_E"] - kinematics.MUON_MASS

    # A frame slimmed to just the recompute inputs lacks stored kinematic
    # columns for the unshifted rows to keep -- rebuild ONLY those missing
    # columns at the nominal BE. Columns the input does carry keep their
    # stored CV values on the unshifted rows (the recompute does not
    # reproduce the stored CV exactly, so prefer the stored values wherever
    # they exist).
    missing = [c for c in _BE_RECOMPUTED_COLS if c not in s.columns]
    if missing:
        for c in missing:
            s[c] = np.nan
        if (~shifted).any():
            cvsub = s.iloc[~shifted].copy()
            recompute_kinematics(cvsub)
            for c in missing:
                s.iloc[~shifted, s.columns.get_loc(c)] = cvsub[c].to_numpy()

    for c in _BE_RECOMPUTED_COLS:
        # the stored columns are float32; upcast so the float64 recompute
        # assigns without pandas' incompatible-dtype warning
        if s[c].dtype != np.float64:
            s[c] = s[c].astype(np.float64)
        s.iloc[shifted, s.columns.get_loc(c)] = sub[c].to_numpy()

    return s

def v_variation(df, setvars):
    df = df[[c for c in df.columns if "univ" not in c]].copy()
    
    for (new, old) in setvars:
        # Raise an error immediately if the 'new' column doesn't exist
        if new not in df.columns:
            raise KeyError(f"Column '{new}' does not exist in the DataFrame.")
            
        # Only overwrite rows where 'old' is NOT NaN
        is_not_nan = df[old].notna()
        df.loc[is_not_nan, new] = df.loc[is_not_nan, old]
            
    return df

def v_chi2alpha(df):
    setvars = [
        ("mu_chi2_of_mu_cand", "mu_chi2alpha_p_of_mu_cand"),
        ("mu_chi2_of_prot_cand",  "mu_chi2alpha_p_of_prot_cand"),
        ("prot_chi2_of_mu_cand", "prot_chi2alpha_p_of_mu_cand"),
        ("prot_chi2_of_prot_cand",  "prot_chi2alpha_p_of_prot_cand"),
    ]
    return v_variation(df, setvars)

def v_chi2beta(df):
    setvars = [
        ("mu_chi2_of_mu_cand", "mu_chi2beta_p_of_mu_cand"),
        ("mu_chi2_of_prot_cand",  "mu_chi2beta_p_of_prot_cand"),
        ("prot_chi2_of_mu_cand", "prot_chi2beta_p_of_mu_cand"),
        ("prot_chi2_of_prot_cand",  "prot_chi2beta_p_of_prot_cand"),
    ]
    return v_variation(df, setvars)

def v_chi2R(df):
    setvars = [
        ("mu_chi2_of_mu_cand", "mu_chi2R_p_of_mu_cand"),
        ("mu_chi2_of_prot_cand",  "mu_chi2R_p_of_prot_cand"),
        ("prot_chi2_of_mu_cand", "prot_chi2R_p_of_mu_cand"),
        ("prot_chi2_of_prot_cand",  "prot_chi2R_p_of_prot_cand"),
    ]
    return v_variation(df, setvars)

def v_chi2dedxbias(df):
    # ICARUS-only dE/dx bias variation (no-op in SBND, where the columns equal CV)
    setvars = [
        ("mu_chi2_of_mu_cand", "mu_chi2dedxbias_of_mu_cand"),
        ("mu_chi2_of_prot_cand",  "mu_chi2dedxbias_of_prot_cand"),
        ("prot_chi2_of_mu_cand", "prot_chi2dedxbias_of_mu_cand"),
        ("prot_chi2_of_prot_cand",  "prot_chi2dedxbias_of_prot_cand"),
    ]
    return v_variation(df, setvars)

def v_flashscale(df, updn):
    TRIG_PE_SCALE     = {1: 0.642, 2: 0.632, 4: 0.358}  # Run -> best-fit s
    TRIG_PE_SCALE_UNC = {1: 0.005, 2: 0.024, 4: 0.017}  # Run -> unc. on s

    # 1. Map the dictionaries to the 'Run' column to get Series for both parts
    unc_series = df["Run"].map(TRIG_PE_SCALE_UNC).astype(float)
    scale_series = df["Run"].map(TRIG_PE_SCALE).astype(float)

    # 2. Divide the two Series element-wise
    f = unc_series / scale_series

    # 3. Perform your vectorized math
    df["flash_maxpe_var"] = df["flash_maxpe"] * (1 - updn * f)

    setvars = [
        ("flash_maxpe", "flash_maxpe_var"),
    ]
    ret = v_variation(df, setvars)

    # Note: Added axis=1 here because drop defaults to dropping rows
    return ret.drop("flash_maxpe_var", axis=1)

def v_chi2smear(df):
    setvars = [
        ("mu_chi2_of_mu_cand", "mu_chi2smear13_of_mu_cand"),
        ("mu_chi2_of_prot_cand",  "mu_chi2smear13_of_prot_cand"),
        ("prot_chi2_of_mu_cand", "prot_chi2smear13_of_mu_cand"),
        ("prot_chi2_of_prot_cand",  "prot_chi2smear13_of_prot_cand"),
    ]
    return v_variation(df, setvars)

def v_chi2hi(df):
    setvars = [
        ("mu_chi2_of_mu_cand", "mu_chi2hi_of_mu_cand"),
        ("mu_chi2_of_prot_cand",  "mu_chi2hi_of_prot_cand"),
        ("prot_chi2_of_mu_cand", "prot_chi2hi_of_mu_cand"),
        ("prot_chi2_of_prot_cand",  "prot_chi2hi_of_prot_cand"),
    ]
    return v_variation(df, setvars)

def v_chi2lo(df):
    setvars = [
        ("mu_chi2_of_mu_cand", "mu_chi2lo_of_mu_cand"),
        ("mu_chi2_of_prot_cand",  "mu_chi2lo_of_prot_cand"),
        ("prot_chi2_of_mu_cand", "prot_chi2lo_of_mu_cand"),
        ("prot_chi2_of_prot_cand",  "prot_chi2lo_of_prot_cand"),
    ]
    return v_variation(df, setvars)

def v_chi22hi(df):
    setvars = [
        ("mu_chi2_of_mu_cand", "mu_chi22hi_of_mu_cand"),
        ("mu_chi2_of_prot_cand",  "mu_chi22hi_of_prot_cand"),
        ("prot_chi2_of_mu_cand", "prot_chi22hi_of_mu_cand"),
        ("prot_chi2_of_prot_cand",  "prot_chi22hi_of_prot_cand"),
    ]
    return v_variation(df, setvars)

def v_chi22lo(df):
    setvars = [
        ("mu_chi2_of_mu_cand", "mu_chi22lo_of_mu_cand"),
        ("mu_chi2_of_prot_cand",  "mu_chi22lo_of_prot_cand"),
        ("prot_chi2_of_mu_cand", "prot_chi22lo_of_mu_cand"),
        ("prot_chi2_of_prot_cand",  "prot_chi22lo_of_prot_cand"),
    ]
    return v_variation(df, setvars)

class SystematicList(object):
    def __init__(self, systs):
        self.systs = systs

    def cov(self, var, cut, bins, NCV, shapeonly=False, fillna=np.nan):
        if len(self.systs) == 0:
            return np.zeros((NCV.size, NCV.size))
        return np.sum([s.cov(var, cut, bins, NCV, shapeonly=shapeonly, fillna=fillna) for s in self.systs], axis=0)

def outern(arrs):
    ret = arrs[0]
    for a in arrs[1:]:
        ret = np.outer(ret, a)

    return ret

class ConcatBins(list):
    """Marker for a per-variable list of bin edges that must be histogrammed
    SEPARATELY and concatenated, instead of forming a joint N-D histogram.

    Pass it as the `bins` argument of Systematic.cov together with a list of
    variables. The universe vector is then [N(var0), N(var1), ...] and the
    covariance comes out in the block form

        [[Cov(x,x), Cov(x,y)],
         [Cov(y,x), Cov(y,y)]]

    which is what conditional_constraint() consumes. This is the two-VARIABLE
    analogue of CorrelatedSystematic, which concatenates two dataframes sharing
    one variable.

    Every systematic class routes its histogramming through histflat(), so the
    concatenated mode works for all of them without per-class handling.
    Area normalization (shapeonly=True) is not meaningful on a concatenated
    vector -- it would normalize the two variables jointly -- and is rejected.
    """

def block_masks(df, cut, nblock):
    """Positional row masks, one per concatenated block.

    `cut` is normally a single column name shared by every block. With
    ConcatBins it may instead be a per-block list of column names, so each block
    selects its own rows -- e.g. one del_p slice per block. Masks are positional
    (.to_numpy()) because the CV frames carry duplicate index labels.
    """
    if isinstance(cut, (list, tuple)):
        if len(cut) != nblock:
            raise ValueError("got %d cut columns for %d blocks" % (len(cut), nblock))
        return [np.asarray(df[c], dtype=bool) for c in cut]

    return [np.asarray(df[cut], dtype=bool)]*nblock


def histflat(df, var, cut, bins, weights, fillna=np.nan):
    """Flat histogram vector for one universe of `df`.

    `weights` is aligned with the whole frame; rows are selected by `cut`.
    Normally this is the joint N-D histogram, flattened. With a ConcatBins each
    variable is histogrammed on its own and the 1-D results are concatenated,
    and `cut` may then be a per-block list (see block_masks).
    """
    w = np.asarray(weights)

    if isinstance(bins, ConcatBins):
        ms = block_masks(df, cut, len(bins))
        return np.concatenate([
            np.histogram(np.asarray(df[var[i]].fillna(fillna))[ms[i]], bins=bins[i],
                         weights=w[ms[i]])[0]
            for i in range(len(bins))])

    m = block_masks(df, cut, 1)[0]
    return np.histogramdd([np.asarray(df[v].fillna(fillna))[m] for v in var],
                          bins=bins, weights=w[m])[0].flatten()

def _binwidths(bins):
    """Per-bin widths of the flattened histogram vector (shapeonly support).

    NB the flatten(): outern() returns the N-D outer product, but the histogram
    it normalizes has been flattened, so the widths must be flattened to match.
    A no-op for a single variable (the only case previously exercised, which is
    why the N-D shapeonly path used to raise a broadcast error).
    """
    if isinstance(bins, ConcatBins):
        raise ValueError("shapeonly area normalization is not defined for "
                         "ConcatBins: it would normalize the concatenated "
                         "variables jointly. Use absolute normalization.")

    return outern([b[1:] - b[:-1] for b in bins]).flatten()

def conditional_constraint(cov, nx):
    """Gaussian conditional constraint: the covariance of the second block of a
    joint covariance, given a measurement of the first.

        cov_cond = Cyy - Cyx Cxx^-1 Cxy

    i.e. the standard Schur-complement / near-detector-style conditional
    constraint (promoted here out of nb/SignalBoxSystematics-ReCAF.ipynb, and
    generalized to unequal block sizes). The conditional covariance depends
    only on the covariance blocks, not on the central values; the returned gain
    matrix K is what shifts the central value,

        mu_y|x = mu_y + K (x_obs - mu_x)

    Parameters
    ----------
    cov : array-like, shape (nx + ny, nx + ny)
        Full joint covariance of (x, y), ordered
        [[Cov(x,x), Cov(x,y)], [Cov(y,x), Cov(y,y)]].
    nx : int
        Size of the constraining (x) block.

    Returns
    -------
    (cov_cond, K, cond_number) : the (ny, ny) conditional covariance, the
        (ny, nx) gain matrix, and the condition number of Cxx. Cxx is inverted
        with np.linalg.pinv when it is ill-conditioned; callers should log
        cond_number.
    """
    cov = np.asarray(cov, dtype=float)
    n = cov.shape[0]
    assert cov.shape == (n, n), "joint covariance must be square"
    assert 0 < nx < n, "nx must split the joint covariance"

    Cxx = cov[:nx, :nx]
    Cxy = cov[:nx, nx:]
    Cyx = cov[nx:, :nx]
    Cyy = cov[nx:, nx:]

    cond_number = np.linalg.cond(Cxx)
    if not np.isfinite(cond_number) or cond_number > 1e12:
        Cxx_inv = np.linalg.pinv(Cxx)
    else:
        Cxx_inv = np.linalg.inv(Cxx)

    K = Cyx @ Cxx_inv

    cov_cond = Cyy - K @ Cxy
    cov_cond = 0.5 * (cov_cond + cov_cond.T)  # enforce symmetry vs. roundoff

    return cov_cond, K, cond_number

class Systematic(object):
    def __init__(self):
        pass
        
    def nuniv(self):
        pass
        
    def univ(self, var, cut, bins, i_univ, fillna=np.nan):
        pass

    # Whether to average the separate universes, or not (i.e. treat them as different uncertainties)
    def avg(self):
        return True # true by default
    
    def cov(self, var, cut, bins, NCV, shapeonly=False, fillna=np.nan):
        if not isinstance(var, list):
            var = [var]
            bins = [bins]

        if shapeonly:
            diff = _binwidths(bins)
            norm = np.sum(NCV*diff)
            if norm > 1e-5:
                NCV = NCV / norm

        N_univ = []
        for i_univ in range(self.nuniv()):
            N = self.univ(var, cut, bins, i_univ, fillna=fillna)

            if shapeonly:
                diff = _binwidths(bins)
                norm = np.sum(N*diff)
                if norm > 1e-5:
                    N = N / norm
            N_univ.append(N)
    
        cov =  np.sum([np.outer(N - NCV, N - NCV) for N in N_univ], axis=0)
        if self.avg():
            cov = cov / self.nuniv()

        return cov

class NormalizationSystematic(Systematic):
    def __init__(self, norm):
       self.norm = norm

    def nuniv(self):
        return 1
        
    def cov(self, var, cut, bins, NCV, shapeonly=False, fillna=np.nan):
        self.CV = NCV
        return super().cov(var, cut, bins, NCV, shapeonly=shapeonly, fillna=fillna)

    def univ(self, var, cut, bins, i_univ, fillna=np.nan):
        assert(i_univ == 0)
        return self.CV*(1 + self.norm)

class SystSampleSystematic(Systematic):
    def __init__(self, df, scale="glob_scale", norm=1.):
        self.df = df
        self.scale = scale
        self.norm = norm
        
    def nuniv(self):
        return 1
        
    def cov(self, var, cut, bins, NCV, shapeonly=False, fillna=np.nan):
        self.CV = NCV
        return super().cov(var, cut, bins, NCV, shapeonly=shapeonly, fillna=fillna)

    def univ(self, var, cut, bins, i_univ, fillna=np.nan):
        assert(i_univ == 0)
        if not isinstance(var, list):
            var = [var]
            bins = [bins]

        return histflat(self.df, var, cut, bins, self.df[self.scale],
                        fillna=fillna)*self.norm + self.CV

class StatSampleSystematic(object):
    def __init__(self, df, scale="glob_scale", norm=1):
        self.df = df
        self.scale = scale
        self.norm = norm
        
    def cov(self, var, cut, bins, NCV, shapeonly=False, fillna=np.nan):
        if not isinstance(var, list):
            var = [var]
            bins = [bins]

        # Poisson variance of weighted events is square of weights
        w2 = np.asarray(self.df[self.scale])**2

        if isinstance(bins, ConcatBins):
            if shapeonly:
                _binwidths(bins)  # raises: not defined for concatenated bins
            # The same MC events fill every block, so the MC statistical
            # uncertainty is CORRELATED across them: element (a, b) of block
            # (i, j) is sum(w^2) over the events falling in bin a of var i AND
            # bin b of var j -- exactly the 2-D weighted histogram, over the
            # events passing BOTH blocks' cuts. For i == j that reduces to the
            # diagonal sum(w^2) of the 1-D case, so the one loop covers both the
            # blocks and their cross terms.
            ms = block_masks(self.df, cut, len(bins))
            vals = [np.asarray(self.df[v].fillna(fillna)) for v in var]
            edges = np.cumsum([0] + [len(b) - 1 for b in bins])
            c = np.zeros((edges[-1], edges[-1]))
            for i in range(len(bins)):
                for j in range(len(bins)):
                    m = ms[i] & ms[j]
                    c[edges[i]:edges[i+1], edges[j]:edges[j+1]] = np.histogramdd(
                        [vals[i][m], vals[j][m]], bins=[bins[i], bins[j]],
                        weights=w2[m])[0]
            return c*self.norm

        v2 = histflat(self.df, var, cut, bins, w2, fillna=fillna)*self.norm

        if shapeonly:
            diff = _binwidths(bins)
            norm = np.sum(NCV*diff)
            return np.diag(v2)/norm**2

        return np.diag(v2)

class CorrelatedSystematic(Systematic):
    def __init__(self, a, b):
        self.systa = a
        self.systb = b

        assert(self.systa.avg() == self.systb.avg())

        if (self.systa.avg() == True and self.systb.avg() == True):
            self._avg = True
        elif (self.systa.avg() == False and self.systb.avg() == False):
            self._avg = False

    def avg(self):
        return self._avg

    def nuniv(self):
        return self.systa.nuniv()

    def cov(self, var, cut, bins, NCV, shapeonly=False, fillna=np.nan):
        NCVa = NCV[:NCV.size//2]
        NCVb = NCV[NCV.size//2:]
        self.systa.cov(var, cut, bins, NCVa, shapeonly=shapeonly)
        self.systb.cov(var, cut, bins, NCVb, shapeonly=shapeonly)
        return super().cov(var, cut, bins, NCV, shapeonly=shapeonly, fillna=fillna)

    def univ(self, var, cut, bins, i_univ, fillna=np.nan):
        Na = self.systa.univ(var, cut, bins, i_univ, fillna=fillna)
        Nb = self.systb.univ(var, cut, bins, i_univ, fillna=fillna)
        N = np.concatenate((Na, Nb))
        return N

class UnCorrelatedSystematic(object):
    def __init__(self, a, b):
        self.systa = a
        self.systb = b

    def cov(self, var, cut, bins, NCV, shapeonly=False, fillna=np.nan):
        NCVa = NCV[:NCV.size//2]
        NCVb = NCV[NCV.size//2:]
        cova = self.systa.cov(var, cut, bins, NCVa, shapeonly=shapeonly, fillna=fillna)
        covb = self.systb.cov(var, cut, bins, NCVb, shapeonly=shapeonly, fillna=fillna)
        cov = np.zeros((cova.shape[0]*2, cova.shape[1]*2))
        cov[:cova.shape[0], :cova.shape[1]] = cova[:]
        cov[cova.shape[0]:, cova.shape[1]:] = covb[:]
        return cov
        
class SampleSystematic(Systematic):
    def __init__(self, dfs, cvdf=None, scale="glob_scale", norm=1):
        if not isinstance(dfs, list):
            dfs = [dfs]
        self.dfs = dfs
        self.scale = scale
        self.cvdf = cvdf
        self.norm = norm
        
    def nuniv(self):
        return len(self.dfs)

    def cov(self, var, cut, bins, NCV, shapeonly=False, fillna=np.nan):
        # compute the CV with __our__ df if configured to
        if self.cvdf is not None:
            if not isinstance(var, list):
                var = [var]
                bins = [bins]
            NCV_lcl = histflat(self.cvdf, var, cut, bins, self.cvdf[self.scale],
                               fillna=fillna)
            c = super().cov(var, cut, bins, NCV_lcl, shapeonly=shapeonly, fillna=fillna)
            # then, scale up the covariance by the ratio of our CV to the _actual_ CV
            scale = NCV/NCV_lcl
            scale[NCV_lcl == 0] = 1
            scale = np.diag(scale)
            c = scale@c@scale
            return c*self.norm**2
        else: # not overwriting the CV, just use the nominal covariance
            return super().cov(var, cut, bins, NCV, shapeonly=shapeonly, fillna=fillna)*self.norm**2
        
    def univ(self, var, cut, bins, i_univ, fillna=np.nan):
        if not isinstance(var, list):
            var = [var]
            bins = [bins]

        d = self.dfs[i_univ]
        return histflat(d, var, cut, bins, d[self.scale], fillna=fillna)

class SelectionSystematic(Systematic):
    """Systematic evaluated by re-histogramming the SAME dataframe with
    alternate boolean selection columns (one universe per column).

    With avg=True (default), an [up, dn] pair of one-sided universes is
    averaged; a single one-sided universe is treated as a symmetrized
    1-sigma variation.
    """
    def __init__(self, df, cuts, scale="glob_scale", avg=True):
        if not isinstance(cuts, list):
            cuts = [cuts]
        self.df = df
        self.cuts = cuts
        self.scale = scale
        self._avg = avg

    def nuniv(self):
        return len(self.cuts)

    def avg(self):
        return self._avg

    def univ(self, var, cut, bins, i_univ, fillna=np.nan):
        # ignores the CV cut name passed in; uses this universe's own cut column
        if not isinstance(var, list):
            var = [var]
            bins = [bins]
        # this universe's own cut column (or per-block list of them)
        c = self.cuts[i_univ]
        return histflat(self.df, var, c, bins, self.df[self.scale], fillna=fillna)

def split_tracks(df, dim, coord, runs=None):
    """Build the split-track universe for muons crossing a detector plane.

    A muon "crosses" the plane at `coord` along dimension `dim` ("x"/"y"/"z")
    when its start and end sit on opposite sides; slc_vtx_* is the established
    stand-in for the track start (the flat GUMP df has no mu_start_*). Crossing
    muons are truncated at the plane: mu_len, mu_end_*, and the muon momentum
    (via muon_range_momentum -- the same length->momentum mapping that filled
    the range-based reco momentum in the CAFs) are recomputed, and the
    downstream kinematics (nu_E_calo, del_p, del_Tp, del_phi, mu_E, mu_T) are
    recalculated with kinematics.py.

    `runs` restricts the split to those run periods (e.g. [4] for the east
    cathode, which lies outside the Run 2 muon fiducial volume).

    Returns (splitdf, crosses): the truncated copy of the crossing rows, and
    the positional boolean mask of those rows in `df`. The caller must
    re-evaluate the selection on splitdf before use.
    """

    vtx = df["slc_vtx_" + dim].to_numpy()
    end = df["mu_end_" + dim].to_numpy()
    crosses = np.isfinite(vtx) & np.isfinite(end) & \
        (np.sign(vtx - coord) != np.sign(end - coord))
    if runs is not None:
        crosses = crosses & df.Run.isin(runs).to_numpy()

    s = df.loc[crosses].copy()

    # fraction of the way along the (straight) vtx->end segment at which the
    # track crosses the plane; 0 < t < 1 by construction for crossers
    t = (coord - s["slc_vtx_" + dim]) / (s["mu_end_" + dim] - s["slc_vtx_" + dim])

    for d in "xyz":
        s["mu_end_" + d] = s["slc_vtx_" + d] + t*(s["mu_end_" + d] - s["slc_vtx_" + d])
    s["mu_len"] = t*s.mu_len

    mu_p = pd.Series(muon_range_momentum(s.mu_len.to_numpy()), index=s.index)
    recompute_kinematics(s, mu_p=mu_p)

    return s, crosses

class TrackSplittingSystematic(Systematic):
    """Systematic where muons crossing a detector plane are split with
    probability `frac`.

    The single universe is the deterministic f-weighted mixture: each crossing
    muon enters at (1 - frac) x its nominal kinematics plus frac x its split
    (truncated) kinematics, with the selection re-evaluated on the split rows
    (the `cut` column must exist in splitdf). One one-sided universe, treated
    as a symmetrized 1-sigma variation (as in SelectionSystematic).

    `crosses` is the positional boolean mask of the crossing rows in `df`
    (NOT index-based -- the notebook CV dfs carry duplicate index labels).
    Build the inputs with split_tracks().
    """
    def __init__(self, df, splitdf, crosses, frac, scale="glob_scale"):
        self.df = df
        self.splitdf = splitdf
        self.crosses = np.asarray(crosses, dtype=bool)
        self.frac = frac
        self.scale = scale

    def nuniv(self):
        return 1

    def univ(self, var, cut, bins, i_univ, fillna=np.nan):
        assert(i_univ == 0)
        if not isinstance(var, list):
            var = [var]
            bins = [bins]

        w = self.df[self.scale].to_numpy()*(1 - self.frac*self.crosses)
        N = histflat(self.df, var, cut, bins, w, fillna=fillna)

        ws = self.splitdf[self.scale].to_numpy()*self.frac
        N += histflat(self.splitdf, var, cut, bins, ws, fillna=fillna)
        return N

class WeightSystematic(Systematic):
    def __init__(self, df, wgts, avg=True, scale="glob_scale"):
        self.df = df
        self.wgts = wgts
        self._nuniv = len(wgts)
        self.scale = scale
        self._avg = avg
        
    def nuniv(self):
        return self._nuniv

    def avg(self):
        return self._avg
        
    def univ(self, var, cut, bins, i_univ, fillna=np.nan):
        if not isinstance(var, list):
            var = [var]
            bins = [bins]

        wgt_v = self.df[self.scale] * self.df[self.wgts[i_univ]].fillna(fillna)
        return histflat(self.df, var, cut, bins, wgt_v, fillna=fillna)

