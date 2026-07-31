import numpy as np

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

def recompute_kinematics(s, mu_p=None, BE=None):
    """Recompute the downstream reco kinematics (nu_E_calo, del_p, del_Tp,
    del_phi, mu_E, mu_T) on the flat df `s`, in place.

    mu_p defaults to sqrt(mu_E^2 - m_mu^2) -- exact, since the reco momentum
    is range-based; BE defaults to kinematics.BE. Pass either to build a
    shifted-kinematics universe from the CV.
    """
    import pandas as pd
    import kinematics

    if BE is None:
        BE = kinematics.BE
    if mu_p is None:
        mu_p = pd.Series(np.sqrt(np.maximum(s.mu_E.to_numpy()**2 - kinematics.MUON_MASS**2, 0)),
                         index=s.index)

    mu_dir = pd.DataFrame({"x": s.mu_dir_x, "y": s.mu_dir_y, "z": s.mu_dir_z})
    p_dir = pd.DataFrame({"x": s.p_dir_x, "y": s.p_dir_y, "z": s.p_dir_z})
    p_p = np.sqrt(np.maximum(s.p_E**2 - kinematics.PROTON_MASS**2, 0))

    tki = kinematics.transverse_kinematics(mu_p, mu_dir, p_p, p_dir, BE=BE)
    s["nu_E_calo"] = kinematics.neutrino_energy(mu_p, mu_dir, p_p, p_dir, BE=BE)
    s["del_p"] = tki["del_p"]
    s["del_Tp"] = tki["del_Tp"]
    s["del_phi"] = tki["del_phi"]
    s["mu_E"] = tki["mu_E"]
    s["mu_T"] = tki["mu_E"] - kinematics.MUON_MASS

    return s

def shift_binding_energy(df, dBE, fraction=1.0, scale="glob_scale"):
    """Universe df for a binding-energy shift: copy of the CV with the reco
    kinematics recomputed under BE -> BE + dBE, scaled per interaction mode:
    the nominal dBE for QE/RES/DIS (genie_mode 0/1/2), dBE*sqrt(2) for MEC
    (genie_mode 10, two-nucleon initial state), and no shift for COH and
    non-neutrino rows (genie_mode 3/NaN -- no struck bound nucleon).

    `fraction` applies the shift to only that fraction of events, as a
    deterministic mixture (as in TrackSplittingSystematic): the returned df
    holds the unshifted rows at (1 - fraction) x `scale` alongside the shifted
    rows at fraction x `scale`.

    The selection cuts use no recomputed column, so 'selected' etc. carry over
    from the CV unchanged.
    """
    import pandas as pd
    import kinematics

    mode = df.genie_mode.to_numpy()
    mode_scl = np.where(np.isin(mode, [0, 1, 2]), 1.,
               np.where(mode == 10, np.sqrt(2.), 0.))

    s = df.copy()
    recompute_kinematics(s, BE=kinematics.BE + dBE*mode_scl)

    if fraction == 1.0:
        return s

    # rebuild the unshifted half with the nominal BE rather than copying the
    # CV columns: the input df may be slimmed to just the recompute inputs
    # (exact -- the recompute round-trips the CV kinematics)
    cv = recompute_kinematics(df.copy())
    cv[scale] = cv[scale]*(1 - fraction)
    s[scale] = s[scale]*fraction
    return pd.concat([cv, s])

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
            diff = outern([b[1:] - b[:-1] for b in bins])
            norm = np.sum(NCV*diff)
            if norm > 1e-5:
                NCV = NCV / norm
        
        N_univ = []
        for i_univ in range(self.nuniv()):
            N = self.univ(var, cut, bins, i_univ, fillna=fillna)
            if shapeonly:
                diff = outern([b[1:] - b[:-1] for b in bins])
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

        return np.histogramdd([self.df.loc[self.df[cut], v].fillna(fillna) for v in var], bins=bins, weights=self.df.loc[self.df[cut], self.scale])[0].flatten()*self.norm + self.CV

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
        w = self.df.loc[self.df[cut], self.scale]**2
        var = np.histogramdd([self.df.loc[self.df[cut], v].fillna(fillna) for v in var], bins=bins, weights=w)[0].flatten()*self.norm
        return np.diag(var)

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
            NCV_lcl = np.histogramdd([self.cvdf.loc[self.cvdf[cut], v].fillna(fillna) for v in var], bins=bins, weights=self.cvdf.loc[self.cvdf[cut], self.scale])[0].flatten()
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

        return np.histogramdd([self.dfs[i_univ].loc[self.dfs[i_univ][cut], v].fillna(fillna) for v in var], bins=bins, weights=self.dfs[i_univ].loc[self.dfs[i_univ][cut], self.scale])[0].flatten()

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
        c = self.cuts[i_univ]
        return np.histogramdd([self.df.loc[self.df[c], v].fillna(fillna) for v in var],
                              bins=bins,
                              weights=self.df.loc[self.df[c], self.scale])[0].flatten()

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
    import pandas as pd

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

        m = self.df[cut].to_numpy()
        w = self.df[self.scale].to_numpy()*(1 - self.frac*self.crosses)
        N = np.histogramdd([self.df[v].fillna(fillna).to_numpy()[m] for v in var],
                           bins=bins, weights=w[m])[0].flatten()

        ms = self.splitdf[cut].to_numpy()
        ws = self.splitdf[self.scale].to_numpy()[ms]*self.frac
        N += np.histogramdd([self.splitdf[v].fillna(fillna).to_numpy()[ms] for v in var],
                            bins=bins, weights=ws)[0].flatten()
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
        return np.histogramdd([self.df.loc[self.df[cut], v].fillna(fillna) for v in var], bins=bins, weights=wgt_v[self.df[cut]])[0].flatten()

