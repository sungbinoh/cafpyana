import numpy as np

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

        wgt_v = self.df[self.scale] * self.df[self.wgts[i_univ]]
        return np.histogramdd([self.df.loc[self.df[cut], v].fillna(fillna) for v in var], bins=bins, weights=wgt_v[self.df[cut]])[0].flatten()

