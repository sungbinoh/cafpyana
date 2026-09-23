"""Exact port of the MAPLE chi2_ALG PID, run on the CAF-stored dE/dx.

Replicates helper_eff_cf...h chi2_ALG() bit-for-bit in algorithm terms:
  - drop the first and last calorimetry point of the track
  - keep points with rr_min < rr < rr_max and non-NaN dedx/rr
  - per point, require 0.5 <= dedx <= 100 (C++: skip if dedx > 100 or < 0.5)
  - template bin lookup a la TProfile::FindBin, with the neighbor-averaging
    patch applied when bin content or error < 1e-6 (using *original* neighbor
    values, as the C++ per-lookup patch does)
  - sigma_dedx = (0.04231 + 0.0001783 dedx^2) * dedx
  - chi2 = sum[ ((dedx - template)/sqrt(err_t^2 + sigma^2))^2 ] / npt
  - npt == 0 -> NaN (C++ 0/0)

The template file is byte-identical between the CAFANA macro's copy
(/exp/icarus/app/users/marterop/dev_areas/dEdxrestemplates.root) and the
larsoft_data copy on cvmfs (verified by md5), so we use whichever is present.
"""
import os
import numpy as np
import pandas as pd
import uproot

TEMPLATE_FILES = [
    "/exp/icarus/app/users/marterop/dev_areas/dEdxrestemplates.root",
    "/cvmfs/larsoft.opensciencegrid.org/products/larsoft_data/v1_02_02/ParticleIdentification/dEdxrestemplates.root",
]

def _load_templates():
    for fname in TEMPLATE_FILES:
        if os.path.exists(fname):
            f = uproot.open(fname)
            out = {}
            for part, hname in [("mu", "dedx_range_mu"), ("pro", "dedx_range_pro"),
                                ("ka", "dedx_range_ka"), ("pi", "dedx_range_pi")]:
                prof = f[hname]
                vals = prof.values()
                errs = prof.errors(error_mode="s")
                edges = prof.axis().edges()
                # neighbor-average patch, using original neighbor values.
                # Pad with 0 (ROOT under/overflow default) for edge bins.
                vpad = np.concatenate([[0.], vals, [0.]])
                epad = np.concatenate([[0.], errs, [0.]])
                vals_patched = np.where(vals < 1e-6, (vpad[:-2] + vpad[2:]) / 2., vals)
                errs_patched = np.where(errs < 1e-6, (epad[:-2] + epad[2:]) / 2., errs)
                out[part] = (edges, vals_patched, errs_patched)
            return out
    raise FileNotFoundError("dEdxrestemplates.root not found in any known location")

_TEMPLATES = _load_templates()
# the C++ code uses the *proton* template binning for FindBin on all particles
_BIN_EDGES = _TEMPLATES["pro"][0]
_NBINS = len(_BIN_EDGES) - 1


def chi2_cafana(hitdf, dedxname="dedx", rr_min=0.0, rr_max=25.0):
    """Compute [chi2_mu, chi2_pro, chi2_ka, chi2_pi] per track from stored dedx.

    hitdf: output of make_trkhitdf (plane-2 calo points) with columns
           `dedxname`, rr, firsthit, lasthit. Index (.., trk, point).
    Returns a DataFrame indexed per-track with columns
           chi2_mu, chi2_pro, chi2_ka, chi2_pi, npt (NaN chi2 when npt == 0).
    """
    dedx = hitdf[dedxname]
    rr = hitdf.rr

    # chi2_ALG hit selection. NaN comparisons are False, as in C++.
    keep = (~hitdf.firsthit) & (~hitdf.lasthit) \
        & (rr < rr_max) & (rr > rr_min) & dedx.notna() & rr.notna() \
        & ~(dedx > 100.) & ~(dedx < 0.5)

    grp_levels = list(range(hitdf.index.nlevels - 1))

    sel = hitdf[keep]
    seldedx = sel[dedxname].values
    selrr = sel.rr.values

    # TProfile::FindBin equivalent (1-based); require 1 <= bin <= nbins
    binidx = np.searchsorted(_BIN_EDGES, selrr, side="right")
    inrange = (binidx >= 1) & (binidx <= _NBINS)

    errdedx = (0.04231 + 0.0001783 * seldedx**2) * seldedx

    out = {}
    for part in ["mu", "pro", "ka", "pi"]:
        _, vals, errs = _TEMPLATES[part]
        binc = np.where(inrange, vals[np.clip(binidx - 1, 0, _NBINS - 1)], np.nan)
        bine = np.where(inrange, errs[np.clip(binidx - 1, 0, _NBINS - 1)], np.nan)
        v = (seldedx - binc)**2 / (bine**2 + errdedx**2)
        s = pd.Series(np.where(inrange, v, np.nan), index=sel.index)
        g = s[inrange].groupby(level=grp_levels)
        out["chi2_" + part] = g.sum() / g.size()

    # npt counts hits passing the dedx window AND the bin-range requirement
    npt = pd.Series(inrange.astype(int), index=sel.index)
    out["npt"] = npt[inrange].groupby(level=grp_levels).size()

    res = pd.DataFrame(out)
    return res
