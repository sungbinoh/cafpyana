"""Track-splitting correction at the three challenging ICARUS planes.

GUMPLE (sbn-rewgted-19+) port of TrackSplittingCorrection.py: cuts come from
analysis_village/gumple/gumple_cuts.py (the production selection vocabulary --
FV via cut_contained/cut_cathode flags, trk cut via cut_np & has_muon &
cut_0shwother, PID at the tuned per-detector thresholds), and file lists come
from mcdata_comparison_gumple. NB the new production also redefines the muon
candidate (longest passing track, was best chi2 ratio), so f is a genuinely
new measurement, not a re-render of the 8-21 result.

ICARUS data shows a sharp excess of reconstructed *muon track endpoints* at
three detector boundaries that MC does not reproduce at the same rate:

    Z = 0                the boundary between the two TPC readout halves in z
    East cathode         x = -210.215 cm  (cryostat 0; Run 4 only -- the
                         East-East TPC is outside the Run 2 fiducial volume)
    West cathode         x = +210.215 cm  (cryostat 1)

The physical cause is track *splitting*: a muon crossing one of these planes is
reconstructed as two tracks, so the muon candidate ends at the plane. Data
splits more often than MC does. See the existing
plots-vertex-recaf-7-27/png/ICARUS-Run{2,4}_simplecosmicrejection_mu_end_*.png
from nb/MCDataComparisonreCAF-Vertex.ipynb, where the excess is visible as one
or two bins sitting on top of an MC peak.

This script quantifies, per plane, the fraction of MC muons *crossing* that
plane which would have to be split to bring MC into agreement with data:

    f = [N_data(W) - a*N_MC(W)] / (a * D)

    W    window |mu_end_d - c| < w around the plane
    a    MC->data normalization from every selected event OUTSIDE a fixed
         +/-15 cm exclusion around all three planes
    D    at-risk MC population = muons whose reco track crosses the plane
         (sign(slc_vtx_d - c) != sign(mu_end_d - c)) plus those already
         ending inside W (i.e. already split in MC)

The uncertainty is derived from data statistics only; MC statistical and
systematic uncertainties are deliberately out of scope.

Everything is measured at two selection stages, identically for MC and data:

    simple   FV & (nu_score > 0.6) & flash cut       -- "simple cosmic rejection"
    gump     FV & CRT veto & two-prong & PID & flash -- the full GUMP selection

and the ratio f(gump)/f(simple) is reported per plane. The two stages are
*overlapping, not nested* (the GUMP chain builds the two-prong cut on the CRT
veto, not on the nu_score cut), so the ratio's uncertainty is computed with the
explicit event-level overlap between them rather than assuming independence.

Angular dependence (ICARUS Run 2). Splitting is a reconstruction effect at a
readout boundary, so how steeply a muon crosses that boundary ought to matter.
For each plane the tracks are binned into equal-population quantiles of the
angle to the plane *normal*,

    theta_n = arccos(|mu_dir_n|),   n = z for the Z=0 gap, x for the cathodes

folded into [0, 90) deg. The quantile edges are POT-weighted quantiles of
theta_n over the *at-risk MC* population, so every bin carries the same
denominator D and the per-bin errors are set purely by the in-window counts.
alpha is refitted per bin -- it varies by ~10% across bins, so using the
inclusive value would bias the end bins. The window half-width is *not* refitted
per bin; the feature width is a property of the plane, not of the angle.
Reported alongside is a flat-line chi2/ndof testing "no angular dependence".

    python analysis_village/gump/TrackSplittingCorrection_GUMPLE.py

Runs both ICARUS periods by default, each in its own subprocess so memory is
reclaimed in between, then writes a combined summary.

ICARUS Run 2 runs on the *full* on-beam stream
(ICARUS_SpringRun2BNB_FullOnBeam.df, 1.994e20 POT). The "unblind" file every
earlier run used is a 1/10 prescale of it, so this is a 10x statistics gain and
nothing else -- and it is what makes the angular measurement above possible at
all. sbn-rewgted-19 has no Run 4 equivalent, so Run 4 stays on the prescaled
stream; that asymmetry is deliberate and is logged at run time. Note OFF_w goes
0.105 -> 1.049 with it, so the w_OFF^2*N_OFF term in var(N_data) is no longer
negligible -- the formulae already carry it, but it is now the term that matters.

NB: `FV` must stay a module-level function in `__main__`. loaddf's load cache is
keyed on the preselection's `__module__ + "." + __qualname__`, and the cache was
written as "__main__.FV".

NB (2026-08-21): the two warnings that used to live here are both obsolete.
loaddf's cache works again -- `_CACHE_VERSION` is defined (loaddf.py:567) and
`load_one` defaults to `cache_dir=CACHE_WITH_INPUT` (loaddf.py:729), writing
<hash>.h5 beside the input .df -- so the first run of this script is a raw parse
and every re-run is near-instant. And the sbn-rewgted-14 .df files *do* carry
all three of loaddf.g4_syst's reinteractions_*_Geant4 knobs (164 knobs total),
so the KeyError that forced `include_syst=False` on sbn-rewgted-10 is gone.
`include_syst=False` is kept anyway, now purely for speed and memory: KEEPCOLS
drops the 100 universes immediately, so loading them is pure waste.
`reweight_aFF=True` is applied before loaddf's include_syst early return, so
`cvwgt` -- and hence `glob_scale` -- is unaffected either way.
"""

import argparse
import json
import os
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (_HERE, os.path.abspath(os.path.join(_HERE, "..", ".."))):
    if _p not in sys.path:
        sys.path.insert(0, _p)

sys.path.insert(0, os.path.join(_HERE, "..", "gumple"))
import gumple_cuts as gc
import loaddf
import mcdata_comparison_gumple as mdc

DETECTORS = ["ICARUS Run2", "ICARUS Run4"]

# Detectors read from the full-statistics on-beam stream instead of the
# 1/10-prescaled *_unblind.df. Run 2 only: sbn-rewgted-19 has no full-stats
# Run 4 file, so that period necessarily stays prescaled. Using the full stream
# here is a sideband measurement of mu_end_* (not a blinded variable), and every
# zoom plot is area-normalized, so no absolute ICARUS rate is exposed.
FULL_ONBEAM_DETECTORS = ["ICARUS Run2"]

DEFAULT_DF_DIR = "/Users/gputnam/Work/osc/sbn-rewgted-20/"
# NB: plots-tracksplit-7-27/, -8-01/ and -8-21/ hold the earlier (reCAF) runs
# and are deliberately not overwritten. Each is a separate *result*, not a
# re-render; this GUMPLE run is a new production (sbn-rewgted-14 -> -19) with
# the new candidate definitions. Do not point a re-run at an old directory.
DEFAULT_PLOTDIR = "/Users/gputnam/Work/osc/cafpyana/plots-gumple-2026-08-26/tracksplit/"

FONTSIZE = 14

# Only these columns are needed. Everything else -- in particular the 100
# systematic universes on the MC -- is dropped immediately after loading.
# p_end_* is here only because FV -> gc.pfv_cut reads it; the selection re-runs
# FV (a no-op, since it was already the load-time preselection) so that the cut
# applied here is verbatim the one in mcdata_comparison.
KEEPCOLS = ["detector", "Run", "nu_score", "flash_maxpe",
            "slc_vtx_x", "slc_vtx_y", "slc_vtx_z",
            "mu_end_x", "mu_end_y", "mu_end_z",
            # the reconstructed muon direction, for the angular binning below.
            # Without these _slim drops them and theta_normal() sees nothing.
            "mu_dir_x", "mu_dir_y", "mu_dir_z",
            "p_end_x", "p_end_y", "p_end_z",
            # FV re-evaluation on the slimmed frame (gumple_cuts):
            # sanity/slcfv need slc_vtx_* (above); the rest are stored flags.
            "cut_contained", "cut_cathode",
            # the rest of the full GUMP chain: CRT veto, trk cut, PID.
            "crthit", "cut_0shwother", "cut_np", "has_muon",
            # gc.pid_cut / get_base_muon_mask inputs
            "mu_chi2_of_mu_cand", "mu_chi2_of_prot_cand",
            "prot_chi2_of_mu_cand", "prot_chi2_of_prot_cand",
            "mu_len", "mu_trackScore", "mu_dist_start",
            "mu_prim_pfp", "mu_contained10"]

# Half-width scan grid, cm. Counting is done directly on the coordinate (no
# histogram binning), so the grid can be as fine as we like.
SCAN_W = np.arange(0.5, 20.001, 0.5)

# Half-width of the region excluded from the normalization sample around every
# plane. Fixed (not tied to the scan) so that `a` is a single number per
# detector and the scan below is not self-referential.
NORM_EXCLUDE_W = 15.0

# ---- angular binning ----
# Number of equal-population bins in the angle to the plane normal.
N_QUANTILE = 4
# Detectors the angular measurement runs for. Run 2 only: it is the only period
# with a full-statistics on-beam stream, and at the prescaled Run 4 statistics
# four angular bins would each hold ~15 in-window data events. --angular-detector
# lifts this without a code change once a Run 4 FullOnBeam file exists.
ANGULAR_DETECTORS = ["ICARUS Run2"]


# ============================================================
# The three planes
# ============================================================
def split_points():
    """(name, key, dim, coordinate, runs-it-applies-to) for each plane.

    Cathode positions are derived from the full-cryostat (Run 4) FV x-ranges,
    as in nb/MCDataComparisonreCAF-Vertex.ipynb, rather than hard-coded.
    """
    fv = gc.ICARUSRun4FVCuts
    east_ctr = 0.5*(fv["C0"]["x"]["min"] + fv["C0"]["x"]["max"])
    west_ctr = 0.5*(fv["C1"]["x"]["min"] + fv["C1"]["x"]["max"])

    return [
        # name             key       dim  coord      detectors
        ("Z = 0",          "z0",     "z", 0.0,       ["ICARUS Run2", "ICARUS Run4"]),
        # The East-East TPC was off in Run 2: ICARUSRun2FVCuts["C0"]["x"]["min"]
        # is -210.22, and trkfv_cut insets that by another 10 cm, so the east
        # cathode lies outside the Run 2 muon-end fiducial volume entirely.
        ("East Cathode",   "east",   "x", east_ctr,  ["ICARUS Run4"]),
        ("West Cathode",   "west",   "x", west_ctr,  ["ICARUS Run2", "ICARUS Run4"]),
    ]


# ============================================================
# Angle to the plane normal
# ============================================================
def theta_normal(frame, dim):
    """Angle [deg] between the reconstructed muon direction and the plane normal.

    `dim` is the plane's normal axis -- "z" for the Z=0 gap, "x" for either
    cathode. mu_dir_* is a unit vector (verified), so the component *is* the
    cosine; taking |.| folds the result into [0, 90) deg, which is what we want:
    a muon crossing the plane in +z and one crossing in -z are the same geometry
    as far as splitting is concerned.

    This is the same construction as `thetadrift` in
    nb/MCDataComparisonPID-reCAF.ipynb (cell 43), generalised to either normal.
    """
    return np.degrees(np.arccos(np.clip(np.abs(frame["mu_dir_" + dim]), 0., 1.)))


def weighted_quantile(x, w, q):
    """POT-weighted quantiles of x. Lifted from nb/MCDataComparisonPID-reCAF.ipynb
    (cell 43); the same function is `_weighted_quantile` in eres_ar23_ar25.py."""
    x = np.asarray(x)
    w = np.asarray(w)
    i = np.argsort(x)
    x, w = x[i], w[i]
    c = np.cumsum(w)
    return np.interp(np.asarray(q)*c[-1], c, x)


def flat_line_chi2(f, ferr):
    """chi2, ndof and p-value for "f is independent of angle".

    The reference is the inverse-variance-weighted mean of the bins themselves,
    so this tests *shape* only -- it says nothing about whether f is non-zero.

    A plain diagonal chi2 is right here: the bins hold disjoint sets of events
    and each one's alpha comes from a disjoint slice of the normalization
    region, so by the same disjoint-Poisson argument stage_overlap relies on,
    cov(f_i, f_j) = 0 for i != j. ndof = nbin - 1 because the mean is fitted.
    """
    f, ferr = np.asarray(f, dtype=float), np.asarray(ferr, dtype=float)
    m = np.isfinite(f) & np.isfinite(ferr) & (ferr > 0)
    f, ferr = f[m], ferr[m]
    if len(f) < 2:
        return np.nan, 0, np.nan
    mean = np.sum(f/ferr**2)/np.sum(1/ferr**2)
    chi2 = float(np.sum((f - mean)**2/ferr**2))
    ndof = len(f) - 1
    return chi2, ndof, float(stats.chi2.sf(chi2, ndof))


# ============================================================
# Preselection
# ============================================================
# NB: keep this at module scope, named FV, with a body identical to
# mcdata_comparison.FV -- loaddf's cache key depends on its __qualname__
# (see the module docstring). The flash cut is applied at selection time.
def FV(df):
    # gumple_cuts equivalent of the old slcfv & mufv & pfv & cathode chain
    # (identical body to mcdata_comparison_gumple.FV -- cache-key contract).
    return gc.sanity_cut(df) & gc.slcfv_cut(df) & df.cut_contained & df.cut_cathode


def simple_cosmic_rej(d):
    """The GUMP "simple cosmic rejection" stage, as in mcdata_comparison."""
    return FV(d) & (d.nu_score > 0.6)


# The full GUMP chain, verbatim from mcdata_comparison.{crtveto,twoprong_cut,pid_cut}.
# NB: this chain does NOT pass through simple_cosmic_rej -- the two-prong cut is
# built on the CRT veto, not on the nu_score cut. The two stages therefore
# overlap without either containing the other, which is why the stage ratio
# below needs an explicit overlap term.
def crtveto(d):
    return FV(d) & ~d.crthit


def twoprong_cut(d):
    # gumple trk_cut (cut_np & has_muon & cut_0shwother) replaces the old
    # isnan(other_shw_length) & isnan(other_trk_length) condition.
    return crtveto(d) & gc.trk_cut(d)


def gump_selection(d):
    return twoprong_cut(d) & gc.pid_cut(d)


# (key, label, cut function). Order matters: the first stage is the one the
# window half-widths are derived from (see resolve_widths), and the stage ratio
# is reported as later/first.
STAGES = [
    ("simple", "Simple Cos. Rej.", simple_cosmic_rej),
    ("gump",   "Full GUMP",        gump_selection),
]
STAGE_LABEL = dict((k, l) for (k, l, _) in STAGES)
WIDTH_STAGE = STAGES[0][0]
RATIO_STAGES = (STAGES[0][0], STAGES[1][0])


# ============================================================
# Loading
# ============================================================
def _slim(df):
    """Drop everything but the columns this analysis needs.

    Called right after each load so the 100 systematic universes on the MC do
    not stay resident. `cvwgt` is kept when present because scale_pot needs it.
    """
    cols = [c for c in KEEPCOLS if c in df.columns]
    if "cvwgt" in df.columns:
        cols = cols + ["cvwgt"]
    return df[cols].copy()


def load_all(detector, files, pot, include_dirt, log):
    """Load MC (+dirt) and ON/OFF data and POT-scale the MC.

    Returns (mcdf, ondf, offdf) after the load-time FV preselection only; the
    per-stage selections are applied by the caller, which needs the stage masks
    side by side to get their overlap. Every load below uses exactly the kwargs
    the vertex notebook used, so all of them hit the loaddf cache.
    """
    mc_files = files["MC_FILES"]
    log("loading MC: %s" % ", ".join(os.path.basename(f) for f in mc_files))
    # include_syst=False: the universes are dropped by _slim anyway, and the
    # uncached path for them crashes on these files (see the module docstring).
    # reweight_aFF is applied before loaddf's include_syst early return, so
    # cvwgt/glob_scale are unaffected. pot_univ is dropped for the same reason
    # it is now moot -- it only builds universes, which live past that return.
    mcdf, match, mcpot = loaddf.loadl(mc_files, njob=min(len(mc_files), 10),
                                      detector=detector, preselection=FV,
                                      reweight_aFF=True, include_syst=False)
    mcdf = _slim(mcdf)
    loaddf.scale_pot(mcdf, mcpot, pot)
    mcdf["dirt"] = False
    log("  MC: %d rows, POT %.4g (scale %.4g)" % (len(mcdf), mcpot, pot/mcpot))

    if include_dirt and files["DIRT_FILES"]:
        dirt_files = files["DIRT_FILES"]
        log("loading dirt: %s" % ", ".join(os.path.basename(f) for f in dirt_files))
        dirt, _, dirtpot = loaddf.loadl(dirt_files, njob=min(len(dirt_files), 10),
                                        detector=detector, preselection=FV,
                                        include_syst=False)
        dirt = _slim(dirt)
        loaddf.scale_pot(dirt, dirtpot, pot)
        dirt["dirt"] = True
        log("  dirt: %d rows, POT %.4g (scale %.4g)" % (len(dirt), dirtpot, pot/dirtpot))
        mcdf = pd.concat([mcdf, dirt], ignore_index=True)

    log("loading ON: %s" % os.path.basename(files["ONBEAM"]))
    ondf, _, _ = loaddf.load(files["ONBEAM"], load_truth=False, include_syst=False,
                             detector=detector, preselection=FV, match_Enu=False)
    ondf = _slim(ondf)

    log("loading OFF: %s" % ", ".join(os.path.basename(f) for f in files["OFFBEAM_FILES"]))
    offs = [loaddf.load(f, load_truth=False, include_syst=False, detector=detector,
                        preselection=FV, offbeampot=True, match_Enu=False)
            for f in files["OFFBEAM_FILES"]]
    offdf = pd.concat([_slim(o[0]) for o in offs], ignore_index=True)

    return (mcdf.reset_index(drop=True), ondf.reset_index(drop=True),
            offdf.reset_index(drop=True))


# ============================================================
# Counting
# ============================================================
def _wsum(df, mask):
    """POT-scaled MC yield under a mask."""
    return float(df.glob_scale[mask].sum())


def _jsonable(o):
    """json default: numpy scalars (notably float32) are not serializable."""
    if hasattr(o, "item"):
        return o.item()
    raise TypeError("not JSON serializable: %r" % (o,))


class PointCounter:
    """All the counting for one plane, for one detector.

    Data counts are ON - OFF*off_w with variance ON + off_w^2*OFF, the standard
    GUMP convention (see mcdata_comparison.make_plot_data).
    """

    def __init__(self, mcdf, ondf, offdf, off_w, alpha, alpha_var, dim, coord):
        self.off_w = off_w
        self.alpha = alpha
        self.alpha_var = alpha_var
        self.coord = coord

        end = "mu_end_" + dim
        vtx = "slc_vtx_" + dim

        self.mc_d = (mcdf[end] - coord).abs().to_numpy()
        self.on_d = (ondf[end] - coord).abs().to_numpy()
        self.off_d = (offdf[end] - coord).abs().to_numpy()
        self.mc_wgt = mcdf.glob_scale.to_numpy()

        # A muon "crosses" the plane when its start and end sit on opposite
        # sides. The flat GUMP df has no mu_start_* column; slc_vtx_* is the
        # established stand-in for the track start (gc.cathode_cut builds its
        # segments the same way). For a straight segment the sign test is
        # exactly the plane-crossing test.
        crosses = np.sign(mcdf[vtx].to_numpy() - coord) != np.sign(mcdf[end].to_numpy() - coord)
        self.n_cross = float(self.mc_wgt[crosses].sum())
        self.n_mc_total = float(self.mc_wgt.sum())

    def counts(self, w):
        """(N_data, var(N_data), N_MC, N_ON, N_OFF) inside |d| < w."""
        n_on = float((self.on_d < w).sum())
        n_off = float((self.off_d < w).sum())
        n_mc = float(self.mc_wgt[self.mc_d < w].sum())
        n_data = n_on - n_off*self.off_w
        var_data = n_on + n_off*self.off_w**2
        return n_data, var_data, n_mc, n_on, n_off

    def result(self, w):
        """Excess, at-risk denominator, split fraction and its data-stat error."""
        n_data, var_data, n_mc, n_on, n_off = self.counts(w)
        a = self.alpha

        excess = n_data - a*n_mc
        # at-risk MC population: crossers + those already ending in the window
        denom_mc = self.n_cross + n_mc
        denom = a*denom_mc

        f = excess/denom if denom > 0 else np.nan
        # f = N_data(W)/(a*D_mc) - n_mc/D_mc  ->  only the first term carries
        # the data-stat error, through N_data(W) and through a. W and the
        # normalization region are disjoint, so the two are uncorrelated.
        if denom > 0:
            var_f = var_data/denom**2 + self.alpha_var*(n_data/(a*denom))**2
        else:
            var_f = np.nan

        # Alternative denominators, reported alongside.
        f_all = excess/(a*self.n_mc_total) if self.n_mc_total > 0 else np.nan
        err_all = (np.sqrt(var_data)/(a*self.n_mc_total)) if self.n_mc_total > 0 else np.nan
        f_peak = excess/(a*n_mc) if n_mc > 0 else np.nan
        err_peak = (np.sqrt(var_data)/(a*n_mc)) if n_mc > 0 else np.nan

        return dict(
            halfwidth=w,
            N_on=n_on, N_off=n_off, N_data=n_data, N_data_err=np.sqrt(var_data),
            N_mc_raw=n_mc, N_mc_scaled=a*n_mc,
            excess=excess, excess_err=np.sqrt(var_data),
            n_cross_mc=self.n_cross, denom_mc=denom_mc, denom_scaled=denom,
            f_split=f, f_split_err=np.sqrt(var_f),
            f_per_selected=f_all, f_per_selected_err=err_all,
            f_peak_growth=f_peak, f_peak_growth_err=err_peak,
        )

    def shells(self, grid):
        """Per-shell excess and its data-stat error, for the width criterion."""
        n_on = np.array([float((self.on_d < w).sum()) for w in grid])
        n_off = np.array([float((self.off_d < w).sum()) for w in grid])
        n_mc = np.array([float(self.mc_wgt[self.mc_d < w].sum()) for w in grid])

        d_on, d_off, d_mc = np.diff(n_on), np.diff(n_off), np.diff(n_mc)
        s = (d_on - d_off*self.off_w) - self.alpha*d_mc
        s_err = np.sqrt(d_on + d_off*self.off_w**2)
        return s, s_err


def pick_from_excess(w_arr, excess, log, cap=10.0, frac=0.90, min_w=1.0, smooth=3):
    """Smallest half-width holding `frac` of the saturated integrated excess.

    See choose_halfwidth for the rationale. Split out from it so the same
    criterion can be applied to a single run's excess curve or to the sum over
    run periods.
    """
    m = w_arr <= cap
    sub, e = np.asarray(w_arr)[m], np.asarray(excess)[m]

    k = smooth
    smoothed = np.convolve(np.pad(e, (k//2, k//2), mode="edge"),
                           np.ones(k)/k, mode="valid")

    if smoothed.max() <= 0:
        w = min(3.0, cap)
        log("  WARNING: no positive excess anywhere in the scan; "
            "reporting a nominal %.1f cm window" % w)
        return w, "no-excess"

    target = frac*smoothed.max()
    for w, ee in zip(sub, e):
        if w >= min_w and ee >= target:
            return float(w), "auto"

    w = float(sub[int(np.argmax(smoothed))])
    log("  WARNING: excess never reached %.0f%% of its smoothed maximum; "
        "falling back to the peak of the smoothed excess, %.1f cm" % (100*frac, w))
    return w, "fallback"


def choose_halfwidth(counter, grid, log, cap=10.0, frac=0.90, min_w=1.0, smooth=3):
    """Smallest half-width holding `frac` of the saturated integrated excess.

    The integrated excess E(w) climbs while the window is still growing into the
    discrepant region and then flattens; the knee is the width of that region.
    Locating the knee from individual 0.5 cm shells does not work here -- with
    only a few thousand selected data events a shell holds a handful of events,
    so shell-to-shell noise swamps the tail and the chosen width ends up driven
    by fluctuations. Instead take the maximum of a 3-point-smoothed E(w) within
    `cap` cm as the saturated value and return the smallest w that reaches
    `frac` of it.

    `cap` bounds the search because track splitting is a reconstruction-scale
    effect: the split endpoint lands within a few cm of the plane, so a
    best-fit width of tens of cm would be an artifact, not a wider feature.
    """
    excess = []
    for w in grid:
        n_data, _, n_mc, _, _ = counter.counts(w)
        excess.append(n_data - counter.alpha*n_mc)
    return pick_from_excess(grid, excess, log, cap=cap, frac=frac,
                            min_w=min_w, smooth=smooth)


def resolve_widths(detectors, plotdir):
    """One half-width per plane, chosen from the excess summed over run periods.

    The width of the splitting feature is set by reconstruction, so it is a
    property of the plane -- not of the run period, and not of the selection
    stage. Letting each run pick its own would just be fitting noise (with ~5k
    selected data events per period, Run 2 and Run 4 picked 7.0 and 3.5 cm for
    the same Z=0 feature). Summing the two excess curves first roughly doubles
    the statistics behind the choice and gives both runs a common window.

    The width is taken from the WIDTH_STAGE ("simple") scan alone and reused for
    every stage. The full GUMP selection keeps only a per-cent of these events,
    so a width fitted there would be pure noise -- and a stage ratio is only
    meaningful if both stages are integrated over the same window.
    """
    curves = {}
    for det in detectors:
        path = os.path.join(plotdir, "%s_%s_scan.csv" % (det.replace(" ", "-"), WIDTH_STAGE))
        if not os.path.exists(path):
            continue
        scan = pd.read_csv(path)
        for name, grp in scan.groupby("point"):
            grp = grp.sort_values("halfwidth")
            w = grp.halfwidth.to_numpy()
            e = grp.excess.to_numpy()
            if name in curves:
                curves[name] = (w, curves[name][1] + e)
            else:
                curves[name] = (w, e)

    widths = {}
    for (name, key, _, _, _) in split_points():
        if name not in curves:
            continue
        w, e = curves[name]
        widths[key], how = pick_from_excess(w, e, lambda m: print(m, flush=True))
        print("common half-width for %-14s : %.1f cm (%s, from %s)"
              % (name, widths[key], how, " + ".join(detectors)), flush=True)
    return widths


# ============================================================
# Plots
# ============================================================
def _style(ax):
    ax.tick_params(axis='both', which='both', direction='in', length=6, width=1.5,
                   labelsize=FONTSIZE, top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)


def zoom_plot(mcdf, ondf, offdf, off_w, alpha, dim, coord, w, name, detector, stage,
              annotate=None):
    """MC (dirt stacked) vs data around one plane, area normalized.

    The analysis window is shaded, but deliberately left out of the title: the
    figure is the MC/data comparison, and the exact half-width belongs on the
    scan figure, which is where it is chosen.

    MC and data are each scaled to unit area over the plotted range, following
    mcdata_comparison.make_plot_data -- the shape comparison is the point, and an
    overall rate offset would otherwise read as a splitting excess.

    NB this is a *display* choice only and does not feed the measurement: the
    quoted split fraction still uses the whole-detector sideband alpha (fitted
    outside +/-NORM_EXCLUDE_W of all three planes). alpha is still applied here
    before normalizing so the stacked dirt fraction stays right, but it cancels
    out of the drawn curve. The peak sits *inside* the normalization range, so
    the excess drawn here is diluted by a few percent relative to the quoted f;
    read the number off the scan figure, not off this one.

    `annotate` is free text drawn top-left under the data legend -- the angular
    bin range, when this is one bin of the angular scan.
    """
    if dim == "x":
        bins = np.linspace(coord - 25, coord + 25, 51)   # 1 cm bins
        xlabel = "Muon End X [cm]"
    else:
        bins = np.linspace(coord - 50, coord + 50, 51)   # 2 cm bins
        xlabel = "Muon End Z [cm]"
    centers = 0.5*(bins[:-1] + bins[1:])
    end = "mu_end_" + dim

    is_dirt = mcdf.dirt.to_numpy()
    n_nu, _ = np.histogram(mcdf[end][~is_dirt], bins=bins,
                           weights=alpha*mcdf.glob_scale[~is_dirt])
    n_dirt, _ = np.histogram(mcdf[end][is_dirt], bins=bins,
                             weights=alpha*mcdf.glob_scale[is_dirt])
    n_mc = n_nu + n_dirt

    n_on, _ = np.histogram(ondf[end], bins=bins)
    n_off, _ = np.histogram(offdf[end], bins=bins)
    n_data = n_on - n_off*off_w
    n_err = np.sqrt(n_on + n_off*off_w**2)

    # ---- area normalization, as in mcdata_comparison.make_plot_data ----
    diff = bins[1:] - bins[:-1]
    mc_norm = float(np.sum(n_mc*diff))
    if mc_norm > 1e-5:
        n_nu, n_dirt, n_mc = n_nu/mc_norm, n_dirt/mc_norm, n_mc/mc_norm
    data_norm = float(np.sum(n_data*diff))
    if data_norm > 1e-5:
        n_data, n_err = n_data/data_norm, n_err/data_norm

    fig, (ax0, ax1) = plt.subplots(2, 1, height_ratios=[3, 1], sharex=True)

    fill = np.array([centers, centers]).T
    ax0.hist(fill, bins=bins, stacked=True, label=["MC $\\nu$ + cosmic", "Dirt"],
             color=["#1e3f54", "#43140b"], weights=np.array([n_nu, n_dirt]).T)
    line = ax0.errorbar(centers, n_data, n_err, color="black", linestyle="none", marker=".")

    for ax in (ax0, ax1):
        ax.axvspan(coord - w, coord + w, color="#d54c28", alpha=0.18, zorder=0)
        ax.axvline(coord, color="#d54c28", linestyle=":", linewidth=1.2, zorder=0)

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(n_mc > 0, n_data/n_mc, np.nan)
        ratio_err = np.where(n_mc > 0, n_err/n_mc, np.nan)
    ax1.errorbar(centers, ratio, ratio_err, color="black", linestyle="none", marker=".")
    ax1.set_ylim([0.5, 1.5])
    ax1.axhline(1, color="red", linestyle="--")
    # both curves are unit-area, so this is a shape ratio, not a rate ratio
    ax1.set_ylabel("Data / MC\n(shape)", fontsize=FONTSIZE-5)

    _style(ax0)
    _style(ax1)
    ax1.set_xlabel(xlabel, fontsize=FONTSIZE, fontweight='bold')
    ax0.set_ylabel("Area Normalized", fontsize=FONTSIZE, fontweight='bold')
    det = detector.replace(" ", "")
    ax0.set_title("$\\bf{%s}$ %s -- %s" % (det, name, STAGE_LABEL[stage]),
                  fontsize=FONTSIZE+2)

    lo, hi = ax0.get_ylim()
    ax0.set_ylim([lo, hi*1.6])
    ld = ax0.legend([line], ["Data\n(ON Beam - OFF)"], frameon=False, loc="upper left",
                    fontsize=10)
    ax0.legend(fontsize=12, loc="upper right", reverse=True)
    ax0.add_artist(ld)
    if annotate:
        ax0.text(0.02, 0.62, annotate, transform=ax0.transAxes, fontsize=11,
                 verticalalignment="top")

    plt.subplots_adjust(hspace=0.05)
    return fig


# NB: scan_plot, angular_plot, summary_plot and ratio_plot below plot fractions
# and ratios, not distributions -- there is no area to normalize. Only zoom_plot
# draws a histogram, and that is the one the area normalization applies to.


def scan_plot(scan, w_star, name, detector, stage):
    """f(w) with its data-stat band, plus the per-shell excess below."""
    fig, (ax0, ax1) = plt.subplots(2, 1, height_ratios=[2, 1], sharex=True)

    w = scan["halfwidth"].to_numpy()
    f = scan["f_split"].to_numpy()
    ferr = scan["f_split_err"].to_numpy()

    ax0.plot(w, 100*f, color="#1e3f54", linewidth=1.8)
    ax0.fill_between(w, 100*(f-ferr), 100*(f+ferr), color="#1e3f54", alpha=0.25, linewidth=0)
    ax0.axvline(w_star, color="#d54c28", linestyle="--", linewidth=1.5,
                label="chosen $w$ = %.1f cm" % w_star)
    ax0.axhline(0, color="gray", linewidth=0.8)
    ax0.set_ylabel("Split fraction [%]", fontsize=FONTSIZE-2, fontweight='bold')
    ax0.legend(fontsize=11, loc="lower right")

    s = scan["shell_excess"].to_numpy()
    serr = scan["shell_excess_err"].to_numpy()
    ax1.errorbar(w, s, serr, color="black", linestyle="none", marker=".", markersize=4)
    ax1.axvline(w_star, color="#d54c28", linestyle="--", linewidth=1.5)
    ax1.axhline(0, color="red", linestyle="--")
    ax1.set_ylabel("Shell\nexcess", fontsize=FONTSIZE-4)
    ax1.set_xlabel("Window half-width [cm]", fontsize=FONTSIZE, fontweight='bold')

    _style(ax0)
    _style(ax1)
    det = detector.replace(" ", "")
    ax0.set_title("$\\bf{%s}$ %s -- %s  half-width scan"
                  % (det, name, STAGE_LABEL[stage]), fontsize=FONTSIZE+2)
    plt.subplots_adjust(hspace=0.05)
    return fig


def angular_plot(rows, incl, name, detector, stage, dim):
    """f +/- sigma against the angle to the plane normal, one point per quantile.

    The inclusive result is drawn as a horizontal band so the bins can be read
    against it, and the flat-line chi2 tests whether the bins are consistent
    with *each other* -- it says nothing about whether f is non-zero.
    """
    fig, ax = plt.subplots(figsize=(7.5, 5))

    lo = np.array([r["theta_lo"] for r in rows])
    hi = np.array([r["theta_hi"] for r in rows])
    f = np.array([100*r["f_split"] for r in rows])
    fe = np.array([100*r["f_split_err"] for r in rows])
    ctr = 0.5*(lo + hi)

    if incl is not None:
        fi, fie = 100*incl["f_split"], 100*incl["f_split_err"]
        ax.axhspan(fi - fie, fi + fie, color="#1e3f54", alpha=0.15, zorder=0,
                   label="inclusive: $(%.2f \\pm %.2f)$ %%" % (fi, fie))
        ax.axhline(fi, color="#1e3f54", linewidth=1.2, linestyle="-", zorder=1)

    ax.errorbar(ctr, f, yerr=fe, xerr=[ctr - lo, hi - ctr], color="#d54c28",
                linestyle="none", marker="o", markersize=5, capsize=0, zorder=3)
    ax.axhline(0, color="gray", linewidth=0.8, zorder=0)

    chi2, ndof, pval = flat_line_chi2([r["f_split"] for r in rows],
                                      [r["f_split_err"] for r in rows])
    if ndof > 0:
        # right-aligned just under the legend: the upper left runs into the
        # legend box for a long label, and the lower left sits on the y=0 line
        ax.text(0.97, 0.84,
                "flat line: $\\chi^2$/ndof = %.2f/%d  (p = %.3f)" % (chi2, ndof, pval),
                transform=ax.transAxes, fontsize=12, horizontalalignment="right",
                verticalalignment="top")

    _style(ax)
    ax.set_xlim([0, 90])
    ax.set_xlabel("$\\theta_{%s}$ to plane normal [deg]" % dim,
                  fontsize=FONTSIZE, fontweight='bold')
    ax.set_ylabel("Split fraction [%]", fontsize=FONTSIZE, fontweight='bold')
    det = detector.replace(" ", "")
    ax.set_title("$\\bf{%s}$ %s -- %s" % (det, name, STAGE_LABEL[stage]),
                 fontsize=FONTSIZE+2)
    ax.legend(fontsize=11, loc="upper right")
    # headroom for the legend and the chi2 line below it
    lo_y, hi_y = ax.get_ylim()
    ax.set_ylim([lo_y, hi_y + 0.45*(hi_y - lo_y)])
    return fig


GROUPS = ["ICARUS Run2", "ICARUS Run4", "Run2+Run4"]
GROUP_COLORS = {"ICARUS Run2": "#1e3f54", "ICARUS Run4": "#d54c28",
                "Run2+Run4": "#315031"}


def summary_plot(rows, stage):
    """f +/- sigma for every point x (Run 2, Run 4, combined), for one stage."""
    points = [p[0] for p in split_points()]

    fig, ax = plt.subplots(figsize=(8, 5))
    for gi, g in enumerate(GROUPS):
        xs, ys, es = [], [], []
        for pi, p in enumerate(points):
            m = [r for r in rows
                 if r["detector"] == g and r["point"] == p and r["stage"] == stage]
            if not m:
                continue
            xs.append(pi + (gi - 1)*0.22)
            ys.append(100*m[0]["f_split"])
            es.append(100*m[0]["f_split_err"])
        if xs:
            ax.errorbar(xs, ys, es, linestyle="none", marker="o", markersize=7,
                        capsize=4, color=GROUP_COLORS[g], label=g)

    ax.set_xticks(range(len(points)))
    ax.set_xticklabels(points, fontsize=FONTSIZE)
    ax.axhline(0, color="gray", linewidth=0.8)
    ax.set_ylabel("Muons to split [%]", fontsize=FONTSIZE, fontweight='bold')
    ax.set_title("$\\bf{ICARUS}$ %s" % STAGE_LABEL[stage], fontsize=FONTSIZE+2)
    ax.legend(fontsize=12)
    _style(ax)
    return fig


def ratio_plot(ratios):
    """f(full GUMP) / f(simple cosmic rejection) per plane."""
    points = [p[0] for p in split_points()]

    fig, ax = plt.subplots(figsize=(8, 5))
    for gi, g in enumerate(GROUPS):
        xs, ys, es = [], [], []
        for pi, p in enumerate(points):
            m = [r for r in ratios if r["detector"] == g and r["point"] == p]
            if not m or not np.isfinite(m[0]["ratio"]):
                continue
            xs.append(pi + (gi - 1)*0.22)
            ys.append(m[0]["ratio"])
            es.append(m[0]["ratio_err"])
        if xs:
            ax.errorbar(xs, ys, es, linestyle="none", marker="o", markersize=7,
                        capsize=4, color=GROUP_COLORS[g], label=g)

    ax.set_xticks(range(len(points)))
    ax.set_xticklabels(points, fontsize=FONTSIZE)
    ax.axhline(1, color="red", linestyle="--", linewidth=1.2)
    ax.axhline(0, color="gray", linewidth=0.8)
    ax.set_ylabel("$f$(Full GUMP) / $f$(Simple Cos. Rej.)",
                  fontsize=FONTSIZE, fontweight='bold')
    ax.set_title("$\\bf{ICARUS}$ splitting correction: stage ratio", fontsize=FONTSIZE+2)
    ax.legend(fontsize=12)
    _style(ax)
    return fig


def _save(fig, plotdir, stem, dosave):
    if not dosave:
        plt.close(fig)
        return
    fig.savefig(os.path.join(plotdir, "png", stem + ".png"), bbox_inches="tight")
    fig.savefig(os.path.join(plotdir, "pdf", stem + ".pdf"), bbox_inches="tight")
    plt.close(fig)


# ============================================================
# Per-detector driver
# ============================================================
def run_angular(smc, son, soff, off_w, norm_sel, dim, coord, w, nq, log):
    """Split fraction in equal-population bins of the angle to the plane normal.

    Returns one row per bin. The quantile edges are POT-weighted quantiles of
    theta_normal over the *at-risk* MC (plane-crossers plus muons already ending
    in the window), so every bin carries the same denominator D and the per-bin
    errors are set purely by the in-window counts.

    alpha is refitted per bin, on the same normalization region as the inclusive
    fit but restricted to that bin. MC/data agreement is itself angle-dependent
    (alpha varies by ~10% across bins in ICARUS Run 2), so reusing the inclusive
    alpha would push that variation into the excess and fake an angular trend.

    The half-width `w` is passed in, not refitted: the width of the splitting
    feature is set by reconstruction, so it is a property of the plane. Fitting
    it per bin would be fitting noise on a quarter of the statistics.
    """
    end, vtx = "mu_end_" + dim, "slc_vtx_" + dim
    wgt = smc.glob_scale.to_numpy()
    mc_d = (smc[end] - coord).abs().to_numpy()
    crosses = (np.sign(smc[vtx].to_numpy() - coord)
               != np.sign(smc[end].to_numpy() - coord))
    inwin = mc_d < w

    if not (crosses | inwin).any():
        log("  no at-risk MC muons -- skipping the angular scan")
        return []

    th_mc = np.asarray(theta_normal(smc, dim))
    th_on = np.asarray(theta_normal(son, dim))
    th_off = np.asarray(theta_normal(soff, dim))

    # mu_dir can be NaN for a candidate with no fitted direction. Such rows fall
    # into no bin (every comparison against NaN is False), so they are dropped
    # from the angular scan while still counting in the inclusive result. Say so
    # rather than letting the bins quietly fail to sum to the inclusive number.
    ok = np.isfinite(th_mc)
    n_bad = int((crosses | inwin).sum() - ((crosses | inwin) & ok).sum())
    if n_bad:
        log("  %d of %d at-risk MC muons have no direction and are dropped from "
            "the angular scan" % (n_bad, int((crosses | inwin).sum())))
    if not ((crosses | inwin) & ok).any():
        log("  no at-risk MC muons with a direction -- skipping the angular scan")
        return []

    # The at-risk denominator is D = n_cross + n_mc(W) (PointCounter.result), in
    # which a muon that both crosses the plane and ends inside the window is
    # counted *twice* -- it is a crosser and an already-split ender at once. So
    # the quantile weights are built the same way, as the concatenation of the
    # two populations rather than their union. That makes sum_j D_j == D exactly
    # and every bin carry the same D/nq, which is the entire point of using
    # quantile edges. (The overlap is well under 1%, so union-vs-multiset barely
    # moves the edges; the multiset is chosen because it makes the identity hold
    # rather than nearly hold.)
    t_ar = np.concatenate([th_mc[crosses & ok], th_mc[inwin & ok]])
    w_ar = np.concatenate([wgt[crosses & ok], wgt[inwin & ok]])
    edges = np.concatenate([[0.0],
                            np.atleast_1d(weighted_quantile(t_ar, w_ar,
                                                            np.arange(1, nq)/nq)),
                            [90.0]])
    log("  angular bins (%d equal-population quantiles of theta_%s over the "
        "at-risk MC): %s deg" % (nq, dim, " / ".join("%.1f" % e for e in edges)))

    rows = []
    for i in range(nq):
        lo, hi = float(edges[i]), float(edges[i+1])
        # the top bin closes at 90 deg inclusive so a track exactly along the
        # plane is not silently dropped
        last = (i == nq - 1)
        def inbin(t):
            return (t >= lo) & ((t <= hi) if last else (t < hi))

        bmc, bon, boff = inbin(th_mc), inbin(th_on), inbin(th_off)

        # ---- per-bin normalization ----
        mc_r = float(wgt[bmc & norm_sel["MC"]].sum())
        on_r = float((bon & norm_sel["ON"]).sum())
        off_r = float((boff & norm_sel["OFF"]).sum())
        if mc_r <= 0 or not bmc.any():
            log("  theta %.1f-%.1f deg: empty normalization region -- skipping"
                % (lo, hi))
            continue
        a = (on_r - off_r*off_w)/mc_r
        a_var = (on_r + off_r*off_w**2)/mc_r**2
        if a <= 0:
            log("  theta %.1f-%.1f deg: non-positive alpha (%.3f) -- skipping"
                % (lo, hi, a))
            continue

        counter = PointCounter(smc[bmc], son[bon], soff[boff], off_w, a, a_var,
                               dim, coord)
        res = counter.result(w)
        res.update(theta_lo=lo, theta_hi=hi, theta_bin=i, theta_dim=dim,
                   alpha=a, alpha_err=float(np.sqrt(a_var)),
                   norm_on=on_r, norm_off=off_r, norm_mc=mc_r)
        rows.append(res)

        # D_mc, not the alpha-scaled denominator: D_mc is the quantity the
        # quantile edges make equal across bins, so printing it is the visible
        # check that the binning did what it claims.
        log("   %5.1f-%5.1f deg  alpha=%.4f  N_on=%4d  D_mc=%7.1f  "
            "f = (%6.2f +/- %5.2f) %%"
            % (lo, hi, a, res["N_on"], res["denom_mc"],
               100*res["f_split"], 100*res["f_split_err"]))

    if rows:
        # sum_j D_j must reproduce the inclusive D up to the direction-less rows
        # dropped above. This is the direct check that the quantile weights were
        # built on the right population -- if it fails, the edges are wrong.
        d_incl = float(wgt[crosses].sum() + wgt[inwin].sum())
        d_bins = float(sum(r["denom_mc"] for r in rows))
        log("  denominator check: sum_j D_j = %.1f vs inclusive D = %.1f (%+.2f%%)"
            % (d_bins, d_incl, 100*(d_bins/d_incl - 1) if d_incl > 0 else np.nan))

        chi2, ndof, pval = flat_line_chi2([r["f_split"] for r in rows],
                                          [r["f_split_err"] for r in rows])
        if ndof > 0:
            log("  flat-line (no angular dependence): chi2/ndof = %.2f/%d, p = %.4f"
                % (chi2, ndof, pval))
    return rows


def run_stage(stage, detector, args, frames, masks, norm_masks, off_w, pot,
              points, widths, log):
    """Normalization, per-plane counting, plots and results for one stage.

    Returns (results, scans, norm_info); norm_info carries the normalization-
    region counts, which the stage-ratio covariance below needs.
    """
    skey, slabel, _ = stage
    mcdf, ondf, offdf = frames
    plotdir = args.plotdir
    dosave = not args.no_save

    log("")
    log("=" * 70)
    log("STAGE: %s  [%s]" % (slabel, skey))
    log("=" * 70)
    for nm, f in (("MC", mcdf), ("ON", ondf), ("OFF", offdf)):
        log("  %s selected: %d / %d" % (nm, int(masks[nm].sum()), len(f)))

    smc = mcdf[masks["MC"]].reset_index(drop=True)
    son = ondf[masks["ON"]].reset_index(drop=True)
    soff = offdf[masks["OFF"]].reset_index(drop=True)

    # ---- normalization: everything outside +/-NORM_EXCLUDE_W of ALL planes ----
    mc_r = _wsum(mcdf, masks["MC"] & norm_masks["MC"])
    on_r = float((masks["ON"] & norm_masks["ON"]).sum())
    off_r = float((masks["OFF"] & norm_masks["OFF"]).sum())
    data_r = on_r - off_r*off_w
    alpha = data_r/mc_r
    alpha_var = (on_r + off_r*off_w**2)/mc_r**2

    log("")
    log("normalization region (outside +/-%.0f cm of all three planes):" % NORM_EXCLUDE_W)
    log("  N_ON = %.0f  N_OFF = %.0f  N_data = %.1f  N_MC = %.1f" % (on_r, off_r, data_r, mc_r))
    log("  alpha = %.4f +/- %.4f" % (alpha, np.sqrt(alpha_var)))

    norm_info = dict(stage=skey, norm_on=on_r, norm_off=off_r, norm_mc=mc_r,
                     alpha=alpha, alpha_err=float(np.sqrt(alpha_var)))

    # the normalization masks arrive on the *unselected* frames (run_detector
    # keeps them there so the stage overlap can be read off directly); restrict
    # them to the selected rows so they line up with smc/son/soff
    norm_sel = {"MC": norm_masks["MC"][masks["MC"]],
                "ON": norm_masks["ON"][masks["ON"]],
                "OFF": norm_masks["OFF"][masks["OFF"]]}

    do_angular = args.angular and detector in args.angular_detector

    results, scans, angular = [], {}, []
    for (name, key, dim, coord, _) in points:
        log("")
        log("-" * 70)
        log("%s  (%s = %.3f cm)  [%s]" % (name, dim, coord, slabel))

        counter = PointCounter(smc, son, soff, off_w, alpha, alpha_var, dim, coord)

        n_widest = counter.counts(SCAN_W[-1])[3]
        if n_widest == 0:
            log("  no selected data events within %.0f cm -- skipping" % SCAN_W[-1])
            continue

        rows = [counter.result(w) for w in SCAN_W]
        scan = pd.DataFrame(rows)
        s, s_err = counter.shells(SCAN_W)
        scan["shell_excess"] = np.append(s, np.nan)
        scan["shell_excess_err"] = np.append(s_err, np.nan)
        scan.insert(0, "point", name)
        scan.insert(0, "stage", skey)
        scans[key] = scan

        if args.scan_only:
            # first pass of the two-pass driver: the scan is all that is needed
            # to choose a width common to both run periods
            log("  scan written; width to be chosen from the combined excess")
            continue

        if args.halfwidth is not None:
            w_star, how = float(args.halfwidth), "user"
        elif key in widths:
            w_star, how = float(widths[key]), "common"
        else:
            w_star, how = choose_halfwidth(counter, SCAN_W, log)
        log("  half-width: %.1f cm (%s)" % (w_star, how))

        res = counter.result(w_star)
        res.update(detector=detector, stage=skey, stage_label=slabel,
                   point=name, key=key, dim=dim, coord=coord,
                   alpha=alpha, alpha_err=float(np.sqrt(alpha_var)),
                   norm_on=on_r, norm_off=off_r, norm_mc=mc_r,
                   halfwidth_how=how, POT=pot, OFF_w=off_w,
                   onbeam=getattr(args, "onbeam_name", ""),
                   onbeam_full=getattr(args, "onbeam_full", False))
        results.append(res)

        log("  N_ON = %.0f  N_OFF = %.0f  ->  N_data = %.1f +/- %.1f"
            % (res["N_on"], res["N_off"], res["N_data"], res["N_data_err"]))
        log("  N_MC (normalized) = %.1f" % res["N_mc_scaled"])
        log("  excess = %.1f +/- %.1f  (%.1f sigma)"
            % (res["excess"], res["excess_err"],
               res["excess"]/res["excess_err"] if res["excess_err"] > 0 else 0))
        log("  MC muons crossing the plane = %.1f ; at-risk total = %.1f"
            % (alpha*res["n_cross_mc"], res["denom_scaled"]))
        log("  SPLIT FRACTION = %.3f%% +/- %.3f%%"
            % (100*res["f_split"], 100*res["f_split_err"]))
        # how much the answer moves if the window choice moves by +/-1 cm
        neigh = ["w=%.1f: %.2f%%" % (w, 100*counter.result(w)["f_split"])
                 for w in (w_star - 1.0, w_star + 1.0) if w >= SCAN_W[0]]
        log("    (stability: %s)" % "; ".join(neigh))
        log("    (per selected muon: %.4f%% +/- %.4f%% ; MC peak growth: %.1f%% +/- %.1f%%)"
            % (100*res["f_per_selected"], 100*res["f_per_selected_err"],
               100*res["f_peak_growth"], 100*res["f_peak_growth_err"]))

        tag = detector.replace(" ", "-")
        fig = zoom_plot(smc, son, soff, off_w, alpha, dim, coord, w_star, name,
                        detector, skey)
        _save(fig, plotdir, "%s_%s_%s_muend_zoom" % (tag, skey, key), dosave)
        fig = scan_plot(scan, w_star, name, detector, skey)
        _save(fig, plotdir, "%s_%s_%s_scan" % (tag, skey, key), dosave)

        # ---- angular scan, at the inclusive half-width ----
        if do_angular:
            log("")
            log("  angular scan (%d bins of theta_%s):" % (args.nquantile, dim))
            arows = run_angular(smc, son, soff, off_w, norm_sel, dim, coord,
                                w_star, args.nquantile, log)
            for r in arows:
                r.update(detector=detector, stage=skey, stage_label=slabel,
                         point=name, key=key, dim=dim, coord=coord,
                         nquantile=args.nquantile, POT=pot, OFF_w=off_w,
                         onbeam=getattr(args, "onbeam_name", ""),
                         onbeam_full=getattr(args, "onbeam_full", False))
            angular += arows

            if arows:
                fig = angular_plot(arows, res, name, detector, skey, dim)
                _save(fig, plotdir, "%s_%s_%s_f_vs_theta" % (tag, skey, key), dosave)
                for r in arows:
                    lo, hi, bi = r["theta_lo"], r["theta_hi"], r["theta_bin"]
                    last = bi == args.nquantile - 1

                    def _sub(fr, lo=lo, hi=hi, last=last):
                        t = np.asarray(theta_normal(fr, dim))
                        return fr[(t >= lo) & ((t <= hi) if last else (t < hi))]

                    fig = zoom_plot(_sub(smc), _sub(son), _sub(soff), off_w,
                                    r["alpha"], dim, coord, w_star, name,
                                    detector, skey,
                                    annotate=r"$\theta_%s$ = %.1f$-$%.1f$^\circ$"
                                             % (dim, lo, hi))
                    _save(fig, plotdir,
                          "%s_%s_%s_muend_zoom_th%d" % (tag, skey, key, bi), dosave)

    return results, scans, norm_info, angular


def stage_overlap(ondf, offdf, masks_a, masks_b, norm_masks, off_w, dim, coord, w):
    """Data-statistical overlap between two selection stages, for one window.

    The two stages share events, so their split fractions are correlated. Both
    f's depend on data only through the in-window data count D and through the
    normalization factor alpha, and the window and the normalization region are
    spatially disjoint, so the only surviving covariances are

        cov(D_a, D_b)         = N_ON(a AND b, window) + w_OFF^2 N_OFF(a AND b, window)
        cov(alpha_a, alpha_b) = [N_ON(a AND b, R) + w_OFF^2 N_OFF(a AND b, R)]
                                / (N_MC_a(R) N_MC_b(R))

    (Poisson counts of disjoint event sets are independent, so the covariance of
    two overlapping counts is the variance of their intersection.) This function
    returns the two numerators; the alpha one is divided by the MC yields in
    combine().
    """
    end = "mu_end_" + dim
    both_on = masks_a["ON"] & masks_b["ON"]
    both_off = masks_a["OFF"] & masks_b["OFF"]

    d_on = (ondf[end] - coord).abs().to_numpy()
    d_off = (offdf[end] - coord).abs().to_numpy()

    cov_data = (float((both_on & (d_on < w)).sum())
                + off_w**2*float((both_off & (d_off < w)).sum()))
    cov_norm = (float((both_on & norm_masks["ON"]).sum())
                + off_w**2*float((both_off & norm_masks["OFF"]).sum()))
    return cov_data, cov_norm


def run_detector(detector, args):
    plotdir = args.plotdir
    dosave = not args.no_save
    widths = json.loads(args.widths) if args.widths else {}
    if dosave or args.scan_only:
        for sub in ("", "png", "pdf"):
            os.makedirs(os.path.join(plotdir, sub), exist_ok=True)

    lines = []

    def log(msg):
        print(msg, flush=True)
        lines.append(str(msg))

    log("=" * 70)
    log("TrackSplittingCorrection: %s" % detector)
    log("=" * 70)

    # ICARUS Run 2 reads the full-statistics on-beam stream; there is no Run 4
    # equivalent in sbn-rewgted-14, so Run 4 stays on the 1/10-prescaled file.
    # Log which one was used -- the asymmetry is deliberate and someone reading
    # two run periods with 10x different errors will want it said out loud.
    use_full = args.full_onbeam and detector in FULL_ONBEAM_DETECTORS
    files = mdc.detector_files(detector, args.df_dir, full_onbeam=use_full)
    # NB "partial stream", not "1/10 prescale": the 1/10 factor is a verified
    # property of ICARUS Run 2's *_unblind.df relative to its FullOnBeam file.
    # No full-stats Run 4 file exists, so Run 4's prescale factor is not known
    # here and must not be asserted.
    log("on-beam stream: %s  [%s]"
        % (os.path.basename(files["ONBEAM"]),
           "FULL STATS" if use_full else "partial stream"))
    args.onbeam_name = os.path.basename(files["ONBEAM"])
    args.onbeam_full = bool(use_full)

    norm = mdc.compute_norm(detector, files, log)
    # cast out of numpy float32 -- gate_delta is single precision, and letting
    # that propagate through every count would both lose precision and make the
    # results unserializable
    off_w, pot = float(norm["OFF_w"]), float(norm["POT"])

    mcdf, ondf, offdf = load_all(detector, files, pot, not args.no_dirt, log)
    frames = (mcdf, ondf, offdf)

    points = [p for p in split_points() if detector in p[4]]
    skipped = [p for p in split_points() if detector not in p[4]]
    for (name, _, _, coord, _) in skipped:
        log("skipping %s (x = %.3f cm): outside the %s fiducial volume"
            % (name, coord, detector))

    # ---- normalization region: outside +/-NORM_EXCLUDE_W of ALL planes ----
    def outside_all(df):
        m = pd.Series(True, index=df.index)
        for (_, _, dim, coord, _) in split_points():
            m &= (df["mu_end_" + dim] - coord).abs() >= NORM_EXCLUDE_W
        return m.to_numpy()

    norm_masks = {"MC": outside_all(mcdf), "ON": outside_all(ondf),
                  "OFF": outside_all(offdf)}

    # ---- per-stage selection masks, kept on the *unselected* frames so the
    # stage overlap (needed for the ratio) can be read off directly ----
    stages = [s for s in STAGES if s[0] in args.stage]
    stage_masks = {}
    for (skey, _, sfn) in stages:
        stage_masks[skey] = {
            nm: (sfn(f) & gc.flash_cut(f)).to_numpy()
            for nm, f in (("MC", mcdf), ("ON", ondf), ("OFF", offdf))}

    results, scans, norms, angular = [], {}, {}, []
    for stage in stages:
        r, s, ni, ar = run_stage(stage, detector, args, frames, stage_masks[stage[0]],
                                 norm_masks, off_w, pot, points, widths, log)
        results += r
        scans[stage[0]] = s
        norms[stage[0]] = ni
        angular += ar

    tag = detector.replace(" ", "-")

    if args.scan_only:
        wrote = []
        for skey, bykey in scans.items():
            if not bykey:
                continue
            path = os.path.join(plotdir, "%s_%s_scan.csv" % (tag, skey))
            pd.concat(bykey.values(), ignore_index=True).to_csv(path, index=False)
            wrote.append(os.path.basename(path))
        if not wrote:
            log("no usable points for %s" % detector)
            return 1
        log("")
        log("wrote %s to %s" % (", ".join(wrote), plotdir))
        return 0

    if not results:
        log("no usable points for %s" % detector)
        return 1

    # ---- stage-overlap terms for the f(gump)/f(simple) ratio ----
    xrows = []
    a, b = RATIO_STAGES
    if a in stage_masks and b in stage_masks:
        log("")
        log("-" * 70)
        log("stage overlap (%s vs %s), for the split-fraction ratio:"
            % (STAGE_LABEL[a], STAGE_LABEL[b]))
        for (name, key, dim, coord, _) in points:
            ra = [r for r in results if r["stage"] == a and r["key"] == key]
            rb = [r for r in results if r["stage"] == b and r["key"] == key]
            if not ra or not rb:
                continue
            w = ra[0]["halfwidth"]
            cov_data, cov_norm = stage_overlap(ondf, offdf, stage_masks[a],
                                               stage_masks[b], norm_masks, off_w,
                                               dim, coord, w)
            cov_alpha = cov_norm/(norms[a]["norm_mc"]*norms[b]["norm_mc"])
            xrows.append(dict(detector=detector, point=name, key=key,
                              stage_a=a, stage_b=b, halfwidth=w,
                              cov_data=cov_data, cov_alpha=cov_alpha))
            log("  %-14s w=%.1f cm : cov(N_data) = %.2f  cov(alpha) = %.3e"
                % (name, w, cov_data, cov_alpha))

    if dosave:
        pd.DataFrame(results).to_csv(
            os.path.join(plotdir, "%s_tracksplit.csv" % tag), index=False)
        for skey, bykey in scans.items():
            if bykey:
                pd.concat(bykey.values(), ignore_index=True).to_csv(
                    os.path.join(plotdir, "%s_%s_scan.csv" % (tag, skey)), index=False)
        with open(os.path.join(plotdir, "%s_tracksplit.json" % tag), "w") as f:
            json.dump(results, f, indent=2, default=_jsonable)
        with open(os.path.join(plotdir, "%s_xstage.json" % tag), "w") as f:
            json.dump(xrows, f, indent=2, default=_jsonable)
        if angular:
            pd.DataFrame(angular).to_csv(
                os.path.join(plotdir, "%s_tracksplit_angular.csv" % tag), index=False)
            with open(os.path.join(plotdir,
                                   "%s_tracksplit_angular.json" % tag), "w") as f:
                json.dump(angular, f, indent=2, default=_jsonable)
        with open(os.path.join(plotdir, "%s_summary.txt" % tag), "w") as f:
            f.write("\n".join(lines) + "\n")
        log("")
        log("wrote %s_{tracksplit.csv,<stage>_scan.csv,tracksplit.json,"
            "xstage.json,summary.txt} to %s" % (tag, plotdir))

    return 0


# ============================================================
# Combination across run periods
# ============================================================
def stage_ratio(parts_a, parts_b, xparts):
    """f_b / f_a and its data-statistical uncertainty, for one plane.

    `parts_*` are the per-run result rows for the two stages (one entry per run
    period, already restricted to the same plane); `xparts` are the matching
    overlap rows. Both stages' f is the same function of data:

        f = sum_r(D_r) / T  -  sum_r(alpha_r M_r) / T ,   T = sum_r alpha_r C_r

    so df/dD_r = 1/T and df/dalpha_r = -D_r/(alpha_r T), and the covariance
    between the stages is the sum over runs of the two overlap terms (different
    run periods are independent data sets).
    """
    def agg(parts):
        excess = sum(r["excess"] for r in parts)
        denom = sum(r["denom_scaled"] for r in parts)
        var = sum(r["excess_err"]**2 for r in parts)
        var_a = sum((r["alpha_err"]*r["N_data"]/r["alpha"])**2 for r in parts)
        return excess/denom, (var + var_a)/denom**2, denom

    fa, va, Ta = agg(parts_a)
    fb, vb, Tb = agg(parts_b)

    cov = 0.0
    for x in xparts:
        ra = [r for r in parts_a if r["detector"] == x["detector"]]
        rb = [r for r in parts_b if r["detector"] == x["detector"]]
        if not ra or not rb:
            continue
        ra, rb = ra[0], rb[0]
        cov += x["cov_data"]/(Ta*Tb)
        cov += (x["cov_alpha"]*(ra["N_data"]/ra["alpha"])*(rb["N_data"]/rb["alpha"])
                / (Ta*Tb))

    if fa == 0 or not np.isfinite(fa):
        return np.nan, np.nan, cov

    ratio = fb/fa
    var_ratio = ratio**2*(vb/fb**2 + va/fa**2 - 2*cov/(fa*fb)) if fb != 0 else np.nan
    return ratio, np.sqrt(var_ratio) if var_ratio >= 0 else np.nan, cov


def combine(detectors, args):
    """Per-run rows plus a combined Run2+Run4 row for each plane and stage,
    then the f(full GUMP)/f(simple cosmic rejection) ratio.

    Each run's MC is POT-scaled to its own data, so the combination is done in
    data units: f = sum(excess) / sum(alpha*D). The East cathode combines over
    Run 4 alone.
    """
    plotdir = args.plotdir
    rows, xrows = [], []
    for det in detectors:
        tag = det.replace(" ", "-")
        path = os.path.join(plotdir, "%s_tracksplit.json" % tag)
        if not os.path.exists(path):
            print("combine: missing %s -- skipping %s" % (path, det), flush=True)
            continue
        with open(path) as f:
            rows += json.load(f)
        xpath = os.path.join(plotdir, "%s_xstage.json" % tag)
        if os.path.exists(xpath):
            with open(xpath) as f:
                xrows += json.load(f)

    if not rows:
        print("combine: nothing to combine", flush=True)
        return 1

    stages = [s for s in STAGES if any(r["stage"] == s[0] for r in rows)]

    combined = []
    for (name, key, dim, coord, _) in split_points():
        for (skey, slabel, _) in stages:
            parts = [r for r in rows if r["point"] == name and r["stage"] == skey]
            if not parts:
                continue
            excess = sum(r["excess"] for r in parts)
            var = sum(r["excess_err"]**2 for r in parts)
            denom = sum(r["denom_scaled"] for r in parts)
            n_data = sum(r["N_data"] for r in parts)
            n_mc = sum(r["N_mc_scaled"] for r in parts)
            # alpha uncertainty: propagate each run's, weighted by its data yield
            var_a = sum((r["alpha_err"]*r["N_data"]/r["alpha"])**2 for r in parts)
            f = excess/denom
            ferr = np.sqrt(var/denom**2 + var_a/denom**2)
            combined.append(dict(
                detector="Run2+Run4", stage=skey, stage_label=slabel,
                point=name, key=key, dim=dim, coord=coord,
                runs=", ".join(r["detector"] for r in parts),
                halfwidth=", ".join("%.1f" % r["halfwidth"] for r in parts),
                N_data=n_data, N_mc_scaled=n_mc, excess=excess, excess_err=np.sqrt(var),
                denom_scaled=denom, f_split=f, f_split_err=ferr))

    allrows = rows + combined
    out = pd.DataFrame(allrows)
    out.to_csv(os.path.join(plotdir, "tracksplit_summary.csv"), index=False)

    # ---- stage ratio, per run period and combined ----
    a, b = RATIO_STAGES
    ratios = []
    if any(r["stage"] == a for r in rows) and any(r["stage"] == b for r in rows):
        for (name, key, dim, coord, _) in split_points():
            groups = [(d, [d]) for d in detectors] + [("Run2+Run4", detectors)]
            for (label, dets) in groups:
                pa = [r for r in rows if r["point"] == name and r["stage"] == a
                      and r["detector"] in dets]
                pb = [r for r in rows if r["point"] == name and r["stage"] == b
                      and r["detector"] in dets]
                if not pa or not pb:
                    continue
                xp = [x for x in xrows if x["point"] == name and x["detector"] in dets]
                ratio, rerr, cov = stage_ratio(pa, pb, xp)
                fa = sum(r["excess"] for r in pa)/sum(r["denom_scaled"] for r in pa)
                fb = sum(r["excess"] for r in pb)/sum(r["denom_scaled"] for r in pb)
                ratios.append(dict(detector=label, point=name, key=key,
                                   f_a=fa, f_b=fb, ratio=ratio, ratio_err=rerr,
                                   cov=cov, n_overlap_rows=len(xp)))
        if ratios:
            pd.DataFrame(ratios).to_csv(
                os.path.join(plotdir, "tracksplit_stage_ratio.csv"), index=False)
            _save(ratio_plot(ratios), plotdir, "tracksplit_stage_ratio", True)

    for (skey, _, _) in stages:
        _save(summary_plot(allrows, skey), plotdir,
              "tracksplit_summary_%s" % skey, True)

    W = 92
    lines = ["", "=" * W,
             "TRACK-SPLITTING CORRECTION  (fraction of plane-crossing MC muons to split)",
             "=" * W]
    for (skey, slabel, _) in stages:
        lines += ["", "STAGE: %s" % slabel, "-" * W,
                  "%-14s %-15s %7s %10s %12s %18s" % (
                      "Point", "Sample", "w [cm]", "excess", "at-risk MC",
                      "split fraction"),
                  "-" * W]
        for (name, _, _, _, _) in split_points():
            for r in allrows:
                if r["point"] != name or r["stage"] != skey:
                    continue
                w = r["halfwidth"] if isinstance(r["halfwidth"], str) \
                    else "%.1f" % r["halfwidth"]
                lines.append("%-14s %-15s %7s %10.1f %12.1f   (%5.2f +/- %4.2f) %%" % (
                    name, r["detector"], w, r["excess"], r["denom_scaled"],
                    100*r["f_split"], 100*r["f_split_err"]))
            lines.append("-" * W)

    if ratios:
        lines += ["", "=" * W,
                  "STAGE RATIO  f(%s) / f(%s)" % (STAGE_LABEL[b], STAGE_LABEL[a]),
                  "(data-stat only; the two stages share events, and the overlap "
                  "is propagated)",
                  "=" * W,
                  "%-14s %-15s %14s %14s %18s" % (
                      "Point", "Sample", "f simple [%]", "f gump [%]", "ratio"),
                  "-" * W]
        for (name, _, _, _, _) in split_points():
            for r in ratios:
                if r["point"] != name:
                    continue
                lines.append("%-14s %-15s %14.2f %14.2f   %8.2f +/- %.2f" % (
                    name, r["detector"], 100*r["f_a"], 100*r["f_b"],
                    r["ratio"], r["ratio_err"]))
            lines.append("-" * W)

    # ---- angular results, merged across detectors ----
    # Concatenated only, never combined across run periods: the quantile edges
    # are fitted per detector/plane/stage, so a Run2+Run4 angular row would be
    # summing over bins that are not the same bins. (Moot today -- only Run 2
    # produces these -- but it is the first thing a Run 4 enablement will trip on.)
    ang = []
    for det in detectors:
        apath = os.path.join(plotdir,
                             "%s_tracksplit_angular.json" % det.replace(" ", "-"))
        if os.path.exists(apath):
            with open(apath) as f:
                ang += json.load(f)
    if ang:
        pd.DataFrame(ang).to_csv(os.path.join(plotdir, "tracksplit_angular.csv"),
                                 index=False)
        W2 = 96
        lines += ["", "=" * W2,
                  "SPLIT FRACTION vs ANGLE TO THE PLANE NORMAL",
                  "(equal-population bins of the at-risk MC; alpha refitted per bin)",
                  "=" * W2]
        for det in detectors:
            for (skey, slabel, _) in STAGES:
                for (name, _, dim, _, _) in split_points():
                    rows_ = [r for r in ang if r["detector"] == det
                             and r["stage"] == skey and r["point"] == name]
                    if not rows_:
                        continue
                    rows_.sort(key=lambda r: r["theta_bin"])
                    chi2, ndof, pval = flat_line_chi2(
                        [r["f_split"] for r in rows_],
                        [r["f_split_err"] for r in rows_])
                    lines.append("-" * W2)
                    lines.append("%s  %s  [%s]   w = %.1f cm"
                                 % (det, name, slabel, rows_[0]["halfwidth"]))
                    lines.append("%18s %8s %10s %12s %20s"
                                 % ("theta_%s [deg]" % dim, "N_ON", "alpha",
                                    "D (MC)", "split fraction"))
                    for r in rows_:
                        lines.append("%8.1f - %7.1f %8.0f %10.4f %12.1f   (%6.2f +/- %5.2f) %%"
                                     % (r["theta_lo"], r["theta_hi"], r["N_on"],
                                        r["alpha"], r["denom_mc"],
                                        100*r["f_split"], 100*r["f_split_err"]))
                    if ndof > 0:
                        lines.append("%18s  chi2/ndof = %.2f/%d,  p = %.4f"
                                     % ("flat line:", chi2, ndof, pval))
        lines.append("-" * W2)

    text = "\n".join(lines)
    print(text, flush=True)
    with open(os.path.join(plotdir, "tracksplit_summary.txt"), "w") as f:
        f.write(text + "\n")
    print("\nwrote tracksplit_summary.{csv,txt}, tracksplit_summary_<stage>.{png,pdf}, "
          "tracksplit_stage_ratio.{csv,png,pdf}%s to %s"
          % (" and tracksplit_angular.csv" if ang else "", plotdir), flush=True)
    return 0


# ============================================================
# Entry point
# ============================================================
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-d", "--detector", action="append", default=None,
                   help="Detector to run: 'ICARUS Run2', 'ICARUS Run4' or 'all'. "
                        "Repeatable. Default: all.")
    p.add_argument("--df-dir", default=DEFAULT_DF_DIR, help="Directory of input .df files")
    p.add_argument("--plotdir", default=DEFAULT_PLOTDIR, help="Output directory")
    p.add_argument("--halfwidth", type=float, default=None,
                   help="Override the automatic window half-width [cm] for every point")
    p.add_argument("--widths", default=None,
                   help="JSON map of point key -> half-width [cm], e.g. "
                        "'{\"z0\": 3.5, \"east\": 2.0}'. Set by the driver's second "
                        "pass so both run periods share a window.")
    p.add_argument("--stage", action="append", default=None,
                   help="Selection stage to run: %s, or 'all'. Repeatable. "
                        "Default: all." % ", ".join(repr(k) for (k, _, _) in STAGES))
    p.add_argument("--scan-only", action="store_true",
                   help="Write only the half-width scan CSV (driver's first pass)")
    p.add_argument("--angular", action="append", default=None,
                   help="Detector to run the angular (theta-to-plane-normal) "
                        "variant for: 'ICARUS Run2', 'ICARUS Run4', 'all' or "
                        "'none'. Repeatable. Default: %s."
                        % ", ".join(ANGULAR_DETECTORS))
    p.add_argument("--nquantile", type=int, default=N_QUANTILE,
                   help="Number of equal-population angular bins [default %d]"
                        % N_QUANTILE)
    p.add_argument("--no-full-onbeam", dest="full_onbeam", action="store_false",
                   help="Use the 1/10-prescaled *_unblind.df on-beam stream for "
                        "ICARUS Run 2 as well, i.e. reproduce the pre-8-21 "
                        "statistics. Useful for isolating a sample change from "
                        "the statistics change.")
    p.add_argument("--no-dirt", action="store_true", help="Skip the dirt sample")
    p.add_argument("--no-save", action="store_true", help="Compute but write nothing")
    p.add_argument("--combine-only", action="store_true",
                   help="Skip the per-detector pass; just rebuild the combined summary")
    args = p.parse_args(argv)

    if not args.detector:
        args.detector = ["all"]

    dets = []
    for d in args.detector:
        if d == "all":
            dets += DETECTORS
        elif d in DETECTORS:
            dets.append(d)
        else:
            p.error("unknown detector %r (choose from %s, or 'all')" % (d, DETECTORS))
    args.detector = dets

    keys = [k for (k, _, _) in STAGES]
    if not args.stage:
        args.stage = ["all"]
    sts = []
    for s in args.stage:
        if s == "all":
            sts += keys
        elif s in keys:
            sts.append(s)
        else:
            p.error("unknown stage %r (choose from %s, or 'all')" % (s, keys))
    # keep STAGES order, drop duplicates
    args.stage = [k for k in keys if k in sts]

    # --angular: same 'all'/'none' resolution as --detector
    if args.angular is None:
        ang = list(ANGULAR_DETECTORS)
    else:
        ang = []
        for d in args.angular:
            if d == "all":
                ang += DETECTORS
            elif d == "none":
                ang = []
                break
            elif d in DETECTORS:
                ang.append(d)
            else:
                p.error("unknown --angular detector %r (choose from %s, 'all' or "
                        "'none')" % (d, DETECTORS))
    # only meaningful for detectors actually being run
    args.angular_detector = [d for d in ang if d in args.detector]
    args.angular = bool(args.angular_detector)
    if ang and not args.angular_detector:
        print("NB: --angular names %s, none of which is being run (%s); the "
              "angular variant is off." % (ang, args.detector), flush=True)

    if args.nquantile < 2:
        p.error("--nquantile must be at least 2 (got %d)" % args.nquantile)

    if WIDTH_STAGE not in args.stage and args.halfwidth is None and not args.widths:
        p.error("the %r stage sets the window half-widths; either include it, or "
                "pass --halfwidth/--widths" % WIDTH_STAGE)

    return args


def main(argv=None):
    args = parse_args(argv)

    if args.combine_only:
        return combine(args.detector, args)

    if len(args.detector) == 1:
        return run_detector(args.detector[0], args)

    # More than one detector. Each period runs in its own subprocess so the
    # multi-GB MC frames are reclaimed in between, and the whole thing runs
    # twice: pass 1 produces the half-width scans, the driver picks one window
    # per plane from the run-summed excess, and pass 2 redoes the numbers and
    # plots with that common window. Loading dominates the runtime and comes
    # straight from the loaddf cache, so the second pass is cheap.
    passthrough = ["--df-dir", args.df_dir, "--plotdir", args.plotdir]
    stage_flags = [x for s in args.stage for x in ("--stage", s)]
    if args.halfwidth is not None:
        passthrough += ["--halfwidth", str(args.halfwidth)]
    for flag in ("no_dirt", "no_save"):
        if getattr(args, flag):
            passthrough.append("--" + flag.replace("_", "-"))
    if not args.full_onbeam:
        passthrough.append("--no-full-onbeam")
    # Safe to pass to both passes: pass 1 is --scan-only, and the angular block
    # sits after run_stage's scan_only early return, so it cannot fire there.
    passthrough += ["--nquantile", str(args.nquantile)]
    passthrough += ([x for d in args.angular_detector for x in ("--angular", d)]
                    if args.angular_detector else ["--angular", "none"])

    def run_pass(extra, label):
        rc = 0
        for det in args.detector:
            cmd = ([sys.executable, os.path.abspath(__file__), "--detector", det]
                   + passthrough + extra)
            print("\n>>> [%s] %s" % (label, " ".join(repr(c) if " " in c else c
                                                     for c in cmd)), flush=True)
            r = subprocess.run(cmd)
            if r.returncode != 0:
                print("!!! %s FAILED (exit %d)" % (det, r.returncode), flush=True)
                rc = r.returncode
        return rc

    if args.halfwidth is not None:
        # an explicit width was given: no need to derive one
        rc = run_pass(stage_flags, "single pass")
    else:
        # pass 1 only needs the width-setting stage; every stage then reuses the
        # windows it picks (see resolve_widths).
        rc = run_pass(["--scan-only", "--stage", WIDTH_STAGE], "pass 1/2: scans")
        if rc != 0:
            return rc
        print("\n" + "=" * 70, flush=True)
        print("choosing one half-width per plane from the run-summed excess "
              "(%s stage)" % STAGE_LABEL[WIDTH_STAGE], flush=True)
        print("=" * 70, flush=True)
        widths = resolve_widths(args.detector, args.plotdir)
        rc = run_pass(["--widths", json.dumps(widths)] + stage_flags,
                      "pass 2/2: results")

    if rc == 0 and not args.no_save:
        rc = combine(args.detector, args)
    return rc


if __name__ == "__main__":
    sys.exit(main())
