"""MC/data comparison plots for the GUMP selection, on the GUMPLE dataframes.

GUMPLE (sbn-rewgted-19+) port of mcdata_comparison.py. Differences from the
reCAF version:
  * cuts come from analysis_village/gumple/gumple_cuts.py -- the PRODUCTION
    selection vocabulary: FV = sanity & slcfv (incl. the ICARUS WW cable
    exclusion) & cut_contained & cut_cathode; cosmic rejection = per-detector
    opening-angle threshold + CRT veto (NO nu_score cut); two-prong =
    gc.trk_cut (cut_np & has_muon & cut_0shwother); PID = gc.pid_cut at the
    tuned per-detector thresholds. Staged distributions therefore differ from
    the reCAF (-14) results by construction.
  * mu_track_score -> mu_trackScore, is_cosmic -> true_iscosmic
  * other_trk_length / other_shw_length are no longer stored; the plots show
    othr_pfp_length / n_shower / n_other instead
  * the BE universe carries the psum_* / n_pfp columns the migrated
    syst.recompute_kinematics requires
  * -19 file layout: split MC CV shards, per-run ICARUS dirt, and no
    *_Overlay_Nom.df detvar CVs -- detvars pair against the rewgt CV instead

Runs the full plot set for one detector ("SBND", "ICARUS Run2" or
"ICARUS Run4"):

  * 5 selection stages x 27 variables, broken down by interaction mode.
    The PID stage is split on proton multiplicity: "PID 1p" (n_pfp == 2,
    the GUMP selection) and "PID Np" (n_pfp > 2, the MAPLE selection).
  * 5 selection stages x 4 chi2 variables, broken down by true PDG

  2026-08-25: the muon-MOMENTUM and neutrino-energy comparisons (mu_p,
  nu_E_calo, nu_E_ccqe, nu_E_frac_diff -- including the del_p-sliced energy
  round) are REMOVED from this script; they live in mcdata_constraint_gumple.py
  (SBND-only, blinded to the 1p selection), which reads their plot specs from
  CONSTRAINT_VARS here. The muon angles (mu_costh, mu_phi), del_p and the
  mu-p opening angle remain. Nothing is absolute-normalized unless --absnorm
  (SBND-only) is passed explicitly.

both with the full systematic covariance, matching the set used by the signal box
(nb/SignalBoxSystematics-ReCAF.ipynb): flux, G4, xsec weights + binding-energy
shift, POT normalization, detector variations, calorimetry/chi2 variations,
trigger (flash-PE scale), SBND cosmic normalization, ICARUS track splitting, MC
stat, and OOAV/dirt.

The one deliberate difference from the signal box: no beam-off stat systematic,
because the off-beam statistical uncertainty is already propagated into the data
error bars in make_plot_data().

Both rounds are area (shape) normalized by default. --absnorm repeats them with
absolute (POT) normalization, writing *_absnorm figures and absolute chi2s.
BLINDING: that flag is restricted to SBND (ABSNORM_DETECTORS) -- an absolute
comparison exposes the data rate, and ICARUS is still blinded.

--blinded draws the MC, its systematic band and the chi2 / ndof / p-value, but
omits the data points from both the main and the ratio panel. It does NOT relax
the --absnorm restriction: the printed chi2 is still a data/MC comparison.

Deliverables per detector, in --plotdir: the pdf/ + png/ figures, one
<det>_chi2_{mode,pdg}_{shape,abs}.tsv table of chi2 / ndof / p-value per
variable, <det>_plotdata.pkl, <det>_chi2.csv and <det>_summary.txt.

    python analysis_village/gump/mcdata_comparison.py --detector all
    python analysis_village/gump/mcdata_comparison.py --detector SBND --absnorm

Passing more than one --detector runs each in its own subprocess so memory is
reclaimed in between.

NB: `FV` must stay a module-level function in `__main__`. loaddf's load cache is
keyed on the preselection's `__module__ + "." + __qualname__`, and the existing
cache was written from the notebook as "__main__.FV". Moving or renaming it
invalidates every cache entry.
"""

import argparse
import os
import pickle
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from multiprocess import Pool
from scipy import stats
from tqdm.auto import tqdm

import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (_HERE, os.path.abspath(os.path.join(_HERE, "..", ".."))):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import pyanalib.pandas_helpers as ph
from makedf.util import *

import kinematics
sys.path.insert(0, os.path.join(_HERE, "..", "gumple"))
import gumple_cuts as gc
import loaddf
import syst

DETECTORS = ["ICARUS Run2", "ICARUS Run4", "SBND"]

DEFAULT_DF_DIR = "/Users/gputnam/Work/osc/sbn-rewgted-20/"
DEFAULT_PLOTDIR = "/Users/gputnam/Work/osc/cafpyana/plots-gumple-2026-08-26/evtsel/"

# ============================================================
# BLINDING: absolute (POT) normalization
# ============================================================
# The area-normalized comparison only ever shows shapes. The absolute-normalized
# one puts data and MC on a common POT scale and so exposes the overall data
# rate -- it must ONLY be run for SBND. ICARUS is still blinded: do not add
# "ICARUS Run2" / "ICARUS Run4" here until ICARUS is unblinded.
ABSNORM_DETECTORS = ["SBND"]

# ============================================================
# BLINDING: per-detector variable blacklist
# ============================================================
# ICARUS is not unblinded in the reconstructed neutrino energy: no data/MC
# comparison of it may be produced -- no figure, no chi2 row, no pickled
# histogram. mu_p is on the list because the CCQE inversion turns the muon
# momentum straight back into a neutrino energy. Everything else is unaffected:
# p_p, del_p, mu_costh, mu_phi, mu_p_opening_angle_deg and the rest still plot.
#
# Hard-coded on purpose -- there is deliberately NO command-line override.
# Empty a detector's tuple only when that detector is unblinded.
# 2026-08-25: the previously blinded variables (nu_E_calo, nu_E_ccqe,
# nu_E_frac_diff, mu_p) are REMOVED from this script's comparisons entirely --
# they live in mcdata_constraint_gumple.py, which is SBND-only and blinded to
# the 1p selection. The per-detector blinding machinery below is kept (empty)
# so a future blinded variable can be added without re-plumbing.
BLIND_NU_VARS = ()
BLIND_VARS = {
    "SBND":        (),
    "ICARUS Run2": BLIND_NU_VARS,
    "ICARUS Run4": BLIND_NU_VARS,
}


def blind_vars(detector):
    """Variables that must not be compared to data for `detector`.

    Fails closed: an unrecognised detector string gets the full blind list, so a
    typo or a newly added sample blinds rather than leaks.
    """
    return tuple(BLIND_VARS.get(detector, BLIND_NU_VARS))


def visible_vars(detector, plotvars):
    """`plotvars` minus what is blinded for `detector`, order preserved."""
    blind = blind_vars(detector)
    return [v for v in plotvars if v not in blind]


def visible_spec(detector, plotvars, bins, labels):
    """Index-aligned (plotvars, bins, labels) with the blinded entries dropped.

    Applied to run_detector's *local* copies. The module-level
    PLOTVARS/LABELS/plot_bins() triple is never mutated: PLOTVARS.index()
    lookups (the del_p binning below, mcdata_constraint._plotvar_spec) and the
    be_cols column list all depend on it, and the systematic universes must
    still carry the blinded columns even when nothing plots them.
    """
    blind = blind_vars(detector)
    keep = [i for i, v in enumerate(plotvars) if v not in blind]
    return ([plotvars[i] for i in keep], [bins[i] for i in keep],
            [labels[i] for i in keep])


# ============================================================
# Globals read by the Pool workers.
#
# make_plot_data() looks these up at module scope, exactly as it did at notebook
# scope. multiprocess uses the fork start method here, so workers inherit them
# without pickling the (multi-GB) frames per task.
# ============================================================
DETECTOR = None
df = None
ONdf = None
OFFdf = None
OFF_w = None
POT = None
systematics_by_cut = None


# ============================================================
# Input files
# ============================================================
# Full-statistics on-beam streams, for the detectors that have one. The ONBEAM
# defaults below are the prescaled/partial streams; detector_files(...,
# full_onbeam=True) swaps in the entry here. ICARUS Run 2's *_unblind.df is
# exactly a 1/10 prescale of its FullOnBeam file (same gate_delta distribution,
# 1.985e19 vs 1.994e20 POT). There is deliberately no ICARUS Run 4 entry --
# sbn-rewgted-19 ships no full-stats Run 4 file.
FULL_ONBEAM = {
    "SBND":        "SBND_SpringBNBData_FullOnBeam.df",
    "ICARUS Run2": "ICARUS_SpringRun2BNB_FullOnBeam.df",
}


def detector_files(detector, df_dir, full_onbeam=False):
    """Input file lists for one detector.

    DETVAR_FILES entries are (detector-to-load-with, (CV file list, variation
    file list, ...)).

    full_onbeam swaps ONBEAM for the full-statistics stream in FULL_ONBEAM, and
    raises for a detector that has none rather than silently handing back the
    prescaled file. Off by default: this module's own results are quoted on the
    prescaled stream and flipping the default would move them without notice.
    TrackSplittingCorrection.py opts in for ICARUS Run 2.

    NB the *(1.-1/200.) Run2 gate correction in compute_norm is a *separate*
    effect (gate loss, not prescale) and still applies either way.
    """
    D = df_dir

    if detector == "SBND":
        onbeam = D + "SBND_SpringBNBData_FixedDev.df"
        offbeam_files = [D + "SBND_SpringBNBOffData.df"]

        mc_files = [D + "SBNDMCCV_%i.df" % i for i in range(3)]
        dirt_files = [D + "SBND_SpringLowEMC.df"]

        detvar_files = [
            ("SBND", ([D + "SBND_SpringMC_Nom.df"],
                      [D + "SBND_SpringMC_0xSCE.df"],
                      [D + "SBND_SpringMC_2xSCE.df"])),
            ("SBND", ([D + "SBND_SpringMC_Nom.df"],
                      [D + "SBND_SpringMC_DENT.df"])),
            ("SBND", ([D + "SBND_SpringMC_WMNom.df"],
                      [D + "SBND_SpringMC_WMXThetaXW.df"])),
            ("SBND", ([D + "SBND_SpringMC_WMNom.df"],
                      [D + "SBND_SpringMC_WMYZ.df"])),
        ]
        detvar_names = ["SCE", "DENT", "WM $X\\theta_{xw}$", "WM $YZ$"]

    elif detector == "ICARUS Run2":
        onbeam = D + "ICARUS_SpringRun2BNB_unblind.df"
        offbeam_files = [D + "ICARUS_SpringRun2BNBOff_unblind.df"]

        mc_files = [D + "ICARUSRun2_SpringMCOverlay_rewgt_%i.df" % i for i in range(3)]
        dirt_files = [D + "ICARUSRun2_Spring_Overlay_Dirt.df"]

        # NB: sbn-rewgted-19 ships no ICARUSRun2_Spring_Overlay_Nom.df, so the
        # detvar CV partner is the rewgt CV sample (same underlying overlay
        # events); match_common_evts restricts both sides to the shared events.
        detvar_files = [
            ("ICARUS Run2", (mc_files,
                             [D + "ICARUSRun2_Spring_Overlay_SCE.df"])),
            ("ICARUS Run2", (mc_files,
                             [D + "ICARUSRun2_Spring_Overlay_WMXThXW.df"])),
        ]
        detvar_names = ["SCE", "WM $X\\theta_{xw}$"]

    elif detector == "ICARUS Run4":
        onbeam = D + "ICARUS_SpringRun4BNB_unblind.df"
        offbeam_files = [D + "ICARUS_SpringRun4BNBOff_unblind.df"]

        mc_files = [D + "ICARUSRun4_SpringMCOverlay_rewgt_%i.df" % i for i in range(4)]
        dirt_files = [D + "ICARUSRun4_Spring_Overlay_Dirt.df"]

        # NB: no ICARUSRun4_Spring_Overlay_Nom.df in sbn-rewgted-19 either --
        # same rewgt-CV pairing as Run 2 above.
        detvar_files = [
            ("ICARUS Run4", (mc_files,
                             [D + "ICARUSRun4_Spring_Overlay_SCE.df"])),
            # NB: the notebook loaded this one with "ICARUS Run2", which applied
            # the Run2 MC flash-PE scale (0.632) instead of Run4's (0.358).
            ("ICARUS Run4", (mc_files,
                             [D + "ICARUSRun4_Spring_Overlay_WMXThXW.df"])),
        ]
        detvar_names = ["SCE", "WM $X\\theta_{xw}$"]

    else:
        raise ValueError("Unknown detector: %r" % detector)

    if full_onbeam:
        if detector not in FULL_ONBEAM:
            raise ValueError(
                "no full-statistics on-beam stream for %r (have: %s)"
                % (detector, ", ".join(sorted(FULL_ONBEAM))))
        onbeam = D + FULL_ONBEAM[detector]

    return dict(ONBEAM=onbeam, OFFBEAM_FILES=offbeam_files, MC_FILES=mc_files,
                DIRT_FILES=dirt_files, DETVAR_FILES=detvar_files,
                DETVAR_NAMES=detvar_names)


def read_dfs(file, key):
    with h5py.File(file, "r") as f:
        keys = [k for k in f.keys() if k.startswith(key)]
    return pd.concat([pd.read_hdf(file, k) for k in keys])


# ============================================================
# Normalization: gate counts, POT, off-beam weight
# ============================================================
def compute_norm(detector, files, log):
    onbeam = files["ONBEAM"]
    offbeam_files = files["OFFBEAM_FILES"]

    if detector == "ICARUS Run2":
        # NB: the old *(1-1/100.) on-beam prescale correction is dropped -- the new
        # "unblind" files are the full (non-prescaled) stream.
        ngates_ON = read_dfs(onbeam, "trig").gate_delta.sum()*(1.-1/200.)
        ngates_OFF = sum(read_dfs(f, "trig").gate_delta.sum() for f in offbeam_files)*(1-1/20.)

        off_w = ngates_ON / ngates_OFF
    elif detector == "ICARUS Run4":
        # NB: the old *(1-1/100.) on-beam prescale correction is dropped -- the new
        # "unblind" files are the full (non-prescaled) stream.
        ngates_ON = read_dfs(onbeam, "trig").gate_delta.sum()*(1.-1/40.)
        ngates_OFF = sum(read_dfs(f, "trig").gate_delta.sum() for f in offbeam_files)*(1-1/20.)

        off_w = ngates_ON / ngates_OFF
    elif "SBND" in detector:
        ngates_ON = read_dfs(onbeam, "bnb").shape[0]
        ngates_OFF = sum(read_dfs(f, "hdr").noffbeambnb.sum() for f in offbeam_files)

        f_factor = 0.0754
        off_w = (1. - f_factor) * (ngates_ON) / (ngates_OFF)

    log("ngates_ON = %r, ngates_OFF = %r, OFF_w = %r" % (ngates_ON, ngates_OFF, off_w))

    if "ICARUS" in detector:
        log("gate_delta mean: ON: %r  OFF: %r" % (
            1/read_dfs(onbeam, "trig").gate_delta.mean(),
            [1/read_dfs(f, "trig").gate_delta.mean() for f in offbeam_files]))

    if "ICARUS" in detector:
        pot = read_dfs(onbeam, "hdr").pot.sum()
    elif "SBND" in detector:
        pot = read_dfs(onbeam, "bnb").TOR875.sum()
        log("TOR875 / 1e19: %r" % (pot / 1e19))

    log("POT = %r" % pot)
    log("N GATES ON / 5e12 POT")
    log("%r" % (5e12*ngates_ON/pot))

    nevt_ON = read_dfs(onbeam, "hdr").shape[0]
    nevt_OFF = sum(read_dfs(f, "hdr").shape[0] for f in offbeam_files)

    nevt = nevt_ON - nevt_OFF*off_w

    log("NEVT_ON = %r, POT = %r, NEVT_ON/1e15 POT = %r" % (nevt_ON, pot, nevt_ON / (pot / 1e15)))
    log("NEVT = %r, POT = %r, NEVT/1e15 POT = %r" % (nevt, pot, nevt / (pot / 1e15)))

    return dict(ngates_ON=ngates_ON, ngates_OFF=ngates_OFF, OFF_w=off_w, POT=pot,
                NEVT_ON=nevt_ON, NEVT_OFF=nevt_OFF, NEVT=nevt)


# ============================================================
# Preselection
# ============================================================
# NB: the flash cut is NOT applied in the load-time preselection so that the
# trigger (flash-PE scale) systematic can vary the threshold in both directions;
# it is applied at selection time (the per-stage cut columns) instead.
#
# NB: keep this function at module scope, named FV -- loaddf's cache key depends
# on its __qualname__ (see the module docstring).
def FV(df):
    # gumple_cuts equivalent of the old slcfv & mufv & pfv & cathode chain:
    # cut_contained is the production's candidate-containment flag (mufv & pfv
    # analog) and cut_cathode its cathode cut; slcfv_cut now also excludes the
    # ICARUS WW dangling-cable region.
    return gc.sanity_cut(df) & gc.slcfv_cut(df) & df.cut_contained & df.cut_cathode


# ============================================================
# Detector variation loading
# ============================================================
def load_detvar(detvar, detector):
    """Load one detvar group and match the frames to their common events."""
    dfs, matches, pots = zip(*[loaddf.loadl(flist, preselection=FV, include_syst=False, detector=detector)
                               for flist in detvar])
    dfs, pots = loaddf.match_common_evts(matches, dfs, pots)
    # The ICARUS detvar CV partner is now the rewgt CV sample (the -19
    # production ships no *_Overlay_Nom.df); an empty match here means the two
    # samples do not actually share events and the pairing is wrong.
    if any(len(d) == 0 for d in dfs):
        raise RuntimeError(
            "detvar matching produced an empty frame -- the CV and variation "
            "samples share no events (FLAG: rewgt-CV pairing invalid; the "
            "*_Overlay_Nom.df detvar CVs need to be produced for sbn-rewgted-19)")
    return dfs, pots


# ============================================================
# PID/calorimetry variations
# ============================================================
def v_variation(d, setvars):
    d = d[[c for c in d.columns if "univ" not in c]].copy()
    for (new, old) in setvars:
        d[new] = d[old]
    return d


def v_calovariation(d, variation):
    setvars = [
        ("mu_chi2_of_mu_cand", "mu_chi2%s_of_mu_cand" % variation),
        ("mu_chi2_of_prot_cand",  "mu_chi2%s_of_prot_cand" % variation),
        ("prot_chi2_of_mu_cand", "prot_chi2%s_of_mu_cand" % variation),
        ("prot_chi2_of_prot_cand",  "prot_chi2%s_of_prot_cand" % variation),
    ]
    return v_variation(d, setvars)


chi2vars = [
    "mu_chi2_of_mu_cand",
    "prot_chi2_of_mu_cand",
    "mu_chi2_of_prot_cand",
    "prot_chi2_of_prot_cand",
]

# Calorimetry variations, per detector, matching the signal box. The tuned
# selection (see the header comment in gump_cuts.py) is "no EMB IC, no sqsmear15,
# double SBND 13": ICARUS takes no EMB (alpha/beta/R) variation, neither detector
# takes sqsmear15, and the SBND dE/dx constant-resolution variation counts double.
# The columns exist for both detectors -- this is a choice, not a data limitation.
CHI2_VARIATIONS = {
    "SBND":   ["smear13", "hi", "alpha_p", "beta_p", "R_p"],
    "ICARUS": ["smear13", "hi"],
}

CHI2_NORM = {
    "SBND":   {"smear13": 2.0},
    "ICARUS": {},
}

CHI2_DETVAR_NAMES = {
    "smear13": "dE/dx Constant Res.",
    "hi":      "Calo. Gain",
    "alpha_p": "EMB $\\alpha$",
    "beta_p":  "EMB $\\beta$",
    "R_p":     "EMB $R$",
}


def chi2_config(detector):
    """(variation list, {variation: covariance norm}) for one detector."""
    key = "SBND" if detector == "SBND" else "ICARUS"
    return CHI2_VARIATIONS[key], CHI2_NORM[key]


# ============================================================
# Truth categorization
# ============================================================
def TrueAV(d):
    vtx = pd.DataFrame({
        "detector": d.detector,
        "Run": d.Run,
        "x": d.true_vtx_x,
        "y": d.true_vtx_y,
        "z": d.true_vtx_z,
    }, index=d.index)
    return gc._fv_cut(vtx, 0, 0, 0, 0)


def OOAV(d):
    return ~np.isnan(d.true_vtx_x) & ~TrueAV(d)


def nuFV(d):
    vtx = pd.DataFrame({
        "detector": d.detector,
        "Run": d.Run,
        "x": d.true_vtx_x,
        "y": d.true_vtx_y,
        "z": d.true_vtx_z,
    }, index=d.index)
    return ~np.isnan(d.true_vtx_x) & gc._fv_cut(vtx, 0, 0, 0, 0)


mode_list = [1, 10, 0]

mode_labels = ["Cosmic", "Dirt", 'Other $\\nu$', '$\\nu_\\mu$ CC RES', '$\\nu_\\mu$ CC MEC', '$\\nu_\\mu$ CC QE']
mode_colors = ["#95af8b", "#43140b", "#c89648", "#1e3f54", "#d54c28", "#315031"]


def breakdown_mode(var, d):
    """Break down variable by interaction mode."""
    numu_cc = (np.abs(d.true_pdg) == 14) & (d.mu_true_p > 0.)
    # NB: the notebook computes `fid = nuFV(d)` here and then immediately
    # overrides it with `fid = ~d.dirt`, so only the latter is ever used. Kept
    # as a comment so the alternative is still visible without paying for a
    # DataFrame build per call that is thrown away.
    # fid = nuFV(d)
    fid = ~d.dirt

    ret = [
        var[np.isnan(d.genie_mode)],
        var[~fid & ~np.isnan(d.genie_mode)],
        var[(~np.any([d.genie_mode == i for i in mode_list], axis=0) | ~numu_cc) & fid & ~np.isnan(d.genie_mode)]
    ] +\
        [var[(d.genie_mode == i) & numu_cc & fid] for i in mode_list]

    return ret


pdg_list = [2212, 13, 211]
pdg_labels = ["$p$", "$\\mu$", "$\\pi^\\pm$", "Other"]
pdg_colors = ["#315031", "#d54c28", "#1e3f54", "#c89648"]


def breakdown_pdg(var, d, particle="p"):
    ret = [var[np.abs(d["%s_true_pdg" % particle] == i)] for i in pdg_list]
    ret.append(var[sum([np.abs(d["%s_true_pdg" % particle] == i) for i in pdg_list]) == 0])
    return ret


def breakdown_pdg_p(var, d):
    return breakdown_pdg(var, d, "p")


def breakdown_pdg_mu(var, d):
    return breakdown_pdg(var, d, "mu")


# ============================================================
# TRIGGER SYSTEMATIC: flash-PE scale-factor uncertainty
# ============================================================
# Best-fit MC PE scale factors from the data/MC fits in FlashMCDataComparison.ipynb
# (applied to flash_maxpe at load time in loaddf.py):
#   SBND: 0.642 +/- 0.005   ICARUS Run2: 0.632 +/- 0.024   ICARUS Run4: 0.358 +/- 0.017
#
# Varying the scale s -> s*(1 -+ f), f = unc/s, is equivalent to varying the
# threshold on flash_maxpe: T -> T/(1 -+ f). The flash cut is applied in the
# per-stage cut columns (NOT the load-time preselection), so both the
# tightened- and loosened-threshold universes are available.
TRIG_PE_SCALE     = {1: 0.642, 2: 0.632, 4: 0.358}  # Run -> best-fit s (matches loaddf.py)
TRIG_PE_SCALE_UNC = {1: 0.005, 2: 0.024, 4: 0.017}  # Run -> unc. on s
# Run -> gump_cuts.flash_cut threshold. NB: SBND used to sit at -1 here, which is
# not the cut gump_cuts applies (2000): both universes then degenerated to "no
# flash cut at all" and blew up the SBND detector systematic.
TRIG_THRESHOLD    = {1: 2000., 2: 5000., 4: 1000.}


def trig_flash_cut_var(d, updn, runs=None):
    """Flash cut with the threshold varied by the +/-1 sigma PE-scale uncertainty.

    updn = +1: scale DOWN 1 sigma -> threshold UP   (fewer events)
    updn = -1: scale UP   1 sigma -> threshold DOWN (more events)
    runs: restrict the variation to these runs (others stay at the CV threshold).
    """
    sel = pd.Series(False, index=d.index)
    for run, T in TRIG_THRESHOLD.items():
        f = TRIG_PE_SCALE_UNC[run] / TRIG_PE_SCALE[run]
        Tvar = T / (1 - updn*f) if (runs is None or run in runs) else T
        sel = sel | ((d.Run == run) & (d.flash_maxpe > Tvar))
    return sel


# ============================================================
# BINDING-ENERGY SYSTEMATIC (SBND + ICARUS)
# ============================================================
# +25 MeV shift of the binding energy in the reco neutrino-energy formula
# (kinematics.BE = 29.5 MeV), built as a shifted-CV universe df and treated as an
# additional cross-section systematic: a single one-sided universe, symmetrized.
# The shift is per interaction mode -- nominal for QE/RES/DIS, x sqrt(2) for MEC,
# zero for COH / non-neutrino rows -- and is applied to half the events
# (fraction=0.5), as in the signal box.
#
# The selection cuts read none of the recomputed columns, so the per-stage cut
# columns carry over from the CV unchanged. del_p is the only *plotted* variable
# the recompute touches; everything else rides along.
BE_SHIFT = 0.025

# Inputs syst.shift_binding_energy / recompute_kinematics need. The universe df is
# slimmed to these plus the plot variables and the cut columns: the CV frame
# carries ~200 weight columns that the universe copy does not need.
# NB: recompute_kinematics is migrated to the GUMPLE psum convention -- it
# requires the stored psum_E / psum_dir_* (and n_pfp) and raises without them.
BE_RECOMPUTE_COLS = ["glob_scale", "genie_mode", "mu_E", "p_E",
                     "mu_dir_x", "mu_dir_y", "mu_dir_z",
                     "p_dir_x", "p_dir_y", "p_dir_z",
                     "psum_E", "psum_p",
                     "psum_dir_x", "psum_dir_y", "psum_dir_z", "n_pfp",
                     # stored kinematics: the unshifted rows of the in-place
                     # BE universe keep these CV values verbatim
                     "nu_E_calo", "nu_E_ccqe", "del_p", "del_Tp", "del_phi"]


# ============================================================
# TRACK-SPLITTING SYSTEMATIC (ICARUS)
# ============================================================
# Data shows an excess of muon-track endpoints at Z=0 and the cathodes that MC
# does not reproduce: muons crossing those planes are split by reconstruction more
# often in data. Split fractions measured by TrackSplittingCorrection.py (see
# TrackSplittingCorrection.md); the Run2+Run4 combined value per plane is applied
# as a one-sided, symmetrized 1-sigma variation, uncorrelated between planes.
# The plane definitions, fractions and run assignments all come from loaddf
# (SPLIT_REGIONS / SPLIT_FRAC / SPLIT_RUNS), so they stay in one place.
#
# NB: the east cathode lies outside the Run 2 muon fiducial volume, so it has no
# Run 2 crossers and contributes nothing to the Run 2 covariance.
TRACKSPLIT_FRAC = loaddf.SPLIT_FRAC

# Cosmic normalization, from the off-beam / intime cosmic MC comparison
# (nb/OffBeamCosmicMCComparison.ipynb). SBND only: the ICARUS MC is an overlay and
# already carries data cosmics, so it takes no separate cosmic normalization.
SBND_COSMIC_NORM = 0.107


# ============================================================
# Selection stages
# ============================================================
# The GUMPLE cosmic rejection: mu-p opening angle at the tuned per-detector
# thresholds (gc.SBND_CUTS / gc.ICARUS_CUTS) with the CRT veto folded in.
# NB: unlike the reCAF chain there is no nu_score cut.
def cosmic_rej(d):
    return FV(d) & gc.cosmic_cut(d)


def twoprong_cut(d):
    # gumple trk_cut = cut_np & has_muon & cut_0shwother
    return cosmic_rej(d) & gc.trk_cut(d)


def pid_cut(d):
    return twoprong_cut(d) & gc.pid_cut(d)


# The PID stage split on proton-candidate multiplicity: exactly one proton
# (the GUMP 1u1p selection) and more than one (the MAPLE 1uNp selection).
# n_pfp counts candidate pfps with the muon included.
def pid_cut_1p(d):
    return pid_cut(d) & (d.n_pfp == 2)


def pid_cut_np(d):
    return pid_cut(d) & (d.n_pfp > 2)


# ============================================================
# Plotting
# ============================================================
FONTSIZE = 14
HAWKS_COLORS = ["#315031", "#d54c28", "#1e3f54", "#c89648", "#43140b", "#95af8b"]


def add_style(ax, xlabel, title="", det="ICARUS"):
    ax.tick_params(axis='both', which='both', direction='in', length=6, width=1.5, labelsize=FONTSIZE, top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    ax.set_xlabel(xlabel, fontsize=FONTSIZE, fontweight='bold')
    ax.set_ylabel('Area Normalized', fontsize=FONTSIZE, fontweight='bold')
    ax.set_title(f"$\\bf{{{det}}}$  {title}", fontsize=FONTSIZE+2)
    ax.legend(fontsize=FONTSIZE)


def syst_band(ax, bins, lo, hi, **kw):
    """Hatched systematic band covering every bin.

    fill_between(bins[:-1], ..., step="post") silently drops the last bin:
    step="post" holds each y from x[i] to x[i+1], so the final value has no
    edge left to extend to. Step over the full edge array with the last value
    repeated instead (the same trick the overlay step line below uses).
    """
    style = dict(facecolor="none", hatch="//", edgecolor="gray", linewidth=0.0)
    style.update(kw)
    return ax.fill_between(bins, np.append(lo, lo[-1]), np.append(hi, hi[-1]),
                           step="post", **style)


def f_chi2(NMC, Ndata, cov):
    # ignore singular entries
    which_bin = NMC > 0

    NMC = NMC[which_bin]
    Ndata = Ndata[which_bin]
    cov = cov[which_bin, :]
    cov = cov[:, which_bin]

    delta = NMC - Ndata
    try:
        cov_inv = np.linalg.inv(cov)
    except np.linalg.LinAlgError as _:
        return -1, which_bin.sum()

    return delta@cov_inv@delta, which_bin.sum()


def make_plot_data(var, bins, cut, mc_weight, breakdown, areanorm, breakdown_labels, breakdown_colors, xlabel, title,
                   det="ICARUS", fillna=np.nan, syst=None, annotate=None):
    # Blinding backstop. Every plotdata in this module is built here, and `det`
    # is always in hand, so this is the one place a blinded variable cannot slip
    # past. run_detector filters upstream so this never fires in normal use --
    # it raises rather than skipping because reaching it means a call site lost
    # its filter, and a silent skip would hide that. NB the "ICARUS" default is
    # not one of DETECTORS, so blind_vars() fails closed on it.
    for _v in (var if isinstance(var, list) else [var]):
        if _v in blind_vars(det):
            raise RuntimeError("%r is blinded for %s -- refusing to build a "
                               "data/MC comparison (see BLIND_VARS)" % (_v, det))

    if syst is None:
        syst = systematics_by_cut[cut]

    pvars = breakdown(df.loc[df[cut], var].fillna(fillna), df[df[cut]])
    weights = breakdown(df.loc[df[cut], mc_weight], df[df[cut]])

    NMC_breakdown = []
    for pvar, w in zip(pvars, weights):
        thisNMC, bins = np.histogram(pvar, bins=bins, weights=w)
        NMC_breakdown.append(thisNMC)

    NMC,_ = np.histogram(df.loc[df[cut], var].fillna(fillna), bins=bins, weights=df.loc[df[cut], mc_weight])
    NMC_abs = NMC
    if areanorm:
        diff = (bins[1:] - bins[:-1])
        norm = np.sum(NMC*diff)
        if norm > 1e-5:
            NMC = NMC / norm
            for i in range(len(NMC_breakdown)):
                NMC_breakdown[i] = NMC_breakdown[i] / norm

    NMC_breakdown = np.array(NMC_breakdown)

    NON,_ = np.histogram(ONdf.loc[ONdf[cut], var].fillna(fillna), bins=bins)
    NOff,_ = np.histogram(OFFdf.loc[OFFdf[cut], var].fillna(fillna), bins=bins)

    N = NON - NOff*OFF_w
    Nerr = np.sqrt(NON + NOff*OFF_w**2)
    if areanorm:
        diff = (bins[1:] - bins[:-1])

        norm = np.sum(N*diff)
        if norm > 1e-5:
            N = N / norm
            Nerr = Nerr / norm

    cov = syst.cov(var, cut, bins, NMC_abs, shapeonly=areanorm, fillna=fillna)
    err = np.sqrt(np.diag(cov))

    cov_w_stat = cov + np.diag(Nerr**2) # add stat uncertainty
    chi2, ndof = f_chi2(NMC, N, cov_w_stat)

    return {
        "det": det,
        "title": title,
        "xlabel": xlabel,
        "bins": bins,
        "areanorm": areanorm,
        "breakdown_labels": breakdown_labels,
        "breakdown_colors": breakdown_colors,
        "NMC_breakdown": NMC_breakdown,
        "NMC_total": NMC,
        "NData": N,
        "NDataErr": Nerr,
        "cov": cov,
        "cov_w_stat": cov_w_stat,
        "chi2": chi2,
        "ndof": ndof,
        "POT": POT,
        # free-text note drawn on the RHS of the plot (used for the del_p range)
        "annotate": annotate,
    }


def ratio_plot(plt, plotdata, blinded=False, overlay=None):
    """Stacked-MC + data ratio plot.

    blinded=True draws the MC and its systematic band, and still prints the
    chi2 / ndof / p-value, but omits the data points from both panels.

    overlay is an optional second MC prediction drawn as a step line over the
    stack (and as its own ratio in the lower panel):
        {"N": array over the bins, "label": str, "color": ..., "linestyle": ...}
    It is drawn before the legends are built, so it picks up a legend entry.
    """
    fig, (ax0, ax1) = plt.subplots(2, 1, height_ratios=[3, 1], sharex=True)
    bins = plotdata["bins"]
    centers = (bins[:-1] + bins[1:])/2

    NMC_breakdown = plotdata["NMC_breakdown"]
    fill = np.array([centers for _ in range(NMC_breakdown.shape[0])]).T
    ax0.hist(fill, bins=bins, stacked=True, label=plotdata["breakdown_labels"],
                    color=plotdata["breakdown_colors"], weights=NMC_breakdown.T)

    NData = plotdata["NData"]
    NDataErr = plotdata["NDataErr"]
    line = None
    if not blinded:
        line = ax0.errorbar(centers, NData, NDataErr, color="black", linestyle="none", marker=".")

    NMC = plotdata["NMC_total"]
    err = np.sqrt(np.diag(plotdata["cov"]))
    syst_band(ax0, bins, NMC-err, NMC+err)

    if "Maximum" in plotdata["xlabel"]:
        ax0.set_yscale("log")

    if not blinded:
        ax1.errorbar(centers, NData/NMC, NDataErr/NMC, color="black", linestyle="none", marker=".")
    ax1.set_ylim([0.5, 1.5])
    ax1.axhline([1], color="red", linestyle="--")
    syst_band(ax1, bins, 1-err/NMC, 1+err/NMC)

    if overlay is not None:
        NO = np.asarray(overlay["N"], dtype=float)
        ostyle = dict(color=overlay.get("color", "black"),
                      linestyle=overlay.get("linestyle", "--"),
                      linewidth=overlay.get("linewidth", 2))
        # step over the full edge array so the last bin is closed off
        ax0.step(bins, np.append(NO, NO[-1]), where="post",
                 label=overlay.get("label", "Overlay"), **ostyle)
        # lower panel: the overlay relative to the same nominal MC the data
        # points are divided by, so the two are directly comparable
        with np.errstate(divide="ignore", invalid="ignore"):
            ax1.step(bins, np.append(NO, NO[-1])/np.append(NMC, NMC[-1]),
                     where="post", **ostyle)

    ax0.tick_params(axis='both', which='both', direction='in', length=6, width=1.5, labelsize=FONTSIZE, top=True, right=True)
    ax1.tick_params(axis='both', which='both', direction='in', length=6, width=1.5, labelsize=FONTSIZE, top=True, right=True)
    for spine in ax0.spines.values():
        spine.set_linewidth(1.5)
    ax1.set_xlabel(plotdata["xlabel"], fontsize=FONTSIZE, fontweight='bold')

    # an explicit ylabel overrides the default (variable-width bins normalized
    # to a fixed width, say); absent one, the old area/POT wording stands
    if plotdata.get("ylabel"):
        ax0.set_ylabel(plotdata["ylabel"], fontsize=FONTSIZE, fontweight='bold')
    elif plotdata["areanorm"]:
        ax0.set_ylabel('Area Normalized', fontsize=FONTSIZE, fontweight='bold')
    else:
        ax0.set_ylabel('Events / %.1f$\\times 10^{19}$ POT' % (plotdata["POT"]/1e19), fontsize=FONTSIZE, fontweight='bold')

    det = plotdata["det"]
    title = plotdata["title"]
    ax0.set_title(f"$\\bf{{{det}}}$ {title}", fontsize=FONTSIZE+2)
    # no data drawn when blinded, so no data legend entry either
    ld = None
    if line is not None:
        ld = ax0.legend([line], ["Data\n(ON Beam - OFF)"], frameon=False, loc="upper left", fontsize=10)

    ax0_l0, ax0_hi = ax0.get_ylim()
    ax0.set_ylim([ax0_l0, ax0_hi*1.7])

    ax0.legend(fontsize=12, loc="upper right", ncol=2, reverse=True)
    if ld is not None:
        ax0.add_artist(ld)

    # RHS note under the legend (.get, so plotdata from older pickles still works)
    annotate = plotdata.get("annotate")
    if annotate:
        ax0.text(0.95, 0.62, annotate, verticalalignment="top", horizontalalignment="right",
                 fontsize=FONTSIZE-2, transform=ax0.transAxes)

    chi2_str = "$\\chi^2_\\mathrm{shape}$" if plotdata["areanorm"] else "$\\chi^2$"
    chi2, ndof, pval = chi2_ndof_pvalue(plotdata)
    # "--" where the p-value is undefined: singular covariance (chi2 = -1) or
    # ndof = 0 (one populated bin, fully consumed by the area normalization)
    ptxt = "%.2g" % pval if np.isfinite(pval) else "--"
    ax0.text(0.05, 0.8, "%s: %.1f / %i\np: %s" % (chi2_str, chi2, ndof, ptxt),
            verticalalignment="top", horizontalalignment="left", fontsize=FONTSIZE-2, transform=ax0.transAxes)

    plt.subplots_adjust(hspace=0.05)
    return fig, ax0, ax1


# ============================================================
# Pool workers (module level so fork'd children resolve them)
# ============================================================
def _inner_mode(dat):
    (v, b, l, cutname) = dat
    return v, make_plot_data(v, b, cutname, "glob_scale", breakdown_mode, True, mode_labels,
                             mode_colors, l, cutname, fillna=-1, det=DETECTOR)


def _inner_pdg(dat):
    (v, b, l, cutname) = dat
    bkdwn = breakdown_pdg_p if "prot_cand" in v else breakdown_pdg_mu
    return v, make_plot_data(v, b, cutname, "glob_scale", bkdwn, True, pdg_labels,
                             pdg_colors, l, cutname, fillna=-1, det=DETECTOR)


# Absolute-(POT-)normalized twins of the two workers above: identical apart from
# areanorm=False, which drops the per-histogram area normalization of MC, data
# and every systematic universe, so the chi2 tests the rate as well as the shape.
# BLINDING: SBND only -- see the guard in run_detector().
def _inner_mode_abs(dat):
    (v, b, l, cutname) = dat
    return v, make_plot_data(v, b, cutname, "glob_scale", breakdown_mode, False, mode_labels,
                             mode_colors, l, cutname, fillna=-1, det=DETECTOR)


def _inner_pdg_abs(dat):
    (v, b, l, cutname) = dat
    bkdwn = breakdown_pdg_p if "prot_cand" in v else breakdown_pdg_mu
    return v, make_plot_data(v, b, cutname, "glob_scale", bkdwn, False, pdg_labels,
                             pdg_colors, l, cutname, fillna=-1, det=DETECTOR)


def _inner_dp(dat):
    """Mode-breakdown worker for a del_p-sliced stage.

    `cutname` is the sliced cut column ("PID Cut lodp"); the plot is titled with
    the base stage and carries the del_p range as an RHS annotation.
    """
    (v, b, l, cutname) = dat
    base, tag, rng = dp_parts(cutname)
    return v, make_plot_data(v, b, cutname, "glob_scale", breakdown_mode, True, mode_labels,
                             mode_colors, l, base, fillna=-1, det=DETECTOR,
                             annotate=dp_annotation(*rng))


def _inner_dp_abs(dat):
    (v, b, l, cutname) = dat
    base, tag, rng = dp_parts(cutname)
    return v, make_plot_data(v, b, cutname, "glob_scale", breakdown_mode, False, mode_labels,
                             mode_colors, l, base, fillna=-1, det=DETECTOR,
                             annotate=dp_annotation(*rng))


def _make_all_plotdata(worker, cutnames, plotvars, bins, labels, nproc):
    out = {}
    with Pool(min(nproc, len(plotvars))) as p:
        for cutname in cutnames:
            done = {}
            inputs = [(v, b, l, cutname) for (v, b, l) in zip(plotvars, bins, labels)]
            for v, plotdata in tqdm(p.imap_unordered(worker, inputs), total=len(inputs),
                                    desc=cutname):
                done[v] = plotdata
            # imap_unordered yields in completion order; reinsert in plotvars
            # order so downstream iteration (chi2 csv, pickle) is reproducible.
            out[cutname] = {v: done[v] for v in plotvars}
    return out


def chi2_ndof_pvalue(p):
    """(chi2, ndof, p-value) for one plotdata entry.

    An area-normalized comparison spends one degree of freedom on the norm, so
    the reported ndof is one below the number of bins f_chi2 actually used;
    absolute normalization spends none. f_chi2 returns chi2 = -1 when the
    covariance is singular -- that has no meaningful p-value, and neither does
    ndof = 0 (a single populated bin, which area normalization fully consumes),
    so both give NaN rather than a made-up number.

    Shared by ratio_plot and chi2_table so the number drawn on the figure and
    the number in the .tsv cannot drift apart.
    """
    chi2 = p["chi2"]
    ndof = p["ndof"] - int(p["areanorm"])
    pval = stats.chi2.sf(chi2, ndof) if (chi2 >= 0 and ndof > 0) else np.nan
    return chi2, ndof, pval


def chi2_table(all_plotdata, cutnames, plotvars, labels):
    """chi2 / ndof / p-value for every (stage, variable) in a plotdata dict."""
    rows = []
    for cname in cutnames:
        for v, l in zip(plotvars, labels):
            p = all_plotdata[cname][v]
            chi2, ndof, pval = chi2_ndof_pvalue(p)
            rows.append(dict(detector=p["det"],
                             norm="shape" if p["areanorm"] else "absolute",
                             cut=cname, var=v, label=l,
                             chi2=chi2, ndof=ndof, pvalue=pval))
    return pd.DataFrame(rows)


def _save_plots(all_plotdata, cutnames, plotvars, plotdir, suffix="", blinded=False):
    n = 0
    for cname in cutnames:
        for v in plotvars:
            pd_ = all_plotdata[cname][v]
            # Second blinding backstop, independent of whoever built the dict:
            # the detector is carried in the plotdata itself.
            assert v not in blind_vars(pd_["det"]), \
                "%s is blinded for %s (see BLIND_VARS)" % (v, pd_["det"])
            fig, _, _ = ratio_plot(plt, pd_, blinded=blinded)
            stem = "%s_%s_%s%s" % (pd_["det"].replace(" ", "-"),
                                   cname.replace(" ", "").replace(".", "").lower(), v, suffix)
            fig.savefig(os.path.join(plotdir, "pdf", stem + ".pdf"), bbox_inches="tight")
            fig.savefig(os.path.join(plotdir, "png", stem + ".png"), bbox_inches="tight")
            plt.close(fig)
            n += 1
    return n


# ============================================================
# Variable definitions
# ============================================================
def nu_E_frac_diff(d):
    """(E_CCQE - E_calo)/E_calo, from the frame's own two estimator columns.

    The per-event disagreement between the muon-only CCQE estimator and the
    calorimetric one: ~0 for genuinely QE-like kinematics, with a tail wherever
    FSI/MEC move energy out of the muon.

    Deliberately reads the two columns off the frame rather than rebuilding
    them from kinematics primitives: every universe that moves either estimator
    (binding energy, track splitting, detector/calorimetry variations) then
    moves this along with it, with no second implementation to keep in step.
    """
    return (d.nu_E_ccqe - d.nu_E_calo)/d.nu_E_calo


extravarnames = [
    "mu_costh",
    "mu_phi",
    "p_costh",
    "p_phi",
    "mu_p",
    "p_p",
    "mu_p_opening_angle_deg",
    "nu_E_ccqe",
    "nu_E_frac_diff",
]

# NB: mu_dir_z / p_dir_z are components of a unit direction vector, i.e. they
# ARE cos(theta). These used to be np.cos(d.mu_dir_z) -- cos(cos(theta)) --
# which folded the whole range into [0.54, 1] (hence the old 0..1 binning).
#
# mu_p / p_p are not stored in the flat df; the same sqrt(E^2 - m^2) inversion
# is what syst.recompute_kinematics uses, so these match the systematic
# universes. The muon reco momentum is range-based, so the inversion is exact.
extravars = [
    lambda d: d.mu_dir_z,
    lambda d: np.arctan2(d.mu_dir_x, d.mu_dir_y)*180/np.pi,
    lambda d: d.p_dir_z,
    lambda d: np.arctan2(d.p_dir_x, d.p_dir_y)*180/np.pi,
    lambda d: np.sqrt(np.maximum(d.mu_E**2 - kinematics.MUON_MASS**2, 0)),
    lambda d: np.sqrt(np.maximum(d.p_E**2 - kinematics.PROTON_MASS**2, 0)),
    # The GUMPLE dfs store mu_p_opening_angle_deg natively (gumple_cuts has no
    # recompute helper); the read-back is idempotent with the
    # "frame[vname] = v(frame)" loop.
    lambda d: d["mu_p_opening_angle_deg"],
    # muon-only CCQE energy estimator, at the nominal binding energy. The BE and
    # track-splitting universes get it from syst.recompute_kinematics instead,
    # which calls the same kinematics function -- one implementation, so the CV
    # and the universes cannot disagree.
    lambda d: kinematics.neutrino_energy_ccqe(
        np.sqrt(np.maximum(d.mu_E**2 - kinematics.MUON_MASS**2, 0)), d.mu_dir_z),
    # must stay AFTER nu_E_ccqe: the loop below assigns in list order, and this
    # reads the column the entry above writes.
    nu_E_frac_diff,
]

CUTS = [
    FV,
    cosmic_rej,
    twoprong_cut,
    pid_cut_1p,
    pid_cut_np,
]

CUTNAMES = [
    "Contained",
    "Cosmic Rej.",
    "Two Prong Cut",
    "PID 1p",
    "PID Np",
]

# NB: other_trk_length / other_shw_length are not stored in the GUMPLE dfs.
# The replacements (lossy -- track and shower lengths are merged): the length of
# the longest non-candidate pfp, plus the shower / other counts behind the
# cut_0shwother two-prong flag.
PLOTVARS = [
    "crlongtrkdiry",
    "nu_score",
    "othr_pfp_length",
    "n_shower",
    "n_other",
    "mu_chi2_of_prot_cand",
    "prot_chi2_of_prot_cand",
    "mu_chi2_of_mu_cand",
    "prot_chi2_of_mu_cand",
    "mu_trackScore",
    "del_p",
    "p_len",

    "slc_vtx_x",
    "slc_vtx_y",
    "slc_vtx_z",
    "mu_end_x",
    "mu_end_y",
    "mu_end_z",
    "p_end_x",
    "p_end_y",
    "p_end_z",

    "mu_costh",
    "mu_phi",
    "p_costh",
    "p_phi",

    "p_p",
    "mu_p_opening_angle_deg",
]

# ============================================================
# Variables REMOVED from this script's comparisons (2026-08-25): the muon
# momentum and the two neutrino-energy estimators (+ their ratio). Their
# MC/data comparisons -- area AND absolute normalized -- live in
# mcdata_constraint_gumple.py, which is SBND-only and 1p-only blinded.
# The columns themselves are still computed (extravars) and carried on every
# systematic universe so the constraint script can build covariances for them;
# their plot specs live here so both scripts stay in step.
# ============================================================
CONSTRAINT_VARS = {
    #  var:            (bins, label)
    "nu_E_calo":      (np.linspace(0.3, 1.2, 10), "$E^{\\nu}_\\mathrm{calo}$ [GeV]"),
    "nu_E_ccqe":      (np.linspace(0.3, 1.2, 10), "$E^{\\nu}_\\mathrm{CCQE}$ [GeV]"),
    "nu_E_frac_diff": (np.linspace(-0.5, 1.0, 16),
                       "$(E^{\\nu}_\\mathrm{CCQE} - E^{\\nu}_\\mathrm{calo}) / E^{\\nu}_\\mathrm{calo}$"),
    "mu_p":           (np.linspace(0.2, 1, 9), "$p_\\mu$ [GeV/c]"),
}

LABELS = [
    "CRLongTrkDirY",
    "$\\nu$ Score",
    "Maximum Other PFP Length [cm]",
    "Number of Other Showers",
    "Number of Other PFPs",
    "Proton Cand. $\\mu$-like PID",
    "Proton Cand. $p$-like PID",
    "Muon Cand. $\\mu$-like PID",
    "Muon Cand. $p$-like PID",
    "Muon Cand. Track Score",
    "$\\delta p$ [GeV/c]",
    "Proton Cand. Length [cm]",

    "Vertex X [cm]",
    "Vertex Y [cm]",
    "Vertex Z [cm]",
    "Muon End X [cm]",
    "Muon End Y [cm]",
    "Muon End Z [cm]",
    "Proton End X [cm]",
    "Proton End Y [cm]",
    "Proton End Z [cm]",

    "Muon $\\cos(\\theta)$",
    "Muon $\\phi$ [$^\\circ$]",
    "Proton $\\cos(\\theta)$",
    "Proton $\\phi$ [$^\\circ$]",

    "$p_p$ [GeV/c]",
    "$\\theta_{\\mu p}$ [$^\\circ$]",
]

# ============================================================
# del_p slices
# ============================================================
# The two neutrino-energy estimators are additionally plotted in three slices of
# the transverse momentum imbalance del_p, which separates QE-like kinematics
# (low del_p) from the FSI/MEC-enriched tail. Each entry: (tag, lo, hi); the tag
# is the plot-file suffix.
#
# NB: unlike every other cut column here, these read del_p -- a column that
# syst.recompute_kinematics *changes*. So the slice must be re-derived on any
# universe that recomputes it (the binding-energy and track-splitting universes
# below) rather than carried over from the CV, or events cannot migrate between
# slices and those systematics are understated.
DP_SLICES = [
    ("lodp",  0.0, 0.2),
    ("middp", 0.2, 0.4),
    ("hidp",  0.4, 0.6),
]

# 2026-08-25: empty -- the del_p-sliced comparisons were of the two
# neutrino-energy estimators, which moved to mcdata_constraint_gumple.py
# (see --dp-bins there). The slice cut columns are still built: the
# constraint script and the BE universe read them.
DP_PLOTVARS = []


def dp_cutname(cname, tag):
    """Cut-column name for stage `cname` restricted to del_p slice `tag`."""
    return "%s %s" % (cname, tag)


def dp_parts(cutname):
    """Inverse of dp_cutname: -> (base cutname, tag, (lo, hi)), or (cutname, None, None)."""
    for tag, lo, hi in DP_SLICES:
        if cutname.endswith(" " + tag):
            return cutname[:-(len(tag) + 1)], tag, (lo, hi)
    return cutname, None, None


def dp_mask(d, lo, hi):
    """Rows whose del_p falls in [lo, hi). NaN del_p is excluded."""
    return (d.del_p >= lo) & (d.del_p < hi)


def dp_annotation(lo, hi):
    return "$%.1f < \\delta p < %.1f$ GeV/c" % (lo, hi)


def dp_cutnames():
    return [dp_cutname(c, tag) for tag, _, _ in DP_SLICES for c in CUTNAMES]


def add_dp_cuts(frame, suffix=""):
    """Add the del_p-sliced cut columns for every stage to `frame`, in place.

    `suffix` picks which base stage column to slice ("" for the nominal stage
    columns, "_trig_up"/"_trig_dn" for the trigger-systematic universes).
    """
    for tag, lo, hi in DP_SLICES:
        m = dp_mask(frame, lo, hi)
        for cname in CUTNAMES:
            frame[dp_cutname(cname, tag) + suffix] = frame[cname + suffix] & m
    return frame


PID_CUTNAMES = [
    "Contained",
    "Cosmic Rej.",
    "Two Prong Cut",
    "PID 1p",
    "PID Np",
]

PID_PLOTVARS = [
    "mu_chi2_of_prot_cand",
    "prot_chi2_of_prot_cand",
    "mu_chi2_of_mu_cand",
    "prot_chi2_of_mu_cand",
]

PID_BINS = [
    np.linspace(0, 60, 13),
    np.linspace(0, 120, 9),
    np.linspace(0, 60, 13),
    np.linspace(0, 300, 21),
]

PID_LABELS = [
    "Proton Cand. $\\chi^2_\\mu$",
    "Proton Cand. $\\chi^2_p$",
    "Muon Cand. $\\chi^2_\\mu$",
    "Muon Cand. $\\chi^2_p$",
]


def plot_bins(detector):
    if detector == "SBND":
        x_range = np.linspace(-200, 200, 21)
        y_range = np.linspace(-200, 200, 21)
        z_range = np.linspace(0, 500, 21)
    else:
        x_range = np.linspace(-360, 260, 21)
        y_range = np.linspace(-182, 135, 21)
        z_range = np.linspace(-895, 895, 21)

    return [
        np.linspace(-1, 0.3, 14),
        np.linspace(0.25, 0.75, 11),
        np.array([-10] + list(np.linspace(0, 100, 11))),  # othr_pfp_length; NaN -> -1 underflow bin
        np.linspace(-0.5, 5.5, 7),                        # n_shower
        np.linspace(-0.5, 5.5, 7),                        # n_other
        np.linspace(0, 60, 13),
        np.linspace(0, 120, 9),
        np.linspace(0, 60, 13),
        np.linspace(0, 300, 21),
        np.linspace(0.4, 1, 7),
        np.linspace(0, 1, 11),
        np.linspace(0, 50, 11),

        x_range,
        y_range,
        z_range,
        x_range,
        y_range,
        z_range,
        x_range,
        y_range,
        z_range,

        np.linspace(-1, 1, 21),
        np.linspace(-180, 180, 21),
        np.linspace(-1, 1, 21),
        np.linspace(-180, 180, 21),

        # (For scale, the 0-50 cm p_len axis corresponds to p_p = 0.2-0.80 GeV/c.)
        np.linspace(0.2, 1.2, 11),
        np.linspace(0, 180, 19),
    ]


# zip() over the three index-aligned lists would silently drop the tail.
assert len(PLOTVARS) == len(LABELS) == len(plot_bins(DETECTORS[0]))
assert len(PID_PLOTVARS) == len(PID_LABELS) == len(PID_BINS)

# Every detector must declare a blinding policy, and every blinded name must be
# a real plot variable -- otherwise a typo in BLIND_VARS silently blinds nothing.
assert set(BLIND_VARS) == set(DETECTORS)
assert all(v in PLOTVARS + PID_PLOTVARS
           for vs in BLIND_VARS.values() for v in vs)


# ============================================================
# Driver
# ============================================================
def prepare_detector(detector, args):
    """Load one detector and build everything the plotting rounds need.

    Sets the module globals (DETECTOR, df, ONdf, OFFdf, OFF_w, POT,
    systematics_by_cut) that make_plot_data and the fork'd Pool workers read at
    module scope, and returns them in a dict together with the normalization,
    the OOAV frame, the output paths and the summary logger.

    Split out of run_detector so other scripts (mcdata_constraint.py) can reuse
    the load half without the plotting half.
    """
    global DETECTOR, df, ONdf, OFFdf, OFF_w, POT, systematics_by_cut

    DETECTOR = detector
    include_dirt = not args.no_dirt
    chi2_variations, chi2_norm = chi2_config(detector)

    # Blinding backstop -- parse_args already rejects this combination, but
    # run_detector is also called directly (tests, notebooks importing it).
    if getattr(args, "absnorm", False) and detector not in ABSNORM_DETECTORS:
        raise ValueError("absolute normalization is only allowed for %s, not %r "
                         "(ICARUS is still blinded)" % (ABSNORM_DETECTORS, detector))

    plotdir = args.plotdir
    os.makedirs(plotdir, exist_ok=True)
    os.makedirs(os.path.join(plotdir, "png"), exist_ok=True)
    os.makedirs(os.path.join(plotdir, "pdf"), exist_ok=True)

    detname = detector.replace(" ", "-")
    summary_lines = []

    def log(msg):
        print(msg, flush=True)
        summary_lines.append(str(msg))

    log("=" * 70)
    log("DETECTOR = %s   INCLUDE_DIRT = %s   ABSNORM = %s   BLINDED = %s"
        % (detector, include_dirt, bool(args.absnorm), bool(getattr(args, "blinded", False))))
    log("calo variations: %s" % ", ".join(
        "%s%s" % (v, "" if chi2_norm.get(v, 1.) == 1. else " (x%g)" % chi2_norm[v])
        for v in chi2_variations))
    log("=" * 70)

    files = detector_files(detector, args.df_dir)

    # ---- normalization ----
    norm = compute_norm(detector, files, log)
    OFF_w = norm["OFF_w"]
    POT = norm["POT"]
    NEVT = norm["NEVT"]

    # ---- MC CV ----
    log("\nLoading MC CV (%d files)..." % len(files["MC_FILES"]))
    df, match, mcpot = loaddf.loadl(files["MC_FILES"], njob=min(len(files["MC_FILES"]), 10),
                                    detector=detector, preselection=FV,
                                    reweight_aFF=True, pot_univ=True)

    # ---- detector variations ----
    log("\nLoading detector variations (%d groups)..." % len(files["DETVAR_FILES"]))
    if len(files["DETVAR_FILES"]) > 0:
        detvars, detvar_pots = zip(*[load_detvar(dv, dv_det)
                                     for (dv_det, dv) in tqdm(files["DETVAR_FILES"])])
    else:
        detvars, detvar_pots = [], []

    # ---- POT scaling ----
    log("MC CV POT scaling: %r" % (loaddf.scale_pot(df, mcpot, POT),))
    for i in range(len(detvars)):
        log(files["DETVAR_NAMES"][i])
        for j in range(len(detvars[i])):
            log("  %r" % (loaddf.scale_pot(detvars[i][j], detvar_pots[i][j], POT),))

    # ---- dirt ----
    if include_dirt and len(files["DIRT_FILES"]) > 0:
        log("\nLoading dirt (%d files)..." % len(files["DIRT_FILES"]))
        dirt, dirtmatch, dirtpot = loaddf.loadl(files["DIRT_FILES"],
                                                njob=min(len(files["DIRT_FILES"]), 10),
                                                detector=detector, preselection=FV,
                                                include_syst=False)
        log("dirt POT scaling: %r" % (loaddf.scale_pot(dirt, dirtpot, POT),))

        # dirt has no systematic-universe weights: set them to 1
        for c in df.columns:
            if "univ" in c:
                dirt[c] = 1

        # dirt has no chi2-variation columns: set them to nominal (no variation)
        for c in chi2vars:
            for chi2_var in chi2_variations:
                dirt[c.replace("chi2", "chi2" + chi2_var)] = dirt[c]

        dirt["dirt"] = True
        df["dirt"] = False
    else:
        dirt = None
        df["dirt"] = False

    # ---- event count bookkeeping ----
    log("")
    log("Number of data triggers/1e15 POT:  %r" % (NEVT*(1e15/POT)))
    log("Number of scaled MC triggers/1e15 POT:  %r" % (match.shape[0]*(1e15/mcpot)))

    if dirt is not None:
        N_dirtevt_per = dirtmatch.shape[0]*(1e15/dirtpot)
        N_dirtevt = dirtmatch.shape[0]*(POT/dirtpot)
    else:
        N_dirtevt_per = 0
        N_dirtevt = 0

    log("Number of scaled MC triggers/1e15 POT (w/dirt):  %r" % (N_dirtevt_per + match.shape[0]*(1e15/mcpot)))
    log("Number of data triggers:  %r" % NEVT)
    log("Number of scaled MC triggers:  %r" % (match.shape[0]*(POT/mcpot)))
    log("Number of scaled MC triggers (w/dirt):  %r" % (match.shape[0]*(POT/mcpot) + N_dirtevt))
    log("")

    # ---- merge dirt into the CV ----
    # it is also included (below) as a 100% syst. unc. via the OOAV sample.
    # NB: df must not be re-bound after this point -- the systematic objects
    # below hold references to it (in-place column edits are fine).
    if dirt is not None:
        df = pd.concat([df[~df.dirt], dirt])
        del dirt

    # Chi2/calorimetry variations built from the full CV df (incl. dirt, whose
    # variation columns were set to nominal); SampleSystematic needs no separate CV.
    chi2_detvars = [v_calovariation(df, V) for V in chi2_variations]

    df_ooav = df[OOAV(df)].copy()

    # ---- data ----
    log("Loading on-beam data...")
    ONdf, _, _ = loaddf.load(files["ONBEAM"], load_truth=False, include_syst=False,
                             detector=detector, preselection=FV)

    log("Loading off-beam data...")
    _offs = [loaddf.load(f, load_truth=False, include_syst=False, detector=detector,
                         preselection=FV, offbeampot=True)
             for f in files["OFFBEAM_FILES"]]
    OFFdf = pd.concat([o[0] for o in _offs])
    OFFPOT = sum(o[2] for o in _offs)  # cross-check only; OFF_w computed above is used for scaling
    log("OFF-beam POT (cross-check only): %r" % OFFPOT)

    # ---- derived variables ----
    for v, vname in zip(extravars, extravarnames):
        for frame in [df, df_ooav, ONdf, OFFdf]:
            frame[vname] = v(frame)

        for i in range(len(detvars)):
            for j in range(len(detvars[i])):
                detvars[i][j][vname] = v(detvars[i][j])

        for i in range(len(chi2_detvars)):
            chi2_detvars[i][vname] = v(chi2_detvars[i])

    # ---- per-stage cut columns ----
    for c, cname in zip(CUTS, CUTNAMES):
        # nominal stage: base cut AND the nominal flash cut
        for frame in [df, df_ooav, ONdf, OFFdf]:
            frame[cname] = c(frame) & gc.flash_cut(frame)

        # trigger-systematic universes (MC CV only; the SelectionSystematic reads these)
        df[cname + "_trig_up"] = c(df) & trig_flash_cut_var(df, +1)
        df[cname + "_trig_dn"] = c(df) & trig_flash_cut_var(df, -1)

        for i in range(len(detvars)):
            for j in range(len(detvars[i])):
                detvars[i][j][cname] = c(detvars[i][j]) & gc.flash_cut(detvars[i][j])

        for i in range(len(chi2_detvars)):
            chi2_detvars[i][cname] = c(chi2_detvars[i]) & gc.flash_cut(chi2_detvars[i])

    # ---- del_p-sliced cut columns ----
    for frame in [df, df_ooav, ONdf, OFFdf]:
        add_dp_cuts(frame)
    add_dp_cuts(df, suffix="_trig_up")
    add_dp_cuts(df, suffix="_trig_dn")
    for i in range(len(detvars)):
        for j in range(len(detvars[i])):
            add_dp_cuts(detvars[i][j])
    for i in range(len(chi2_detvars)):
        add_dp_cuts(chi2_detvars[i])

    # ---- binding-energy universe ----
    # Built here, after the per-stage cut columns exist. The *stage* cuts read
    # none of the recomputed columns, so they carry over from the CV unchanged --
    # but the del_p slices do, so they are re-derived from the universe's own
    # (shifted) del_p below.
    be_cols = sorted(set(BE_RECOMPUTE_COLS + PLOTVARS + PID_PLOTVARS
                         + list(CONSTRAINT_VARS)
                         + list(CUTNAMES) + dp_cutnames()))
    be_df = syst.shift_binding_energy(df[be_cols], BE_SHIFT, fraction=0.5)
    add_dp_cuts(be_df)
    # shift_binding_energy recomputes both nu_E_calo and nu_E_ccqe under the
    # shifted BE, but their ratio came along as a plain CV copy -- re-derive it,
    # or the shift (which moves the two estimators by different amounts) cancels
    # out of this variable and its systematic is understated.
    # NB: re-running the whole extravars loop here would be WRONG -- the
    # nu_E_ccqe lambda rebuilds that column at the nominal kinematics.BE and
    # would overwrite the shifted value.
    be_df["nu_E_frac_diff"] = nu_E_frac_diff(be_df)
    log("Binding-energy universe: %d rows (dE_miss = %.0f MeV on 50%% of events)"
        % (len(be_df), BE_SHIFT*1000))
    # be_df has the same rows and weights as the CV: shift_binding_energy
    # recomputes the shifted kinematics in place on a deterministic half of
    # the rows and leaves the rest at their stored CV values.
    def _yield(frame, cut):
        return float(frame.loc[frame[cut], "glob_scale"].sum())
    log("  del_p slice migration under the BE shift (%s): %s"
        % (CUTNAMES[-1],
           ", ".join("%s %+.1f" % (tag,
                                   _yield(be_df, dp_cutname(CUTNAMES[-1], tag))
                                   - _yield(df, dp_cutname(CUTNAMES[-1], tag)))
                     for tag, _, _ in DP_SLICES)))

    # ---- ICARUS track-splitting universes ----
    tracksplit_systs = []
    if "ICARUS" in detector:
        run = 2 if detector == "ICARUS Run2" else 4
        for name, (dim, coord) in loaddf.SPLIT_REGIONS.items():
            if run not in loaddf.SPLIT_RUNS[name]:
                # plane does not apply to this run period (the east cathode is
                # outside the Run 2 muon FV) -> no systematic, skipped for the
                # same reason as the no-crossers case below
                log("Track split %-14s: does not apply to Run %d, skipped" % (name, run))
                continue
            sdf, crosses = syst.split_tracks(df, dim, coord, runs=[run])
            if not crosses.any():
                # No crossers -> a null systematic. Skipped rather than added as a
                # zero term because the cuts do not survive an empty frame
                # (gc.cathode_cut reads df.detector.iloc[0]).
                log("Track split %-14s: no crossing muons in Run %d, skipped" % (name, run))
                continue
            # split_tracks recomputes the stored kinematics (nu_E_calo, del_p,
            # mu_E, ...) but the derived plot variables are plain copies carried
            # over from the CV -- mu_p is a function of the truncated mu_E, so it
            # has to be rebuilt or it stays at its CV value in this universe.
            for v, vname in zip(extravars, extravarnames):
                sdf[vname] = v(sdf)
            # mu_end_* moved, so every stage cut must be re-evaluated on the split
            # rows (FV includes gc.mufv_cut, and TrackSplittingSystematic.univ
            # indexes splitdf by whichever stage is being plotted)
            for c, cname in zip(CUTS, CUTNAMES):
                sdf[cname] = c(sdf) & gc.flash_cut(sdf)
            # split_tracks recomputes del_p, so re-slice rather than carry over
            add_dp_cuts(sdf)
            tracksplit_systs.append(
                syst.TrackSplittingSystematic(df, sdf, crosses, TRACKSPLIT_FRAC[name]))
            log("Track split %-14s: f = %5.2f%%, %6d crossing muons "
                "(%5d still pass the PID cut after splitting)"
                % (name, 100*TRACKSPLIT_FRAC[name], int(crosses.sum()),
                   int(sdf[CUTNAMES[-1]].sum())))

    # ---- cosmic normalization ----
    # SBND only; the ICARUS overlay MC already carries data cosmics.
    cosmic_systs = []
    if detector == "SBND":
        cosmic_systs.append(syst.SystSampleSystematic(df[df.true_iscosmic], norm=SBND_COSMIC_NORM))

    # ---- systematics ----
    # NB: no beam-off stat systematic here -- the off-beam statistical uncertainty
    # is already propagated into the data error bars in make_plot_data.
    base_systematics = syst.SystematicList([
        loaddf.FluxSystematic(df),
        loaddf.G4Systematic(df),
        # XSec: the weight-based knobs plus the binding-energy shift
        syst.SystematicList([loaddf.XSecSystematic(df), syst.SampleSystematic(be_df)]),
        syst.NormalizationSystematic(0.005),                                       # POT norm
        syst.SystematicList(
            [syst.SampleSystematic(list(dv[1:]), cvdf=dv[0]) for dv in detvars] +  # detector variations
            [syst.SampleSystematic(d, norm=chi2_norm.get(v, 1.))                   # calo/chi2 variations
             for d, v in zip(chi2_detvars, chi2_variations)] +
            cosmic_systs +                                                         # cosmic norm (SBND only)
            tracksplit_systs +                                                     # track splitting (ICARUS only)
            [syst.StatSampleSystematic(df)]                                        # MC stat
        ),
        syst.SystSampleSystematic(df_ooav),                                        # OOAV/dirt, 100% unc.
    ])

    # Per-cut: add the trigger (flash-PE scale) SelectionSystematic for that stage.
    systematics_by_cut = {
        cname: syst.SystematicList([
            base_systematics,
            syst.SelectionSystematic(df, ["%s_trig_up" % cname, "%s_trig_dn" % cname]),
        ])
        for cname in list(CUTNAMES) + dp_cutnames()
    }

    return dict(detector=detector, detname=detname, plotdir=plotdir,
                norm=norm, df=df, df_ooav=df_ooav, ONdf=ONdf, OFFdf=OFFdf,
                OFF_w=OFF_w, POT=POT, systematics_by_cut=systematics_by_cut,
                log=log, summary_lines=summary_lines)


def run_detector(detector, args):
    ctx = prepare_detector(detector, args)
    plotdir = ctx["plotdir"]
    detname = ctx["detname"]
    norm = ctx["norm"]
    log = ctx["log"]
    summary_lines = ctx["summary_lines"]

    # ---- blinding ----
    # Filter here, before anything is histogrammed, rather than at save time:
    # the blinded numbers then never exist, so the figures, the chi2 tables, the
    # legacy csv and the pickle are all covered by construction. It also skips
    # the per-variable systematic covariance, which is the expensive part.
    bins = plot_bins(detector)
    vis_vars, vis_bins, vis_labels = visible_spec(detector, PLOTVARS, bins, LABELS)
    pid_vars, pid_bins, pid_labels = visible_spec(detector, PID_PLOTVARS, PID_BINS,
                                                  PID_LABELS)
    if blind_vars(detector):
        log("BLINDED for %s: %s -- no plots, no chi2, not pickled"
            % (detector, ", ".join(blind_vars(detector))))

    # ---- plots: interaction-mode breakdown ----
    log("\nBuilding mode-breakdown plot data (%d cuts x %d vars)..." % (len(CUTNAMES), len(vis_vars)))
    all_plotdata = _make_all_plotdata(_inner_mode, CUTNAMES, vis_vars, vis_bins, vis_labels, args.nproc)

    n_saved = 0
    if not args.no_save:
        n_saved += _save_plots(all_plotdata, CUTNAMES, vis_vars, plotdir,
                               blinded=args.blinded)

    # ---- plots: PDG breakdown ----
    all_plotdata_pdg = {}
    if not args.skip_pdg:
        log("\nBuilding PDG-breakdown plot data (%d cuts x %d vars)..." % (len(PID_CUTNAMES), len(pid_vars)))
        all_plotdata_pdg = _make_all_plotdata(_inner_pdg, PID_CUTNAMES, pid_vars,
                                              pid_bins, pid_labels, args.nproc)
        if not args.no_save:
            n_saved += _save_plots(all_plotdata_pdg, PID_CUTNAMES, pid_vars, plotdir,
                                   suffix="_bkdwnpdg", blinded=args.blinded)

    # ---- plots: absolute (POT) normalization ----
    # The same two rounds with areanorm=False: data and MC on a common POT scale,
    # so the chi2 is absolute and tests the rate as well as the shape.
    all_plotdata_abs = {}
    all_plotdata_pdg_abs = {}
    if args.absnorm:
        log("\nBuilding absolute-normalized mode-breakdown plot data...")
        all_plotdata_abs = _make_all_plotdata(_inner_mode_abs, CUTNAMES, vis_vars, vis_bins,
                                              vis_labels, args.nproc)
        if not args.no_save:
            n_saved += _save_plots(all_plotdata_abs, CUTNAMES, vis_vars, plotdir,
                                   suffix="_absnorm", blinded=args.blinded)

        if not args.skip_pdg:
            log("\nBuilding absolute-normalized PDG-breakdown plot data...")
            all_plotdata_pdg_abs = _make_all_plotdata(_inner_pdg_abs, PID_CUTNAMES,
                                                      pid_vars, pid_bins, pid_labels,
                                                      args.nproc)
            if not args.no_save:
                n_saved += _save_plots(all_plotdata_pdg_abs, PID_CUTNAMES, pid_vars,
                                       plotdir, suffix="_bkdwnpdg_absnorm",
                                       blinded=args.blinded)

    # ---- plots: energy estimators in del_p slices ----
    # Same stages and mode breakdown, restricted to each del_p slice. Saved with
    # the slice tag as the filename suffix (..._lodp / _middp / _hidp).
    # NB the index lookups go against the unfiltered PLOTVARS/LABELS on purpose --
    # that is why visible_spec() never mutates them.
    dp_vars = visible_vars(detector, DP_PLOTVARS)
    dp_bins = [bins[PLOTVARS.index(v)] for v in dp_vars]
    dp_labels = [LABELS[PLOTVARS.index(v)] for v in dp_vars]

    all_plotdata_dp = {}
    all_plotdata_dp_abs = {}
    if not dp_vars:
        # DP_PLOTVARS is empty: the del_p-sliced energy comparisons moved to
        # mcdata_constraint_gumple.py. The empty stores make the dp chi2
        # tables below skip themselves.
        log("\ndel_p-sliced energy round: moved to mcdata_constraint_gumple.py, skipped")
    else:
        log("\nBuilding del_p-sliced plot data (%d slices x %d cuts x %d vars)..."
            % (len(DP_SLICES), len(CUTNAMES), len(dp_vars)))
        for tag, lo, hi in DP_SLICES:
            sliced = [dp_cutname(c, tag) for c in CUTNAMES]
            for worker, store, sfx, on in ((_inner_dp, all_plotdata_dp, "", True),
                                           (_inner_dp_abs, all_plotdata_dp_abs, "_absnorm",
                                            bool(args.absnorm))):
                if not on:
                    continue
                byname = _make_all_plotdata(worker, sliced, dp_vars, dp_bins,
                                            dp_labels, args.nproc)
                # re-key "PID Cut lodp" -> "PID Cut" so the saved stem keeps the plain
                # stage name and the slice shows up only as the filename suffix
                store[tag] = {dp_parts(k)[0]: val for k, val in byname.items()}
                if not args.no_save:
                    n_saved += _save_plots(store[tag], CUTNAMES, dp_vars, plotdir,
                                           suffix="_" + tag + sfx, blinded=args.blinded)

    log("\nSaved %d figures (pdf+png) to %s" % (n_saved, plotdir))

    # ---- chi2 tables ----
    # One TSV per (breakdown, normalization): chi2, ndof and p-value per variable.
    tables = [("mode_shape", all_plotdata, CUTNAMES, vis_vars, vis_labels),
              ("pdg_shape", all_plotdata_pdg, PID_CUTNAMES, pid_vars, pid_labels),
              ("mode_abs", all_plotdata_abs, CUTNAMES, vis_vars, vis_labels),
              ("pdg_abs", all_plotdata_pdg_abs, PID_CUTNAMES, pid_vars, pid_labels)]

    chi2dfs = {}
    for tag, apd, cnames, pvars, labs in tables:
        if not apd:
            continue
        chi2dfs[tag] = chi2_table(apd, cnames, pvars, labs)

    # del_p slices: one table per normalization, with the slice as a column
    for name, store in (("dp_shape", all_plotdata_dp), ("dp_abs", all_plotdata_dp_abs)):
        if not store:
            continue
        parts = []
        for tag, lo, hi in DP_SLICES:
            if tag not in store:
                continue
            t = chi2_table(store[tag], CUTNAMES, dp_vars, dp_labels)
            t.insert(2, "dp", tag)
            t.insert(3, "dp_lo", lo)
            t.insert(4, "dp_hi", hi)
            parts.append(t)
        if parts:
            chi2dfs[name] = pd.concat(parts, ignore_index=True)

    if not args.no_save:
        for tag, tdf in chi2dfs.items():
            tsv_path = os.path.join(plotdir, "%s_chi2_%s.tsv" % (detname, tag))
            tdf.to_csv(tsv_path, sep="\t", index=False, float_format="%.6g")
            log("Wrote %s" % tsv_path)

    # legacy combined csv (shape-normalized rounds only), kept for anything
    # already reading it
    rows = []
    for tag, apd in (("mode", all_plotdata), ("pdg", all_plotdata_pdg)):
        for cname, byvar in apd.items():
            for v, p in byvar.items():
                rows.append(dict(breakdown=tag, cut=cname, var=v,
                                 chi2=p["chi2"], ndof=p["ndof"] - int(p["areanorm"])))
    chi2df = pd.DataFrame(rows)

    if not args.no_save:
        chi2_path = os.path.join(plotdir, "%s_chi2.csv" % detname)
        chi2df.to_csv(chi2_path, index=False)
        log("Wrote %s" % chi2_path)

        # The pickle carries NData/NDataErr/cov verbatim and is unaffected by
        # --blinded, so it is the costliest thing to get wrong. Check the keys
        # actually about to be written rather than trusting the upstream filter.
        _bad = set(blind_vars(detector))
        for _store in (all_plotdata, all_plotdata_pdg, all_plotdata_abs,
                       all_plotdata_pdg_abs):
            for _byvar in _store.values():
                assert not (set(_byvar) & _bad), "blinded var in plotdata pickle"
        for _store in (all_plotdata_dp, all_plotdata_dp_abs):
            for _bytag in _store.values():
                for _byvar in _bytag.values():
                    assert not (set(_byvar) & _bad), "blinded var in dp pickle"

        pkl_path = os.path.join(plotdir, "%s_plotdata.pkl" % detname)
        with open(pkl_path, "wb") as f:
            pickle.dump({"detector": detector, "norm": norm,
                         "mode": all_plotdata, "pdg": all_plotdata_pdg,
                         "mode_abs": all_plotdata_abs, "pdg_abs": all_plotdata_pdg_abs,
                         "dp": all_plotdata_dp, "dp_abs": all_plotdata_dp_abs,
                         "dp_slices": DP_SLICES}, f)
        log("Wrote %s" % pkl_path)

    for tag, tdf in chi2dfs.items():
        nbad = int((tdf.chi2 < 0).sum())
        log("chi2 %-10s: %d entries (%d singular)" % (tag, len(tdf), nbad))

    if not args.no_save:
        summary_path = os.path.join(plotdir, "%s_summary.txt" % detname)
        with open(summary_path, "w") as f:
            f.write("\n".join(summary_lines) + "\n")
        print("Wrote %s" % summary_path, flush=True)

    return 0


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-d", "--detector", action="append", default=None,
                   help="Detector to run: 'SBND', 'ICARUS Run2', 'ICARUS Run4', or "
                        "'all'. Repeatable.")
    p.add_argument("--df-dir", default=DEFAULT_DF_DIR, help="Directory of input .df files")
    p.add_argument("--plotdir", default=DEFAULT_PLOTDIR, help="Output directory")
    p.add_argument("--nproc", type=int, default=10, help="Pool size for plot building")
    p.add_argument("--no-dirt", action="store_true", help="Skip the dirt sample")
    p.add_argument("--no-save", action="store_true", help="Build plots but write nothing")
    p.add_argument("--skip-pdg", action="store_true", help="Skip the PDG-breakdown round")
    p.add_argument("--absnorm", action="store_true",
                   help="Also run the absolute-(POT-)normalized comparison. SBND ONLY: "
                        "ICARUS is still blinded (see ABSNORM_DETECTORS)")
    p.add_argument("--blinded", action="store_true",
                   help="Draw MC, its systematic band and the chi2/ndof/p-value, but no "
                        "data points in either panel. NB: this does NOT relax the "
                        "--absnorm detector restriction -- the printed chi2 is still a "
                        "data/MC comparison")
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

    # Fail before anything is loaded rather than part-way through a multi-detector run.
    bad = [d for d in dets if d not in ABSNORM_DETECTORS]
    if args.absnorm and bad:
        p.error("--absnorm is only allowed for %s (requested %s). Absolute normalization "
                "exposes the data rate and ICARUS is still blinded." % (ABSNORM_DETECTORS, bad))

    return args


def main(argv=None):
    args = parse_args(argv)

    if len(args.detector) == 1:
        return run_detector(args.detector[0], args)

    # More than one detector: run each in its own subprocess so memory is
    # reclaimed in between.
    passthrough = ["--df-dir", args.df_dir, "--plotdir", args.plotdir,
                   "--nproc", str(args.nproc)]
    for flag in ("no_dirt", "no_save", "skip_pdg", "absnorm", "blinded"):
        if getattr(args, flag):
            passthrough.append("--" + flag.replace("_", "-"))

    rc = 0
    for det in args.detector:
        cmd = [sys.executable, os.path.abspath(__file__), "--detector", det] + passthrough
        print("\n>>> %s" % " ".join(repr(c) if " " in c else c for c in cmd), flush=True)
        r = subprocess.run(cmd)
        if r.returncode != 0:
            print("!!! %s FAILED (exit %d)" % (det, r.returncode), flush=True)
            rc = r.returncode
    return rc


if __name__ == "__main__":
    sys.exit(main())
