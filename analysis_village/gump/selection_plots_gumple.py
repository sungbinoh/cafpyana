"""CV-only GUMP selection plots for the technote (GUMPLE / sbn-rewgted-2x).

Script refactor of the cut-optimization notebook
SignalBoxSystematics_rewgt10_OGopt_muTrackShared_PIDperDet_customCalo_muMin40_OGeval.ipynb
(cafpyana root), keeping only its figure output and porting it from the reCAF
rewgt10 dataframes + gump_cuts to the GUMPLE production + gumple_cuts. The
optimizer, the figure-of-merit machinery and the 1D cut scans are NOT ported:
the technote's 1D_scan.png documents the original optimization and is kept.

Three figure groups (all central value, no systematics; MC weighted by
glob_scale = POT scale x cvwgt, dirt and beam-off data included as
"Non-Fid." / "In-Time Cosmic"):

  cuts       staged cut-variable breakdowns (SBND | ICARUS columns), one page
             per production cut group, stacked by true final state, with the
             gumple_cuts thresholds drawn:
               cut_breakdown_cosmic       nu_score (NO cut in production -- it
                                          is shown for reference only), mu-p
                                          opening angle < 160 deg, at the
                                          preselection + flash stage
               cut_breakdown_track_score  mu_trackScore > 0.5 after cosmic rej.
               cut_breakdown_mu_len       40 < L_mu < 400 cm (the lower bound
                                          is already part of the production
                                          preselection, so nothing below it
                                          survives into the dataframes)
               cut_breakdown_mu_cand      muon-candidate chi2_mu / chi2_p
               cut_breakdown_p_cand       proton-candidate chi2_mu / chi2_p
  eff        efficiency.png: ABSOLUTE selection efficiency -- selected events
             over ALL generated neutrino interactions with true vertex in the
             fiducial volume -- vs true muon / proton / charged-pion momentum
             and, for pi0 events, true neutrino energy; SBND and ICARUS
             overlaid. The denominator is built from the mcnu table (every
             generated neutrino, reconstructed or not) joined to the GENIE
             event record (evtrec, status-1 final-state particles) for the
             leading-particle momenta; those reproduce the production's stored
             true_*_p exactly, so numerator and denominator share one
             definition. Both carry the POT scale and the aFF CV reweight.
             Six panels: muon / proton / charged-pion (vs momentum), pi0
             events, numu CC QE events, and numu CC 1u1p events (exactly one
             proton with T_p > 50 MeV, no charged pion, no pi0) vs E_nu.
  signalbox  near/far signal-box distributions: SBND and ICARUS stacked by
             interaction mode + a totals panel with the SBND/ICARUS ratio, for
             delta p, E_reco (nu_E_calo) and E_reco in three delta p slices

Outputs (in --plotdir/{png,pdf}/):
  cut_breakdown_{cosmic,track_score,mu_len,mu_cand,p_cand}
  efficiency
  signalbox_{del_p,nu_E_calo}_near_far, signalbox_nu_E_calo_{lodp,middp,hidp}_near_far

Normalization follows the notebook: SBND to 1e20 POT; ICARUS Run 2 to 2e20 and
Run 4 to 3e20, then concatenated (5e20). Input files come from
mcdata_comparison_gumple.detector_files (3 SBND CV shards, 3+4 ICARUS).

    python selection_plots_gumple.py --df-dir ../sbn-rewgted-21/ --plotdir <dir>
                                     [--only cuts|eff|signalbox] [--no-dirt] [--no-offbeam]

Run with cwd = analysis_village/gump. NB: FV must stay a module-scope function
named FV -- loaddf's cache key hashes the preselection's qualname.
"""

import argparse
import os
import sys
import warnings

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
warnings.simplefilter(action="ignore", category=FutureWarning)

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (_HERE, os.path.abspath(os.path.join(_HERE, "..", ".."))):
    #if _p not in sys.path:
    sys.path.insert(0, _p)
sys.path.insert(0, os.path.join(_HERE, "..", "gumple"))

import gumple_cuts as gc  # noqa: E402
import loaddf  # noqa: E402
import mcdata_comparison_gumple as mdc  # noqa: E402

DEFAULT_DF_DIR = "/Users/gputnam/Work/osc/sbn-rewgted-21/"
DEFAULT_PLOTDIR = "/Users/gputnam/Work/osc/cafpyana/plots-gumple-2026-09-02-rewgted21/selection/"

# POT normalization (as in the notebook)
SBND_POT = 1e20
ICARUS_RUN2_POT = 2e20
ICARUS_RUN4_POT = 3e20

# 1u1p signal definition: proton kinetic-energy threshold and the momentum it
# corresponds to (KE = sqrt(p^2 + m^2) - m)
P_KE_THRESHOLD_KE = 0.050          # GeV
_M_P = 0.938272
P_KE_THRESHOLD = np.sqrt((P_KE_THRESHOLD_KE + _M_P)**2 - _M_P**2)   # 0.310 GeV/c

FONTSIZE = 12
plt.rcParams["figure.max_open_warning"] = 0


# ----------------------------------------------------------------------------
# Preselection -- keep at module scope, named FV (loaddf cache key)
# ----------------------------------------------------------------------------
def FV(df):
    return gc.sanity_cut(df) & gc.slcfv_cut(df) & df.cut_contained & df.cut_cathode


# ----------------------------------------------------------------------------
# Loading
# ----------------------------------------------------------------------------
TRUTH_COLS = ["genie_mode", "nu_E", "true_pdg", "true_iscc", "true_isnc", "true_iscosmic",
              "true_vtx_x", "true_vtx_y", "true_vtx_z", "true_mu_p", "true_p_p", "true_p2_p",
              "true_cpi_p", "true_pi0_p"]


def _load_group(files, detector, pot, log, **kw):
    df, match, mcpot = loaddf.loadl(files, njob=min(len(files), 10), detector=detector,
                                    preselection=FV, include_syst=False, **kw)
    log("  %s: %d rows, %.3g POT -> scale %r" % (detector, len(df), mcpot,
                                                  loaddf.scale_pot(df, mcpot, pot)))
    return df, mcpot


# ----------------------------------------------------------------------------
# Truth denominator: every generated neutrino (mcnu) + final-state kinematics
# from the GENIE event record (evtrec)
# ----------------------------------------------------------------------------
def _n_splits(fname, key):
    with h5py.File(fname, "r") as f:
        return len([k for k in f.keys() if k.startswith(key + "_")])


def _evtrec_final_state(er):
    """Per GENIE event (__ntuple, evtrec entry): leading final-state (GHEP
    status 1) muon / proton / charged-pion momentum [GeV/c] and the pi0 count.
    Verified on sbn-rewgted-21 to reproduce the production's true_mu_p /
    true_p_p / true_cpi_p / true_pi0_p (from rec.mc.nu.prim) bit-for-bit."""
    fs = er[er.status == 1]
    p = np.sqrt(fs.px**2 + fs.py**2 + fs.pz**2)
    pdg = fs.pdg
    lv = ["__ntuple", "entry"]

    def lead(mask):
        return p[mask].groupby(level=lv).max()

    is_p = pdg == 2212
    out = pd.DataFrame({
        "mu_p": lead(pdg.abs() == 13),
        "p_p": lead(is_p),
        "cpi_p": lead(pdg.abs() == 211),
        "n_pi0": fs[pdg == 111].groupby(level=lv).size(),
        "n_cpi": fs[pdg.abs() == 211].groupby(level=lv).size(),
        # protons above the 1u1p kinetic-energy threshold
        "n_p50": fs[is_p & (p > P_KE_THRESHOLD)].groupby(level=lv).size(),
    })
    out = out.reindex(fs.groupby(level=lv).size().index)
    for c in ("n_pi0", "n_cpi", "n_p50"):
        out[c] = out[c].fillna(0).astype(int)
    return out


def _aff_cvwgt(fname, idf):
    """The axial-form-factor CV reweight per mcnu row, exactly as loaddf applies
    it to the reconstructed rows (product of the cv/morph column of every
    xsec_cv_rwgt knob); 1 where the weight table is absent."""
    try:
        w = pd.read_hdf(fname, "wgt_%i" % idf)[loaddf.xsec_cv_rwgt]
    except (KeyError, ValueError):
        return None
    cv = pd.Series(1., index=w.index)
    for k in loaddf.xsec_cv_rwgt:
        cvcol = "cv" if "cv" in w[k].columns else "morph"
        cv = cv * w[k][cvcol]
    return cv


def load_truth(files, detector, log):
    """All generated neutrinos of the MC files: one row per mcnu entry with
    detector/Run, nu_E, pdg, is_cc/is_nc, genie_mode, true vertex, the
    final-state leading momenta, in_fv (true vertex inside the analysis FV,
    gumple_cuts._fv_cut with the slice-vertex insets) and cvwgt.

    NB: not de-duplicated the way loaddf de-duplicates the reconstructed rows
    (a ~0.1% effect on sbn-rewgted-21)."""
    frames = []
    for fname in files:
        n = _n_splits(fname, "mcnu")
        for idf in range(n):
            m = pd.read_hdf(fname, "mcnu_%i" % idf)
            try:
                er = pd.read_hdf(fname, "evtrec_%i" % idf)
            except KeyError:
                raise RuntimeError("%s split %d has no evtrec table -- cannot build the "
                                   "truth denominator" % (fname, idf))
            kin = _evtrec_final_state(er)
            key = pd.MultiIndex.from_arrays(
                [m.index.get_level_values("__ntuple"), m.genie_evtrec_idx.astype(int)],
                names=["__ntuple", "entry"])
            kin = kin.reindex(key)
            kin.index = m.index
            t = pd.DataFrame({
                "detector": m.detector.astype(str), "Run": m.Run,
                "nu_E": m.nu_E, "pdg": m.pdg, "is_cc": m.is_cc.astype(bool),
                "is_nc": m.is_nc.astype(bool), "genie_mode": m.genie_mode,
                "pos_x": m.pos_x, "pos_y": m.pos_y, "pos_z": m.pos_z,
                "mu_p": kin.mu_p, "p_p": kin.p_p, "cpi_p": kin.cpi_p,
                "n_pi0": kin.n_pi0.fillna(0).astype(int),
                "n_cpi": kin.n_cpi.fillna(0).astype(int),
                "n_p50": kin.n_p50.fillna(0).astype(int),
            }, index=m.index)
            cv = _aff_cvwgt(fname, idf)
            t["cvwgt"] = 1. if cv is None else cv.reindex(m.index).fillna(1.).to_numpy()
            frames.append(t.reset_index(drop=True))
        log("  truth: %s -> %d splits, %d neutrinos" % (os.path.basename(fname), n,
                                                       sum(len(x) for x in frames[-n:])))
    t = pd.concat(frames, ignore_index=True)
    # the loaddf-style detector label ("ICARUS Run2") drives the FV choice
    t["detector"] = detector
    vtx = pd.DataFrame({"detector": t.detector, "Run": t.Run, "x": t.pos_x, "y": t.pos_y,
                        "z": t.pos_z})
    t["in_fv"] = gc._fv_cut(vtx).to_numpy()
    return t


def load_detector(det, df_dir, log, include_dirt=True, include_offbeam=True, with_truth=True):
    """One combined frame per detector ('SBND' or 'ICARUS'), columns
    sample in {mc, dirt, offbeam} and glob_scale at the target POT; plus the
    generated-neutrino truth frame (load_truth) at the same POT scale, or None."""
    frames = []
    truths = []
    runs = [("SBND", SBND_POT)] if det == "SBND" else \
           [("ICARUS Run2", ICARUS_RUN2_POT), ("ICARUS Run4", ICARUS_RUN4_POT)]
    for run_det, pot in runs:
        files = mdc.detector_files(run_det, df_dir)
        log("loading %s MC CV (%d files)" % (run_det, len(files["MC_FILES"])))
        mc, mcpot = _load_group(files["MC_FILES"], run_det, pot, log, reweight_aFF=True)
        mc["sample"] = "mc"
        frames.append(mc)
        if with_truth:
            log("loading %s generated-neutrino truth (mcnu + evtrec)" % run_det)
            t = load_truth(files["MC_FILES"], run_det, log)
            t["glob_scale"] = (pot / mcpot) * t.cvwgt   # same scale as the MC rows
            log("  %s: %d generated neutrinos, %.3f in FV, sum w (FV) = %.1f" % (
                run_det, len(t), t.in_fv.mean(), t.glob_scale[t.in_fv].sum()))
            truths.append(t)
        if include_dirt and files["DIRT_FILES"]:
            log("loading %s dirt" % run_det)
            dirt, _ = _load_group(files["DIRT_FILES"], run_det, pot, log)
            dirt["sample"] = "dirt"
            frames.append(dirt)
        if include_offbeam and files["OFFBEAM_FILES"]:
            log("loading %s beam-off data" % run_det)
            off, _ = _load_group(files["OFFBEAM_FILES"], run_det, pot, log,
                                 offbeampot=True, load_truth=False)
            off["sample"] = "offbeam"
            for c in TRUTH_COLS:
                if c not in off.columns:
                    off[c] = np.nan
            frames.append(off)
    df = pd.concat(frames, ignore_index=True, copy=False)
    for c in TRUTH_COLS:
        if c not in df.columns:
            raise KeyError("missing truth column %r after load" % c)
    df["dirt"] = df["sample"] == "dirt"
    df["offbeam"] = df["sample"] == "offbeam"
    # opening angle is stored by the production; recompute only if absent
    if "mu_p_opening_angle_deg" not in df.columns:
        cos = (df.mu_dir_x * df.p_dir_x + df.mu_dir_y * df.p_dir_y + df.mu_dir_z * df.p_dir_z)
        df["mu_p_opening_angle_deg"] = np.degrees(np.arccos(np.clip(cos, -1, 1)))
    log("  combined %s: %d rows (mc %d, dirt %d, offbeam %d), sum w = %.1f" % (
        det, len(df), (df["sample"] == "mc").sum(), df.dirt.sum(), df.offbeam.sum(),
        df.glob_scale.sum()))
    truth = pd.concat(truths, ignore_index=True) if truths else None
    return df, truth


# ----------------------------------------------------------------------------
# Truth categorization
# ----------------------------------------------------------------------------
def _pos(s):
    return s.notna() & (s > 0)


def _cosmic_like(d):
    gm = d.genie_mode
    return d.offbeam | gm.isna() | (gm < 0) | d.true_iscosmic.fillna(False).astype(bool)


def _nonfid(d):
    # dirt sample, or an MC neutrino whose true vertex lies outside the active volume
    with np.errstate(invalid="ignore"):
        ooav = mdc.OOAV(d)
    return d.dirt | (ooav & ~_cosmic_like(d))


# final-state categories: (label, color), stacked bottom -> top
FS_CATS = [
    ("$\\nu_\\mu$ CC 1$\\mu$1p0$\\pi$", "#315031"),
    ("Cosmic / Non-Fid.", "#5f8aa3"),
    ("NC", "#c89648"),
    ("Other $\\nu$", "#e6dcb8"),
    ("$\\nu_\\mu$ CC 1$\\mu$0p0$\\pi$", "#1e3f54"),
    ("$\\nu_\\mu$ CC 1$\\mu$Np0$\\pi$", "#95af8b"),
    ("$\\nu_\\mu$ CC 1$\\mu$1p1$\\pi$", "#7a6a2a"),
    ("$\\nu_\\mu$ CC Other", "#b5541c"),
]


def final_state(d):
    """Integer category per row, indexing FS_CATS."""
    cat = np.full(len(d), 7, dtype=int)
    cosmic = _cosmic_like(d) | _nonfid(d)
    nc = d.true_isnc.fillna(0).astype(bool) & ~cosmic
    iscc = d.true_iscc.fillna(0).astype(bool) & ~cosmic
    numucc = iscc & (np.abs(d.true_pdg) == 14)
    othernu = iscc & ~numucc
    n_p = _pos(d.true_p_p).astype(int) + _pos(d.true_p2_p).astype(int)
    has_cpi = _pos(d.true_cpi_p)
    has_pi0 = _pos(d.true_pi0_p)
    has_pi = has_cpi | has_pi0
    cat[(numucc & (n_p == 1) & ~has_pi).to_numpy()] = 0
    cat[cosmic.to_numpy()] = 1
    cat[nc.to_numpy()] = 2
    cat[othernu.to_numpy()] = 3
    cat[(numucc & (n_p == 0) & ~has_pi).to_numpy()] = 4
    cat[(numucc & (n_p >= 2) & ~has_pi).to_numpy()] = 5
    cat[(numucc & (n_p == 1) & has_cpi & ~has_pi0).to_numpy()] = 6
    return cat


# interaction-mode categories for the signal-box plots (bottom -> top)
MODE_CATS = [
    ("$\\nu_\\mu$ CC QE", "#315031"),
    ("In-Time Cosmic", "#6b4c7a"),
    ("Out-of-Time Cosmic", "#95af8b"),
    ("Non-Fid.", "#d54c28"),
    ("NC", "#1e3f54"),
    ("Other $\\nu$", "#c89648"),
    ("$\\nu_\\mu$ CC MEC", "#43140b"),
    ("$\\nu_\\mu$ CC RES", "#4c6b7a"),
]


def mode_category(d):
    cat = np.full(len(d), 5, dtype=int)   # default: other nu
    nonfid = _nonfid(d)
    cosmic = _cosmic_like(d) & ~d.offbeam & ~nonfid
    numucc = d.true_iscc.fillna(0).astype(bool) & (np.abs(d.true_pdg) == 14) & ~nonfid & ~cosmic & ~d.offbeam
    nc = d.true_isnc.fillna(0).astype(bool) & ~nonfid & ~cosmic & ~d.offbeam
    cat[(numucc & (d.genie_mode == 0)).to_numpy()] = 0
    cat[d.offbeam.to_numpy()] = 1
    cat[cosmic.to_numpy()] = 2
    cat[nonfid.to_numpy()] = 3
    cat[nc.to_numpy()] = 4
    cat[(numucc & (d.genie_mode == 10)).to_numpy()] = 6
    cat[(numucc & (d.genie_mode == 1)).to_numpy()] = 7
    return cat


# ----------------------------------------------------------------------------
# Selection stages (production cuts, gumple_cuts thresholds)
# ----------------------------------------------------------------------------
def det_th(d, key):
    return pd.Series(gc._det_cut_th(d.detector, key), index=d.index)


def stage_masks(d, log=None):
    """Cumulative masks of the production selection, in the order the cut
    breakdown pages show them. Every mask is a boolean numpy array."""
    presel = gc.presel_cut(d) & gc.flash_cut(d)
    cosmic = gc.cosmic_cut(d)
    ts = d.mu_trackScore >= det_th(d, "musel_track_score_min")
    lenc = d.mu_len.between(det_th(d, "musel_len_th_min"), det_th(d, "musel_len_th_max"))
    twoprong = d.n_pfp == 2
    mupid = (d.mu_chi2_of_mu_cand < det_th(d, "musel_muscore_th")) & \
            (d.prot_chi2_of_mu_cand > det_th(d, "musel_pscore_th"))
    ppid = (d.mu_chi2_of_prot_cand > det_th(d, "psel_muscore_th")) & \
           (d.prot_chi2_of_prot_cand < det_th(d, "psel_pscore_th"))
    m = {}
    m["presel"] = presel.to_numpy()
    m["cosmic"] = (presel & cosmic).to_numpy()
    m["track_score"] = (presel & cosmic & ts).to_numpy()
    m["mu_len"] = (presel & cosmic & ts & lenc).to_numpy()
    m["twoprong"] = (presel & cosmic & ts & lenc & twoprong).to_numpy()
    m["mu_cand"] = (presel & cosmic & ts & lenc & twoprong & mupid).to_numpy()
    m["final"] = (presel & cosmic & ts & lenc & twoprong & mupid & ppid).to_numpy()
    if log is not None and "gump_sel" in d.columns:
        # the production's own flag (presel & cosmic & flash & trk & pid & n_pfp==2)
        prod = d.gump_sel.fillna(False).astype(bool).to_numpy()
        mc = (d["sample"] == "mc").to_numpy()
        agree = (m["final"] == prod)[mc].mean()
        log("  final selection vs stored gump_sel (MC rows): agree %.4f; "
            "ours %d, stored %d" % (agree, m["final"][mc].sum(), prod[mc].sum()))
    return m


# ----------------------------------------------------------------------------
# Plot helpers
# ----------------------------------------------------------------------------
def _style(ax, xlabel, ylabel, title):
    ax.tick_params(axis="both", which="both", direction="in", length=5, width=1.2,
                   labelsize=FONTSIZE - 1, top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
    ax.set_xlabel(xlabel, fontsize=FONTSIZE)
    ax.set_ylabel(ylabel, fontsize=FONTSIZE)
    ax.set_title(title, fontsize=FONTSIZE + 1)
    ax.grid(alpha=0.25)


def stacked_hist(ax, x, w, cat, cats, bins, order=None):
    """Stack the weighted histogram of x by category. Returns the total."""
    order = order if order is not None else range(len(cats))
    xs, ws, labels, colors = [], [], [], []
    for i in order:
        sel = cat == i
        xs.append(x[sel])
        ws.append(w[sel])
        labels.append(cats[i][0])
        colors.append(cats[i][1])
    ax.hist(xs, bins=bins, weights=ws, stacked=True, color=colors, label=labels,
            edgecolor="white", linewidth=0.3)
    tot, _ = np.histogram(x, bins=bins, weights=w)
    return tot


def pot_label(det):
    return "Events / $10^{20}$ POT" if det == "SBND" else "Events / $5\\times10^{20}$ POT"


def cut_lines(ax, lines, keep):
    """lines: [(x, label)] dashed verticals; keep: (lo, hi) with None for open
    ends -> a 'KEEP' arrow between them (or from the cut to the axis edge)."""
    for x, lab in lines:
        ax.axvline(x, color="tab:blue", linestyle="--", linewidth=1.3, label=lab)
    if keep is None:
        return
    xlo, xhi = ax.get_xlim()
    lo, hi = keep
    y = ax.get_ylim()[1] * 0.92
    a = xlo + 0.02 * (xhi - xlo) if lo is None else lo
    b = xhi - 0.02 * (xhi - xlo) if hi is None else hi
    if a >= b:
        return
    ax.annotate("", xy=(a, y), xytext=(b, y),
                arrowprops=dict(arrowstyle="<->" if (lo is not None and hi is not None)
                                else ("<-" if lo is not None else "->"), color="black", lw=1.2))
    ax.text((a + b) / 2, y, "KEEP", ha="center", va="bottom", fontsize=FONTSIZE - 2,
            fontweight="bold")


def cuts_for(det):
    return gc.SBND_CUTS if det == "SBND" else gc.ICARUS_CUTS


# ----------------------------------------------------------------------------
# Figure group 1: staged cut breakdowns
# ----------------------------------------------------------------------------
def cut_pages():
    """name -> list of panel specs (var, xlabel, bins, stage, cut spec).
    cut spec is a callable(cutdict) -> (lines, keep)."""
    inf = np.inf

    def gt(key, fmt):
        def f(c):
            v = c[key]
            if not np.isfinite(v):
                return [], None
            return [(v, fmt % v)], (v, None)
        return f

    def lt(key, fmt):
        def f(c):
            v = c[key]
            if not np.isfinite(v):
                return [], None
            return [(v, fmt % v)], (None, v)
        return f

    def between(klo, khi, fmtlo, fmthi):
        def f(c):
            lo, hi = c[klo], c[khi]
            lines = [(lo, fmtlo % lo), (hi, fmthi % hi)]
            return lines, (lo, hi if np.isfinite(hi) else None)
        return f

    def nocut(c):
        return [], None

    return {
        "cosmic": [
            ("nu_score", "Neutrino score", np.linspace(0.0, 1.0, 41), "presel", nocut,
             "no cut in production"),
            ("mu_p_opening_angle_deg", "Muon-proton opening angle [deg]",
             np.linspace(0, 180, 37), "presel",
             lt("max_opening_angle", "Cut: angle $< %g^\\circ$"), None),
        ],
        "track_score": [
            ("mu_trackScore", "Muon-candidate track score", np.linspace(0.0, 1.0, 41),
             "cosmic", gt("musel_track_score_min", "Cut: score $> %g$"), None),
        ],
        "mu_len": [
            ("mu_len", "Muon-candidate length [cm]", np.linspace(0, 500, 26), "track_score",
             between("musel_len_th_min", "musel_len_th_max",
                     "Min: $L_\\mu > %g$ cm (production presel.)", "Max: $L_\\mu < %g$ cm"), None),
        ],
        "mu_cand": [
            ("mu_chi2_of_mu_cand", "Muon-candidate $\\chi^2_\\mu$", np.linspace(0, 115, 24),
             "twoprong", lt("musel_muscore_th", "Cut: $\\chi^2_\\mu < %g$"), "no cut in ICARUS"),
            ("prot_chi2_of_mu_cand", "Muon-candidate $\\chi^2_p$", np.linspace(0, 300, 31),
             "twoprong", gt("musel_pscore_th", "Cut: $\\chi^2_p > %g$"), None),
        ],
        "p_cand": [
            ("mu_chi2_of_prot_cand", "Proton-candidate $\\chi^2_\\mu$", np.linspace(0, 100, 26),
             "mu_cand", nocut, "no cut"),
            ("prot_chi2_of_prot_cand", "Proton-candidate $\\chi^2_p$", np.linspace(0, 300, 31),
             "mu_cand", lt("psel_pscore_th", "Cut: $\\chi^2_p < %g$"), None),
        ],
    }


def plot_cut_pages(dfs, masks, cats, plotdir, save, log):
    pages = cut_pages()
    for name, specs in pages.items():
        nrow = len(specs)
        fig, axes = plt.subplots(nrow, 2, figsize=(13, 3.6 * nrow), squeeze=False)
        for r, (var, xlabel, bins, stage, cutspec, note) in enumerate(specs):
            for c, det in enumerate(["SBND", "ICARUS"]):
                ax = axes[r, c]
                d = dfs[det]
                m = masks[det][stage]
                x = d[var].to_numpy(dtype=float)[m]
                w = d.glob_scale.to_numpy(dtype=float)[m]
                ok = np.isfinite(x)
                stacked_hist(ax, x[ok], w[ok], cats[det][m][ok], FS_CATS, bins)
                _style(ax, xlabel, pot_label(det), "%s: %s" % (det, xlabel.split(" [")[0]))
                lines, keep = cutspec(cuts_for(det))
                ax.set_ylim(0, ax.get_ylim()[1] * 1.15)
                cut_lines(ax, lines, keep)
                if note and (not lines or (det == "ICARUS" and "ICARUS" in note)):
                    ax.text(0.98, 0.97, note, transform=ax.transAxes, ha="right", va="top",
                            fontsize=FONTSIZE - 2, style="italic")
                ax.legend(fontsize=FONTSIZE - 4, ncol=2, loc="upper right",
                          bbox_to_anchor=(0.98, 0.88), frameon=True)
        fig.tight_layout()
        if save:
            for ext in ("png", "pdf"):
                fig.savefig(os.path.join(plotdir, ext, "cut_breakdown_%s.%s" % (name, ext)),
                            dpi=200 if ext == "png" else None, bbox_inches="tight")
        plt.close(fig)
        log("  wrote cut_breakdown_%s" % name)


# ----------------------------------------------------------------------------
# Figure group 2: efficiencies
# ----------------------------------------------------------------------------
def _numucc(d, pdgcol):
    return d.is_cc.astype(bool).to_numpy() & (np.abs(d[pdgcol].to_numpy(dtype=float)) == 14) \
        if "is_cc" in d.columns else \
        d.true_iscc.fillna(0).astype(bool).to_numpy() & (np.abs(d[pdgcol].to_numpy(dtype=float)) == 14)


def _num_1u1p(d):
    """Reconstructed rows: numu CC, exactly one proton above the KE threshold
    (leading above, second-leading below or absent), no charged pion, no pi0."""
    p1 = d.true_p_p.to_numpy(dtype=float)
    p2 = d.true_p2_p.to_numpy(dtype=float)
    one_p = (np.nan_to_num(p1) > P_KE_THRESHOLD) & ~(np.nan_to_num(p2) > P_KE_THRESHOLD)
    no_pi = d.true_cpi_p.isna().to_numpy() & d.true_pi0_p.isna().to_numpy()
    return _numucc(d, "true_pdg") & one_p & no_pi


def _den_1u1p(t):
    return _numucc(t, "pdg") & (t.n_p50.to_numpy() == 1) & (t.n_cpi.to_numpy() == 0) \
        & (t.n_pi0.to_numpy() == 0)


# (title, numerator class fn [evt frame], denominator class fn [truth frame],
#  numerator x, denominator x, xlabel, bins). The particle classes mean "a
# final-state particle of that kind exists"; QE and 1u1p are event classes,
# binned in true neutrino energy.
EFF_SPECS = [
    ("Muon efficiency",
     lambda d: _pos(d.true_mu_p).to_numpy(), lambda t: _pos(t.mu_p).to_numpy(),
     "true_mu_p", "mu_p", "True muon momentum [GeV/$c$]", np.linspace(0.0, 2.5, 26)),
    ("Proton efficiency",
     lambda d: _pos(d.true_p_p).to_numpy(), lambda t: _pos(t.p_p).to_numpy(),
     "true_p_p", "p_p", "True leading-proton momentum [GeV/$c$]", np.linspace(0.0, 2.0, 21)),
    ("Charged-pion efficiency",
     lambda d: _pos(d.true_cpi_p).to_numpy(), lambda t: _pos(t.cpi_p).to_numpy(),
     "true_cpi_p", "cpi_p", "True charged-pion momentum [GeV/$c$]", np.linspace(0.0, 1.5, 21)),
    ("$\\pi^0$-event efficiency",
     lambda d: _pos(d.true_pi0_p).to_numpy(), lambda t: (t.n_pi0.to_numpy() > 0),
     "nu_E", "nu_E", "True neutrino energy [GeV]", np.linspace(0.0, 3.0, 21)),
    ("$\\nu_\\mu$ CC QE efficiency",
     lambda d: _numucc(d, "true_pdg") & (d.genie_mode.to_numpy(dtype=float) == 0),
     lambda t: _numucc(t, "pdg") & (t.genie_mode.to_numpy(dtype=float) == 0),
     "nu_E", "nu_E", "True neutrino energy [GeV]", np.linspace(0.0, 3.0, 21)),
    ("$\\nu_\\mu$ CC 1$\\mu$1p efficiency ($T_p > 50$ MeV, no $\\pi$)",
     _num_1u1p, _den_1u1p,
     "nu_E", "nu_E", "True neutrino energy [GeV]", np.linspace(0.0, 3.0, 21)),
]


def true_in_fv(d):
    """Reconstructed rows whose TRUE vertex lies in the analysis FV (same
    gumple_cuts._fv_cut as the truth denominator's in_fv)."""
    vtx = pd.DataFrame({"detector": d.detector, "Run": d.Run, "x": d.true_vtx_x,
                        "y": d.true_vtx_y, "z": d.true_vtx_z}, index=d.index)
    with np.errstate(invalid="ignore"):
        return gc._fv_cut(vtx).to_numpy() & d.true_vtx_x.notna().to_numpy()


def plot_efficiency(dfs, truths, masks, plotdir, save, log):
    """Absolute selection efficiency: selected MC events with true vertex in the
    FV, over ALL generated neutrinos with true vertex in the FV (mcnu), both
    POT- and aFF-weighted, per final-state class."""
    fig, axes = plt.subplots(3, 2, figsize=(11, 12))
    ymax = 0.0
    for ax, (title, ncls, dcls, nx, dx, xlabel, bins) in zip(axes.flat, EFF_SPECS):
        for det, color in [("SBND", "tab:blue"), ("ICARUS", "tab:orange")]:
            d, t = dfs[det], truths[det]
            # numerator: selected MC rows, true vertex in FV, class present
            num = (d["sample"] == "mc").to_numpy() & masks[det]["final"] & true_in_fv(d) \
                & ncls(d)
            xn = d[nx].to_numpy(dtype=float)
            wn = d.glob_scale.to_numpy(dtype=float)
            # denominator: every generated neutrino, true vertex in FV, class present
            den = t.in_fv.to_numpy() & dcls(t)
            xd = t[dx].to_numpy(dtype=float)
            wd = t.glob_scale.to_numpy(dtype=float)
            nnum, _ = np.histogram(xn[num & np.isfinite(xn)], bins=bins, weights=wn[num & np.isfinite(xn)])
            nden, _ = np.histogram(xd[den & np.isfinite(xd)], bins=bins, weights=wd[den & np.isfinite(xd)])
            with np.errstate(divide="ignore", invalid="ignore"):
                eff = np.where(nden > 0, nnum / nden, np.nan)
            centers = 0.5 * (bins[1:] + bins[:-1])
            ax.plot(centers, eff, marker="o", ms=4, color=color, label=det)
            ymax = max(ymax, np.nanmax(eff) if np.isfinite(eff).any() else 0)
            log("  %-26s %-6s selected %.1f / generated-in-FV %.1f = %.4f" % (
                title.replace("$", ""), det, nnum.sum(), nden.sum(),
                nnum.sum() / nden.sum() if nden.sum() > 0 else np.nan))
        _style(ax, xlabel, "Selected / all generated $\\nu$ in FV", title)
        ax.legend(fontsize=FONTSIZE - 1)
    for ax in axes.flat:
        top = max(l.get_ydata()[np.isfinite(l.get_ydata())].max() if np.isfinite(l.get_ydata()).any()
                  else 0 for l in ax.get_lines())
        ax.set_ylim(0, min(1.0, top * 1.2) if top > 0 else 1.0)
    fig.suptitle("GUMP 1$\\mu$1p selection efficiency (denominator: every generated "
                 "neutrino interaction with true vertex in the fiducial volume)",
                 fontsize=FONTSIZE - 1, y=1.0)
    fig.tight_layout()
    if save:
        for ext in ("png", "pdf"):
            fig.savefig(os.path.join(plotdir, ext, "efficiency.%s" % ext),
                        dpi=200 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    log("  wrote efficiency")


# ----------------------------------------------------------------------------
# Figure group 3: near/far signal-box distributions
# ----------------------------------------------------------------------------
def plot_near_far(dfs, masks, mcats, var, bins, xlabel, short, plotdir, save, log,
                  slice_var=None, slice_range=None, tag=""):
    fig = plt.figure(figsize=(16, 4.4))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.05], height_ratios=[3, 1.2], hspace=0.05,
                          wspace=0.3)
    ax_s = fig.add_subplot(gs[:, 0])
    ax_i = fig.add_subplot(gs[:, 1])
    ax_t = fig.add_subplot(gs[0, 2])
    ax_r = fig.add_subplot(gs[1, 2], sharex=ax_t)
    order = list(range(len(MODE_CATS)))
    totals = {}
    for ax, det in [(ax_s, "SBND"), (ax_i, "ICARUS")]:
        d = dfs[det]
        m = masks[det]["final"].copy()
        if slice_var is not None:
            sv = d[slice_var].to_numpy(dtype=float)
            m &= (sv >= slice_range[0]) & (sv < slice_range[1])
        x = d[var].to_numpy(dtype=float)[m]
        w = d.glob_scale.to_numpy(dtype=float)[m]
        ok = np.isfinite(x)
        totals[det] = stacked_hist(ax, x[ok], w[ok], mcats[det][m][ok], MODE_CATS, bins, order)
        _style(ax, xlabel, pot_label(det), "$\\bf{%s}$  Signal box: %s" % (det, short))
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], fontsize=FONTSIZE - 4, loc="upper right")
        if slice_var is not None:
            ax.text(0.03, 0.95, "$%g \\leq \\delta p < %g$ GeV/c" % slice_range,
                    transform=ax.transAxes, ha="left", va="top", fontsize=FONTSIZE - 2)
    centers = 0.5 * (bins[1:] + bins[:-1])
    for det, color in [("SBND", "tab:blue"), ("ICARUS", "tab:orange")]:
        ax_t.step(bins, np.append(totals[det], totals[det][-1]), where="post", color=color,
                  label=det, lw=1.8)
    _style(ax_t, "", "Weighted events", "Signal box: %s: totals" % short)
    ax_t.legend(fontsize=FONTSIZE - 2)
    plt.setp(ax_t.get_xticklabels(), visible=False)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(totals["ICARUS"] > 0, totals["SBND"] / totals["ICARUS"], np.nan)
    ax_r.step(bins, np.append(ratio, ratio[-1]), where="post", color="tab:blue", lw=1.8)
    ax_r.axhline(1.0, color="gray", linestyle=":", lw=1)
    _style(ax_r, xlabel, "SBND / ICARUS", "")
    ax_r.set_ylim(0, max(2.5, np.nanmax(ratio) * 1.1) if np.isfinite(ratio).any() else 2.5)
    stem = "signalbox_%s%s_near_far" % (var, tag)
    if save:
        for ext in ("png", "pdf"):
            fig.savefig(os.path.join(plotdir, ext, "%s.%s" % (stem, ext)),
                        dpi=200 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    log("  wrote %s" % stem)


def plot_signalbox(dfs, masks, mcats, plotdir, save, log):
    delp_bins = np.linspace(0.0, 0.6, 13)
    ereco_bins = np.linspace(0.2, 1.5, 16)
    plot_near_far(dfs, masks, mcats, "del_p", delp_bins, "$\\delta p$ [GeV/c]", "$\\delta p$",
                  plotdir, save, log)
    plot_near_far(dfs, masks, mcats, "nu_E_calo", ereco_bins, "$E_{\\rm reco}$ [GeV]",
                  "$E_{\\rm reco}$", plotdir, save, log)
    for tag, lo, hi in mdc.DP_SLICES:
        plot_near_far(dfs, masks, mcats, "nu_E_calo", ereco_bins, "$E_{\\rm reco}$ [GeV]",
                      "$E_{\\rm reco}$", plotdir, save, log,
                      slice_var="del_p", slice_range=(lo, hi), tag="_" + tag)


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--df-dir", default=DEFAULT_DF_DIR)
    p.add_argument("--plotdir", default=DEFAULT_PLOTDIR)
    p.add_argument("--only", action="append", choices=["cuts", "eff", "signalbox"], default=None,
                   help="figure groups to make (default: all). Repeatable.")
    p.add_argument("--no-dirt", action="store_true")
    p.add_argument("--no-offbeam", action="store_true")
    p.add_argument("--no-save", action="store_true")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    df_dir = args.df_dir if args.df_dir.endswith("/") else args.df_dir + "/"
    groups = set(args.only or ["cuts", "eff", "signalbox"])
    save = not args.no_save
    plotdir = args.plotdir
    if save:
        for sub in ("", "png", "pdf"):
            os.makedirs(os.path.join(plotdir, sub), exist_ok=True)
    lines = []

    def log(msg):
        print(msg, flush=True)
        lines.append(str(msg))

    log("=" * 70)
    log("selection_plots_gumple: df-dir %s  groups %s  dirt %s  offbeam %s" % (
        df_dir, ",".join(sorted(groups)), not args.no_dirt, not args.no_offbeam))
    log("=" * 70)

    dfs, truths, masks, fcats, mcats = {}, {}, {}, {}, {}
    for det in ["SBND", "ICARUS"]:
        dfs[det], truths[det] = load_detector(det, df_dir, log, include_dirt=not args.no_dirt,
                                             include_offbeam=not args.no_offbeam,
                                             with_truth="eff" in groups)
        masks[det] = stage_masks(dfs[det], log)
        fcats[det] = final_state(dfs[det])
        mcats[det] = mode_category(dfs[det])
        d = dfs[det]
        w = d.glob_scale.to_numpy(dtype=float)
        log("  %s stage yields (weighted): %s" % (
            det, ", ".join("%s %.1f" % (k, w[v].sum()) for k, v in masks[det].items())))
        log("  %s final-state fractions after selection: %s" % (det, ", ".join(
            "%s %.3f" % (FS_CATS[i][0].replace("$", "").replace("\\", ""),
                         w[masks[det]["final"] & (fcats[det] == i)].sum() /
                         max(w[masks[det]["final"]].sum(), 1e-9))
            for i in range(len(FS_CATS)))))

    if "cuts" in groups:
        log("\ncut breakdown pages")
        plot_cut_pages(dfs, masks, fcats, plotdir, save, log)
    if "eff" in groups:
        log("\nefficiencies")
        plot_efficiency(dfs, truths, masks, plotdir, save, log)
    if "signalbox" in groups:
        log("\nnear/far signal-box distributions")
        plot_signalbox(dfs, masks, mcats, plotdir, save, log)

    if save:
        with open(os.path.join(plotdir, "selection_summary.txt"), "w") as f:
            f.write("\n".join(lines) + "\n")
    log("\ndone -> %s" % plotdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
