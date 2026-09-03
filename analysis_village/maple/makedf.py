"""MAPLE (1mu + N>1 p) dataframe builders for cafpyana.
Detector support:
  ICARUS -- the original selection, bit-compatible with the CAFANA port.
  SBND   -- generalization: SBND geometry for FV (high-YZ volume dropped,
            same 10 cm insets / 50 cm zback) and containment; the
            ICARUS-only cosmic rejection (CRT top veto, cryo-light,
            bar-flash) passes trivially; the GUMP cathode-crossing veto
            (no track in the slice may cross the x=0 cathode) joins the
            containment stage; PID uses SBND-calibrated chi2 with the
            same MAPLE thresholds as ICARUS.

Post-hoc candidate scheme:
  Candidates are FIXED at df production; the chi2 cuts can be re-applied
  after df production under any calorimetric chi2 variation on the fixed
  candidates (the GUMP SignalBoxSystematics pattern; see
  maple_sel.maple_selection).
    - Muon candidate: legacy-GUMP rule -- the track (trackScore > 0)
      starting within 10 cm of the vertex, longer than 40 cm, with the
      smallest NOMINAL chi2u/chi2p ratio; the chi2 cuts on it are applied
      post-hoc (calo variations re-evaluate the cuts, not the choice).
    - Candidate protons: every other pfp passing the id_pfp gates with
      dist_start < 10 -- exactly the set the old id_pfp split into
      pion/proton by chi2.  Kinematics (recoE, ...) assume all of them are
      protons; the evt df stores the worst-case cut variable over them
      (mu_chi2{v}_of_prot_cand = min chi2u, prot_chi2{v}_of_prot_cand =
      max chi2p, min_proton_ke, min_proton_track_score) so the post-hoc
      selection can require that every candidate is a proton.
  Column naming follows GUMP (slc_vtx_*, nu_E, mu_*/p_* candidate
  variables, {mu,prot}_chi2{v}_of_{mu,prot}_cand for the nominal chi2, all
  calorimetric variations, and the CAFANA-compat chi2 under v="cafana").

  Intentional differences from the old (chi2-in-candidate) selection:
    - A vertex-attached track with ke_proton < 40 MeV was previously
      classified by its shower energy (possibly UNKNOWN -> ignored); now it
      is a candidate proton and fails the min_proton_ke cut.
    - A slice whose longest chi2-free muon candidate fails the chi2 cuts is
      rejected instead of falling back to a shorter track.
    - The pid_mode/_alt machinery is gone: cafana-PID selection is
      maple_sel.maple_selection(df, "cafana").

Selection option:
  selection="none"   -- keep all slices, with all cut booleans (default)
  selection="presel" -- keep slices passing the PID-free MAPLE preselection
                        (sanity + FV + all-track containment + cathode + a
                        muon and >=1 candidate proton; no CRT veto or
                        cryo-light)
  selection="full"   -- keep slices passing the full MAPLE selection
                        (evaluated with the nominal chi2)
"""
import numpy as np
import pandas as pd

from pyanalib.pandas_helpers import *
from makedf.util import *
from makedf.makedf import (
    loadbranches, make_slcdf, make_trkdf, make_trkhitdf, make_crthitdf, make_opflashdf,
    make_hdrdf, make_triggerdf, make_potdf_bnb, make_mcnudf,
    make_genie_evtrec_df, _build_genie_evtrec_df,
)
from makedf import chi2pid

from analysis_village.gump.kinematics import *
import analysis_village.gumple.gumple_cuts as gmpl
from analysis_village.maple import chi2pid_cafana
from makedf.branches import (
    crtpmtbranches, shwbranches,
    mcbranches, mcprimbranches, mcprimvisEbranches, trueparticlebranches,
)

# PFPID enum from the helper header
PID_UNKNOWN, PID_PROTON, PID_PION, PID_SHOWER, PID_OTHER = 0, 1, 2, 3, 4
# TruthClass: 0=1mu1p, 1=1muNp, 2=Other, 3=Cosmic, 4=Invalid
CLS_1MU1P, CLS_1MUNP, CLS_OTHER, CLS_COSMIC, CLS_INVALID = 0, 1, 2, 3, 4

# All chi2 variation suffixes stored in the evt df (when do_calo_syst=True),
# plus the CAFANA-compat chi2 stored under the same scheme.
SCALE_SMEAR_VARIATIONS = ["lo", "hi", "2lo", "2hi", "smear5", "smear13", "sqsmear15"]
CALO_VARIATIONS = ["cv", "alpha_p", "alpha_m", "beta_p", "beta_m", "R_p", "R_m", "R_p25", "dedxbias"]
CHI2_VARIATIONS = SCALE_SMEAR_VARIATIONS + CALO_VARIATIONS + ["cafana"]

# The four per-slice chi2 candidate columns, GUMP naming ("%s" is the
# variation suffix; "" is nominal).
CHI2_CAND_COLS = [
    "mu_chi2%s_of_mu_cand",
    "prot_chi2%s_of_mu_cand",
    "mu_chi2%s_of_prot_cand",
    "prot_chi2%s_of_prot_cand",
]


def _flatcols(df):
    """Flatten a multiindex-column dataframe to underscore-joined names."""
    df = df.copy()
    df.columns = ["_".join([str(c) for c in (col if isinstance(col, tuple) else (col,)) if str(c) != ""])
                  for col in df.columns]
    return df

def _reindex(series, index, fill):
    return series.reindex(index).fillna(fill)


# =====================================================================
# Truth classification (classification_type_MC port)
# =====================================================================
def maple_truth_classdf(f, det, run):
    """Per-(entry, mc.nu index) MAPLE truth classification.

    Returns a DataFrame with columns:
      maple_class (0=1mu1p, 1=1muNp, 2=Other, 4=Invalid -- Cosmic is a
      slice-level concept, handled in the evt df), n_mu, mu_length,
      veto, uncontained, true_visible_Enu
    """
    mc = _flatcols(loadbranches(f["recTree"], mcbranches).rec.mc.nu)
    mc["detector"] = det
    mc["Run"] = run
    if mc.empty:
        out = pd.DataFrame(columns=["maple_class", "n_mu", "mu_length",
                                    "veto", "uncontained", "true_visible_Enu"])
        return out

    nuidx = mc.index

    prim = _flatcols(loadbranches(f["recTree"], mcprimbranches + mcprimvisEbranches).rec.mc.nu.prim)
    prim["detector"] = det
    prim["Run"] = run
    prim.index.names = ["entry", "inu", "iprim"]

    tp = _flatcols(loadbranches(f["recTree"], trueparticlebranches).rec.true_particles)
    tp.index.names = ["entry", "itp"]

    # Only primaries with G4ID >= 0 and cryostat >= 0 participate (C++ `continue`)
    prim = prim[(prim.G4ID >= 0) & (prim.cryostat >= 0)]

    apdg = np.abs(prim.pdg)
    is_mu = apdg == 13
    is_p = apdg == 2212
    is_cpi = apdg == 211
    is_pi0 = apdg == 111
    is_n = apdg == 2112 
    is_gamma = apdg == 22
    is_e = apdg == 11
    is_charged = is_mu | is_p | is_cpi | is_e

    # own deposited energy on plane 2 in the primary's cryostat [MeV]
    visE_own = np.where(prim.cryostat == 0, prim.plane_I0_I2_visE, prim.plane_I1_I2_visE)
    prim = prim.assign(visE_own=visE_own)

    # ---- daughters: true_particles with parent == prim.G4ID (same entry) ----
    prim_r = prim.reset_index()
    tp_r = tp.reset_index()
    m = prim_r[["entry", "inu", "iprim", "G4ID", "cryostat", "detector", "Run"]].merge(
        tp_r, left_on=["entry", "G4ID"], right_on=["entry", "parent"],
        suffixes=("", "_d"))
    if len(m):
        # daughter visE evaluated at the PRIMARY's cryostat (as in C++)
        m["d_visE"] = np.where(m.cryostat == 0, m.plane_I0_I2_visE, m.plane_I1_I2_visE)
        m["d_apdg"] = np.abs(m.pdg)
        m["d_is_gamma"] = m.d_apdg == 22
        m["d_charged"] = m.d_apdg.isin([13, 2212, 211, 11])
        # daughter containment check: skip cryostat<0 or end == -9999
        d_valid_cont = (m.cryostat_d >= 0) & m.d_charged & \
            ~((m.end_x == -9999) | (m.end_y == -9999) | (m.end_z == -9999))
        m["d_uncont"] = d_valid_cont & ~gmpl.prefix_fv_cut(m, "end")
        # note: for the pi0-gamma veto and visE sums, C++ has no cryostat check
        # on the daughter itself
        g = m.groupby(["entry", "inu", "iprim"])
        dep_daughters = g.d_visE.sum()
        pi0_gamma_veto = g.apply(lambda x: bool(((x.d_is_gamma) & (x.d_visE > 0.0)).any())) # use PION_KE_MIN
        d_uncont_any = g.d_uncont.any()
    else:
        dep_daughters = pd.Series(dtype=float)
        pi0_gamma_veto = pd.Series(dtype=bool)
        d_uncont_any = pd.Series(dtype=bool)

    dep_daughters = _reindex(dep_daughters, prim.index, 0.0) # use PION_KE_MIN
    pi0_gamma_veto = _reindex(pi0_gamma_veto, prim.index, False).astype(bool)
    d_uncont_any = _reindex(d_uncont_any, prim.index, False).astype(bool)

    # ---- per-primary vetoes ----
    cpi_veto = is_cpi & (prim.visE_own > 0.0) # use PION_KE_MIN
    pi0_veto = is_pi0 & pi0_gamma_veto
    gamma_veto = is_gamma & ((prim.visE_own + dep_daughters) > 0.0) # use PION_KE_MIN
    p_depE = prim.visE_own + dep_daughters
    prim_veto = cpi_veto | pi0_veto | gamma_veto

    # ---- containment (all_contained_mc): charged primaries, no -9999 skip ----
    prim_uncont = is_charged & ~gmpl.prefix_fv_cut(prim, "end")

    # ---- per-nu aggregation ----
    grp = prim.groupby(level=["entry", "inu"])

    def agg(series, fill=False):
        s = series.groupby(level=["entry", "inu"]).any() if series.dtype == bool else \
            series.groupby(level=["entry", "inu"]).sum()
        s.index.names = nuidx.names
        return _reindex(s, nuidx, fill)

    n_mu = agg(is_mu.astype(int), 0)
    n_p = agg(is_p.astype(int), 0)
    n_pi = agg(is_cpi.astype(int), 0)
    n_pi0 = agg(is_pi0.astype(int), 0)
    n_n = agg(is_n.astype(int), 0)
    veto_any = agg(prim_veto, False).astype(bool)
    uncont_any = (agg(prim_uncont, False).astype(bool) |
                  agg(d_uncont_any, False).astype(bool))

    # last muon in loop order sets muon_length / E_mu (C++ overwrites per muon)
    mu_rows = prim[is_mu]
    mu_last = mu_rows.groupby(level=["entry", "inu"]).last()
    mu_last.index.names = nuidx.names
    mu_length = _reindex(mu_last.length, nuidx, np.nan)
    mu_p_GeV = np.sqrt(mu_last.genp_x**2 + mu_last.genp_y**2 + mu_last.genp_z**2)
    E_mu_vis = np.sqrt(mu_p_GeV**2 + MUON_MASS**2)  # MeV
    E_mu_vis = _reindex(E_mu_vis, nuidx, 0.0)

    # visible energy: protons 
    p_rows = prim[is_p]
    p_ke = kinetic_energy(PROTON_MASS, np.sqrt(p_rows.genp_x**2 + p_rows.genp_y**2 + p_rows.genp_z**2))
    E_p_sum = (p_ke + BE).groupby(level=["entry", "inu"]).sum()
    E_p_sum.index.names = nuidx.names
    E_p_sum = _reindex(E_p_sum, nuidx, 0.0)

    true_visible_Enu = (E_p_sum + E_mu_vis)  # GeV

    # ---- classification ----
    pos_nan = mc.position_x.isna() | mc.position_y.isna() | mc.position_z.isna()
    not_numucc = (np.abs(mc.pdg) != 14) | (mc.iscc == 0)
    not_fv = ~gmpl.prefix_fv_cut(mc, "position")
    #good_mu = (n_mu == 1) & (mu_length > MIN_MUON_LENGTH) & (mu_length < MAX_MUON_LENGTH)

    maple_class = np.select(
        [pos_nan,
         not_numucc | not_fv | veto_any | uncont_any,
         (n_mu == 1) & (n_p == 1),
         (n_mu == 1) & (n_p > 1)],
         #good_mu & (n_p == 1),
         #good_mu & (n_p > 1)],
        [CLS_INVALID, CLS_OTHER, CLS_1MU1P, CLS_1MUNP],
        default=CLS_OTHER)

    return pd.DataFrame({
        "maple_class": maple_class,
        "nmu": n_mu,
        "np": n_p,
        "npi": n_pi,
        "npi0": n_pi0,
        "nn": n_n,
        "mu_length": mu_length,
        "veto": veto_any,
        "uncontained": uncont_any,
        "true_visible_Enu": true_visible_Enu,
    }, index=nuidx)


# =====================================================================
# Per-pfp candidate machinery (chi2-free; chi2 cuts live in maple_sel)
# =====================================================================
def _find_candidates(P, use_chi2=True):
    """Muon + candidate-pfp finding and fixed pfp counting.

    P: flat per-pfp frame (index entry, slc, pfp), with the nominal chi2
    columns already attached (PID_calcs runs first).
    Returns (mu_ilocs [per-slice Index of muon rows],
             is_prot_cand [bool Series over P],
             counts DataFrame per slice [n_pfp, n_pfp_no_calo, n_shower, n_other]).

    Muon candidate, selected by `use_chi2`:
      use_chi2=True (default) -- legacy-GUMP rule: among tracks with any
        track score (trackScore > 0) starting within 10 cm of the vertex --
        with the addendum that the track must be longer than 40 cm -- take
        the smallest NOMINAL chi2u/chi2p ratio.  NaN ratios (no-calo
        tracks) sort last, so a slice whose only qualifying tracks lack
        calo still yields a muon, as in legacy GUMP.  The candidate is
        FIXED under calorimetric variations (legacy semantics: variations
        re-evaluate the chi2 cuts on the fixed candidate, never the
        choice).
      use_chi2=False -- chi2-free MAPLE rule: the longest track passing
        get_base_muon_mask (trackScore >= 0.5, dist_start <= 10,
        len in [40, 400], primary, contained, same-side).

    A candidate pfp is any pandora pfp that is primary, has calo points,
    and starts within 10 cm of the vertex; n_pfp counts the candidate pfps
    in the slice (the muon included).  Candidate protons are the non-muon
    candidate pfps -- the set whose worst-case chi2 aggregates the proton
    PID cut is applied to.  The remaining pfps keep the chi2-independent
    shower/other/unknown classification, so the counts are fixed under
    calorimetric variations.
    """

    if use_chi2:
        keep = (P.trackScore > 0.0) & (P.dist_start < 10.0) & (P.len > 40.0)
    else:
        keep = gmpl.get_base_muon_mask(P, level="trk")

    cand = P[keep]
    if len(cand) and use_chi2:
        ratio = cand.chi2u_cafpyana / cand.chi2p_cafpyana
        first = cand.assign(_mu_ratio=ratio).sort_values(
            "_mu_ratio", na_position="last", kind="stable").groupby(level=[0, 1]).head(1)
        mu_ilocs = pd.Series(
            list(first.index),
            index=pd.MultiIndex.from_arrays(
                [first.index.get_level_values(0), first.index.get_level_values(1)],
                names=P.index.names[:2]))
    elif len(cand):
        mu_ilocs = cand.len.groupby(level=[0, 1]).idxmax()
    else:
        mu_ilocs = pd.Series(dtype=object)

    is_mu = pd.Series(False, index=P.index)
    if len(mu_ilocs):
        is_mu.loc[pd.Index(mu_ilocs.values)] = True

    # ---- candidate pfps: primary, has calo points, starts within 10 cm ----
    # NaN semantics: NaN ncalo/dist_start fail the comparisons -> not a candidate.
    is_cand = P.prim_pfp.astype(bool) & (P.ncalo > 0) & (P.dist_start < 10.0)
    is_prot_cand = is_cand & ~is_mu

    # would-be candidates without calo: primary, starts within 10 cm, but no
    # calo points -- invisible to n_pfp/PID (classified "unknown" below)
    is_cand_no_calo = P.prim_pfp.astype(bool) & (P.ncalo == 0) & (P.dist_start < 10.0)

    # ---- id_pfp gates (the chi2-free skip conditions), kept for the
    # shower/other/unknown classification of the remaining pfps ----
    unknown0 = (~P.prim_pfp) | P.start_x.isna() | P.end_x.isna() | P.len.isna()
    unknown1 = P.min_dist > 50.0 # VTX_MAX_DIST
    no_calo = P.ncalo == 0
    gate = ~unknown0 & ~unknown1 & ~no_calo

    # remaining pfps: chi2-independent shower/other/unknown classification
    shw_unknown = P.shw_energy2.isna()
    is_shower = P.shw_energy2 > 0.0 # PION_KE_MIN
    pid_rest = pd.Series(np.select(
        [~gate, shw_unknown, is_shower],
        [PID_UNKNOWN, PID_UNKNOWN, PID_SHOWER],
        default=PID_OTHER), index=P.index)
    pid_rest[is_mu | is_prot_cand] = -1

    grp = lambda s: s.groupby(level=[0, 1]).sum()
    counts = pd.DataFrame({
        # muon included via the union: the muon counts toward n_pfp even in
        # the (rare) case it misses the candidate definition itself
        "n_pfp": grp((is_cand | is_mu).astype(int)),
        "n_pfp_no_calo": grp(is_cand_no_calo.astype(int)),
        "n_shower": grp((pid_rest == PID_SHOWER).astype(int)),
        "n_other": grp((pid_rest == PID_OTHER).astype(int)),
    })

    return mu_ilocs, is_prot_cand, counts


def _proton_aggregates(P, is_prot_cand, mincols, maxcols):
    """Per-slice worst-case cut variables over the candidate protons.

    mincols are aggregated with min (for lower-bound cuts), maxcols with
    max (for upper-bound cuts).  skipna=False semantics: a NaN value on any
    candidate makes the aggregate NaN, so the downstream cut fails (the old
    per-pfp selection required chi2.notna()).

    Returns a DataFrame per slice with columns "min_<col>" / "max_<col>".
    """
    g = P[is_prot_cand].groupby(level=[0, 1])
    out = {}
    for cols, fn in ((mincols, "min"), (maxcols, "max")):
        for col in cols:
            cg = g[col]
            agg = getattr(cg, fn)()
            agg[cg.count() < cg.size()] = np.nan
            out["%s_%s" % (fn, col)] = agg
    return pd.DataFrame(out)

def add_crt(f, S):
    # Add in crt hit matching for ICARUS
    crt = make_crthitdf(f)
    S = S.join(((crt.time > -1) & (crt.time < 1.8) & (crt.plane != 50)).groupby(level=[0]).any().rename("crthit"))
    S = S.join(((crt.time > -1) & (crt.time < 1.8) & (crt.plane != 50) & (crt.truth.bestmatch_id != -1)).groupby(level=[0]).any().rename("crthit_ismc"))
    S["crthit"] = S.crthit.fillna(False).astype(bool)
    S["crthit_ismc"] = S.crthit_ismc.fillna(False).astype(bool)
    return S

def add_pmt(f, S, DETECTOR):
    # Flash value for trigger emulation. Note: these need to be scaled per-detector, per-Run
    flashes = make_opflashdf(f)
    if DETECTOR == "ICARUS":
        timename = "firsttime"
    elif DETECTOR == "SBND":
        timename = "time"
    intime = (flashes[timename] > -5) & (flashes[timename] < 5)
    maxpe = (flashes.totalpe*intime).groupby(level=[0]).max().rename("flash_maxpe")
    S = S.join(maxpe)
    sumpe = (flashes.totalpe*intime).groupby(level=[0]).sum().rename("flash_sumpe")
    S = S.join(sumpe)

    flash_cryo0 = (flashes.totalpe * intime * (flashes.cryo == 0)).groupby(level=[0]).max().rename("flash_maxpe_cryo0")
    flash_cryo1 = (flashes.totalpe * intime * (flashes.cryo == 1)).groupby(level=[0]).max().rename("flash_maxpe_cryo1")
    S = S.join(flash_cryo0)
    S = S.join(flash_cryo1)
    return S

def fetch_metadata(f):
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det

    if det.empty:
        return pd.DataFrame()
    if 1 == det.unique():
        DETECTOR = "SBND"
    elif 2 == det.unique():
        DETECTOR = "ICARUS"
    else:
        raise ValueError("MAPLE needs rec.hdr.det == 1 (SBND) or 2 (ICARUS); got %s" % det.unique())
    run = loadbranches(f["recTree"], ["rec.hdr.run"]).rec.hdr.run
    RUN = 1 if DETECTOR == "SBND" else (2 if run.iloc[0] < 12960 else 4)
    ismc = bool(loadbranches(f["recTree"], ["rec.hdr.ismc"]).rec.hdr.ismc.iloc[0])
    return DETECTOR, RUN, ismc

def fetch_info(f):
    DETECTOR, RUN, ismc = fetch_metadata(f)

    # ------------------------------------------------------------------
    # slice frame
    # ------------------------------------------------------------------
    slcdf = make_slcdf(f)
    S = pd.DataFrame({
        "slc_vtx_x": slcdf.slc.vertex.x,
        "slc_vtx_y": slcdf.slc.vertex.y,
        "slc_vtx_z": slcdf.slc.vertex.z,
        "charge_center_z": slcdf.slc.charge_center.z,
        "nu_score": slcdf.slc.nu_score,
        "crlongtrkdiry": slcdf.slc.nuid.crlongtrkdiry,
        "tmatch_idx": slcdf.slc.tmatch.idx,
        "tmatch_eff": slcdf.slc.tmatch.eff,
        "tmatch_pur": slcdf.slc.tmatch.pur,
        "nu_E": slcdf.slc.truth.E,
        "true_pdg": slcdf.slc.truth.pdg,
        "true_iscc": slcdf.slc.truth.iscc,
        "true_isnc": slcdf.slc.truth.isnc,
        "true_iscosmic": (slcdf.slc.truth.pdg == -1),
        "true_genie_mode": slcdf.slc.truth.genie_mode,
        "true_vtx_x": slcdf.slc.truth.position.x,
        "true_vtx_y": slcdf.slc.truth.position.y,
        "true_vtx_z": slcdf.slc.truth.position.z,
        # slice-level true signal particles (GUMP true_mu_*/true_p_*/true_p2_*)
        "true_mu_p": slcdf.slc.truth.mu.totp,
        "true_mu_dir_x": slcdf.slc.truth.mu.dir.x,
        "true_mu_dir_y": slcdf.slc.truth.mu.dir.y,
        "true_mu_dir_z": slcdf.slc.truth.mu.dir.z,
        "true_mu_end_x": slcdf.slc.truth.mu.end.x,
        "true_mu_end_y": slcdf.slc.truth.mu.end.y,
        "true_mu_end_z": slcdf.slc.truth.mu.end.z,
        "true_p_p": slcdf.slc.truth.p.totp,
        "true_p_dir_x": slcdf.slc.truth.p.dir.x,
        "true_p_dir_y": slcdf.slc.truth.p.dir.y,
        "true_p_dir_z": slcdf.slc.truth.p.dir.z,
        "true_p_end_x": slcdf.slc.truth.p.end.x,
        "true_p_end_y": slcdf.slc.truth.p.end.y,
        "true_p_end_z": slcdf.slc.truth.p.end.z,
        "true_p2_p": slcdf.slc.truth.p2.totp,
        "true_p2_dir_x": slcdf.slc.truth.p2.dir.x,
        "true_p2_dir_y": slcdf.slc.truth.p2.dir.y,
        "true_p2_dir_z": slcdf.slc.truth.p2.dir.z,
        "true_p2_end_x": slcdf.slc.truth.p2.end.x,
        "true_p2_end_y": slcdf.slc.truth.p2.end.y,
        "true_p2_end_z": slcdf.slc.truth.p2.end.z,
        # slice-level true leading charged pion / photon / neutral pion
        "true_cpi_p": slcdf.slc.truth.cpi.totp,
        "true_cpi_dir_x": slcdf.slc.truth.cpi.dir.x,
        "true_cpi_dir_y": slcdf.slc.truth.cpi.dir.y,
        "true_cpi_dir_z": slcdf.slc.truth.cpi.dir.z,
        "true_cpi_end_x": slcdf.slc.truth.cpi.end.x,
        "true_cpi_end_y": slcdf.slc.truth.cpi.end.y,
        "true_cpi_end_z": slcdf.slc.truth.cpi.end.z,
        "true_g_p": slcdf.slc.truth.g.totp,
        "true_g_dir_x": slcdf.slc.truth.g.dir.x,
        "true_g_dir_y": slcdf.slc.truth.g.dir.y,
        "true_g_dir_z": slcdf.slc.truth.g.dir.z,
        "true_g_end_x": slcdf.slc.truth.g.end.x,
        "true_g_end_y": slcdf.slc.truth.g.end.y,
        "true_g_end_z": slcdf.slc.truth.g.end.z,
        "true_pi0_p": slcdf.slc.truth.pi0.totp,
        "true_pi0_dir_x": slcdf.slc.truth.pi0.dir.x,
        "true_pi0_dir_y": slcdf.slc.truth.pi0.dir.y,
        "true_pi0_dir_z": slcdf.slc.truth.pi0.dir.z,
        "true_pi0_end_x": slcdf.slc.truth.pi0.end.x,
        "true_pi0_end_y": slcdf.slc.truth.pi0.end.y,
        "true_pi0_end_z": slcdf.slc.truth.pi0.end.z,
        "ismc" : ismc,
    })
    S["slice_index"] = S.index.get_level_values(1)
    S["detector"] = DETECTOR
    S["Run"] = RUN
    S["ismc"] = ismc

    # ------------------------------------------------------------------
    # per-pfp frame
    # ------------------------------------------------------------------
    trkdf = make_trkdf(f, False)  # NO vertex-distance pre-filter: MAPLE needs all pfps
    keys = set(f["recTree"].keys())
    shwdf = loadbranches(f["recTree"], [b for b in shwbranches if b in keys]).rec.slc.reco.pfp.shw
    shwdf.index.names = trkdf.index.names

    P = pd.DataFrame({
        "start_x": trkdf.pfp.trk.start.x,
        "start_y": trkdf.pfp.trk.start.y,
        "start_z": trkdf.pfp.trk.start.z,
        "end_x": trkdf.pfp.trk.end.x,
        "end_y": trkdf.pfp.trk.end.y,
        "end_z": trkdf.pfp.trk.end.z,
        "len": trkdf.pfp.trk.len,
        "dir_x": trkdf.pfp.trk.dir.x,
        "dir_y": trkdf.pfp.trk.dir.y,
        "dir_z": trkdf.pfp.trk.dir.z,
        "p_muon": trkdf.pfp.trk.rangeP.p_muon,
        "p_pion": trkdf.pfp.trk.rangeP.p_pion,
        "p_proton": trkdf.pfp.trk.rangeP.p_proton,
        "trackScore": trkdf.pfp.trackScore,
        "true_genp_x": trkdf.pfp.trk.truth.p.genp.x,
        "true_genp_y": trkdf.pfp.trk.truth.p.genp.y,
        "true_genp_z": trkdf.pfp.trk.truth.p.genp.z,
        "true_pdg": trkdf.pfp.trk.truth.p.pdg,
        "true_end_x": trkdf.pfp.trk.truth.p.end.x,
        "true_end_y": trkdf.pfp.trk.truth.p.end.y,
        "true_end_z": trkdf.pfp.trk.truth.p.end.z,
        "shw_energy2": shwdf.plane.I2.energy,
    })

    P["prim_pfp"] = trkdf.pfp.parent_is_primary.fillna(False).astype(bool)
    P["detector"] = DETECTOR
    P["Run"] = RUN
    P["ismc"] = ismc

    # broadcast slice vertex onto pfps
    P = P.join(S[["slc_vtx_x", "slc_vtx_y", "slc_vtx_z"]])
    P["dist_start"] = np.sqrt((P.start_x - P.slc_vtx_x)**2 + (P.start_y - P.slc_vtx_y)**2 + (P.start_z - P.slc_vtx_z)**2)
    dist_end = np.sqrt((P.end_x - P.slc_vtx_x)**2 + (P.end_y - P.slc_vtx_y)**2 + (P.end_z - P.slc_vtx_z)**2)

    # std::min(a, b) semantics: b if b < a else a  (NaN b -> a; NaN a -> NaN)
    P["min_dist"] = np.where(np.isnan(dist_end), P.dist_start, np.minimum(P.dist_start, dist_end))
    # track-endpoint containment: 10 cm insets on every face including z-back
    # (is_trk=True), matching the legacy GUMP mufv/pfv candidate-endpoint cut
    P["contained10"] = gmpl.prefix_fv_cut(P, "end", is_trk=True)
    P["ke_pion"] = kinetic_energy(PION_MASS, np.sqrt((P.dir_x * P.p_pion)**2 + (P.dir_y * P.p_pion)**2 + (P.dir_z * P.p_pion)**2))
    P["ke_proton"] = kinetic_energy(PROTON_MASS, np.sqrt((P.dir_x * P.p_proton)**2 + (P.dir_y * P.p_proton)**2 + (P.dir_z * P.p_proton)**2))
  
    if "ICARUS" in DETECTOR:
        S = add_crt(f, S)
    else:
        S["crthit"] = False 
        S["crthit_ismc"] = False 

    S = add_pmt(f, S, DETECTOR)

    return S, P, DETECTOR, RUN, ismc

# plane-flavor tags for the alternative chi2 calculations (do_alt_chi2)
#   p2trim/p2trim2/p2trim3 -- collection plane with the last 1.5/2/3 cm of the
#                             track removed (always built under do_alt_chi2)
#   p0        -- front induction plane (plane 0)      (do_ind_chi2 only)
#   p1        -- middle induction plane (plane 1)     (do_ind_chi2 only)
#   bestplane -- per-pfp, the plane (0/1/2) with the most calo hits (do_ind_chi2 only)
#
# collection-plane trim variations: tag -> rr_min (cm of track-end hits dropped)
TRIM_CHI2_TAGS = {"p2trim": 1.5, "p2trim2": 2.0, "p2trim3": 3.0}
# induction-plane + best-plane tags (only built when do_ind_chi2 is on;
# bestplane needs the p0/p1 chi2 to exist, so it rides with them)
IND_PLANE_TAGS = ["p0", "p1", "bestplane"]

# Default proton PID: the collection plane for proton candidates at a large
# angle to the drift (X) direction (theta_x > 47 deg), and the last-2cm-trimmed
# collection plane (p2trim2) for candidates nearly along the drift (theta_x <=
# 47 deg), where the reconstructed track end -- which the collection-plane chi2
# samples heavily -- is most distorted.  The muon candidate always uses the
# untrimmed collection plane.
PROT_TRIM_ANGLE_X_DEG = 47.0
DEFAULT_PROT_TRIM_TAG = "p2trim2"   # must be a key of TRIM_CHI2_TAGS

def _alt_plane_tags(do_ind_chi2):
    tags = list(TRIM_CHI2_TAGS.keys())
    if do_ind_chi2:
        tags = IND_PLANE_TAGS + tags
    return tags

def _chi2_flavors(do_calo_syst, do_cafana_chi2):
    """The per-plane chi2 flavors: nominal cafpyana (+ cafana when do_cafana_chi2
    is on), plus the full calorimetric variation suite when do_calo_syst is on."""
    flavors = ["cafpyana"]
    if do_cafana_chi2:
        flavors.append("cafana")
    if do_calo_syst:
        flavors = flavors + SCALE_SMEAR_VARIATIONS + CALO_VARIATIONS
    return flavors

def _add_dedx_variations(trkhitdf, DETECTOR, ismc, do_calo_syst):
    """Add the recomputed dE/dx column (dedx_redo) and, when do_calo_syst is on,
    all scale/smear and calorimetric-variation dedx columns to a plane's
    trkhitdf (in place). Ported from gump make_pandora_no_cuts_df, including
    the gump detector-specific scale sizes."""
    # gump-style dE/dx on recomputed dE/dx (detector gains + calibration)
    trkhitdf["dedx_redo"] = chi2pid.dedx(trkhitdf, gain=DETECTOR, calibrate=DETECTOR, isMC=ismc)
    if not do_calo_syst:
        return

    if DETECTOR == "ICARUS":
        scale_lo, scale_hi, scale_2lo, scale_2hi = 0.99, 1.01, 0.98, 1.02
        calo_var_params = chi2pid.ICARUS_CALO_VARIATIONS
    else:
        scale_lo, scale_hi, scale_2lo, scale_2hi = 0.98, 1.02, 0.96, 1.04
        calo_var_params = chi2pid.SBND_CALO_VARIATIONS
    trkhitdf["dedx_lo"] = chi2pid.dedx(trkhitdf, gain=DETECTOR, calibrate=DETECTOR, isMC=ismc, scale=scale_lo)
    trkhitdf["dedx_hi"] = chi2pid.dedx(trkhitdf, gain=DETECTOR, calibrate=DETECTOR, isMC=ismc, scale=scale_hi)
    trkhitdf["dedx_2lo"] = chi2pid.dedx(trkhitdf, gain=DETECTOR, calibrate=DETECTOR, isMC=ismc, scale=scale_2lo)
    trkhitdf["dedx_2hi"] = chi2pid.dedx(trkhitdf, gain=DETECTOR, calibrate=DETECTOR, isMC=ismc, scale=scale_2hi)
    trkhitdf["dedx_smear5"] = chi2pid.dedx(trkhitdf, gain=DETECTOR, calibrate=DETECTOR, isMC=ismc, smear=0.05)
    trkhitdf["dedx_smear13"] = chi2pid.dedx(trkhitdf, gain=DETECTOR, calibrate=DETECTOR, isMC=ismc, smear=0.13)
    trkhitdf["dedx_sqsmear15"] = chi2pid.dedx(trkhitdf, gain=DETECTOR, calibrate=DETECTOR, isMC=ismc, sqrt_smear=0.15)
    for c_var in CALO_VARIATIONS:
        if c_var == "dedxbias":
            # dedxbias is ICARUS-only: scale corrected dE/dx up by the dE/dx
            # spline. In SBND it is a no-op equal to CV.
            if DETECTOR == "ICARUS":
                trkhitdf["dedx_dedxbias"] = chi2pid.dedx(
                    trkhitdf, gain=DETECTOR, calibrate=DETECTOR, isMC=ismc,
                    dedx_bias=True)
            else:
                trkhitdf["dedx_dedxbias"] = trkhitdf["dedx_cv"]
        elif c_var == "R_p25":
            # R_p25 is SBND-only (large R+0.25 recombination test); no-op = CV on ICARUS
            if DETECTOR == "ICARUS":
                trkhitdf["dedx_R_p25"] = trkhitdf["dedx_cv"]
            else:
                trkhitdf["dedx_R_p25"] = chi2pid.dedx(
                    trkhitdf, gain=DETECTOR, calibrate=DETECTOR, isMC=ismc,
                    new_calo_params=calo_var_params["R_p25"])
        else:
            trkhitdf["dedx_%s" % c_var] = chi2pid.dedx(
                trkhitdf, gain=DETECTOR, calibrate=DETECTOR, isMC=ismc,
                new_calo_params=calo_var_params[c_var])

def _write_plane_chi2(trkhitdf, P, tag, do_calo_syst, cosmic, do_cafana_chi2, rr_min=None):
    """Write chi2u_{tag}<flavor> / chi2p_{tag}<flavor> columns into P for one
    plane's trkhitdf (dedx variation columns must already be attached via
    _add_dedx_variations). tag="" reproduces the collection-plane column names;
    otherwise use "p0_", "p1_", "p2trim_", etc. rr_min, if set, drops hits with
    rr < rr_min from the chi2 hit selection."""
    # gump-style chi2 on recomputed dE/dx
    P["chi2u_%scafpyana" % tag] = chi2pid.chi2u(trkhitdf, dedxname="dedx_redo", rr_min=rr_min)[0]
    P["chi2p_%scafpyana" % tag] = chi2pid.chi2p(trkhitdf, dedxname="dedx_redo", rr_min=rr_min)[0]

    # CAFANA-compat chi2 on stored dedx (opt-in)
    if do_cafana_chi2:
        cafana = chi2pid_cafana.chi2_cafana(trkhitdf, rr_min=(rr_min if rr_min is not None else 0.0))
        P["chi2u_%scafana" % tag] = cafana.chi2_mu
        P["chi2p_%scafana" % tag] = cafana.chi2_pro

    if do_calo_syst:
        for var in SCALE_SMEAR_VARIATIONS + CALO_VARIATIONS:
            P["chi2u_%s%s" % (tag, var)] = chi2pid.chi2u(trkhitdf, dedxname="dedx_%s" % var, rr_min=rr_min)[0]
            P["chi2p_%s%s" % (tag, var)] = chi2pid.chi2p(trkhitdf, dedxname="dedx_%s" % var, rr_min=rr_min)[0]

        # Don't apply variations to (Overlay) cosmics
        for var in SCALE_SMEAR_VARIATIONS + CALO_VARIATIONS:
            P.loc[cosmic, "chi2u_%s%s" % (tag, var)] = P.loc[cosmic, "chi2u_%scafpyana" % tag]
            P.loc[cosmic, "chi2p_%s%s" % (tag, var)] = P.loc[cosmic, "chi2p_%scafpyana" % tag]

def PID_calcs(f, P, DETECTOR, ismc, do_calo_syst=True, do_alt_chi2=False,
              do_cafana_chi2=False, do_ind_chi2=False):
    # ------------------------------------------------------------------
    # PID on the collection plane (plane 2)
    # ------------------------------------------------------------------
    trkhitdf = make_trkhitdf(f)

    # number of plane-2 calo points (compute_chi2 returns {} when empty)
    ncalo = trkhitdf.groupby(level=[0, 1, 2]).size()
    P["ncalo"] = _reindex(ncalo, P.index, 0).astype(int)

    # cosmic (Overlay) mask -- variations are pinned to nominal for these
    cosmic = P.true_genp_x.isna()

    # dedx + chi2 on the collection plane (tag="" keeps the legacy column names)
    _add_dedx_variations(trkhitdf, DETECTOR, ismc, do_calo_syst)
    _write_plane_chi2(trkhitdf, P, "", do_calo_syst, cosmic, do_cafana_chi2)

    # p2trim2 (collection plane, last 2 cm removed) is half of the default
    # proton PID (the angle blend), so it is ALWAYS computed, independent of
    # do_alt_chi2.  It reuses the plane-2 dedx columns attached just above.
    _write_plane_chi2(trkhitdf, P, DEFAULT_PROT_TRIM_TAG + "_", do_calo_syst,
                      cosmic, do_cafana_chi2, rr_min=TRIM_CHI2_TAGS[DEFAULT_PROT_TRIM_TAG])

    # ------------------------------------------------------------------
    # Alternative chi2 flavors: the collection plane with the last 1.5/3 cm of
    # the track removed and, when do_ind_chi2 is on, the front/middle induction
    # planes and the best plane. Each carries the same full variation suite as
    # the collection plane. (p2trim2 is always built just above.)
    # ------------------------------------------------------------------
    if do_alt_chi2:
        flavors = _chi2_flavors(do_calo_syst, do_cafana_chi2)

        # collection plane with the last 1.5/3 cm removed (reuse plane-2 dedx cols)
        for tag, rr_min in TRIM_CHI2_TAGS.items():
            if tag == DEFAULT_PROT_TRIM_TAG:
                continue  # already built above (always)
            _write_plane_chi2(trkhitdf, P, tag + "_", do_calo_syst, cosmic,
                              do_cafana_chi2, rr_min=rr_min)

        # induction planes + best-plane selection (opt-in)
        if do_ind_chi2:
            # front (p0) and middle (p1) induction planes: fresh hitdfs.
            # make_trkhitdf sets hd["plane"]=plane so dedx() uses the right gains.
            ncalo_by_plane = {2: P["ncalo"]}
            for tag, plane in [("p0_", 0), ("p1_", 1)]:
                hd = make_trkhitdf(f, plane)
                ncalo_p = hd.groupby(level=[0, 1, 2]).size()
                ncol = "ncalo_%s" % tag.rstrip("_")
                P[ncol] = _reindex(ncalo_p, P.index, 0).astype(int)
                ncalo_by_plane[plane] = P[ncol]
                _add_dedx_variations(hd, DETECTOR, ismc, do_calo_syst)
                _write_plane_chi2(hd, P, tag, do_calo_syst, cosmic, do_cafana_chi2)

            # bestplane: per-pfp select the plane with the most calo hits, preferring
            # the collection plane on ties (stack order [p2, p0, p1] -> argmax).
            ncalo_stack = np.vstack([ncalo_by_plane[2].values,
                                     ncalo_by_plane[0].values,
                                     ncalo_by_plane[1].values])
            best = np.argmax(ncalo_stack, axis=0)  # 0 -> plane2, 1 -> plane0, 2 -> plane1
            conds = [best == 0, best == 1, best == 2]
            bestcols = {}
            for flavor in flavors:
                for uorp in ["chi2u", "chi2p"]:
                    bestcols["%s_bestplane_%s" % (uorp, flavor)] = np.select(
                        conds,
                        [P["%s_%s" % (uorp, flavor)],       # collection plane (untagged)
                         P["%s_p0_%s" % (uorp, flavor)],
                         P["%s_p1_%s" % (uorp, flavor)]],
                        default=np.nan)
            # concat all bestplane columns at once (avoids DataFrame fragmentation)
            P = pd.concat([P, pd.DataFrame(bestcols, index=P.index)], axis=1)

    return P

def fetch_candidates(S, P, do_calo_syst, use_chi2=True, do_alt_chi2=False,
                     do_cafana_chi2=False, do_ind_chi2=False):
    # ------------------------------------------------------------------
    # fixed candidates; chi2 cuts are applied post-hoc (maple_sel)
    # ------------------------------------------------------------------
    # evt-df chi2 suffix -> per-pfp chi2 flavor ("" = nominal cafpyana)
    chi2_suffixes = {"": "cafpyana"}
    if do_cafana_chi2:
        chi2_suffixes["cafana"] = "cafana"
    if do_calo_syst:
        chi2_suffixes.update({v: v for v in SCALE_SMEAR_VARIATIONS + CALO_VARIATIONS})
    if do_alt_chi2:
        # alternative-plane flavors: per-pfp col name == evt-df suffix
        for tag in _alt_plane_tags(do_ind_chi2):
            for flavor in _chi2_flavors(do_calo_syst, do_cafana_chi2):
                key = "%s_%s" % (tag, flavor)
                chi2_suffixes[key] = key

    mu_ilocs, is_prot_cand, counts = _find_candidates(P, use_chi2=use_chi2)
    counts = counts.reindex(S.index).fillna(0).astype(int)
    has_mu = pd.Series(False, index=S.index)
    if len(mu_ilocs):
        has_mu.loc[mu_ilocs.index] = True

    for c in counts.columns:
        S[c] = counts[c]
    S["has_muon"] = has_mu
    # muon + at least one other candidate pfp (same semantics as the old
    # n_proton > 0 once combined with has_muon)
    S["cut_np"] = counts.n_pfp >= 2
    S["cut_0shwother"] = (counts.n_shower == 0) & (counts.n_other == 0)

    # Default proton PID is angle-dependent: the untrimmed collection plane for
    # candidates at a large angle to the drift (theta_x > 47 deg) and p2trim2 for
    # candidates nearly along the drift.  Build the per-pfp blend for every
    # collection-plane flavor (nominal + calorimetric variations), then map each
    # evt-df chi2 suffix to its (chi2u, chi2p) per-pfp columns: default-plane
    # flavors use the blend, alternate-plane tags stay as computed, and (under
    # do_alt_chi2) a "p2" tag exposes the pure untrimmed collection plane.
    default_flavors = _chi2_flavors(do_calo_syst, do_cafana_chi2)
    costh = np.cos(np.radians(PROT_TRIM_ANGLE_X_DEG))
    use_trim = P.dir_x.abs() >= costh  # theta_x <= 47 deg (NaN dir -> False -> collection)
    blend = {}
    for fl in default_flavors:
        for uorp in ("chi2u", "chi2p"):
            base = P["%s_%s" % (uorp, fl)]
            trim = P["%s_%s_%s" % (uorp, DEFAULT_PROT_TRIM_TAG, fl)]
            blend["%s_protblend_%s" % (uorp, fl)] = base.where(~use_trim, trim)
    P = pd.concat([P, pd.DataFrame(blend, index=P.index)], axis=1)

    # evt-df chi2 suffix -> (per-pfp chi2u col, per-pfp chi2p col) for the proton
    prot_chi2_cols = {}
    for suff, fl in chi2_suffixes.items():
        if fl in default_flavors:  # collection-plane flavor -> angle blend
            prot_chi2_cols[suff] = ("chi2u_protblend_%s" % fl, "chi2p_protblend_%s" % fl)
        else:                      # alternate-plane tag -> pure per-plane column
            prot_chi2_cols[suff] = ("chi2u_%s" % fl, "chi2p_%s" % fl)
    if do_alt_chi2:                # p2-only: the pure untrimmed collection plane
        for fl in default_flavors:
            prot_chi2_cols["p2_%s" % fl] = ("chi2u_%s" % fl, "chi2p_%s" % fl)

    # worst-case cut variables over the candidate protons (min for cuts
    # with direction >, max for cuts with direction <): min chi2u / max chi2p
    # over ALL non-muon candidate pfps, so the proton PID cut on the max chi2p
    # only passes when every candidate is proton-like
    u_cols = list(dict.fromkeys(c for c, _ in prot_chi2_cols.values()))
    p_cols = list(dict.fromkeys(c for _, c in prot_chi2_cols.values()))
    aggs = _proton_aggregates(
        P, is_prot_cand,
        mincols=["trackScore", "ke_proton"] + u_cols,
        maxcols=["dist_start"] + p_cols)
    aggs = aggs.reindex(S.index)
    for suff, (ucol, pcol) in prot_chi2_cols.items():
        S["mu_chi2%s_of_prot_cand" % suff] = aggs["min_%s" % ucol]
        S["prot_chi2%s_of_prot_cand" % suff] = aggs["max_%s" % pcol]
    S["min_proton_track_score"] = aggs.min_trackScore
    S["min_proton_ke"] = aggs.min_ke_proton
    S["max_proton_dist_start"] = aggs.max_dist_start

    # ------------------------------------------------------------------
    # candidate variables (muon + leading candidate proton)
    # ------------------------------------------------------------------
    truthcols = ["true_pdg", "true_genp_x", "true_genp_y", "true_genp_z",
                 "true_end_x", "true_end_y", "true_end_z"]
    mucols = ["len", "start_x", "start_y", "start_z", "end_x", "end_y", "end_z", "dir_x", "dir_y", "dir_z",
              "p_muon", "trackScore", "dist_start", "prim_pfp", "contained10"] + truthcols + \
        ["chi2u_%s" % fl for fl in chi2_suffixes.values()] + \
        ["chi2p_%s" % fl for fl in chi2_suffixes.values()]

    if len(mu_ilocs):
        mu = P.loc[pd.Index(mu_ilocs.values), mucols].copy()
        mu.index = mu_ilocs.index
    else:
        mu = pd.DataFrame(columns=mucols, dtype=float)
    mu = mu.reindex(S.index)

    # leading proton: the longest-length non-muon candidate pfp (with len > 0);
    # all p_* / true_pcand_* output columns come from this same pfp
    prodf = P[is_prot_cand & (P.len > 0)]
    pcols = ["len", "start_x", "start_y", "start_z", "end_x", "end_y", "end_z",
             "dir_x", "dir_y", "dir_z", "dist_start",
             "p_proton", "ke_proton", "trackScore", "chi2u_cafpyana", "chi2p_cafpyana",
             "chi2u_%s_cafpyana" % DEFAULT_PROT_TRIM_TAG,
             "chi2p_%s_cafpyana" % DEFAULT_PROT_TRIM_TAG] + truthcols
    if len(prodf):
        p_ilocs = prodf.len.groupby(level=[0, 1]).idxmax()
        pro = P.loc[pd.Index(p_ilocs.values), pcols].copy()
        pro.index = p_ilocs.index
    else:
        pro = pd.DataFrame(columns=pcols, dtype=float)
    pro = pro.reindex(S.index)

    # muon 4-momentum (MeV) and proton KE sum for recoE
    p_mu_x = mu.p_muon * mu.dir_x
    p_mu_y = mu.p_muon * mu.dir_y
    p_mu_z = mu.p_muon * mu.dir_z
    p_mu_mag = np.sqrt(p_mu_x**2 + p_mu_y**2 + p_mu_z**2)
    E_mu = np.sqrt(p_mu_mag**2 + MUON_MASS**2)

    proton_ke_sum = (P.ke_proton[is_prot_cand] + BE).groupby(level=[0, 1]).sum()
    proton_ke_sum = _reindex(proton_ke_sum, S.index, np.nan)
    found_proton = _reindex(is_prot_cand.groupby(level=[0, 1]).any(), S.index, False).astype(bool)

    recoE = np.where(has_mu & found_proton, E_mu + proton_ke_sum, -999.0)

    # ------------------------------------------------------------------
    # psum: the summed proton-candidate system -- vector-sum momentum,
    # summed KE and total energy, and unit direction over ALL candidate
    # protons (the chi2-free non-muon candidate pfps, before PID cuts).
    # NaN semantics follow _proton_aggregates: a NaN momentum/direction on
    # any candidate (or no candidates at all) makes the psum NaN.
    # ------------------------------------------------------------------
    pc = P[is_prot_cand]
    pcg = pd.DataFrame({
        "px": pc.p_proton * pc.dir_x,
        "py": pc.p_proton * pc.dir_y,
        "pz": pc.p_proton * pc.dir_z,
        "ke": pc.ke_proton,
        "E": pc.ke_proton + PROTON_MASS,
    }).groupby(level=[0, 1])
    psum = pcg.sum()
    cnt = pcg.count()
    sz = pcg.size()
    for c in psum.columns:
        psum.loc[cnt[c] < sz, c] = np.nan
    psum = psum.reindex(S.index)
    psum_p = np.sqrt(psum.px**2 + psum.py**2 + psum.pz**2)
    S["psum_p"] = psum_p
    S["psum_ke"] = psum.ke
    S["psum_E"] = psum.E
    S["psum_dir_x"] = psum.px / psum_p
    S["psum_dir_y"] = psum.py / psum_p
    S["psum_dir_z"] = psum.pz / psum_p

    # transverse / angular variables: the opening angle uses the LEADING
    # proton candidate; the TKI (del_*) uses the summed proton system, so
    # the variables generalize to the Np case
    p_mu = mu.p_muon
    dir_mu = pd.DataFrame({"x":mu.dir_x, "y":mu.dir_y, "z":mu.dir_z})

    p_p = pro.p_proton
    dir_p = pd.DataFrame({"x":pro.dir_x, "y":pro.dir_y, "z":pro.dir_z})

    # opening angle between the muon and the leading-length proton candidate:
    # cosine = dot of the direction vectors normalized by their magnitudes
    # (dirs are unit vectors, but normalize explicitly for safety)
    mu_dirn = np.sqrt(np.einsum("ij,ij->i", dir_mu, dir_mu))
    p_dirn = np.sqrt(np.einsum("ij,ij->i", dir_p, dir_p))
    valid = np.isfinite(mu_dirn) & np.isfinite(p_dirn) & (mu_dirn > 0) & (p_dirn > 0)
    cosang = np.full(len(S), np.nan, dtype=float)
    dot = np.einsum("ij,ij->i", dir_mu, dir_p)
    cosang[valid] = dot[valid] / (mu_dirn[valid] * p_dirn[valid])
    cosang = np.clip(cosang, -1.0, 1.0)
    ang = np.degrees(np.arccos(cosang))
    S["mu_p_opening_angle_deg"] = ang

    dir_psum = pd.DataFrame({"x": S.psum_dir_x, "y": S.psum_dir_y, "z": S.psum_dir_z})
    # n_proton generalizes the TKI residual mass to the Np case (the summed
    # p_E carries each proton's rest mass); sz is the per-slice candidate
    # count from the psum groupby above.
    tki_psum = transverse_kinematics(p_mu, dir_mu, S.psum_p, dir_psum, p_E=S.psum_E,
                                     n_proton=sz.reindex(S.index))

    del_p = tki_psum['del_p']
    del_Tp = tki_psum['del_Tp']
    del_phi = tki_psum['del_phi']
    del_alpha = tki_psum['del_alpha']
    mu_E = tki_psum['mu_E']
    # p_E stays the LEADING proton candidate's on-shell energy (the summed
    # system energy is stored as psum_E)
    p_E = np.sqrt(pro.p_proton**2 + PROTON_MASS**2)

    S["nu_E_calo"] = recoE
    S["mu_len"] = mu.len
    S["mu_dist_start"] = mu.dist_start
    S["mu_prim_pfp"] = mu.prim_pfp
    S["mu_contained10"] = mu.contained10
    S["mu_end_x"] = mu.end_x
    S["mu_end_y"] = mu.end_y
    S["mu_end_z"] = mu.end_z
    S["mu_start_x"] = mu.start_x
    S["mu_start_y"] = mu.start_y
    S["mu_start_z"] = mu.start_z
    S["mu_dir_x"] = mu.dir_x
    S["mu_dir_y"] = mu.dir_y
    S["mu_dir_z"] = mu.dir_z
    S["mu_trackScore"] = mu.trackScore
    for suff, fl in chi2_suffixes.items():
        S["mu_chi2%s_of_mu_cand" % suff] = mu["chi2u_%s" % fl]
        S["prot_chi2%s_of_mu_cand" % suff] = mu["chi2p_%s" % fl]

    # muon-candidate truth (truth of the reco track matched to the muon cand);
    # GUMP kept both the mu_true_* alias and the true_mucand_* canonical name
    mu_true_p = magdf(pd.DataFrame({"x": mu.true_genp_x, "y": mu.true_genp_y, "z": mu.true_genp_z}))
    S["mu_true_p"] = mu_true_p
    S["true_mucand_p"] = mu_true_p
    S["mu_true_pdg"] = mu.true_pdg
    S["true_mucand_pdg"] = mu.true_pdg
    S["true_mucand_dir_x"] = mu.true_genp_x / mu_true_p
    S["true_mucand_dir_y"] = mu.true_genp_y / mu_true_p
    S["true_mucand_dir_z"] = mu.true_genp_z / mu_true_p
    S["true_mucand_end_x"] = mu.true_end_x
    S["true_mucand_end_y"] = mu.true_end_y
    S["true_mucand_end_z"] = mu.true_end_z

    S["p_len"] = pro.len
    S["p_ke"] = pro.ke_proton  # GeV, range-based
    S["p_T"] = pro.ke_proton   # legacy GUMP alias (p_E - PROTON_MASS on-shell)
    S["p_start_x"] = pro.start_x
    S["p_start_y"] = pro.start_y
    S["p_start_z"] = pro.start_z
    S["p_end_x"] = pro.end_x
    S["p_end_y"] = pro.end_y
    S["p_end_z"] = pro.end_z
    S["p_dir_x"] = pro.dir_x
    S["p_dir_y"] = pro.dir_y
    S["p_dir_z"] = pro.dir_z
    S["p_dist_to_vertex"] = pro.dist_start  # legacy GUMP name
    S["p_trackScore"] = pro.trackScore
    S["p_track_score"] = pro.trackScore  # legacy GUMP alias
    # leading-proton default chi2 = the same angle blend as *_of_prot_cand
    use_trim_lead = pro.dir_x.abs() >= np.cos(np.radians(PROT_TRIM_ANGLE_X_DEG))
    S["mu_chi2_of_lead_prot"] = pro.chi2u_cafpyana.where(
        ~use_trim_lead, pro["chi2u_%s_cafpyana" % DEFAULT_PROT_TRIM_TAG])
    S["prot_chi2_of_lead_prot"] = pro.chi2p_cafpyana.where(
        ~use_trim_lead, pro["chi2p_%s_cafpyana" % DEFAULT_PROT_TRIM_TAG])

    # leading-proton-candidate truth (truth of the reco track matched to the
    # leading proton cand); GUMP kept both p_true_* and true_pcand_* names
    p_true_p = magdf(pd.DataFrame({"x": pro.true_genp_x, "y": pro.true_genp_y, "z": pro.true_genp_z}))
    S["p_true_p"] = p_true_p
    S["true_pcand_p"] = p_true_p
    S["p_true_pdg"] = pro.true_pdg
    S["true_pcand_pdg"] = pro.true_pdg
    S["true_pcand_dir_x"] = pro.true_genp_x / p_true_p
    S["true_pcand_dir_y"] = pro.true_genp_y / p_true_p
    S["true_pcand_dir_z"] = pro.true_genp_z / p_true_p
    S["true_pcand_end_x"] = pro.true_end_x
    S["true_pcand_end_y"] = pro.true_end_y
    S["true_pcand_end_z"] = pro.true_end_z

    # sub-leading proton length: 2nd-longest proton candidate (same is_prot_cand
    # set as the leading proton -> no track-score gate). Sort ascending by len,
    # then nth(-2) -> 2nd largest (legacy GUMP idiom); slices with <2 candidates
    # drop out and become NaN after reindex. nth keeps the per-pfp index, so
    # drop the pfp level before aligning onto the slice frame.
    subldf = P[is_prot_cand & (P.len > 0)]
    if len(subldf):
        subl = subldf.sort_values("len").len.groupby(level=[0, 1]).nth(-2).droplevel(2)
    else:
        subl = pd.Series(dtype=float)
    S["subl_proton_length"] = subl.reindex(S.index)

    # longest other pfp: length of the longest pfp in the slice that is
    # neither the muon candidate nor the leading-length proton candidate
    # (ANY pfp counts here, candidate or not); NaN if no such pfp
    excl = pd.Series(False, index=P.index)
    if len(mu_ilocs):
        excl.loc[pd.Index(mu_ilocs.values)] = True
    if len(prodf):
        excl.loc[pd.Index(p_ilocs.values)] = True
    othr = P.len[~excl & P.len.notna()].groupby(level=[0, 1]).max()
    S["othr_pfp_length"] = othr.reindex(S.index)

    # longest "far" primary shower: max length over primary pfps classified as
    # showers (trackScore < 0.5) starting 10-50 cm from the slice vertex. Feeds
    # the MAPLE-side far-shower veto (gumple_cuts). NaN if no such pfp.
    far_shw = P.prim_pfp & (P.trackScore < 0.5) & P.dist_start.between(10.0, 50.0) & P.len.notna()
    far = P.len[far_shw].groupby(level=[0, 1]).max()
    S["max_far_shw_len"] = far.reindex(S.index)

    S["del_p"] = del_p
    S["del_Tp"] = del_Tp
    S["del_phi"] = del_phi
    S["del_alpha"] = del_alpha
    S["mu_E"] = mu_E
    S["p_E"] = p_E

    return S

# =====================================================================
# Main evt builder
# =====================================================================
def make_maple_evt_df(f, selection="none", do_calo_syst=True, use_chi2=True, do_alt_chi2=False,
                      do_cafana_chi2=False, do_ind_chi2=False):
    # use_chi2: muon-candidate selection behavior -- True (default) picks the
    # legacy-GUMP min-chi2u/chi2p-ratio track (len > 40 cm), False the
    # longest chi2-free base-mask track (see _find_candidates).
    # Get a slice level df (S) a pfp level df (P) and some meta-data
    S, P, DETECTOR, RUN, ismc = fetch_info(f)

    # After df prod, limited trk info access. So cut here.
    # Check that all tracks are contained
    valid_c = P.start_x.notna() & P.end_x.notna() & P.len.notna()
    bad_contain = valid_c & ((P.end_x * P.slc_vtx_x < 0) | ~P.contained10)
    any_bad = bad_contain.groupby(level=[0, 1]).any()
    S["cut_contained"] = ~_reindex(any_bad, S.index, False).astype(bool)

    # SBND cathode-crossing veto 
    if DETECTOR == "SBND":
        cross = gmpl.sbnd_cathode_crossing(
            P.slc_vtx_x[valid_c], P.slc_vtx_y[valid_c], P.slc_vtx_z[valid_c],
            P.end_x[valid_c], P.end_y[valid_c], P.end_z[valid_c])
        any_cross = pd.Series(cross, index=P.index[valid_c]).groupby(level=[0, 1]).any()
        S["cut_cathode"] = ~_reindex(any_cross, S.index, False).astype(bool)
    else:
        S["cut_cathode"] = True

    # Again, since we have limited trk info access after df prod
    # we want to grab information about candidate mu and p tracks here.
    # Start with PID calculations and then do some basic candidate ID.
    P = PID_calcs(f, P, DETECTOR, ismc, do_calo_syst=do_calo_syst, do_alt_chi2=do_alt_chi2,
                  do_cafana_chi2=do_cafana_chi2, do_ind_chi2=do_ind_chi2)
    S = fetch_candidates(S, P, do_calo_syst, use_chi2=use_chi2, do_alt_chi2=do_alt_chi2,
                         do_cafana_chi2=do_cafana_chi2, do_ind_chi2=do_ind_chi2)

    # ------------------------------------------------------------------
    # Reco_class (classification_type_debug port, via mcnu-level classification)
    # ------------------------------------------------------------------
    clsdf = maple_truth_classdf(f, det=DETECTOR, run=RUN)
    cls_lookup = clsdf.maple_class if len(clsdf) else pd.Series(dtype=float)
    key = pd.MultiIndex.from_arrays([S.index.get_level_values(0), S.tmatch_idx.fillna(-1).astype(int)])
    nu_class = pd.Series(cls_lookup.reindex(key).values, index=S.index)

    S["Reco_class"] = np.select(
        [S.tmatch_idx < 0,
         nu_class == CLS_1MU1P,
         nu_class == CLS_1MUNP,
         ~gmpl.prefix_fv_cut(S, "true_vtx"),
         np.abs(S.true_pdg) == 12,
         S.true_iscc == 0,
         (S.true_iscc == 1) & (S.true_genie_mode == 0),
         (S.true_iscc == 1) & (S.true_genie_mode == 1),
         (S.true_iscc == 1) & (S.true_genie_mode == 2),
         (S.true_iscc == 1) & ((S.true_genie_mode == 3) | (S.true_genie_mode == 4)),
         (S.true_iscc == 1) & (S.true_genie_mode == 10)],
        [3, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11], default=12)

    # ------------------------------------------------------------------
    # assemble sBruce columns
    S["flash_maxpe"] = S.flash_maxpe
    S["flash_maxpe_cryo0"] = S.flash_maxpe_cryo0
    S["flash_maxpe_cryo1"] = S.flash_maxpe_cryo1

    chain = gmpl.maple_cut_chain(S)
    S["cut_presel"] = chain.cut_presel
    S["cut_cosmic"] = chain.cut_cosmic
    S["cut_flash"] = chain.cut_flash
    S["cut_trk"] = chain.cut_trk
    S["cut_muon"] = chain.cut_muon
    S["cut_protons"] = chain.cut_protons
    S["cut_far_shw"] = chain.cut_far_shw
    # exclusive selections split on n_pfp: gump = base chain & n_pfp==2 (1u1p),
    # maple = base chain & n_pfp>2 (1uN>1p). NB maple_sel used to be the
    # inclusive base chain; 1u1p analyses must use gump_sel.
    S["gump_sel"] = chain.gump_sel
    S["maple_sel"] = chain.maple_sel

    # ------------------------------------------------------------------
    # selection option
    # ------------------------------------------------------------------
    if selection == "none":
        pass
    elif selection == "presel":
        S = S[S.cut_presel]
    elif selection == "full":
        S = S[S.gump_sel | S.maple_sel]
    else:
        raise ValueError("selection must be 'none', 'presel', or 'full'")

    return S


# thin wrappers for configs -----------------------------------------------
def make_maple_evt_nosel_df(f):
    return make_maple_evt_df(f, selection="none", do_calo_syst=True)

def make_maple_evt_presel_df(f):
    return make_maple_evt_df(f, selection="presel", do_calo_syst=True)

# non-CV / detector-variation MC: preselection only, calo variations skipped
# (those are only needed for the CV reweighting envelope, and computing them
# substantially slows down processing).
def make_maple_evt_presel_nocalo_df(f):
    return make_maple_evt_df(f, selection="presel", do_calo_syst=False)

# data: preselection only, no truth (mcnu built separately/omitted) and no calo
# variations.
def make_maple_evt_presel_data_df(f):
    return make_maple_evt_df(f, selection="presel", do_calo_syst=False)

def make_maple_evt_fullsel_df(f):
    return make_maple_evt_df(f, selection="full", do_calo_syst=True)

def make_maple_evt_fullsel_data_df(f):
    return make_maple_evt_df(f, selection="full", do_calo_syst=False)


# =====================================================================
# mcnu builder
# =====================================================================
def make_maple_nudf(f):
    mc = _flatcols(loadbranches(f["recTree"], mcbranches).rec.mc.nu)
    if mc.empty:
        return pd.DataFrame()

    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det
    DETECTOR = "SBND" if 1 == det.unique() else "ICARUS"
    run = loadbranches(f["recTree"], ["rec.hdr.run"]).rec.hdr.run
    RUN = 1 if DETECTOR == "SBND" else (2 if run.iloc[0] < 12960 else 4)

    cls = maple_truth_classdf(f, det=DETECTOR, run=RUN)

    nudf = pd.DataFrame({
        "nu_E": mc.E,
        "detector": DETECTOR,
        "Run": RUN,
        "pdg": mc.pdg,
        "is_cc": mc.iscc,
        "is_nc": mc.isnc,
        "genie_mode": mc.genie_mode,
        "pos_x": mc.position_x,
        "pos_y": mc.position_y,
        "pos_z": mc.position_z,
        "baseline": mc.baseline,
        "time": mc.time,
        # link to the GENIE event record (evtrec table entry index)
        "genie_evtrec_idx": mc.genie_evtrec_idx,
        "maple_class": cls.maple_class,
        "is_1mu1p_maple": cls.maple_class == CLS_1MU1P,
        "is_1muNp_maple": cls.maple_class == CLS_1MUNP,
        "nmu": cls.nmu,
        "np": cls.np,
        "npi": cls.npi,
        "npi0": cls.npi0,
        "nn": cls.nn,
        "mu_length": cls.mu_length,
        "veto_particles": cls.veto,
        "uncontained_truth": cls.uncontained,
        "true_visible_Enu": cls.true_visible_Enu,
    })
    nudf["is_sig"] = (cls.maple_class == CLS_1MU1P) | (cls.maple_class == CLS_1MUNP)
    nudf["is_other_numucc"] = cls.maple_class == CLS_1MUNP
    nudf["is_fv"] = gmpl.prefix_fv_cut(nudf, "pos")
    nudf["ind"] = nudf.index.get_level_values(1)
    nudf["detector"] = DETECTOR
    nudf["Run"] = RUN

    nudf.columns = pd.MultiIndex.from_tuples([(col, '') for col in nudf.columns])

    return nudf

gump_ar23_weights = [
    # CCQE
    "GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape",
    'GENIEReWeight_SBN_v1_multisigma_CoulombCCQE',

    # MEC
    'GENIEReWeight_SBN_v1_multisigma_NormCCMEC',
    'GENIEReWeight_SBN_v1_multisigma_NormNCMEC',
    "GENIEReWeight_SBN_v1_multisigma_DecayAngMEC",

    # RES
    "GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi",
    "GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad",
    "GENIEReWeight_SBN_v1_multisigma_MaCCRES",
    "GENIEReWeight_SBN_v1_multisigma_MaNCRES",
    "GENIEReWeight_SBN_v1_multisigma_MvCCRES",
    "GENIEReWeight_SBN_v1_multisigma_MvNCRES",
    "GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma",
    "GENIEReWeight_SBN_v1_multisigma_RDecBR1eta",

    # Non-Res
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi',

    # DIS
    # "GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse",
    'GENIEReWeight_SBN_v1_multisigma_AhtBY',
    'GENIEReWeight_SBN_v1_multisigma_BhtBY',
    'GENIEReWeight_SBN_v1_multisigma_CV1uBY',
    'GENIEReWeight_SBN_v1_multisigma_CV2uBY',

    # COH
    "GENIEReWeight_SBN_v1_multisigma_NormCCCOH",
    "GENIEReWeight_SBN_v1_multisigma_NormNCCOH",

    # FSI
    # "GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse",
    # "GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse",
    'GENIEReWeight_SBN_v1_multisigma_MFP_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrCEx_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrInel_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrAbs_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi',

    # NCEL
    'GENIEReWeight_SBN_v1_multisigma_MaNCEL',
    'GENIEReWeight_SBN_v1_multisigma_EtaNCEL',
]

# Systematics introduced by Ar23+
gump_ar23p_weights = [
    "CCQETemplateReweight_SBN_v3_LFGToSF_q0bin0",
    "CCQETemplateReweight_SBN_v3_LFGToSF_q0bin1",
    "CCQETemplateReweight_SBN_v3_LFGToSF_q0bin2",
    "CCQETemplateReweight_SBN_v3_LFGToSF_q0bin3",
    "CCQETemplateReweight_SBN_v3_LFGToSF_q0bin4",

    "CCQETemplateReweight_SBN_v3_LFGToHF_q0bin0",
    "CCQETemplateReweight_SBN_v3_LFGToHF_q0bin1",
    "CCQETemplateReweight_SBN_v3_LFGToHF_q0bin2",
    "CCQETemplateReweight_SBN_v3_LFGToHF_q0bin3",
    "CCQETemplateReweight_SBN_v3_LFGToHF_q0bin4",

    "CCQETemplateReweight_SBN_v3_HFToCRPA_q0bin0",
    "CCQETemplateReweight_SBN_v3_HFToCRPA_q0bin1",
    "CCQETemplateReweight_SBN_v3_HFToCRPA_q0bin2",
    "CCQETemplateReweight_SBN_v3_HFToCRPA_q0bin3",
    "CCQETemplateReweight_SBN_v3_HFToCRPA_q0bin4",

    "QEInterference_SBN_v3_QEIntf_dial_0",
    "QEInterference_SBN_v3_QEIntf_dial_1",
    "QEInterference_SBN_v3_QEIntf_dial_2",
    "QEInterference_SBN_v3_QEIntf_dial_3",
    "QEInterference_SBN_v3_QEIntf_dial_4",
    "QEInterference_SBN_v3_QEIntf_dial_5",

    "GENIEReWeight_SBN_v3_FrG4LoE_N",
    "GENIEReWeight_SBN_v3_FrG4M1E_N",
    "GENIEReWeight_SBN_v3_FrG4M2E_N",
    "GENIEReWeight_SBN_v3_FrG4HiE_N",
    "GENIEReWeight_SBN_v3_FrINCLLoE_N",
    "GENIEReWeight_SBN_v3_FrINCLM1E_N",
    "GENIEReWeight_SBN_v3_FrINCLM2E_N",
    "GENIEReWeight_SBN_v3_FrINCLHiE_N",
    "GENIEReWeight_SBN_v3_MFPLoE_N",
    "GENIEReWeight_SBN_v3_MFPM1E_N",
    "GENIEReWeight_SBN_v3_MFPM2E_N",
    "GENIEReWeight_SBN_v3_MFPHiE_N",

    "ZExpPCAWeighter_SBN_v3_MvA_b1",
    "ZExpPCAWeighter_SBN_v3_MvA_b2",
    "ZExpPCAWeighter_SBN_v3_MvA_b3",
    "ZExpPCAWeighter_SBN_v3_MvA_b4",

    "MECq0q3InterpWeighting_SBN_v3_SuSAToVal_MECResponse_q0bin0",
    "MECq0q3InterpWeighting_SBN_v3_SuSAToVal_MECResponse_q0bin1",
    "MECq0q3InterpWeighting_SBN_v3_SuSAToVal_MECResponse_q0bin2",
    "MECq0q3InterpWeighting_SBN_v3_SuSAToVal_MECResponse_q0bin3",
    "MECq0q3InterpWeighting_SBN_v3_SuSAToMar_MECResponse_q0bin0",
    "MECq0q3InterpWeighting_SBN_v3_SuSAToMar_MECResponse_q0bin1",
    "MECq0q3InterpWeighting_SBN_v3_SuSAToMar_MECResponse_q0bin2",
    "MECq0q3InterpWeighting_SBN_v3_SuSAToMar_MECResponse_q0bin3",

    "CCQEXSecCorr_SBN_v3_CCQEXSecCorr",
    "GENIEReWeight_SBN_v3_FrKin_PiProFix_N",
]

# Other systematics we keep for extra info
extra_weights = [
    'GENIEReWeight_SBN_v1_multisigma_RPA_CCQE',
    'GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE',
    'GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE',
    'GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE',
    'GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE',
    'GENIEReWeight_SBN_v1_multisigma_MFP_N',
    'GENIEReWeight_SBN_v1_multisigma_FrCEx_N',
    'GENIEReWeight_SBN_v1_multisigma_FrInel_N',
    'GENIEReWeight_SBN_v1_multisigma_FrAbs_N',
    'GENIEReWeight_SBN_v1_multisigma_FrPiProd_N',

    "PionAbsWeighter_SBN_v3_QuasiDeuteronFraction",
    "GENIEReWeight_SBN_v3_FrG4_N",
    "GENIEReWeight_SBN_v3_FrINCL_N",

    "ZExpPCAWeighter_SBN_v3_Deut_b1",
    "ZExpPCAWeighter_SBN_v3_Deut_b2",
    "ZExpPCAWeighter_SBN_v3_Deut_b3",
    "ZExpPCAWeighter_SBN_v3_Deut_b4",
    "GENIEReWeight_SBN_v3_FrKin_PiProBias_N",
]

g4_weights = [
    "reinteractions_piminus_Geant4",
    "reinteractions_piplus_Geant4",
    "reinteractions_proton_Geant4"
]

gump_genie_reknob_systematics = gump_ar23_weights + gump_ar23p_weights + extra_weights + g4_weights

# =====================================================================
# wgt builder
# =====================================================================
def make_maple_wgtdf(f):
    """Systematic-weight dataframe (mcnu + multisim/multisigma weights).

    Requests the standard GENIE reweight set, restricted to the psets
    actually present in the input file (samples like ReCAF2026 do not carry
    every pset in the default list).
    """
    from makedf import geniesyst
    if "globalTree" in f:
        avail = list(f["globalTree"]["global/wgts/wgts.name"].arrays(library="np")["wgts.name"][0])
    else:
        avail = []
    systs = [s for s in geniesyst.regen_systematics if s in avail]

    systs.extend(g4_weights)

    missing = len(geniesyst.regen_systematics) + len(g4_weights) - len(systs)
    if missing:
        print("make_maple_wgtdf: %d requested GENIE systematics absent in file, using %d" % (missing, len(systs)))
    return make_mcnudf(f, include_weights=True, multisim_nuniv=100, genie_systematics=systs)


def make_maple_rewgtdf(f):

    """Systematic-weight dataframe with the GUMP CV reweight knob set.

    Same weight request as gump make_gump_nurewgtdf, so the wgt table is
    column-compatible with the GUMP sbn-rewgted CV productions.
    """

    syst = gump_genie_reknob_systematics + g4_weights
    ret = make_mcnudf(f, include_weights=True, slim=False,
                       genie_systematics=list(set(syst)))

    return ret
# =====================================================================
# GENIE event record (evtrec) builder
# =====================================================================
def _read_genie_evtrec_subprocess(path, timeout=900):
    """Run genie_evtrec.read_genie_evtrec in a fresh python interpreter.

    pyROOT deadlocks when first used inside a forked multiprocessing Pool
    worker (as spawned by NTupleGlob.dataframes): the worker either hangs at
    recycling (maxtasksperchild=1) or dies without delivering its result,
    stalling the pool forever.  A fresh exec'd interpreter has none of the
    inherited fork state and exits cleanly, so the raw-object read is done
    there and the numpy arrays are shipped back via pickle.

    Raises on subprocess failure or if the read exceeds `timeout` seconds
    (a stuck read must fail loudly rather than hang the production).
    """
    import os
    import pickle
    import subprocess
    import sys
    import tempfile

    with tempfile.NamedTemporaryFile(suffix=".pkl") as tf:
        code = (
            "import pickle\n"
            "from makedf import genie_evtrec\n"
            "d = genie_evtrec.read_genie_evtrec(%r)\n"
            "with open(%r, 'wb') as f:\n"
            "    pickle.dump(d, f)\n" % (str(path), tf.name)
        )
        # repo root on sys.path so `makedf` is importable regardless of cwd
        repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        env = dict(os.environ)
        env["PYTHONPATH"] = repo + os.pathsep + env.get("PYTHONPATH", "")
        res = subprocess.run(
            [sys.executable, "-c", code],
            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
            timeout=timeout, env=env)
        if res.returncode != 0:
            raise IOError(
                "GENIE event-record subprocess failed (exit %d) for %s:\n%s"
                % (res.returncode, path, res.stderr.decode(errors="replace")[-2000:]))
        with open(tf.name, "rb") as f:
            return pickle.load(f)

def make_maple_evtrec_df(f):
    """GENIE event record (evtrec) table; same schema as make_genie_evtrec_df.

    The flat-StdHep path is pure uproot and pool-safe, so it goes through the
    core builder unchanged.  The raw genie::NtpMCEventRecord path (used e.g.
    by the ICARUS ReCAF2026 files) needs pyROOT + the GENIE libraries, which
    deadlock inside forked Pool workers -- that read is isolated in a fresh
    interpreter via _read_genie_evtrec_subprocess.
    """
    if "GenieEvtRecTree" not in f:
        return pd.DataFrame([])
    if "GenieEvtRec.StdHepPdg" in f["GenieEvtRecTree"].keys():
        return make_genie_evtrec_df(f)
    path = getattr(f, "file_path", None)
    if path is None:
        path = f._file.file_path
    d = _read_genie_evtrec_subprocess(path)
    ret = _build_genie_evtrec_df(d)
    return ret
