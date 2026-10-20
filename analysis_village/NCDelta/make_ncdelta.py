"""Dataframe makers for the Pandora NC Delta -> N gamma (1g1p / 1g0p) selection.

One row per slice with the leading primary shower and track as flat columns.
By default the topology and pre-selection from cuts.py are applied here to
keep the dataframes small; the *_nopresel makers skip them.
"""

import numpy as np
import pandas as pd

from makedf.makedf import (
    loadbranches, make_slcdf, make_pfpdf, make_trkhitdf, make_mcdf,
    make_hdrdf, make_potdf_bnb, barycenterFMbranches,
)
from makedf.branches import trueparticlebranches
from makedf.util import mag

from analysis_village.NCDelta import cuts

PROTON_MASS = 938.272  # MeV
SLC_LEVELS = [0, 1]    # (entry, rec.slc..index)


def fetch_metadata(f):
    """(DETECTOR, RUN, ismc). Same as analysis_village/mcs_pid."""
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det

    if det.empty:
        raise ValueError("fetch_metadata: rec.hdr.det is empty; cannot identify detector.")

    dets = pd.unique(det)
    if len(dets) != 1:
        raise ValueError("fetch_metadata: file mixes detectors (rec.hdr.det = %s)." % dets)

    if dets[0] == 1:
        DETECTOR = "SBND"
    elif dets[0] == 2:
        DETECTOR = "ICARUS"
    else:
        raise ValueError("df maker needs rec.hdr.det == 1 (SBND) or 2 (ICARUS); got %s" % dets[0])
    run = loadbranches(f["recTree"], ["rec.hdr.run"]).rec.hdr.run
    RUN = 1 if DETECTOR == "SBND" else (2 if run.iloc[0] < 12960 else 4)
    ismc = bool(loadbranches(f["recTree"], ["rec.hdr.ismc"]).rec.hdr.ismc.iloc[0])
    return DETECTOR, RUN, ismc


# ---------------------------------------------------------------------------
# Track calorimetry from collection-plane hits
# ---------------------------------------------------------------------------

def _good_hits(hitdf, dedxname="dedx"):
    return (hitdf[dedxname] > 0) & (hitdf[dedxname] < 1000.)


def truncated_mean_dedx(hitdf, nsigma=1.0, dedxname="dedx"):
    """Mean dE/dx of hits within nsigma RMS of the track median.

    Single-pass version of LArSoft's TruncMean (MicroBooNE iterates).
    """
    d = hitdf.loc[_good_hits(hitdf, dedxname), dedxname]
    lv = list(range(d.index.nlevels - 1))
    g = d.groupby(level=lv)
    med = g.transform("median")
    rms = g.transform("std").fillna(0.)  # one-hit tracks keep their hit
    keep = (d - med).abs() <= nsigma * rms
    return d[keep].groupby(level=lv).mean()


def calo_energy(hitdf, dedxname="dedx"):
    """Per-track sum of dE/dx * pitch [MeV]."""
    good = _good_hits(hitdf, dedxname)
    e = (hitdf[dedxname] * hitdf.pitch)[good]
    return e.groupby(level=list(range(e.index.nlevels - 1))).sum()


# ---------------------------------------------------------------------------
# Leading primary shower / track per slice
# ---------------------------------------------------------------------------

def _leading(df, key, prefix):
    """Row with the largest `key` per slice; NaN keys never win."""
    lead = df.sort_values(key, na_position="first").groupby(level=SLC_LEVELS).tail(1)
    lead = lead.reset_index(level=-1)
    lead = lead.rename(columns={lead.columns[0]: "pfp_idx"})
    lead.columns = [prefix + c for c in lead.columns]
    return lead


def _shower_table(pfp):
    s = pfp.pfp.shw
    e = pd.concat([s.plane.I0.energy, s.plane.I1.energy, s.plane.I2.energy], axis=1) * 1000.
    e_max = e.max(axis=1)
    e2 = e.iloc[:, 2]
    # MicroBooNE "#GE#": collection-plane energy, falling back to the max plane.
    ge = e2.where(e2 > 0, e_max)
    return pd.DataFrame({
        "energy": ge,
        "energy_max": e_max,
        "energy_best": s.bestplane_energy * 1000.,
        "energy_I0": e.iloc[:, 0], "energy_I1": e.iloc[:, 1], "energy_I2": e2,
        "dedx_best": s.bestplane_dEdx,
        "dedx_I2": s.plane.I2.dEdx,
        "nhits_I0": s.plane.I0.nHits, "nhits_I1": s.plane.I1.nHits, "nhits_I2": s.plane.I2.nHits,
        "start_x": s.start.x, "start_y": s.start.y, "start_z": s.start.z,
        "end_x": s.end.x, "end_y": s.end.y, "end_z": s.end.z,
        "dir_x": s.dir.x, "dir_y": s.dir.y, "dir_z": s.dir.z,
        "conv_gap": s.conversion_gap,
        "density": s.density,
        "open_angle": s.open_angle,
        "len": s.len,
        "trackScore": pfp.pfp.trackScore,
        "true_pdg": s.truth.p.pdg,
        "true_startE": s.truth.p.startE * 1000.,
    }, index=pfp.index)


def _track_table(pfp):
    t = pfp.pfp.trk
    return pd.DataFrame({
        "start_x": t.start.x, "start_y": t.start.y, "start_z": t.start.z,
        "end_x": t.end.x, "end_y": t.end.y, "end_z": t.end.z,
        "dir_x": t.dir.x, "dir_y": t.dir.y, "dir_z": t.dir.z,
        "len": t.len,
        "p_proton": t.rangeP.p_proton * 1000.,
        "chi2_proton": t.chi2pid.I2.chi2_proton,
        "chi2_muon": t.chi2pid.I2.chi2_muon,
        "pida": t.chi2pid.I2.pida,
        "trackScore": pfp.pfp.trackScore,
        "true_pdg": t.truth.p.pdg,
        "true_startE": t.truth.p.startE * 1000.,
    }, index=pfp.index)


# ---------------------------------------------------------------------------
# Truth
# ---------------------------------------------------------------------------

def _truth_columns(t):
    """Flat truth columns from a make_mcdf-style frame (slc.truth or mc.nu)."""
    return pd.DataFrame({
        "true_E": t.E * 1000.,
        "true_pdg": t.pdg,
        "true_iscc": t.iscc,
        "true_isnc": t.isnc,
        "true_genie_mode": t.genie_mode,
        "true_w": t.w,
        "true_vtx_x": t.position.x,
        "true_vtx_y": t.position.y,
        "true_vtx_z": t.position.z,
        "true_npi0": t.npi0,
        "true_ngamma": t.ng,
        "true_np": t.np,
        "true_np_20MeV": t.np_20MeV,
        "true_np_50MeV": t.np_50MeV,
        "true_npi_30MeV": t.npi_30MeV,
        "true_nmu_27MeV": t.nmu_27MeV,
    })


def _leading_photon_E(f, prefix):
    """Energy [MeV] of the highest-energy primary photon per interaction."""
    prim = loadbranches(f["recTree"], [prefix + ".prim.pdg", prefix + ".prim.genE"])
    for lvl in prefix.split("."):
        prim = prim[lvl]
    prim = prim.prim
    photons = prim.genE[prim.pdg == 22]
    return (photons.groupby(level=SLC_LEVELS).max() * 1000.).rename("true_leading_gamma_E")


# ---------------------------------------------------------------------------
# pi0 decay-photon truth
# ---------------------------------------------------------------------------
# The dominant background is NC 1pi0 in which only one of the two decay photons
# is reconstructed. Those photons are GEANT4 secondaries, so they appear in
# neither true_ngamma nor true_leading_gamma_E (both count GENIE primaries, and
# the pi0 is itself the primary). The columns below describe them, so the
# failure mode can be identified per event:
#
#   escaped        gamma2_contained == 0, or conv_dist running out of the TPC
#   too soft       gamma2_E below the clustering threshold
#   merged         small opening_angle and small gamma_sep
#   not reco'd     real gamma2_visE inside the TPC but no second shower PFP
#
# Only the last two are addressable by reconstruction work; the first two are
# physics. All energies MeV, lengths cm, angles rad.

PI0_PDG = 111
PHOTON_PDG = 22


def _tp_visE(tp):
    """Collection-plane visible energy per true particle [MeV], or NaN.

    The branch is rec.true_particles.plane.<cryo>.2.visE; the sub-level naming
    differs between uproot/awkward versions, so fall back rather than fail.
    """
    for cryo in ("0", 0):
        try:
            return tp.plane[cryo][2].visE * 1000.
        except Exception:
            pass
    try:
        cols = [c for c in tp.columns if "visE" in [str(x) for x in c]]
        if cols:
            return tp[cols[0]] * 1000.
    except Exception:
        pass
    return None


PI0_TRUTH_COLS = ["true_pi0_ngamma", "true_pi0_gamma1_E", "true_pi0_gamma2_E",
                  "true_pi0_gamma2_visE", "true_pi0_gamma2_contained",
                  "true_pi0_gamma2_conv_dist", "true_pi0_gamma_sep",
                  "true_pi0_opening_angle"]


def _pi0_photon_summary(flat, entry):
    """Per-(entry, interaction) pi0 photon summary from a flat particle table.

    `flat` has one row per true particle with columns pdg, G4ID, parent, genE,
    cont_tpc, visE, s{x,y,z} (production), e{x,y,z} (conversion/stop),
    p{x,y,z} (momentum) and iid (interaction id). Split out from the CAF
    loading so it can be unit-tested without a CAF.
    """
    pi0 = flat[flat.pdg == PI0_PDG][[entry, "G4ID", "iid", "genE"]]
    gam = flat[flat.pdg == PHOTON_PDG]
    if pi0.empty or gam.empty:
        return None

    g = gam.merge(pi0.rename(columns={"G4ID": "parent", "iid": "pi0_iid",
                                      "genE": "pi0_E"}),
                  on=[entry, "parent"], how="inner")
    if g.empty:
        return None

    key = [entry, "pi0_iid"]
    # most energetic pi0 per interaction, then its photons by energy
    g = g.sort_values(key + ["pi0_E", "genE"], ascending=[True, True, False, False])
    g = g[g.parent == g.groupby(key).parent.transform("first")]
    g["rank"] = g.groupby(key + ["parent"]).cumcount()

    g1 = g[g["rank"] == 0].set_index(key)
    g2 = g[g["rank"] == 1].set_index(key)

    out = pd.DataFrame({"true_pi0_ngamma": g.groupby(key).size()})
    out["true_pi0_gamma1_E"] = g1.genE
    out["true_pi0_gamma2_E"] = g2.genE
    out["true_pi0_gamma2_visE"] = g2.visE
    out["true_pi0_gamma2_contained"] = g2.cont_tpc
    out["true_pi0_gamma2_conv_dist"] = mag(g2.ex - g2.sx, g2.ey - g2.sy, g2.ez - g2.sz)
    out["true_pi0_gamma_sep"] = mag(g1.ex - g2.ex, g1.ey - g2.ey, g1.ez - g2.ez)
    d1 = mag(g1.px, g1.py, g1.pz)
    d2 = mag(g2.px, g2.py, g2.pz)
    cosang = (g1.px * g2.px + g1.py * g2.py + g1.pz * g2.pz) / (d1 * d2)
    out["true_pi0_opening_angle"] = np.arccos(np.clip(cosang, -1., 1.))
    return out


def _pi0_photon_truth(f, index):
    """Per-slice summary of the pi0 decay photons in the matched interaction.

    Returns a DataFrame reindexed onto `index`; all-NaN if the CAF does not
    carry the needed collections.
    """
    empty = pd.DataFrame(np.nan, index=index, columns=PI0_TRUTH_COLS)
    try:
        tp = loadbranches(f["recTree"], trueparticlebranches).rec.true_particles
    except Exception as e:
        print("Warning: could not load true_particles (%s); pi0 photon truth "
              "columns will be NaN." % type(e).__name__)
        return empty

    visE = _tp_visE(tp)
    flat = pd.DataFrame({
        "pdg": tp.pdg, "G4ID": tp.G4ID, "parent": tp.parent,
        "genE": tp.genE * 1000.,
        "cont_tpc": tp.cont_tpc,
        "sx": tp.start.x, "sy": tp.start.y, "sz": tp.start.z,
        "ex": tp.end.x, "ey": tp.end.y, "ez": tp.end.z,
        "px": tp.genp.x, "py": tp.genp.y, "pz": tp.genp.z,
        "iid": tp.interaction_id,
    })
    flat["visE"] = np.nan if visE is None else visE
    flat = flat.reset_index()
    entry = flat.columns[0]

    out = _pi0_photon_summary(flat, entry)
    if out is None:
        return empty

    try:
        sid = loadbranches(f["recTree"], ["rec.slc.truth.prim.interaction_id"])
        sid = sid.rec.slc.truth.prim.interaction_id
        sid = sid.groupby(level=SLC_LEVELS).first().rename("iid")
    except Exception as e:
        print("Warning: no rec.slc.truth.prim.interaction_id (%s); pi0 photon "
              "truth cannot be matched to slices." % type(e).__name__)
        return empty

    sid = sid.reset_index()
    idx_names = list(sid.columns[:2])
    joined = sid.merge(out.reset_index().rename(columns={"pi0_iid": "iid"}),
                       left_on=[sid.columns[0], "iid"], right_on=[entry, "iid"],
                       how="left").set_index(idx_names)
    return joined.reindex(index)[PI0_TRUTH_COLS]


# ---------------------------------------------------------------------------
# Slice dataframe
# ---------------------------------------------------------------------------

def make_ncdelta_slcdf(f, apply_preselection=True, include_truth=None):
    DETECTOR, RUN, ismc = fetch_metadata(f)
    if include_truth is None:
        include_truth = ismc

    # ---- slices ----
    slcdf = make_slcdf(f)
    slc = slcdf.slc
    S = pd.DataFrame({
        "slc_vtx_x": slc.vertex.x,
        "slc_vtx_y": slc.vertex.y,
        "slc_vtx_z": slc.vertex.z,
        "is_clear_cosmic": slc.is_clear_cosmic.astype(bool),
        "nu_score": slc.nu_score,
        "crlongtrkdiry": slc.nuid.crlongtrkdiry,
        "opt0_score": slc.opt0.score,
        "opt0_time": slc.opt0.time,
        "opt0_measPE": slc.opt0.measPE,
        "opt0_hypoPE": slc.opt0.hypoPE,
    }, index=slcdf.index)

    if include_truth:
        S["tmatch_idx"] = slc.tmatch.idx
        S["tmatch_eff"] = slc.tmatch.eff
        S["tmatch_pur"] = slc.tmatch.pur

    if "rec.slc.barycenterFM.score" in f["recTree"].keys():
        bary = loadbranches(f["recTree"], barycenterFMbranches).rec.slc.barycenterFM
        S["bary_score"] = bary.score
        S["bary_deltaZ"] = bary.deltaZ

    # MicroBooNE only considered Pandora's neutrino candidate.
    S = S[~S.is_clear_cosmic].drop(columns="is_clear_cosmic")

    # nu_score rank within the event (1 = highest), before pre-selection.
    S["nu_score_rank"] = S.nu_score.groupby(level=0).rank(ascending=False, method="first")

    S["detector"] = DETECTOR
    S["Run"] = RUN
    S["ismc"] = ismc
    S["vtx_dist_to_wall"] = cuts.dist_to_active_boundary(S.slc_vtx_x, S.slc_vtx_y, S.slc_vtx_z, DETECTOR)

    # ---- primary PFPs in the kept slices ----
    pfp = make_pfpdf(f, update_shw=False)
    pfp = pfp[pfp.index.droplevel(-1).isin(S.index)]
    primary = pfp.pfp.parent_is_primary.fillna(False).astype(bool)
    score = pfp.pfp.trackScore
    is_trk = primary & (score >= cuts.TRK_SCORE_CUT)
    is_shw = primary & (score >= 0) & (score < cuts.TRK_SCORE_CUT)
    is_unc = primary & ~(score >= 0)  # NaN / negative score
    # Shower-like PFPs deeper in the hierarchy, e.g. a pi0's second photon
    # attached to another particle. Used by the second-shower veto in cuts.py.
    is_nonprim_shw = ~primary & (score >= 0) & (score < cuts.TRK_SCORE_CUT)

    for name, m in (("n_trk", is_trk), ("n_shw", is_shw), ("n_unclassified", is_unc),
                    ("n_nonprim_shw", is_nonprim_shw)):
        S[name] = m.groupby(level=SLC_LEVELS).sum().reindex(S.index, fill_value=0).astype(int)

    if apply_preselection:
        S = S[cuts.topology_cut(S)]

    # ---- leading shower (highest #GE#) and track (longest) ----
    shw = _shower_table(pfp[is_shw])
    shw = _leading(shw[shw.index.droplevel(-1).isin(S.index)], "energy", "shw_")
    trk = _track_table(pfp[is_trk])
    trk = _leading(trk[trk.index.droplevel(-1).isin(S.index)], "len", "trk_")
    del pfp

    D = S.join(shw, how="left").join(trk, how="left")

    # ---- track calorimetry, only for the kept leading tracks ----
    D["trk_dedx_trunc"] = np.nan
    D["trk_calo_E"] = np.nan
    has_trk = D.trk_pfp_idx.notna()
    if has_trk.any():
        sub = D.index[has_trk]
        lead_trk_idx = pd.MultiIndex.from_arrays([
            sub.get_level_values(0),
            sub.get_level_values(1),
            D.trk_pfp_idx[has_trk].astype(int).values,
        ])
        hits = make_trkhitdf(f, plane=2)
        hits = hits[hits.index.droplevel(-1).isin(lead_trk_idx)]
        D.loc[has_trk, "trk_dedx_trunc"] = truncated_mean_dedx(hits).reindex(lead_trk_idx).values
        D.loc[has_trk, "trk_calo_E"] = calo_energy(hits).reindex(lead_trk_idx).values
        del hits

    # ---- derived quantities ----
    dx = D.shw_start_x - D.slc_vtx_x
    dy = D.shw_start_y - D.slc_vtx_y
    dz = D.shw_start_z - D.slc_vtx_z
    D["shw_conv_dist"] = mag(dx, dy, dz)
    D["shw_implied_dir_x"] = dx / D.shw_conv_dist
    D["shw_implied_dir_y"] = dy / D.shw_conv_dist
    D["shw_implied_dir_z"] = dz / D.shw_conv_dist

    D["shw_start_dist_to_wall"] = cuts.dist_to_active_boundary(D.shw_start_x, D.shw_start_y, D.shw_start_z, DETECTOR)
    D["shw_end_dist_to_wall"] = cuts.dist_to_active_boundary(D.shw_end_x, D.shw_end_y, D.shw_end_z, DETECTOR)
    D["trk_start_dist_to_wall"] = cuts.dist_to_active_boundary(D.trk_start_x, D.trk_start_y, D.trk_start_z, DETECTOR)
    D["trk_end_dist_to_wall"] = cuts.dist_to_active_boundary(D.trk_end_x, D.trk_end_y, D.trk_end_z, DETECTOR)

    D["trk_start_dist_to_vtx"] = mag(D.trk_start_x - D.slc_vtx_x, D.trk_start_y - D.slc_vtx_y, D.trk_start_z - D.slc_vtx_z)
    D["trk_shw_start_dist"] = mag(D.trk_start_x - D.shw_start_x, D.trk_start_y - D.shw_start_y, D.trk_start_z - D.shw_start_z)

    D["cos_trk_shw_implied"] = (D.trk_dir_x * D.shw_implied_dir_x + D.trk_dir_y * D.shw_implied_dir_y
                                + D.trk_dir_z * D.shw_implied_dir_z)
    D["cos_trk_shw"] = D.trk_dir_x * D.shw_dir_x + D.trk_dir_y * D.shw_dir_y + D.trk_dir_z * D.shw_dir_z

    D["trk_ke_p"] = np.sqrt(D.trk_p_proton**2 + PROTON_MASS**2) - PROTON_MASS
    # Proton + photon invariant mass; peaks near the Delta(1232) for signal.
    e_p = D.trk_ke_p + PROTON_MASS
    D["delta_mass"] = np.sqrt(PROTON_MASS**2 + 2. * D.shw_energy * (e_p - D.trk_p_proton * D.cos_trk_shw_implied))

    if apply_preselection:
        D = D[cuts.preselection(D)]

    # ---- slice truth ----
    if include_truth:
        T = _truth_columns(slcdf.slc.truth).reindex(D.index)
        T["true_leading_gamma_E"] = _leading_photon_E(f, "rec.slc.truth").reindex(D.index)
        D = D.join(T)
        # pi0 decay-photon truth, for studying the missed second shower
        D = D.join(_pi0_photon_truth(f, D.index))
    else:
        D = D.drop(columns=[c for c in D.columns if "true_" in c])

    return D


def make_ncdelta_slcdf_mc(f):
    return make_ncdelta_slcdf(f, apply_preselection=True, include_truth=True)


def make_ncdelta_slcdf_data(f):
    return make_ncdelta_slcdf(f, apply_preselection=True, include_truth=False)


def make_ncdelta_slcdf_mc_nopresel(f):
    return make_ncdelta_slcdf(f, apply_preselection=False, include_truth=True)


def make_ncdelta_slcdf_data_nopresel(f):
    return make_ncdelta_slcdf(f, apply_preselection=False, include_truth=False)


# ---------------------------------------------------------------------------
# True neutrinos (efficiency denominators)
# ---------------------------------------------------------------------------

def make_ncdelta_nudf(f):
    DETECTOR, RUN, ismc = fetch_metadata(f)
    N = _truth_columns(make_mcdf(f))
    N["true_leading_gamma_E"] = _leading_photon_E(f, "rec.mc.nu").reindex(N.index)
    N["detector"] = DETECTOR
    N["Run"] = RUN
    return N
