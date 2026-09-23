import numpy as np
import pandas as pd
from makedf.makedf import (
    loadbranches, make_slcdf, make_trkdf, make_trkhitdf, make_crthitdf, make_opflashdf,
    make_hdrdf, make_triggerdf, make_potdf_bnb, make_mcnudf, trkmcsbranches, barycenterFMbranches
)

from makedf import chi2pid
from pyanalib.pandas_helpers import *

def fetch_metadata(f):
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det

    if det.empty:
        return pd.DataFrame()
    if 1 == det.unique():
        DETECTOR = "SBND"
    elif 2 == det.unique():
        DETECTOR = "ICARUS"
    else:
        raise ValueError("df maker needs rec.hdr.det == 1 (SBND) or 2 (ICARUS); got %s" % det.unique())
    run = loadbranches(f["recTree"], ["rec.hdr.run"]).rec.hdr.run
    RUN = 1 if DETECTOR == "SBND" else (2 if run.iloc[0] < 12960 else 4)
    ismc = bool(loadbranches(f["recTree"], ["rec.hdr.ismc"]).rec.hdr.ismc.iloc[0])
    return DETECTOR, RUN, ismc

def make_mcs_slcdf(f):
    DETECTOR, RUN, ismc = fetch_metadata(f)

    slcdf = make_slcdf(f)
    S = pd.DataFrame({
        "slc_vtx_x": slcdf.slc.vertex.x,
        "slc_vtx_y": slcdf.slc.vertex.y,
        "slc_vtx_z": slcdf.slc.vertex.z,
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
        "ismc" : ismc,
    })
    S["slice_index"] = S.index.get_level_values(1)
    S["detector"] = DETECTOR
    S["Run"] = RUN
    S["ismc"] = ismc
    return S

def make_geom1u1p_trkdf(f):
    DETECTOR, RUN, ismc = fetch_metadata(f)
    trkdf = make_trkdf(f, False, mcs=True)

    slcdf = make_slcdf(f)

    barycenterFM_df = loadbranches(f["recTree"], barycenterFMbranches).rec

    S = pd.DataFrame({
        "slc_vtx_x": slcdf.slc.vertex.x,
        "slc_vtx_y": slcdf.slc.vertex.y,
        "slc_vtx_z": slcdf.slc.vertex.z,
        "nu_score": slcdf.slc.nu_score,
        "baryscore":barycenterFM_df.slc.barycenterFM.score})

    trkhitdf = make_trkhitdf(f)
    chi2u = chi2pid.chi2u(trkhitdf)[0]
    chi2p = chi2pid.chi2p(trkhitdf)[0]
    chi2pi = chi2pid.chi2pi(trkhitdf)[0]

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
        "trackScore": trkdf.pfp.trackScore,
        # PID
        "chi2u": chi2u,
        "chi2p": chi2p,
        "chi2pi": chi2pi,

    })

    P = P.join(S[["slc_vtx_x", "slc_vtx_y", "slc_vtx_z", "nu_score", "baryscore"]])

    dist_start = np.sqrt((P.start_x - P.slc_vtx_x)**2 + (P.start_y - P.slc_vtx_y)**2 + (P.start_z - P.slc_vtx_z)**2)
    P["dist_start"] = dist_start
    P["prim_pfp"] = trkdf.pfp.parent_is_primary.fillna(False).astype(bool)
    P["detector"] = DETECTOR
    P["Run"] = RUN
    P["ismc"] = ismc
    return P

def make_geom1u1p_trkdf_mc(f):
    P = make_geom1u1p_trkdf(f)
    trkdf = make_trkdf(f, False, mcs=False)

    TrueP = pd.DataFrame({
        "true_pdg": trkdf.pfp.trk.truth.p.pdg,
        "true_end_x": trkdf.pfp.trk.truth.p.start.x,
        "true_end_x": trkdf.pfp.trk.truth.p.start.y,
        "true_end_x": trkdf.pfp.trk.truth.p.start.z,
        "true_end_x": trkdf.pfp.trk.truth.p.end.x,
        "true_end_y": trkdf.pfp.trk.truth.p.end.y,
        "true_end_z": trkdf.pfp.trk.truth.p.end.z,
    })

    return P.join(TrueP)

def make_mcs_trkdf(f):
    DETECTOR, RUN, ismc = fetch_metadata(f)
    trkdf = make_trkdf(f, False, mcs=True)
    keys = set(f["recTree"].keys())

    trkhitdf = make_trkhitdf(f)
    chi2u = chi2pid.chi2u(trkhitdf)[0]
    chi2p = chi2pid.chi2p(trkhitdf)[0]
    chi2pi = chi2pid.chi2pi(trkhitdf)[0]

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
        "trackScore": trkdf.pfp.trackScore,
        "true_pdg": trkdf.pfp.trk.truth.p.pdg,
        "true_end_x": trkdf.pfp.trk.truth.p.start.x,
        "true_end_x": trkdf.pfp.trk.truth.p.start.y,
        "true_end_x": trkdf.pfp.trk.truth.p.start.z,
        "true_end_x": trkdf.pfp.trk.truth.p.end.x,
        "true_end_y": trkdf.pfp.trk.truth.p.end.y,
        "true_end_z": trkdf.pfp.trk.truth.p.end.z,
        # PID
        "chi2u": chi2u,
        "chi2p": chi2p,
        "chi2pi": chi2pi,
        # MCS
        "maxlen": trkdf.pfp.trk.mcsP.maxlen,
        "maxang": trkdf.pfp.trk.mcsP.maxang,
        "minlen": trkdf.pfp.trk.mcsP.minlen,
        "minang": trkdf.pfp.trk.mcsP.minang,
        "avglen": trkdf.pfp.trk.mcsP.avglen,
        "avgang": trkdf.pfp.trk.mcsP.avgang,
        "stdlen": trkdf.pfp.trk.mcsP.stdlen,
        "stdang": trkdf.pfp.trk.mcsP.stdang,
    })

    P["prim_pfp"] = trkdf.pfp.parent_is_primary.fillna(False).astype(bool)
    P["detector"] = DETECTOR
    P["Run"] = RUN
    P["ismc"] = ismc
    return P
