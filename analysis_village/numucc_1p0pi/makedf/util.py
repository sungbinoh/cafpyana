import pandas as pd
import numpy as np
from pyanalib.variable_calculator import *
from pyanalib.pandas_helpers import *
from makedf.constants import *
from makedf.util import *


def cut_clear_cosmic(df):
    return df[df.slc.is_clear_cosmic == 0]

def cut_vertex_in_fv(df, det="SBND"):
    return df[InFV(df.slc.vertex, det=det)]

def cut_nu_score(df, th=0.5):
    return df[df.slc.nu_score > th]

def get_valid_trks(df):
    return df[df.pfp.trk.producer != 4294967295]

def cut_good_trks(trkdf):
    mask = (trkdf.pfp.trk.len > 0) &\
         (trkdf.pfp.pfochar.vtxdist < 100) #&\
    return trkdf[mask]

def get_trk_info(evtdf, trkdf, save_ntrks=3):
    nlevels = len(trkdf.index.names)
    ntrks = trkdf.pfp.id.groupby(level=list(range(nlevels-1))).count()
    ntrks.reindex(evtdf.index, fill_value=0)
    evtdf["n_trks"] = ntrks

    good_trks = cut_good_trks(trkdf).copy()
    ntrks = good_trks.pfp.id.groupby(level=list(range(nlevels-1))).count()
    ntrks.reindex(evtdf.index, fill_value=0)
    evtdf["n_good_trks"] = ntrks

    trks_sorted = trkdf.sort_values(by=('pfp','trk','len'), ascending=False)
    good_trks_sorted = good_trks.sort_values(by=('pfp','trk','len'), ascending=False)
    # get 'ntrks' longest tracks
    for i in range(save_ntrks):
        trk_i = good_trks_sorted.groupby(level=list(range(nlevels-1))).nth(i)
        trk_i.columns = pd.MultiIndex.from_tuples([tuple(["nocut_trk" + str(i+1)] + list(c)) for c in trk_i.columns])
        evtdf = multicol_merge(evtdf, trk_i.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")

        good_trk_i = good_trks_sorted.groupby(level=list(range(nlevels-1))).nth(i)
        good_trk_i.columns = pd.MultiIndex.from_tuples([tuple(["trk" + str(i+1)] + list(c)) for c in good_trk_i.columns])
        evtdf = multicol_merge(evtdf, good_trk_i.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")

    return evtdf


def cut_2prong(df):
    # return df[(df.n_good_trks == 2) & (df.n_trks <= 3)]
    return df[(df.n_good_trks == 2)]

def cut_2prong_contained(df, det="SBND"):
    return df[InFV(df.trk1.pfp.trk.start, det=det) & InFV(df.trk1.pfp.trk.end, det=det) \
        & InFV(df.trk2.pfp.trk.start, det=det) & InFV(df.trk2.pfp.trk.end, det=det)]

def cut_2prong_trackscore(df, trackscore_th=0.5):
    return df[(df.trk1.pfp.trackScore > trackscore_th) & (df.trk2.pfp.trackScore > trackscore_th)]

def cut_2prong_vtxdist(df, vtxdist_th=1.5):
    return df[(df.trk1.pfp.pfochar.vtxdist < vtxdist_th) & (df.trk2.pfp.pfochar.vtxdist < vtxdist_th)]

def get_mu_p_candidate(df, 
                       mu_chi2mu_th=30, mu_chi2p_th=100, mu_len_th=50, qual_th=0.25,
                       p_chi2mu_th=30, p_chi2p_th=90, p_len_th=0):

    nlevels = len(df.index.names)

    trks = pd.concat([df.trk1, df.trk2])

    chimu_avg = avg_chi2(trks, "chi2_muon")
    chip_avg = avg_chi2(trks, "chi2_proton")

    mcs_range_diff = np.abs((trks.pfp.trk.rangeP.p_muon - trks.pfp.trk.mcsP.fwdP_muon) / trks.pfp.trk.rangeP.p_muon)

    mu_cut = (chimu_avg > 0) & (chimu_avg < mu_chi2mu_th) & \
            (chip_avg > mu_chi2p_th) & \
            (trks.pfp.trk.len > mu_len_th) & \
            (mcs_range_diff < qual_th)

    mu_candidate = trks[mu_cut]
    mu_candidate = mu_candidate.groupby(level=list(range(nlevels))).nth(0)

    mu_candidate.columns = pd.MultiIndex.from_tuples([tuple(["mu"] + list(c)) for c in mu_candidate.columns])
    df = multicol_merge(df, mu_candidate, left_index=True, right_index=True, how="left", validate="one_to_one")

    # TODO: keep & use original trk index?
    not_mu_candidate = pd.concat([trks[~mu_cut], trks[mu_cut].groupby(level=list(range(nlevels))).nth(1)])
    chip_avg = avg_chi2(not_mu_candidate, "chi2_proton")
    p_candidate = not_mu_candidate[(chip_avg > 0) & (chip_avg < p_chi2p_th) & (not_mu_candidate.pfp.trk.len > p_len_th)]
    p_candidate = p_candidate.groupby(level=list(range(nlevels))).nth(0)

    p_candidate.columns = pd.MultiIndex.from_tuples([tuple(["p"] + list(c)) for c in p_candidate.columns])
    df = multicol_merge(df, p_candidate, left_index=True, right_index=True, how="left", validate="one_to_one")

    return df

def cut_has_mu(df):
    return df[~np.isnan(df.mu.pfp.trk.producer)]

def cut_has_p(df):
    return df[~np.isnan(df.p.pfp.trk.producer)]

def cut_mu_kinematics(df, mu_Plo_th=0.22, mu_Phi_th= 1):
    return df[(df.mu.pfp.trk.rangeP.p_muon > mu_Plo_th) & (df.mu.pfp.trk.rangeP.p_muon < mu_Phi_th)]

def cut_p_kinematics(df, p_Plo_th=0.3, p_Phi_th= 1):
    return df[(df.p.pfp.trk.rangeP.p_proton > p_Plo_th) & (df.p.pfp.trk.rangeP.p_proton < p_Phi_th)]