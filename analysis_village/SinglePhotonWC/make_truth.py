"""Truth table for the inclusive single-photon search (see cuts.py for the
signal definition).

One row per true neutrino (rec.mc.nu), index (entry, nu index), so Wire-Cell
candidates join on (entry, tmatch_idx). Columns, all prefixed true_:

  * E, pdg, iscc, isnc, genie_mode, w, resnum, vtx_x/y/z: rec.mc.nu
  * primary counts: npi0, neta, npipm, np_35MeV (protons above 35 MeV KE),
    nn, nmu, nmu_100MeV, mu_KE (leading muon, GeV), ne (e+-)
  * prompt photons (GENIE primaries and decay photons of primary hadrons):
      ngamma_20MeV  photons with E >= 20 MeV
      ngamma_rec    of those, reconstructable (shower starts >= 3 cm inside
                    the TPC)
      ngamma_out    of those, not reconstructable: converts outside the TPC
                    or < 3 cm from its boundary
      ndecayg_soft  decay photons below 20 MeV
      nshower       reconstructable showers after merging photons < 20 deg
                    apart; nmerged = ngamma_rec - nshower
  * shw1_* / shw2_*: the two most energetic showers. E [GeV] (sum of its
    photons), gap [cm] (vertex to the nearest conversion point), costh (of
    its leading photon, to the beam), nphot, origin (22 for a GENIE-primary
    photon, else the pdg of the primary hadron it comes from: 111, 221, ...)

The shower start is the photon trajectory end: SBND's Geant4 does not keep
EM shower daughters, and the CAF stores only in-TPC points (-9999 if none),
so an escaping photon ends on the boundary and fails the 3 cm inset.
"""

import numpy as np
import pandas as pd

from pyanalib.pandas_helpers import loadbranches

from analysis_village.SinglePhotonWC import cuts

PROTON_MASS = 0.938272   # GeV
MUON_MASS = 0.105658     # GeV
G4_DECAY = 3             # caf::g4_process_ kG4Decay
LEPTONS = [11, 12, 13, 14, 15, 16]
MAX_DECAY_DEPTH = 4      # eta -> 3pi0 -> 6gamma needs 2


def _load(f, prefix, leaves):
    df = loadbranches(f["recTree"], [prefix + l for l in leaves])
    for level in prefix.strip(".").split("."):
        df = df[level]
    df.columns = [("_".join(c for c in col if c) if isinstance(col, tuple) else col).replace(".", "_")
                  for col in df.columns]
    return df


def _prompt_photons(prim, tp):
    """Photons from the neutrino interaction: primaries, and the decay chains
    of primary hadrons. Columns entry, nu, origin plus the tp columns."""
    p = prim.reset_index()
    p.columns = ["entry", "nu", "iprim"] + list(p.columns[3:])
    t = tp.reset_index()
    t.columns = ["entry", "itp"] + list(t.columns[2:])
    key = ["entry", "G4ID"]

    out = [t.merge(p.loc[p.pdg == 22, key + ["nu"]], on=key).assign(origin=cuts.ORIGIN_PRIMARY)]
    front = p.loc[~p.pdg.abs().isin(LEPTONS + [22]), key + ["nu", "pdg"]].rename(columns={"pdg": "origin"})
    dec = t[t.start_process == G4_DECAY]
    for _ in range(MAX_DECAY_DEPTH):
        if not len(front):
            break
        kids = dec.merge(front.rename(columns={"G4ID": "parent"}), on=["entry", "parent"])
        out.append(kids[kids.pdg == 22])
        front = kids.loc[~kids.pdg.abs().isin(LEPTONS + [22]), key + ["nu", "origin"]]
    g = pd.concat(out, ignore_index=True)
    p_ = np.sqrt(g.genp_x**2 + g.genp_y**2 + g.genp_z**2)
    for c in "xyz":
        g["d" + c] = g["genp_" + c] / p_
    return g


def _merge_showers(r):
    """Shower label per reconstructable photon: photons closer than
    MERGE_ANGLE to each other (single linkage) form one shower."""
    r = r.sort_values(["entry", "nu", "genE"], ascending=[True, True, False]).reset_index(drop=True)
    r["shower"] = r.groupby(["entry", "nu"]).cumcount()
    multi = r[r.groupby(["entry", "nu"]).genE.transform("size") > 1]
    for _, g in multi.groupby(["entry", "nu"]):
        d = g[["dx", "dy", "dz"]].to_numpy()
        ang = np.degrees(np.arccos(np.clip(d @ d.T, -1, 1)))
        lab = np.arange(len(g))
        for i in range(len(g)):
            for j in range(i):
                if ang[i, j] < cuts.MERGE_ANGLE:
                    lab[lab == lab[i]] = lab[j]
        r.loc[g.index, "shower"] = lab
    return r


def make_sp_truthdf(f):
    nu = _load(f, "rec.mc.nu.", ["E", "pdg", "iscc", "isnc", "genie_mode", "w", "resnum",
                                 "position.x", "position.y", "position.z"])
    N = nu.rename(columns={"position_x": "vtx_x", "position_y": "vtx_y", "position_z": "vtx_z"})
    idx = N.index
    if not len(N):
        return N

    prim = _load(f, "rec.mc.nu.prim.", ["pdg", "G4ID", "genE"])
    pdg = prim.pdg.abs()
    cnt = lambda m: m.groupby(level=[0, 1]).sum().reindex(idx).fillna(0).astype(int)
    muKE = (prim.genE - MUON_MASS).where(pdg == 13)
    N["npi0"] = cnt(pdg == 111)
    N["neta"] = cnt(pdg == 221)
    N["npipm"] = cnt(pdg == 211)
    N["np_35MeV"] = cnt((pdg == 2212) & (prim.genE - PROTON_MASS > 0.035))
    N["nn"] = cnt(pdg == 2112)
    N["nmu"] = cnt(pdg == 13)
    N["nmu_100MeV"] = cnt(muKE >= cuts.MUON_KE_MAX)
    N["mu_KE"] = muKE.groupby(level=[0, 1]).max().reindex(idx)
    N["ne"] = cnt(pdg == 11)

    tp = _load(f, "rec.true_particles.", ["pdg", "G4ID", "parent", "start_process", "genE",
                                          "genp.x", "genp.y", "genp.z", "end.x", "end.y", "end.z"])
    g = _prompt_photons(prim, tp)
    vtx = N[["vtx_x", "vtx_y", "vtx_z"]].copy()
    vtx.index = vtx.index.set_names(["entry", "nu"])
    g = g.join(vtx, on=["entry", "nu"])
    g["gap"] = np.sqrt((g.end_x - g.vtx_x)**2 + (g.end_y - g.vtx_y)**2 + (g.end_z - g.vtx_z)**2)
    hard = g.genE >= cuts.PHOTON_E_MIN
    g["rec"] = hard & (cuts.dist_to_active_boundary(g.end_x, g.end_y, g.end_z) >= cuts.SHOWER_INSET)

    per_nu = lambda m: m.groupby([g.entry, g.nu]).sum().reindex(idx).fillna(0).astype(int).values
    N["ngamma_20MeV"] = per_nu(hard)
    N["ngamma_rec"] = per_nu(g.rec)
    N["ngamma_out"] = per_nu(hard & ~g.rec)
    N["ndecayg_soft"] = per_nu(~hard & (g.origin != cuts.ORIGIN_PRIMARY))

    r = _merge_showers(g[g.rec])
    s = r.groupby(["entry", "nu", "shower"]).agg(
        E=("genE", "sum"), gap=("gap", "min"), costh=("dz", "first"),
        nphot=("genE", "size"), origin=("origin", "first")).reset_index()
    s = s.sort_values(["entry", "nu", "E"], ascending=[True, True, False])
    s["rank"] = s.groupby(["entry", "nu"]).cumcount()
    s = s.set_index(["entry", "nu"])
    N["nshower"] = s.groupby(level=[0, 1]).size().reindex(idx.set_names(["entry", "nu"])).fillna(0).astype(int).values
    N["nmerged"] = N.ngamma_rec - N.nshower
    for k, cols in [(0, ["E", "gap", "costh", "nphot", "origin"]), (1, ["E", "gap"])]:
        sk = s[s["rank"] == k][cols].add_prefix("shw%d_" % (k + 1))
        sk = sk.reindex(idx.set_names(["entry", "nu"]))
        for c in sk:
            N[c] = sk[c].values
    return N.add_prefix("true_")
