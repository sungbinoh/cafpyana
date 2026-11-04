"""True final-state kinematics for the NC Delta -> N gamma truth study.

One row per true neutrino (rec.mc.nu), index (entry, nu index), so it joins
to candidates on (entry, tmatch_idx) like make_ncdeltawc_truthdf, whose
columns it extends with:

  * true_np_35MeV, true_nn: proton count above the WC proton threshold,
    neutron count
  * g1_* / g2_*: the two most energetic GENIE-primary photons. E [GeV],
    costh (to the beam, z), gap [cm], conv (end process is pair
    conversion), lost (born in the TPC, but it leaves before it ends)
  * p1_*: the most energetic primary proton. KE [GeV], len [cm] (in-TPC
    track length, NaN if it starts outside the TPC), costh,
    and cos_gp, the angle to g1
  * pi0_*: the most energetic primary pi0 and its two decay photons
    (Geant4 daughters, so from rec.true_particles). pi0_E, pi0_g1_* /
    pi0_g2_* (leading / sub-leading photon, same leaves as g1_*),
    pi0_open (opening angle, deg), pi0_asym = (E1 - E2) / (E1 + E2)

gap = distance from where the photon is made (the vertex, or the pi0 decay
point, ~1e-6 cm from it) to where its trajectory ends: the pair conversion
(or photoabsorption, conv = 0). SBND's Geant4 does not keep EM shower
daughters, so an earlier Compton scatter, which does not end the photon, is
not seen; above ~50 MeV conversion dominates. The CAF stores only the in-TPC
part of a trajectory (start/end, -9999 if none), so gap is NaN unless the
photon is made and ends inside the TPC.
"""

import numpy as np
import pandas as pd

from pyanalib.pandas_helpers import loadbranches

from analysis_village.NCDeltaWC.make_ncdeltawc import make_ncdeltawc_truthdf

PROTON_MASS = 0.938272   # GeV
PHOTON_E_MIN = 0.020     # GeV; above the ~3-9 MeV nuclear de-excitation gammas
G4_CONV = 36             # caf::g4_process_ kG4conv
NOPOS = -9000.           # CAF positions are -9999 outside the TPC

NU_LEVELS = [0, 1]       # (entry, rec.mc.nu index)
PART_LEAVES = ["pdg", "G4ID", "parent", "genE", "length", "end_process", "contained",
               "gen.x", "gen.y", "gen.z", "start.x", "start.y", "start.z",
               "end.x", "end.y", "end.z", "genp.x", "genp.y", "genp.z"]


def _flat(df):
    df.columns = ["_".join(c for c in col if c).replace(".", "_") for col in df.columns]
    return df


def _particles(f, prefix):
    df = loadbranches(f["recTree"], [prefix + l for l in PART_LEAVES])
    for level in prefix.strip(".").split("."):
        df = df[level]
    return _flat(df)


def _with_dir(df):
    p = np.sqrt(df.genp_x**2 + df.genp_y**2 + df.genp_z**2)
    for c in "xyz":
        df["dir_" + c] = df["genp_" + c] / p
    return df


def _photon_cols(g, prefix):
    """Columns <prefix>_E/costh/gap/conv/lost for one photon row per neutrino."""
    born_in = (g.start_x > NOPOS) & (np.abs(g.start_x - g.gen_x) + np.abs(g.start_y - g.gen_y)
                                     + np.abs(g.start_z - g.gen_z) < 0.01)
    ends_in = born_in & (g.contained == 1)
    gap = np.sqrt((g.end_x - g.gen_x)**2 + (g.end_y - g.gen_y)**2 + (g.end_z - g.gen_z)**2)
    out = pd.DataFrame(index=g.index)
    out[prefix + "_E"] = g.genE
    out[prefix + "_costh"] = g.dir_z
    out[prefix + "_gap"] = gap.where(ends_in)
    out[prefix + "_conv"] = (g.end_process == G4_CONV).astype(float)
    out[prefix + "_lost"] = (born_in & ~ends_in).astype(float)
    return out


def _nth(df, n):
    """n-th most energetic row per neutrino (n = 0 is the leading one)."""
    df = df.sort_values("genE", ascending=False)
    return df.groupby(level=NU_LEVELS).nth(n).droplevel(-1) if len(df) else df.droplevel(-1)


def make_ncdeltawc_truthkinedf(f):
    N = make_ncdeltawc_truthdf(f)

    prim = _with_dir(_particles(f, "rec.mc.nu.prim."))
    tp = _particles(f, "rec.true_particles.")
    tp = _with_dir(tp[tp.pdg == 22])

    pdg = prim.pdg.abs()
    pKE = prim.genE - PROTON_MASS
    cnt = lambda m: m.groupby(level=NU_LEVELS).sum().reindex(N.index).fillna(0).astype(int)
    N["true_np_35MeV"] = cnt((pdg == 2212) & (pKE > 0.035))
    N["true_nn"] = cnt(pdg == 2112)

    # GENIE-primary photons
    gam = prim[(prim.pdg == 22) & (prim.genE > PHOTON_E_MIN)]
    g1, g2 = _nth(gam, 0), _nth(gam, 1)
    N = N.join(_photon_cols(g1, "g1")).join(_photon_cols(g2, "g2"))

    # leading proton
    p1 = _nth(prim[pdg == 2212], 0)
    P = pd.DataFrame({"p1_KE": p1.genE - PROTON_MASS, "p1_len": p1.length.where(p1.start_x > NOPOS), "p1_costh": p1.dir_z})
    N = N.join(P)
    gp = (g1[["dir_x", "dir_y", "dir_z"]] * p1[["dir_x", "dir_y", "dir_z"]].reindex(g1.index)).sum(axis=1, min_count=3)
    N["p1_cos_gp"] = gp.reindex(N.index)

    # leading pi0 and its decay photons (Geant4 daughters)
    pi0 = _nth(prim[prim.pdg == 111], 0)
    N["pi0_E"] = pi0.genE
    d = tp.reset_index()
    d.columns = ["entry", "i"] + list(d.columns[2:])
    lead = pd.DataFrame({"entry": pi0.index.get_level_values(0), "nu": pi0.index.get_level_values(1),
                         "pi0_G4ID": pi0.G4ID.values})
    d = d.merge(lead, left_on=["entry", "parent"], right_on=["entry", "pi0_G4ID"])
    d = d.set_index(["entry", "nu", "i"])
    d.index.names = N.index.names + ["i"]
    pg1, pg2 = _nth(d, 0), _nth(d, 1)
    N = N.join(_photon_cols(pg1, "pi0_g1")).join(_photon_cols(pg2, "pi0_g2"))
    cosop = (pg1[["dir_x", "dir_y", "dir_z"]] * pg2[["dir_x", "dir_y", "dir_z"]].reindex(pg1.index)).sum(axis=1, min_count=3)
    N["pi0_open"] = np.degrees(np.arccos(cosop.clip(-1, 1))).reindex(N.index)
    N["pi0_asym"] = (N.pi0_g1_E - N.pi0_g2_E) / (N.pi0_g1_E + N.pi0_g2_E)
    return N
