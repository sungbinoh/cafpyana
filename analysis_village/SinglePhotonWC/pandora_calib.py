"""Pandora per-wire charge (pandwire table) through the cafpyana dE/dx
calibration chain (makedf/chi2pid.dqdx(gain="SBND", calibrate="SBND") and
dedx), without importing chi2pid (it reads CVMFS at import):

  dQ/dx [e/cm] = hit integral / pitch x exp(t_drift / tau) x YZ(y, z) / gain
  dE/dx        = the ellipsoid modified box inverted (makedf/calo.py)

  YZ maps      sbnd_data YZmaps/yz_correction_map_data1e20.root (data),
               yz_correction_map_mcp2025b5e18.root (MC), per plane and TPC
  lifetime     data: CalibrationDatabase/tpc_elifetime.db per run (IOV), east /
               west; MC: 35 ms (SBND_CALO_PARAMS["etau"], the simulated value).
               t_drift = t / 2000 - 0.2 ms with t the window centre in ticks
  gain         chi2pid SBND_CALO_PARAMS["gains"] (wc_calo.PANDORA_GAIN)
  recombination  wc_calo.RECOMB: "mc" (ICARUS, the simulation) or "data"
               (the cafpyana SBND data fit); 0.5 kV/cm, 1.38434 g/cm3

The reconstruction's own CAF dQ/dx differs: summed snippet integrals, the
sbndcode YZ maps, and no lifetime correction of the MC (its NormTools are YZ
only although the simulation uses 35 ms): ~2-2.5% on data / MC.

The table has one row per wire without the hit's y, z: they are taken on the
straight line from the track end (pandtrk end_*) towards its start at the
wire's residual range. sbnd_data: --sbnd-data, else $SBND_DATA_DIR, else the
CVMFS version chi2pid uses.
"""

import os
import sqlite3

import numpy as np
import pandas as pd
import uproot

from analysis_village.SinglePhotonWC import wc_calo

SBND_DATA_CVMFS = "/cvmfs/sbnd.opensciencegrid.org/products/sbnd/sbnd_data/v01_41_00"   # chi2pid's sbnd_data_v
YZ_MAPS = {"data": "YZmaps/yz_correction_map_data1e20.root", "mc": "YZmaps/yz_correction_map_mcp2025b5e18.root"}
ETAU_MC_MS = 35.0
EFIELD, RHO, WION = 0.5, 1.38434, 23.6e-6
TRK_LEVELS = ["file", "__ntuple", "entry", "rec.slc..index", "rec.slc.reco.pfp..index"]


def sbnd_data_dir(path=None):
    d = path or os.environ.get("SBND_DATA_DIR") or SBND_DATA_CVMFS
    if not os.path.exists(os.path.join(d, YZ_MAPS["data"])):
        raise SystemExit("no sbnd_data YZ maps under %s: pass --sbnd-data <sbnd_data directory>" % d)
    return d


class Calibration:
    def __init__(self, path=None):
        d = sbnd_data_dir(path)
        self.maps = {}
        for kind, fn in YZ_MAPS.items():
            f = uproot.open(os.path.join(d, fn))
            for tpc in (0, 1):
                for plane in (0, 1, 2):
                    h = f["CzyHist_%d_%d" % (plane, tpc)]
                    self.maps[(kind, plane, tpc)] = (h.axis(0).edges(), h.axis(1).edges(), h.values())
        con = sqlite3.connect(os.path.join(d, "CalibrationDatabase/tpc_elifetime.db"))
        et = pd.read_sql("select * from tpc_elifetime_data", con).set_index("__iov_id")
        iov = pd.read_sql("select * from tpc_elifetime_iovs where active=1", con)
        con.close()
        iov["run"] = iov.begin_time % 1000000000
        self.iov = iov.sort_values("run")
        self.etau = et[["etau_sce_spatial_east", "etau_sce_spatial_west"]]

    def yz(self, kind, plane, tpc, y, z):
        out = np.ones(len(y))
        for t in (0, 1):
            s = tpc == t
            if not s.any():
                continue
            ze, ye, v = self.maps[(kind, plane, t)]
            zi = np.clip(np.searchsorted(ze, z[s]) - 1, 0, len(ze) - 2)
            yi = np.clip(np.searchsorted(ye, y[s]) - 1, 0, len(ye) - 2)
            val = v[zi, yi]
            out[s] = np.where(val < 1e-3, 1.0, val)
        return out

    def lifetime(self, runs, tpc):
        """Data electron lifetime [ms] per row, east (tpc 0) / west (tpc 1), from the run's IOV."""
        k = np.searchsorted(self.iov.run.to_numpy(), runs, side="right") - 1
        ids = np.where(k >= 0, self.iov.iov_id.to_numpy()[np.clip(k, 0, None)], -1)
        e = self.etau.reindex(ids)
        return np.where(tpc == 0, e.etau_sce_spatial_east.to_numpy(), e.etau_sce_spatial_west.to_numpy())


def dedx(dqdx_e, phi_deg, recomb):
    """dE/dx [MeV/cm] from dQ/dx [e/cm], ellipsoid modified box (makedf/calo.recombination_cor)."""
    a, b90, r = recomb
    phi = np.radians(phi_deg)
    beta = b90 / np.sqrt(np.sin(phi) ** 2 + np.cos(phi) ** 2 / r ** 2) / (RHO * EFIELD)
    return (np.exp(dqdx_e * WION * beta) - a) / beta


def calibrate(wires, trk, cal, isdata, runs=None):
    """Per-wire columns dqdx_e (e/cm), dedx (the sample's own recombination: "data" for data, "mc" for MC)
    and dedx_mcrecomb / dedx_datarecomb (one model on both). wires: pandwire rows; trk: pandtrk (any plane,
    for the track end points); runs: per-row run numbers (data)."""
    ends = trk[["start_x", "start_y", "start_z", "end_x", "end_y", "end_z"]]
    ends = ends[~ends.index.droplevel("plane").duplicated()]
    ends.index = ends.index.droplevel("plane")
    w = wires.join(ends, on=TRK_LEVELS, how="inner")
    rr = 0.5 * (w.rr_min + w.rr_max)
    L = np.sqrt((w.start_x - w.end_x) ** 2 + (w.start_y - w.end_y) ** 2 + (w.start_z - w.end_z) ** 2).replace(0, np.nan)
    y = (w.end_y + rr * (w.start_y - w.end_y) / L).to_numpy(float)
    z = (w.end_z + rr * (w.start_z - w.end_z) / L).to_numpy(float)
    tpc = (w.x > 0).to_numpy().astype(int)
    plane = w.index.get_level_values("plane").to_numpy()
    tdrift = (0.5 * (w.t_start + w.t_end) / 2000.0 - 0.2).to_numpy(float)
    kind = "data" if isdata else "mc"
    s = np.ones(len(w))
    for p in range(3):
        m = plane == p
        s[m] = cal.yz(kind, p, tpc[m], y[m], z[m])
    if isdata:
        r = runs.reindex(w.index.droplevel(["rec.slc..index", "rec.slc.reco.pfp..index", "plane", "wire"])).to_numpy()
        tau = cal.lifetime(np.nan_to_num(r, nan=-1).astype(int), tpc)
    else:
        tau = np.full(len(w), ETAU_MC_MS)
    gain = np.asarray(wc_calo.PANDORA_GAIN[kind])[plane]
    q = (w.integral / w.pitch).to_numpy(float) * np.exp(tdrift / tau) * s / gain
    phi = w.phi.to_numpy(float)
    return w.assign(dqdx_e=q, dedx=dedx(q, phi, wc_calo.RECOMB[kind]), dedx_mcrecomb=dedx(q, phi, wc_calo.RECOMB["mc"]),
                    dedx_datarecomb=dedx(q, phi, wc_calo.RECOMB["data"]))
