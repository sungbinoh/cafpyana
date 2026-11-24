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
wire's residual range.

No CVMFS needed: the YZ maps are copied into data/sbnd_yz_maps_cafpyana.npz
(sbnd_data v01_41_00, identical at v01_42_00; float16, <= 0.2% from the
originals) and the lifetime database into ETAU_TABLES below (the active
IOVs: first run, east and west lifetime): v01_41_00, the version chi2pid
reads (default), and v01_42_00, the one the v10_14_02 production used.
Calibration(path) reads an sbnd_data directory instead.
"""

import os
import sqlite3

import numpy as np
import pandas as pd
import uproot

from analysis_village.SinglePhotonWC import wc_calo

SBND_DATA_CVMFS = "/cvmfs/sbnd.opensciencegrid.org/products/sbnd/sbnd_data/v01_41_00"   # chi2pid's sbnd_data_v
YZ_NPZ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "sbnd_yz_maps_cafpyana.npz")
# sbnd_data CalibrationDatabase/tpc_elifetime.db, active IOVs (etau_sce_spatial_east / _west)
ETAU_TABLES = {
    "v01_41_00": [   # (first run of the IOV, east, west) [ms], 38 active IOVs
        (18250, 60.9700, 36.8590), (18255, 52.0762, 32.3816), (18262, 116.1950, 47.5632), (18265, 47.7290, 38.8096),
        (18270, 57.1716, 36.0394), (18271, 125.5190, 54.6021), (18273, 76.8861, 40.2454), (18280, 67.7690, 61.8949),
        (18307, 111.8910, 77.1901), (18314, 127.2600, 65.0363), (18328, 49.8161, 34.3523), (18335, 72.2975, 39.3286),
        (18348, 104.3680, 32.7435), (18355, 37.5529, 41.6395), (18359, 65.4133, 35.1604), (18365, 146.4100, 55.9247),
        (18386, 107.8150, 52.5417), (18391, 73.7853, 38.2213), (18410, 51.1139, 43.6863), (18415, 68.1969, 37.3505),
        (18423, 43.6382, 37.6252), (18432, 56.9198, 35.9546), (18441, 90.4286, 61.7683), (18444, 55.2165, 48.9228),
        (18445, 43.3810, 32.3672), (18451, 197.2360, 101.2660), (18461, 58.7540, 40.9139), (18467, 80.8610, 37.8951),
        (18483, 49.0176, 40.3590), (18486, 67.5196, 182.2560), (18488, 65.7848, 56.1764), (18489, 73.2001, 49.3025),
        (18491, 145.7290, 41.9435), (18529, 64.6974, 33.2740), (18552, 56.0874, 41.0037), (18556, 50.2334, 40.6079),
        (18585, 24.5972, 51.5207), (18591, 78.1933, 36.6851),
    ],
    "v01_42_00": [   # (first run of the IOV, east, west) [ms], 49 active IOVs
        (18250, 42.6767, 30.0608), (18255, 48.0464, 32.8810), (18259, 77.2473, 39.8263), (18262, 49.4648, 43.0885),
        (18265, 74.1096, 39.3267), (18270, 80.8600, 37.8958), (18271, 47.7015, 42.7502), (18273, 50.1613, 41.6532),
        (18280, 57.3194, 41.2698), (18307, 57.0281, 46.3913), (18314, 68.4934, 51.6435), (18328, 125.3410, 53.9303),
        (18333, 109.9010, 52.1802), (18335, 102.2900, 50.2390), (18345, 84.5556, 51.2750), (18348, 65.8269, 40.1955),
        (18355, 57.5052, 55.0617), (18359, 91.1667, 46.1552), (18365, 53.3254, 37.8356), (18386, 66.7814, 49.9537),
        (18391, 89.2891, 46.6613), (18410, 49.7560, 41.3011), (18415, 72.3951, 40.5428), (18423, 53.4514, 34.9721),
        (18432, 48.4315, 35.3465), (18441, 50.6979, 41.3105), (18444, 52.0775, 32.3818), (18445, 57.1739, 36.0398),
        (18451, 86.0501, 35.5156), (18461, 62.8224, 39.8326), (18467, 64.6993, 33.2740), (18468, 49.0175, 40.3595),
        (18469, 62.5104, 40.8349), (18483, 114.1440, 54.0143), (18486, 152.0550, 56.4619), (18488, 125.5260, 54.6022),
        (18489, 97.5651, 46.6381), (18503, 66.0484, 55.3477), (18504, 45.3696, 33.5664), (18505, 10.1738, 9.3057),
        (18507, 14.7585, 13.7243), (18513, 31.4551, 22.2020), (18516, 57.3127, 30.7900), (18522, 83.0489, 46.5719),
        (18524, 136.5810, 62.8018), (18538, 103.6520, 68.3525), (18552, 218.8860, 69.6839), (18556, 112.9880, 56.2053),
        (18585, 177.2880, 96.9980),
    ],
}
YZ_MAPS = {"data": "YZmaps/yz_correction_map_data1e20.root", "mc": "YZmaps/yz_correction_map_mcp2025b5e18.root"}
ETAU_MC_MS = 35.0
EFIELD, RHO, WION = 0.5, 1.38434, 23.6e-6
TRK_LEVELS = ["file", "__ntuple", "entry", "rec.slc..index", "rec.slc.reco.pfp..index"]


class Calibration:
    def __init__(self, path=None, lifetime_db="v01_41_00"):
        """path: an sbnd_data directory to read the maps and lifetime database from; None: the embedded copies
        (lifetime_db: "v01_41_00" or "v01_42_00")."""
        if path is None:
            z = np.load(YZ_NPZ)
            ze, ye = z["z_edges"], z["y_edges"]
            self.maps = {(k, p, t): (ze, ye, z["%s_p%d_t%d" % (k, p, t)].astype(float))
                         for k in ("data", "mc") for p in (0, 1, 2) for t in (0, 1)}
            tab = pd.DataFrame(ETAU_TABLES[lifetime_db], columns=["run", "etau_sce_spatial_east", "etau_sce_spatial_west"])
            tab["iov_id"] = np.arange(len(tab))
            self.iov = tab[["iov_id", "run"]]
            self.etau = tab.set_index("iov_id")[["etau_sce_spatial_east", "etau_sce_spatial_west"]]
            return
        d = path
        if not os.path.exists(os.path.join(d, YZ_MAPS["data"])):
            raise SystemExit("no sbnd_data YZ maps under %s" % d)
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
