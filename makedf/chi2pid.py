from . import calo
import pandas as pd
import numpy as np
import sqlite3
import uproot

larsoft_data_v = "v1_02_02"
icarus_data_v = "v10_06_00"
sbnd_data_v = "v01_35_00"

rr_max_cut_chi2 = 25. ## for resolving MC's hit RR cut, after fixing the issue, put this value to 26.

#### == use pandora_df_calo_update to apply these changes
ICARUS_CALO_PARAMS = {
    "alpha_emb": 0.904,
    "beta_90": 0.204,
    "R_emb": 1.25,
    "gains": [0.016751, 0.012755, 0.012513],
    "c_cal_frac": [1., 1., 1.],
}

#### Default
#### == For each element, the first entry is for MC and the second entry is for the data
SBND_CALO_PARAMS = {
    "alpha_emb": [0.904, 0.904],
    "beta_90": [0.204, 0.204],
    "R_emb": [1.25, 1.25],
    "gains": [
        [0.0203521, 0.0202351, 0.0200727], ## MC
        [0.0223037, 0.0219534, 0.0215156]], ## Data
    "c_cal_frac": [1., 1., 1.],
    "etau": [100., 35.], ## first value for MC and second value for data
}


def chi2(hitdf, exprr, expdedx, experr, dedxname="dedx"):
    dedx_exp = pd.cut(hitdf.rr, exprr, labels=expdedx).astype(float)
    dedx_err = pd.cut(hitdf.rr, exprr, labels=experr).astype(float)

    dedx_res = (0.04231 + 0.0001783*hitdf[dedxname]**2)*hitdf[dedxname]

    v_chi2 = (hitdf[dedxname] - dedx_exp)**2 / (dedx_err**2 + dedx_res**2)

    when_chi2 = (hitdf.rr < rr_max_cut_chi2) & ~hitdf.firsthit & ~hitdf.lasthit & (hitdf[dedxname] < 1000.)

    chi2_group = v_chi2[when_chi2].groupby(level=list(range(hitdf.index.nlevels-1)))

    return chi2_group.sum() / chi2_group.size(), chi2_group.size()

def chi2u(hitdf, dedxname="dedx"):
    return chi2(hitdf, muon_rr, muon_dedx, muon_yerr, dedxname)

def chi2p(hitdf, dedxname="dedx"):
    return chi2(hitdf, proton_rr, proton_dedx, proton_yerr, dedxname)

def chi2par(hitdf, dedxname="dedx", par=""):
    if par == "muon":
        return chi2u(hitdf, dedxname)
    elif par == "proton":
        return chi2p(hitdf, dedxname)
    else:
        raise ValueError(f"Invalid par={par!r}. Expected 'muon' or 'proton'.")

def chi2_ndof(hitdf):
    when_chi2 = (hitdf.rr < rr_max_cut_chi2) & ~hitdf.firsthit & ~hitdf.lasthit & (hitdf.dedx < 1000.)
    chi2_group = v_chi2[when_chi2].groupby(level=list(range(hitdf.index.nlevels-1)))

    return chi2_group.size()

def dqdx(dqdxdf, gain=None, calibrate=None, isMC=False):
    if calibrate == "ICARUS": 
        # get raw dqdx
        dqdx = dqdxdf.integral / dqdxdf.pitch

        # compute y-scale
        ybin = _yz_ybin(dqdxdf.y, yz_ybin)
        zbin = _yz_zbin(dqdxdf.z, yz_zbin)
        iov = _yz_iov(dqdxdf.run)
        itpc = dqdxdf.tpc // 2 + dqdxdf.cryo*2
        plane = dqdxdf.plane

        yzdf = pd.DataFrame({"ybin": ybin, "zbin": zbin, "itpc": itpc, "plane": plane, "iov": iov})
        yz_scale = yzdf.merge(IC_yz_cal_df, on=["iov", "itpc", "plane", "ybin", "zbin"], how="left", validate="many_to_one").scale
        yz_scale[yz_scale == -999.000000] = 1.
        yz_scale = np.clip(yz_scale, 0.7, 1.3).fillna(1)
        yz_scale.index = dqdxdf.index

        # compute lifetime correction
        iov = _etau_iov(dqdxdf.run)
        etaudf = pd.DataFrame({"itpc": itpc, "iov": iov})
        etau = etaudf.merge(IC_etau_cal_df, on=["iov", "itpc"], how="left", validate="many_to_one").etau
        etau = etau.fillna(np.inf)
        etau.index = dqdxdf.index
        etau[dqdxdf.run == 1] = 3.5e3 # set MC lifetime to MC default

        # compute TPC scale
        iov = _tpc_iov(dqdxdf.run)
        tpcdf = pd.DataFrame({"itpc": itpc, "plane": plane, "iov": iov})
        tpc_scale = tpcdf.merge(IC_tpc_cal_df, on=["iov", "itpc", "plane"], how="left", validate="many_to_one").scale
        tpc_scale.index = dqdxdf.index

        # apply the corrections
        t0 = 0 # assume in time
        tick_period = 0.4 # us
        tanode = 850 # ticks
        tdrift = dqdxdf.t*tick_period - t0 - tanode*tick_period

        dqdx = dqdx * tpc_scale * np.exp(tdrift / etau) / yz_scale
    elif calibrate == "SBND": # TODO: add calibrations?
        # get raw dqdx
        dqdx  = dqdxdf.integral / dqdxdf.pitch

        this_yz_cal_df = SBND_yz_cal_mc_df   if isMC else SBND_yz_cal_data_df
        this_yz_zbin = yz_zbin_sbnd_mc if isMC else yz_zbin_sbnd_data
        this_yz_ybin = yz_ybin_sbnd_mc if isMC else yz_ybin_sbnd_data
        this_etau_df   = SBND_etau_cal_mc_df if isMC else SBND_etau_cal_data_df

        # compute y-scale
        ybin = _yz_ybin(dqdxdf.y, this_yz_ybin)
        zbin = _yz_zbin(dqdxdf.z, this_yz_zbin)
        bin_test_df = pd.DataFrame({"y":dqdxdf.y, "z": dqdxdf.z, "ybin": ybin, "zbin": zbin, "dqdx": dqdx})
        dqdxdf['iov'] = 0 ## FIXME: once SBND has time dep. calo, it should be updated
        iov = dqdxdf.iov ## FIXME: once SBND has time dep. calo, it should be updated
        itpc = dqdxdf.tpc
        plane = dqdxdf.plane

        yzdf = pd.DataFrame({"ybin": ybin, "zbin": zbin, "itpc": itpc, "plane": plane, "iov": iov})
        yz_scale = yzdf.merge(this_yz_cal_df, on=["iov", "itpc", "plane", "ybin", "zbin"], how="left", validate="many_to_one").scale
        yz_scale[yz_scale < 1e-6] = 1.
        yz_scale = yz_scale.fillna(1)
        yz_scale.index = dqdxdf.index

        #yzdf['rr'] = dqdxdf.rr
        #yzdf['scale'] = yz_scale
        #print(yzdf[yzdf.rr < 26.].head(50))
        # compute lifetime correction
        iov = dqdxdf.iov ## FIXME: once SBND has time dep. calo, it should be updated
        etaudf = pd.DataFrame({"itpc": itpc, "iov": iov})
        etau = etaudf.merge(this_etau_df, on=["iov", "itpc"], how="left", validate="many_to_one").etau
        etau = etau.fillna(np.inf)
        etau.index = dqdxdf.index

        # apply the corrections
        t0 = 0 # assume in time
        tdrift = dqdxdf.t / 2000. - 0.2
        dqdx = dqdx * np.exp(tdrift / etau) * yz_scale
        #dqdx = dqdx / yz_scale

    else: # if not specified, rely on input calibration
        dqdx = dqdxdf.dqdx

    # apply gain
    if gain == "ICARUS":
        gains = ICARUS_CALO_PARAMS["gains"]
        gain_perhit = pd.Series(1.0, dqdxdf.index)
        for iplane, g in enumerate(gains):
            gain_perhit[dqdxdf.plane == iplane] = 1.0/g
    elif gain == "SBND": # TODO
        gains = SBND_CALO_PARAMS["gains"][0] if isMC else SBND_CALO_PARAMS["gains"][1]
        gain_perhit = pd.Series(1.0, dqdxdf.index)
        for iplane, g in enumerate(gains):
            gain_perhit[dqdxdf.plane == iplane] = 1.0/g
    else:
        gain_perhit = 1

    return dqdx*gain_perhit

def dedx(dqdxdf, gain=None, calibrate=None, plane=2, isMC=False, smear=-1, scale=1):
    dqdx_v = dqdx(dqdxdf, gain=gain, calibrate=calibrate, isMC=isMC)
    if gain == "ICARUS":
        scalegain = ICARUS_CALO_PARAMS['c_cal_frac'][plane]
    elif gain == "SBND":
        scalegain = SBND_CALO_PARAMS['c_cal_frac'][plane]
    else:
        scalegain = 1.

    if gain == "SBND":
        this_alpha_emb = SBND_CALO_PARAMS["alpha_emb"][0] if isMC else SBND_CALO_PARAMS["alpha_emb"][1]
        this_beta_90 = SBND_CALO_PARAMS["beta_90"][0] if isMC else SBND_CALO_PARAMS["beta_90"][1]
        this_R_emb = SBND_CALO_PARAMS["R_emb"][0] if isMC else SBND_CALO_PARAMS["R_emb"][1]
        dedx = calo.recombination_cor(scale*dqdx_v/scalegain, dqdxdf.phi, dqdxdf.efield, dqdxdf.rho, this_alpha_emb, this_beta_90, this_R_emb)
    else:
        dedx = calo.recombination_cor(scale*dqdx_v/scalegain, dqdxdf.phi, dqdxdf.efield, dqdxdf.rho)

    if smear > 0:
        dedx = dedx*np.random.normal(loc=1., scale=smear, size=dedx.size)

    return dedx

def _yz_ybin(y, yz_ybin):
    return np.searchsorted(yz_ybin, y) - 1

def _yz_zbin(z, yz_zbin):
    return np.searchsorted(yz_zbin, z) - 1

def _yz_iov(run): 
    iov = __iov(run, IC_yz_cal_iovdf)
    iov[run == 1] = 4 # MC default to Run 4
    return iov

def _etau_iov(run):
    iov = __iov(run, IC_etau_cal_iovdf)
    iov[run == 1] = -1 # MC default to no run
    return iov

def _tpc_iov(run):
    iov = __iov(run, IC_tpc_cal_iovdf)
    iov[run == 1] = 3 # MC default to Run 4
    return iov
    

def __iov(run, df):
    return pd.cut(run, list(df.run) + [np.inf], labels=df.iov).astype(float).fillna(-1).astype(int)

##############################
# EXPECTED dE/dx FILES
##############################
datadir = "/cvmfs/larsoft.opensciencegrid.org/products/larsoft_data/" + larsoft_data_v + "/ParticleIdentification/"
fhist = datadir + "dEdxrestemplates.root"

profp = uproot.open(fhist)["dedx_range_pro"]
profmu = uproot.open(fhist)["dedx_range_mu"]

proton_dedx = profp.values()
proton_rr = profp.axis().edges()
proton_yerr = profp.errors(error_mode="s")
for i in range(len(proton_yerr)):
    if proton_yerr[i] < 1e-6:
        proton_yerr[i] = (proton_yerr[i-1] + proton_yerr[i+1]) / 2
    if proton_dedx[i] < 1e-6:
        proton_dedx[i] = (proton_dedx[i-1] + proton_dedx[i+1]) / 2

muon_dedx = profmu.values()
muon_rr = profmu.axis().edges()
muon_rr_center = profmu.axis().centers()
muon_yerr = profmu.errors(error_mode="s")

##############################
# ICARUS TPC calo files
##############################
# ICARUS CALIBRATION DATABASES
IC_yz_cal_f = "/cvmfs/icarus.opensciencegrid.org/products/icarus/icarus_data/" + icarus_data_v + "/icarus_data/database/tpc_yz_correction_allplanes_data.db"
IC_yz_cal_db = "tpc_yz_correction_allplanes_data_data"
IC_yz_cal_iov = "tpc_yz_correction_allplanes_data_iovs"

IC_etau_cal_f = "/cvmfs/icarus.opensciencegrid.org/products/icarus/icarus_data/" + icarus_data_v + "/icarus_data/database/tpc_elifetime_data.db"
IC_etau_cal_db = "tpc_elifetime_data_data"
IC_etau_cal_iov = "tpc_elifetime_data_iovs"

IC_tpc_cal_f = "/cvmfs/icarus.opensciencegrid.org/products/icarus/icarus_data/" + icarus_data_v + "/icarus_data/database/tpc_dqdxcalibration_allplanes_data.db"
IC_tpc_cal_db = "tpc_dqdxcalibration_allplanes_data_data"
IC_tpc_cal_iov = "tpc_dqdxcalibration_allplanes_data_iovs"

# LOAD THE YZ CALIBRATION
yz_ybin = np.linspace(-180, 130, 32)
yz_ylos = yz_ybin[:-1]
yz_yhis = yz_ybin[1:]
yz_ys = (yz_ylos + yz_yhis) / 2.

yz_zbin = np.linspace(-900, 900, 181)
yz_zlos = yz_zbin[:-1]
yz_zhis = yz_zbin[1:]
yz_zs = (yz_zlos + yz_zhis) / 2.

conn = sqlite3.connect(IC_yz_cal_f)
cursor = conn.cursor()
cursor.execute("SELECT * FROM %s" % IC_yz_cal_db)
rows = cursor.fetchall()
data = list(zip(*rows))
IC_yz_cal_df = pd.DataFrame({
  "iov": data[0],
  "plane": data[2],
  "tpc": data[3],
  "ybin": data[4],
  "zbin": data[5],
  "scale": data[6]
})
IC_yz_cal_df["itpc"] = 0 
IC_yz_cal_df.loc[IC_yz_cal_df.tpc == "EE", "itpc"] = 0
IC_yz_cal_df.loc[IC_yz_cal_df.tpc == "EW", "itpc"] = 1
IC_yz_cal_df.loc[IC_yz_cal_df.tpc == "WE", "itpc"] = 2
IC_yz_cal_df.loc[IC_yz_cal_df.tpc == "WW", "itpc"] = 3
del IC_yz_cal_df["tpc"]

cursor.execute("SELECT * FROM %s WHERE ACTIVE=1" % IC_yz_cal_iov)
rows = cursor.fetchall()
data = list(zip(*rows))
IC_yz_cal_iovdf = pd.DataFrame({
  "iov": data[0],
  "begin_time": data[1],
})
IC_yz_cal_iovdf["run"] = IC_yz_cal_iovdf.begin_time % 1000000000
IC_yz_cal_iovdf.sort_values(by="run", inplace=True)
conn.close()

# LOAD THE LIFETIME CALIBRATION
conn = sqlite3.connect(IC_etau_cal_f)
cursor = conn.cursor()
cursor.execute("SELECT * FROM %s" % IC_etau_cal_db)
rows = cursor.fetchall()
data = list(zip(*rows))
IC_etau_cal_df = pd.DataFrame({
  "iov": data[0],
  "itpc": data[1],
  "etau": data[2]
})

cursor.execute("SELECT * FROM %s WHERE ACTIVE=1" % IC_etau_cal_iov)
rows = cursor.fetchall()
data = list(zip(*rows))
IC_etau_cal_iovdf = pd.DataFrame({
  "iov": data[0],
  "begin_time": data[1],
})
IC_etau_cal_iovdf["run"] = IC_etau_cal_iovdf.begin_time % 1000000000
IC_etau_cal_iovdf.sort_values(by="run", inplace=True)
conn.close()

# LOAD THE TPC SCALE
conn = sqlite3.connect(IC_tpc_cal_f)
cursor = conn.cursor()
cursor.execute("SELECT * FROM %s" % IC_tpc_cal_db)
rows = cursor.fetchall()
data = list(zip(*rows))
IC_tpc_cal_df = pd.DataFrame({
  "iov": data[0],
  "plane": data[2],
  "tpc": data[3],
  "scale": data[4]
})
IC_tpc_cal_df["itpc"] = 0 
IC_tpc_cal_df.loc[IC_tpc_cal_df.tpc == "EE", "itpc"] = 0
IC_tpc_cal_df.loc[IC_tpc_cal_df.tpc == "EW", "itpc"] = 1
IC_tpc_cal_df.loc[IC_tpc_cal_df.tpc == "WE", "itpc"] = 2
IC_tpc_cal_df.loc[IC_tpc_cal_df.tpc == "WW", "itpc"] = 3
del IC_tpc_cal_df["tpc"]

cursor.execute("SELECT * FROM %s WHERE ACTIVE=1" % IC_tpc_cal_iov)
rows = cursor.fetchall()
data = list(zip(*rows))
IC_tpc_cal_iovdf = pd.DataFrame({
  "iov": data[0],
  "begin_time": data[1],
})
IC_tpc_cal_iovdf["run"] = IC_tpc_cal_iovdf.begin_time % 1000000000
IC_tpc_cal_iovdf.sort_values(by="run", inplace=True)
conn.close()

##############################
# SBND TPC calo files
##############################
SBND_yz_cal_mc_f = "/cvmfs/sbnd.opensciencegrid.org/products/sbnd/sbnd_data/" + sbnd_data_v + "/YZmaps/yz_correction_map_mcp2025b5e18.root"
SBND_yz_cal_data_f = "/cvmfs/sbnd.opensciencegrid.org/products/sbnd/sbnd_data/" + sbnd_data_v + "/YZmaps/yz_correction_map_data1e20.root"

yz_zbin_sbnd_mc = []
yz_ybin_sbnd_mc = []
yz_zbin_sbnd_data = []
yz_ybin_sbnd_data = []

def call_sbnd_yz_corr(map_f):
    maps = []
    z_edges = y_edges = None
    for tpc, plane in [(t, p) for t in range(2) for p in range(3)]:
        this_hist = uproot.open(map_f)["CzyHist_" + str(plane) + "_" + str(tpc)]
        this_yz_zbin = this_hist.axis(0).edges()
        this_yz_ybin = this_hist.axis(1).edges()

        if z_edges is None:
            z_edges, y_edges = this_yz_zbin, this_yz_ybin

        corr = this_hist.values()

        corr_mi = pd.MultiIndex.from_product(
            [range(corr.shape[0]), range(corr.shape[1])],
            names=["zbin", "ybin"]
        )

        corr_df = pd.Series(corr.ravel(), index=corr_mi, name="scale").to_frame()
        corr_df = corr_df.reset_index()

        corr_df['itpc'] = tpc
        corr_df['plane'] = plane
        corr_df['iov'] = 0
        new_order = ["iov", "plane", "ybin", "zbin", "scale", 'itpc']
        corr_df = corr_df[new_order]
        maps.append(corr_df)

    out_df = pd.concat(maps, ignore_index=True)
    return out_df, z_edges, y_edges

SBND_yz_cal_mc_df, yz_zbin_sbnd_mc, yz_ybin_sbnd_mc = call_sbnd_yz_corr(SBND_yz_cal_mc_f)
SBND_yz_cal_data_df, yz_zbin_sbnd_data, yz_ybin_sbnd_data = call_sbnd_yz_corr(SBND_yz_cal_data_f)

SBND_etau_cal_mc_df = pd.DataFrame( {'iov': [0, 0], 'itpc': [0, 1], 'etau': [SBND_CALO_PARAMS["etau"][0], SBND_CALO_PARAMS["etau"][0]]})
SBND_etau_cal_data_df = pd.DataFrame( {'iov': [0, 0], 'itpc': [0, 1], 'etau': [SBND_CALO_PARAMS["etau"][1], SBND_CALO_PARAMS["etau"][1]]})
