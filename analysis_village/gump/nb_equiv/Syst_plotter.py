import gump_cuts as gc
import pandas as pd
import os
import h5py
import sys
import loaddf
import syst as syst

from cycler import *
#from plot_tools import *
workspace_root = os.getcwd()
sys.path.insert(0, workspace_root + "/../../")
import pyanalib.pandas_helpers as ph
from pyanalib.split_df_helpers import *
from pyanalib.ntuple_glob import *
from makedf.util import *
from rwt_map import *

RECO = "PANDORA"
FONTSIZE=14

# import importlib
def scale_pot(df, pot, desired_pot):
    """Scale DataFrame by desired POT."""
    print(f"POT: {pot}\nScaling to: {desired_pot}")
    scale = desired_pot / pot
    df['glob_scale'] = scale
    return pot, scale

def add_style(ax, xlabel, title="", det="ICARUS", ylabel='Events / $10^{20}$ POT'):
    ax.tick_params(axis='both', which='both', direction='in', length=6, width=1.5, labelsize=FONTSIZE, top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    ax.set_xlabel(xlabel, fontsize=FONTSIZE, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=FONTSIZE, fontweight='bold')
    ax.set_title(f"$\\bf{{{det}}}$  {title}", fontsize=FONTSIZE+2)
    ax.legend(fontsize=FONTSIZE)

def v_variation(df, setvars):
    df = df[[c for c in df.columns if "univ" not in c]].copy()
    for (new, old) in setvars:
        df[new] = df[old]
    return df

def v_chi2smear(df):
    setvars = [
        ("mu_chi2_of_mu_cand", "mu_chi2smear_of_mu_cand"),
        ("mu_chi2_of_p_cand",  "mu_chi2smear_of_prot_cand"),
        ("prot_chi2_of_mu_cand", "prot_chi2smear_of_mu_cand"),
        ("prot_chi2_of_p_cand",  "prot_chi2smear_of_prot_cand"),
    ]
    return v_variation(df, setvars)


def v_chi2hi(df):
    setvars = [
        ("mu_chi2_of_mu_cand", "mu_chi2hi_of_mu_cand"),
        ("mu_chi2_of_p_cand",  "mu_chi2hi_of_prot_cand"),
        ("prot_chi2_of_mu_cand", "prot_chi2hi_of_mu_cand"),
        ("prot_chi2_of_p_cand",  "prot_chi2hi_of_prot_cand"),
    ]
    return v_variation(df, setvars)

def FV(df, det):    
    is_spine = "SPINE" in RECO
    
    ret = gc.slcfv_cut(df, det) & gc.mufv_cut(df, det) & gc.pfv_cut(df, det) 
    
    if is_spine:
        ret = ret & (df.is_time_contained)
    
    return ret
    
def simple_cosmic_rej(df):
    is_spine = "SPINE" in RECO
    return FV(df) & (df.crlongtrkdiry > -0.3)

def cosmic_cut(df):
    if "SPINE" not in RECO:
        return df.nu_score > 0.4
    else:
        return df.is_time_contained
        
def crtveto(df):
    return ~df.crthit

def twoprong_cut(df):
    return np.isnan(df.other_shw_length) & np.isnan(df.other_trk_length)

def pid_cut(df):
    is_spine = "SPINE" in RECO
    if not is_spine:
        return twoprong_cut(df) & gc.pid_cut_df(df)
    else:
        return twoprong_cut(df) & (df.prot_chi2_of_prot_cand > 0.6) & (df.mu_chi2_of_mu_cand > 0.6)

def signalbox(df, det):
    return FV(df, det) & cosmic_cut(df) & twoprong_cut(df) & pid_cut(df) # & (~df.crthit)

def ICARUS_dircut(df):
    vtx = pd.DataFrame({
        "x": df.true_vtx_x,
        "y": df.true_vtx_y,
        "z": df.true_vtx_z,
    })
    return ~np.isnan(df.true_vtx_x) & ~gc._fv_cut(vtx, "ICARUS", 0, 0, 0, 0)

def main():
    PLOTDIR = "/exp/sbnd/app/users/nrowe/cafpyana/analysis_village/gump/signal-box-syst-plots/"

    DOSAVE = True
    os.makedirs(PLOTDIR, exist_ok=True)
    os.makedirs(PLOTDIR + "/png", exist_ok=True)
    os.makedirs(PLOTDIR + "/pdf", exist_ok=True)

    FONTSIZE = 14
    HAWKS_COLORS = ["#315031", "#d54c28", "#1e3f54", "#c89648", "#43140b", "#95af8b"]

    #RECO = "PANDORA"

    DF_DIR = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-4/"
    
    # SCV_FILES = [DF_DIR + "SBND_SpringMC_rewgt_%i.df" % i for i in range(5)]
    SCV_FILES = [DF_DIR + "rewegt-try1/SBND_SpringMC_rewgt_%i.df" % i for i in range(5)]
    SDIRT_FILES = [DF_DIR + "rewegt-try1/SBND_SpringLowEMC_rewgt_%i.df" % i for i in range(10)]
    SBEAMOFF_FILE = DF_DIR + "SBND_SpringBNBOffData_5000.df"
    SDETVAR_FILES = [
        DF_DIR + "SBND_SpringMC_WMXThetaXW.df",
        DF_DIR + "SBND_SpringMC_WMYZ.df",
        DF_DIR + "SBND_SpringMC_0xSCE.df",
        DF_DIR + "SBND_SpringMC_2xSCE.df",
    ]
    SDETVAR_NAMES = ["WM $X\\theta_{xw}$", "WM $YZ$", "0x SCE", "2x SCE"]
    
    ICV_FILES = [DF_DIR + "ICARUS_SpringMCOverlay_rewgt.df"]
    IDIRT_FILES = [DF_DIR + "ICARUS_SpringMCDirt_slimwgt.df"]
    IBEAMOFF_FILE = DF_DIR + "ICARUS_SpringRun2BNBOff_unblind_prescaled.df"
    IDETVAR_FILES = []
    IDETVAR_NAMES = []

    #F = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-4/SBND_SpringBNBOffData_5000.df"

    #with h5py.File(F, "r") as f:
    #    print(f.keys())

    IGOAL_POT = 5e20
    SGOAL_POT = 1e20
    
    # importlib.reload(loaddf)
    # importlib.reload(gc)
    # importlib.reload(syst)
    
    Idf, Ipot = loaddf.loadl(ICV_FILES, njob=min(len(ICV_FILES), 10))
    Idirt, Idirtpot = loaddf.loadl(IDIRT_FILES, njob=min(len(IDIRT_FILES), 10), include_syst=False)
    Ioffbeam, Ioffbeampot = loaddf.load(IBEAMOFF_FILE, offbeampot_ICARUS=True,  include_syst=False, load_truth=False, hdrname="trig_%i")
    # Idetvars, Idetvar_pots = zip(*tqdm([loaddf.load(f, include_syst=False, evtname="evt_%i", hdrname="hdr_%i") for f in IDETVAR_FILES]))
    
    
    Sdf, Spot = loaddf.loadl(SCV_FILES, njob=min(len(SCV_FILES), 10))
    Sdirt, Sdirtpot = loaddf.loadl(SDIRT_FILES, njob=min(len(SDIRT_FILES), 10), include_syst=False)
    Soffbeam, Soffbeampot = loaddf.load(SBEAMOFF_FILE, offbeampot_SBND=True, include_syst=False, load_truth=False)
    Sdetvars, Sdetvar_pots = zip(*tqdm([loaddf.load(f, include_syst=False) for f in SDETVAR_FILES]))
    
    Sdetvars = list(Sdetvars)
    Idetvars = []
    Idetvar_pots = 0
    
    print("SBND CV")
    scale_pot(Sdf, Spot, SGOAL_POT)
    print("SBND Dirt")
    scale_pot(Sdirt, Sdirtpot, SGOAL_POT)
    print("SBND Beam OFF")
    scale_pot(Soffbeam, Soffbeampot, SGOAL_POT)
    for i in range(len(Sdetvars)):
        print("SBND", SDETVAR_NAMES[i])
        scale_pot(Sdetvars[i], Sdetvar_pots[i], SGOAL_POT)
    
    print("ICARUS CV")
    scale_pot(Idf, Ipot, IGOAL_POT)
    print("ICARUS Dirt")
    scale_pot(Idirt, Idirtpot, IGOAL_POT)
    print("ICARUS Beam OFF")
    scale_pot(Ioffbeam, Ioffbeampot, IGOAL_POT)
    for i in range(len(Idetvars)):
        print("ICARUS", IDETVAR_NAMES[i])
        scale_pot(Idetvars[i], Idetvar_pots[i], IGOAL_POT)
    
    bins = np.linspace(0, 60, 21)
    
    _ = plt.hist(Sdf.mu_chi2_of_mu_cand, bins=bins, density=True, histtype="step", linewidth=2, label="CV")
    _ = plt.hist(Sdf.mu_chi2smear_of_mu_cand, bins=bins, density=True, histtype="step", linewidth=2, label="Smeared dE/dx")
    _ = plt.hist(Sdf.mu_chi2hi_of_mu_cand, bins=bins, density=True, histtype="step", linewidth=2, label="Gain Hi")
    add_style(plt.gca(), "$\\chi^2_\\mu$ of Muon Candidate", det="SBND", ylabel="Area Normalized")
    plt.savefig(PLOTDIR+"png/SBNDPID.png")
    plt.clf()   

    bins = np.linspace(0, 60, 21)
    
    _ = plt.hist(Idf.mu_chi2_of_mu_cand, bins=bins, density=True, histtype="step", linewidth=2, label="CV")
    _ = plt.hist(Idf.mu_chi2smear_of_mu_cand, bins=bins, density=True, histtype="step", linewidth=2, label="Smeared dE/dx")
    _ = plt.hist(Idf.mu_chi2hi_of_mu_cand, bins=bins, density=True, histtype="step", linewidth=2, label="Gain High")
    add_style(plt.gca(), "$\\chi^2_\\mu$ of Muon Candidate", det="ICARUS", ylabel="Area Normalized")    
    plt.savefig(PLOTDIR+"png/ICARUSPID.png")
    plt.clf()   
   
    
    Sdetvars += [v_chi2smear(Sdf), v_chi2hi(Sdf)]
    Idetvars += [v_chi2smear(Idf), v_chi2hi(Idf)]
    
    SDETVAR_NAMES += ["Smeared dE/dx", "Gain Hi"]
    IDETVAR_NAMES += ["Smeared dE/dx", "Gain Hi"]
        
    Sdf["selected"] = signalbox(Sdf, "SBND")
    Sdirt["selected"] = signalbox(Sdirt, "SBND")
    Soffbeam["selected"] = signalbox(Soffbeam, "SBND")
    
    for i in range(len(Sdetvars)):
        Sdetvars[i]["selected"] = signalbox(Sdetvars[i], "SBND")
    
    Idf["selected"] = signalbox(Idf, "ICARUS")
    Idirt["selected"] = signalbox(Idirt, "ICARUS") & ICARUS_dircut(Idirt)
    Ioffbeam["selected"] = signalbox(Ioffbeam, "ICARUS")
    
    for i in range(len(Idetvars)):
        Idetvars[i]["selected"] = signalbox(Idetvars[i], "ICARUS")
    
    Sdf["true_E"] = Sdf.true_E.fillna(-1)
    Sdirt["true_E"] = Sdirt.true_E.fillna(-1)
    Soffbeam["true_E"] = -1
    
    for i in range(len(Sdetvars)):
        Sdetvars[i]["true_E"] = Sdetvars[i].true_E.fillna(-1)
    
    Idf["true_E"] = Idf.true_E.fillna(-1)
    Idirt["true_E"] = Idirt.true_E.fillna(-1)
    Ioffbeam["true_E"] = -1
    
    for i in range(len(Idetvars)):
        Idetvars[i]["true_E"] = Idetvars[i].true_E.fillna(-1)
    
    Ssystematics = [
        loaddf.FluxSystematic(Sdf),
        loaddf.XSecSystematic(Sdf),
        syst.SystematicList([syst.SampleSystematic(d) for d in Sdetvars]),
        syst.SystSampleSystematic(Sdirt),
        syst.StatSampleSystematic(Soffbeam, norm=0.1) # TODO: change after unblinding. Simulate scaling up stats by 10x.
    ]
    
    Isystematics = [
        loaddf.FluxSystematic(Idf),
        loaddf.XSecSystematic(Idf),
        syst.SystematicList([syst.SampleSystematic(d) for d in Idetvars]),
        syst.SystSampleSystematic(Idirt),
        syst.StatSampleSystematic(Ioffbeam, norm=0.1) # TODO: change after unblinding. Simulate scaling up stats by 10x.
    ]
    
    systematics = [
        syst.CorrelatedSystematic(loaddf.FluxSystematic(Sdf), loaddf.FluxSystematic(Idf)),
        syst.CorrelatedSystematic(loaddf.XSecSystematic(Sdf), loaddf.XSecSystematic(Idf)),
        syst.UnCorrelatedSystematic(
            syst.SystematicList([syst.SampleSystematic(d) for d in Sdetvars] +\
                                [syst.SystSampleSystematic(Sdirt), syst.StatSampleSystematic(Soffbeam, norm=0.1) # TODO: change after unblinding. Simulate scaling up stats by 10x.
                                ]),
            syst.SystematicList([syst.SampleSystematic(d) for d in Idetvars] +\
                                [syst.SystSampleSystematic(Idirt), syst.StatSampleSystematic(Ioffbeam, norm=0.1) # TODO: change after unblinding. Simulate scaling up stats by 10x.
                                ]))
    ]
    
    
    labels = [
        "Flux",
        "XSec",
        "Detector",
        "Dirt",
        "Beam Off",
        "Stat",
    ]
    
    var = "nu_E_calo"
    wgt = "glob_scale"
    cut = "selected"
    # bins = np.linspace(0, 1.5, 11)[2:]
    bins = np.array([0.0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5])
    centers = (bins[1:] + bins[:-1]) / 2
    
    _ = plt.hist(Sdf.loc[Sdf[cut], var], bins=bins, weights=Sdf.loc[Sdf[cut], wgt], histtype="step", linewidth=2, label="CV")
    _ = plt.hist(Sdirt.loc[Sdirt[cut], var], bins=bins, weights=Sdirt.loc[Sdirt[cut], wgt], histtype="step", linewidth=2, label="Dirt")
    _ = plt.hist(Soffbeam.loc[Soffbeam[cut], var], bins=bins, weights=Soffbeam.loc[Soffbeam[cut], wgt], histtype="step", linewidth=2, label="Beam Off")
    
    add_style(plt.gca(), "Reco. Neutrino Energy [GeV]", det="SBND")
    plt.yscale("log")
    plt.savefig(PLOTDIR+"png/SBND.png")
    plt.clf()

    _ = plt.hist(Idf.loc[Idf[cut], var], bins=bins, weights=Idf.loc[Idf[cut], wgt], histtype="step", linewidth=2, label="CV")
    _ = plt.hist(Idirt.loc[Idirt[cut], var], bins=bins, weights=Idirt.loc[Idirt[cut], wgt], histtype="step", linewidth=2, label="Dirt")
    _ = plt.hist(Ioffbeam.loc[Ioffbeam[cut], var], bins=bins, weights=Ioffbeam.loc[Ioffbeam[cut], wgt], histtype="step", linewidth=2, label="Beam Off")
    
    add_style(plt.gca(), "Reco. Neutrino Energy [GeV]", det="ICARUS")
    plt.yscale("log")
    plt.savefig(PLOTDIR+"png/ICARUS.png")
    plt.clf()
    
    _ = plt.hist(Sdf.loc[Sdf[cut], var], bins=bins, weights=Sdf.loc[Sdf[cut], wgt], histtype="step", linewidth=2, label="CV")
    
    for i in range(20):
        w = Sdf[wgt]*Sdf["flux_univ%i" % i]
        label = "Flux Univ's" if i == 0 else None
        _ = plt.hist(Sdf.loc[Sdf[cut], var], bins=bins, weights=w[Sdf[cut]], histtype="step", linewidth=1, color="gray", label=label)
    
    
    add_style(plt.gca(), "Reco. Neutrino Energy [GeV]", det="SBND")
    plt.savefig(PLOTDIR+"png/FluxSBND.png")
    plt.clf()

    
    _ = plt.hist(Idf.loc[Idf[cut], var], bins=bins, weights=Idf.loc[Idf[cut], wgt], histtype="step", linewidth=2, label="CV")
    
    for i in range(20):
        w = Idf[wgt]*Idf["flux_univ%i" % i]
        label = "Flux Univ's" if i == 0 else None
        _ = plt.hist(Idf.loc[Idf[cut], var], bins=bins, weights=w[Idf[cut]], histtype="step", linewidth=1, color="gray", label=label)
    
    add_style(plt.gca(), "Reco. Neutrino Energy [GeV]", det="ICARUS")
    plt.savefig(PLOTDIR+"png/FluxICARUS.png")
    plt.clf()   
    
    
    _ = plt.hist(Sdf.loc[Sdf[cut], var], bins=bins, weights=Sdf.loc[Sdf[cut], wgt], histtype="step", linewidth=2, label="CV")
    
    for i,s in enumerate(loaddf.xsec_syst):
        w = Sdf[wgt]*Sdf["%s_univ" % s]
        label = "XSec Syst's" if i == 0 else None
        _ = plt.hist(Sdf.loc[Sdf[cut], var], bins=bins, weights=w[Sdf[cut]], histtype="step", linewidth=1, color="gray", label=label)
    
    
    add_style(plt.gca(), "Reco. Neutrino Energy [GeV]", det="SBND")
    plt.savefig(PLOTDIR+"png/SBNDXSecCV.png")
    plt.clf()   
  
    xsec_systs_toplot = loaddf.xsec_syst
    
    CV = np.histogram(Sdf.loc[Sdf[cut], var], bins=bins, weights=Sdf.loc[Sdf[cut], wgt])[0]
    
    vs = []
    ls = []
    
    for i,s in enumerate(xsec_systs_toplot):
        w = Sdf[wgt]*Sdf["%s_univ" % s]
        label = "XSec Syst's Univ" if i == 0 else None
        v = np.histogram(Sdf.loc[Sdf[cut], var], bins=bins, weights=w[Sdf[cut]])[0]
        _ = plt.hist(centers, bins=bins, weights=v/CV, histtype="step", linewidth=2, label=label)
    
    plt.ylim([0.8, 1.2])
    plt.axhline([1], color="gray", linestyle="--")   
    add_style(plt.gca(), "Reco. Neutrino Energy [GeV]", det="SBND")
    plt.savefig(PLOTDIR+"png/SBNDXSecUniv.png")
    plt.clf()   
    
    CV = np.histogram(Idf.loc[Idf[cut], var], bins=bins, weights=Idf.loc[Idf[cut], wgt])[0]
    
    for i,s in enumerate(xsec_systs_toplot):
        w = Idf[wgt]*Idf["%s_univ" % s]
        label = "XSec Syst's Univ" if i == 0 else None
        v = np.histogram(Idf.loc[Idf[cut], var], bins=bins, weights=w[Idf[cut]])[0]
        _ = plt.hist(centers, bins=bins, weights=v/CV, histtype="step", linewidth=2, label=label)
    
    plt.ylim([0.8, 1.2])
    plt.axhline([1], color="gray", linestyle="--")
    plt.savefig(PLOTDIR+"png/ICARUSXSecUniv.png")
    plt.clf()   
    
    _ = plt.hist(Idf.loc[Idf[cut], var], bins=bins, weights=Idf.loc[Idf[cut], wgt], histtype="step", linewidth=2, label="CV")
    
    for i,s in enumerate(loaddf.xsec_syst):
        w = Idf[wgt]*Idf["%s_univ" % s]
        label = "XSec Syts's" if i == 0 else None
        _ = plt.hist(Idf.loc[Idf[cut], var], bins=bins, weights=w[Idf[cut]], histtype="step", linewidth=1, color="gray", label=label)
    
    
    add_style(plt.gca(), "Reco. Neutrino Energy [GeV]", det="ICARUS")
    plt.savefig(PLOTDIR+"png/ICARUSXSecUnivCV.png")
    plt.clf()      
    _ = plt.hist(Sdf.loc[Sdf[cut], var], bins=bins, weights=Sdf.loc[Sdf[cut], wgt], histtype="step", linewidth=2, label="CV")
    
    for i in range(len(Sdetvars)):
        _ = plt.hist(Sdetvars[i].loc[Sdetvars[i][cut], var], bins=bins, weights=Sdetvars[i].loc[Sdetvars[i][cut], wgt],
                histtype="step", linewidth=2, label=SDETVAR_NAMES[i])
    
    add_style(plt.gca(), "Reco. Neutrino Energy [GeV]", det="SBND")
    plt.savefig(PLOTDIR+"png/SBNDDet.png")
    plt.clf()      
   
    _ = plt.hist(Idf.loc[Idf[cut], var], bins=bins, weights=Idf.loc[Idf[cut], wgt], histtype="step", linewidth=2, label="CV")
    
    for i in range(len(Idetvars)):
        _ = plt.hist(Idetvars[i].loc[Idetvars[i][cut], var], bins=bins, weights=Idetvars[i].loc[Idetvars[i][cut], wgt],
                histtype="step", linewidth=2, label=IDETVAR_NAMES[i])
    
    add_style(plt.gca(), "Reco. Neutrino Energy [GeV]", det="ICARUS")
    plt.savefig(PLOTDIR+"png/ICARUSDet.png")
    plt.clf()      
   
    SCV = np.histogram(Sdf.loc[Sdf[cut], var], bins=bins, weights=Sdf.loc[Sdf[cut], wgt])[0]
    Scovs = [s.cov(var, cut, bins, SCV) for s in Ssystematics]
    Scovs.append(np.diag(SCV))
    
    ICV = np.histogram(Idf.loc[Idf[cut], var], bins=bins, weights=Idf.loc[Idf[cut], wgt])[0]
    Icovs = [s.cov(var, cut, bins, ICV) for s in Isystematics]
    Icovs.append(np.diag(ICV))
    
    for c, l in zip(Scovs, labels):
        _ = plt.hist(centers, bins=bins, weights=(np.sqrt(np.diag(c))/SCV), label=l, histtype="step", linewidth=2)
    
    plt.legend()
    add_style(plt.gca(), "Reco. Neutrino Energy [GeV]", det="SBND", ylabel="Fractional Uncertainty")
    plt.savefig(PLOTDIR+"png/SBNDCov.png")
    plt.clf()      
   
    for c, l in zip(Icovs, labels):
        _ = plt.hist(centers[ICV>0], bins=bins, weights=(np.sqrt(np.diag(c))/ICV)[ICV>0], label=l, histtype="step", linewidth=2)
    
    plt.legend()
    add_style(plt.gca(), "Reco. Neutrino Energy [GeV]", det="ICARUS", ylabel="Fractional Uncertainty")
    plt.savefig(PLOTDIR+"png/ICARUSCov.png")
    plt.clf()      
   
    covs = [s.cov(var, cut, bins, np.concatenate((SCV, ICV))) for s in systematics]
    covs.append(np.diag(np.concatenate((SCV, ICV))))
    
    def corr_f(cov):
        err = np.sqrt(np.diag(cov))
        err[err==0] = 1
        err_inv = np.diag(1/err)
        return err_inv@cov@err_inv
    
    plt.title("Flux Systematics")
    plt.imshow(corr_f(covs[0]), origin="lower", cmap="bwr", vmin=-1, vmax=1)
    plt.colorbar()
    plt.savefig(PLOTDIR+"png/FluxCorr.png")
    plt.clf()      
    
    plt.title("XSec Systematics")
    
    plt.imshow(corr_f(covs[1]), origin="lower", cmap="bwr", vmin=-1, vmax=1)
    plt.colorbar()
    plt.savefig(PLOTDIR+"png/XSecCorr.png")
    plt.clf()        

    plt.title("Detector Systematics")
    
    plt.imshow(corr_f(covs[2]), origin="lower", cmap="bwr", vmin=-1, vmax=1)
    plt.colorbar()
    plt.savefig(PLOTDIR+"png/DetCorr.png")
    plt.clf()        
   
    def ratio_cov_full(x, y, cov):
        """
        Covariance of r = x / y given the full covariance of (x, y).
    
        Parameters
        ----------
        x, y : array-like, shape (n,)
            Central values
        cov : array-like, shape (2n, 2n)
            Full covariance matrix of (x, y)
    
            Ordering must be:
            cov = [[Cov(x,x), Cov(x,y)],
                   [Cov(y,x), Cov(y,y)]]
    
        Returns
        -------
        cov_r : ndarray, shape (n, n)
            Covariance matrix of r
        """
        n = len(x)
        assert cov.shape == (2*n, 2*n)
    
        # Protect against division by zero
        eps = 1e-12
        y_safe = np.where(np.abs(y) < eps, eps, y)
    
        Dx = np.diag(1.0 / y_safe)
        Dy = np.diag(-x / y_safe**2)
    
        # Full Jacobian: shape (n, 2n)
        J = np.hstack([Dx, Dy])
    
        return J @ cov @ J.T
    
    ratio = SCV / ICV
    ratio = ratio / ratio.mean()
    
    ratio_cov_flux = ratio_cov_full(SCV, ICV, covs[0])
    ratio_cov_xsec = ratio_cov_full(SCV, ICV, covs[1])
    ratio_cov_det = ratio_cov_full(SCV, ICV, covs[2])
    ratio_cov_stat = ratio_cov_full(SCV, ICV, covs[3])
    ratio_cov_all = ratio_cov_full(SCV, ICV, np.sum(covs, axis=0))
    
    _ = plt.hist(centers, bins=bins, weights=ratio, histtype="step", linewidth=2)
    add_style(plt.gca(), "Reco. Neutrino Energy [GeV]", det="SBND / ICARUS", ylabel="SBND/ICARUS Ratio")
    plt.axhline([1], color="black", linestyle="--")
    plt.savefig(PLOTDIR+"png/Ratio.png")
    plt.clf()        



    _ = plt.hist(centers, bins=bins, weights=np.sqrt(np.diag(ratio_cov_flux))/ratio, histtype="step", linewidth=2, label="Flux Only")
    _ = plt.hist(centers, bins=bins, weights=np.sqrt(np.diag(ratio_cov_xsec))/ratio, histtype="step", linewidth=2, label="XSec Only")
    _ = plt.hist(centers, bins=bins, weights=np.sqrt(np.diag(ratio_cov_det))/ratio, histtype="step", linewidth=2, label="Det. Only")
    _ = plt.hist(centers, bins=bins, weights=np.sqrt(np.diag(ratio_cov_stat))/ratio, histtype="step", linewidth=2, label="Stat. Only")
    _ = plt.hist(centers, bins=bins, weights=np.sqrt(np.diag(ratio_cov_all))/ratio, histtype="step", linewidth=2, label="All")
    # plt.ylim([0, 0.25])
    
    add_style(plt.gca(), "Reco. Neutrino Energy [GeV]", det="SBND / ICARUS", ylabel="Uncertainty on SBND/ICARUS Ratio")
    plt.savefig(PLOTDIR+"png/RatioUncert.png")
    plt.clf()        
   
    reco_bins = bins
    true_bins = np.array([-2] + list(reco_bins) + [np.inf])
    true_bins_centers = (true_bins[:-1] + true_bins[1:])/2
    true_bins_centers, true_bins
    
    SCV2 = np.histogramdd([Sdf.loc[Sdf[cut], "nu_E_calo"], Sdf.loc[Sdf[cut], "true_E"]],
                          bins=[reco_bins, true_bins], weights=Sdf.loc[Sdf[cut], wgt])[0]
    
    ICV2 = np.histogramdd([Idf.loc[Idf[cut], "nu_E_calo"], Idf.loc[Idf[cut], "true_E"]],
                          bins=[reco_bins, true_bins], weights=Idf.loc[Idf[cut], wgt])[0]
    
    SCVT = np.histogram(Sdf.loc[Sdf[cut], "true_E"], bins=true_bins, weights=Sdf.loc[Sdf[cut], wgt])[0]
    ICVT = np.histogram(Idf.loc[Idf[cut], "true_E"], bins=true_bins, weights=Idf.loc[Idf[cut], wgt])[0]
    
    CV2D = np.concatenate((SCV2.T, ICV2.T), axis=0) # .flatten()
    CV2D = np.concatenate((CV2D, CV2D), axis=1)
    
    Imig = ICV2.T / ICV2.sum(axis=1)
    Smig = SCV2.T / SCV2.sum(axis=1)
    
    migration = np.block([
        [Smig, np.zeros(Smig.shape)],
        [np.zeros(Imig.shape), Imig]
    ])
    
    def unfold(C):
        return migration@C@migration.T
    
    ratioT = SCVT / ICVT
    ratioT = ratioT / ratioT.mean()
    
    ratioT_cov_flux = ratio_cov_full(SCVT, ICVT, unfold(covs[0]))
    ratioT_cov_xsec = ratio_cov_full(SCVT, ICVT, unfold(covs[1]))
    ratioT_cov_det = ratio_cov_full(SCVT, ICVT, unfold(covs[2]))
    ratioT_cov_stat = ratio_cov_full(SCVT, ICVT, unfold(covs[3]))
    ratioT_cov_all = ratio_cov_full(SCVT, ICVT, unfold(np.sum(covs, axis=0)))
    
    plot_bins = true_bins[1:-1]
    
    _ = plt.hist(true_bins_centers[1:-1], bins=plot_bins, weights=(np.sqrt(np.diag(ratioT_cov_flux))/ratioT)[1:-1],
                 histtype="step", linewidth=2, label="Flux Only")
    _ = plt.hist(true_bins_centers[1:-1], bins=plot_bins, weights=(np.sqrt(np.diag(ratioT_cov_xsec))/ratioT)[1:-1],
                 histtype="step", linewidth=2, label="XSec Only")
    _ = plt.hist(true_bins_centers[1:-1], bins=plot_bins, weights=(np.sqrt(np.diag(ratioT_cov_det))/ratioT)[1:-1],
                 histtype="step", linewidth=2, label="Det. Only")
    _ = plt.hist(true_bins_centers[1:-1], bins=plot_bins, weights=(np.sqrt(np.diag(ratioT_cov_stat))/ratioT)[1:-1],
                 histtype="step", linewidth=2, label="Stat. Only")
    _ = plt.hist(true_bins_centers[1:-1], bins=plot_bins, weights=(np.sqrt(np.diag(ratioT_cov_all))/ratioT)[1:-1],
                 histtype="step", linewidth=2, label="All")
    # plt.ylim([0, 0.25])
    
    add_style(plt.gca(), "True Neutrino Energy [GeV]", det="SBND / ICARUS", ylabel="Unfolded Unc. on Ratio")
    plt.savefig(PLOTDIR+"png/Unfolded.png")
    plt.clf()        
if __name__ == "__main__":
    main()
