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
#import loaddf

RECO = "PANDORA"
FONTSIZE=14

# import importlib
def scale_pot(df, pot, desired_pot):
    """Scale DataFrame by desired POT."""
    print(f"POT: {pot}\nScaling to: {desired_pot}")
    scale = desired_pot / pot
    df['glob_scale'] = scale
    return pot, scale

def add_style(ax, xlabel, title="", det="ICARUS", ylabel='Events / $10^{20}$ POT', legend_loc=None, legend_ncol=1, legend_title=None):
    ax.tick_params(axis='both', which='both', direction='in', length=6, width=1.5, labelsize=FONTSIZE, top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    ax.set_xlabel(xlabel, fontsize=FONTSIZE, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=FONTSIZE, fontweight='bold')
    ax.set_title(f"$\\bf{{{det}}}$  {title}", fontsize=FONTSIZE+2)
    ax.legend(fontsize=FONTSIZE-1, loc=legend_loc, ncol=legend_ncol, title=legend_title, title_fontsize=FONTSIZE)
def v_variation(df, setvars):
    df = df[[c for c in df.columns if "univ" not in c]].copy()
    for (new, old) in setvars:
        df[new] = df[old]
    return df

def v_chi2smear(df):
    setvars = [
        ("mu_chi2_of_mu_cand", "mu_chi2smear13_of_mu_cand"),
        ("mu_chi2_of_p_cand",  "mu_chi2smear13_of_prot_cand"),
        ("prot_chi2_of_mu_cand", "prot_chi2smear13_of_mu_cand"),
        ("prot_chi2_of_p_cand",  "prot_chi2smear13_of_prot_cand"),
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

def FVSBND(df):
    return FV(df, "SBND")

def TrueAV(df, det):
    vtx = pd.DataFrame({'x': df.true_vtx_x,
                       'y': df.true_vtx_y,
                       'z': df.true_vtx_z}, index=df.index)
    return gc._fv_cut(vtx, det, 0, 0, 0, 0)

def TrueAVSBND(df):
    return TrueAV(df, "SBND")

def TrueAVICARUS(df):
    return TrueAV(df, "ICARUS")

def OOAVSBND(df):
    return ~np.isnan(df.true_vtx_x) & ~TrueAVSBND(df)

def OOAVICARUS(df):
    return ~np.isnan(df.true_vtx_x) & ~TrueAVICARUS(df)

def ICARUS_dirtcut(df):
    # mbox from CV sample
    xlo = -378.49
    ylo = -191.86
    zlo = -904.950652270838

    xhi = 378.49
    yhi = 144.96
    zhi = 904.950652270838
    vtx = pd.DataFrame({'x': df.true_vtx_x,
                       'y': df.true_vtx_y,
                       'z': df.true_vtx_z}, index=df.index)

    return ~((vtx.x > xlo) & (vtx.x < xhi) & (vtx.y > ylo) & (vtx.y < yhi) & (vtx.z > zlo) & (vtx.z < zhi))

def main():
    PLOTDIR = "/exp/sbnd/app/users/nrowe/cafpyana/analysis_village/gump/debug/"

    DOSAVE = True
    os.makedirs(PLOTDIR, exist_ok=True)

    FONTSIZE = 14
    HAWKS_COLORS = ["#315031", "#d54c28", "#1e3f54", "#c89648", "#43140b", "#95af8b"]

    #RECO = "PANDORA"

    DF_DIR = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/"
    
    # SCV_FILES = [DF_DIR + "SBND_SpringMC_rewgt_%i.df" % i for i in range(5)]
    SCV_FILES = [DF_DIR + "SBND_SpringMC_rewgt_%i.df" % i for i in range(20)]
    SDIRT_FILE = DF_DIR + "SBND_SpringLowEMC.df"
    #SDIRT_FILES = [DF_DIR + "SBND_SpringLowEMC_rewgt_%i.df" % i for i in range(10)]
    SBEAMOFF_FILE = DF_DIR + "SBND_SpringBNBOffData_5000.df"

    SDETVAR_FILES = [
        DF_DIR + "SBND_SpringMC_Nom.df",
        DF_DIR + "SBND_SpringMC_WMXThetaXW.df",
        DF_DIR + "SBND_SpringMC_WMYZ.df",
        # DF_DIR + "SBND_SpringMC_0xSCE.df",
        DF_DIR + "SBND_SpringMC_2xSCE.df",
    ]
    SDETVAR_NAMES = ["Nominal",
                     "WM $X\\theta_{xw}$", "WM $YZ$", 
                     # "0x SCE", 
                     "2x SCE"]

    SGOAL_POT = 1e20

    Sdf, Smatch, Spot = loaddf.loadl(SCV_FILES, njob=min(len(SCV_FILES), 10), preselection=FVSBND, reweight_aFF=True)
    Sdirt, Sdirtmatch, Sdirtpot = loaddf.load(SDIRT_FILE, preselection=FVSBND, include_syst=False)
    Soffbeam, Soffbeammatch, Soffbeampot = loaddf.load(SBEAMOFF_FILE, preselection=FVSBND, offbeampot_SBND=True, include_syst=False, load_truth=False)
    Sdetvars, Sdetvarsmatch, Sdetvar_pots = zip(*tqdm([loaddf.load(f, preselection=FVSBND, include_syst=False) for f in SDETVAR_FILES]))


    for c in Sdf.columns:
        if "_univ" in c:
            Sdirt[c] = 1

    if "dirt" not in Sdf.columns:
        Sdf["dirt"] = False
        Sdirt["dirt"] = True

    Sdetvars, Sdetvar_pots = loaddf.match_common_evts(Sdetvarsmatch, Sdetvars, Sdetvar_pots)
    Sdetvars = list(Sdetvars)
    
    print("SBND CV")
    scale_pot(Sdf, Spot, SGOAL_POT)
    print("SBND Dirt")
    scale_pot(Sdirt, Sdirtpot, SGOAL_POT)
    print("SBND Beam OFF")
    scale_pot(Soffbeam, Soffbeampot, SGOAL_POT)
    for i in range(len(Sdetvars)):
        loaddf.scale_pot(Sdetvars[i], Sdetvar_pots[i], SGOAL_POT)

    SBND_OOT_COSMIC_SCALE = 0.9961317495469143

    Sdf.loc[Sdf.is_cosmic, "glob_scale"] *= SBND_OOT_COSMIC_SCALE
    Sdirt.loc[Sdirt.is_cosmic, "glob_scale"] *= SBND_OOT_COSMIC_SCALE
    for i in range(len(Sdetvars)):
        Sdetvars[i].loc[Sdetvars[i].is_cosmic, "glob_scale"] *= SBND_OOT_COSMIC_SCALE

    Sdf = pd.concat([Sdf[~Sdf.dirt], Sdirt])
    Schi2_detvars = [v_chi2smear(Sdf), v_chi2hi(Sdf)]

    #bins = np.array([0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.25, 1.5, 3.0]) 
    bins = np.array([0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.25, 1.5]) 
    Sdetvars += [v_chi2smear(Sdf), v_chi2hi(Sdf)]
    
    #SDETVAR_NAMES += ["Smeared dE/dx", "Gain Hi"]
    Sdf["selected"] = signalbox(Sdf, "SBND")
    Sdirt["selected"] = signalbox(Sdirt, "SBND")
    Soffbeam["selected"] = signalbox(Soffbeam, "SBND")
   
    for i in range(len(Sdetvars)):
        Sdetvars[i]["selected"] = signalbox(Sdetvars[i], "SBND")

    for i in range(len(Schi2_detvars)):
        Schi2_detvars[i]["selected"] = signalbox(Schi2_detvars[i], "SBND")
   
    Sdf["true_E"] = Sdf.true_E.fillna(-1)
    Sdirt["true_E"] = Sdirt.true_E.fillna(-1)
    Soffbeam["true_E"] = -1
    
    for i in range(len(Sdetvars)):
        Sdetvars[i]["true_E"] = Sdetvars[i].true_E.fillna(-1)
   
    Ssystematics = [
        loaddf.FluxSystematic(Sdf),
        loaddf.XSecSystematic(Sdf),
        syst.NormalizationSystematic(0.02),
        syst.SystematicList(
            [syst.SampleSystematic(d, cvdf=Sdetvars[0]) for d in Sdetvars[1:]] + # detector variation samples
            [syst.SampleSystematic(d) for d in Schi2_detvars] # chi2 variations
        ),
        syst.SystSampleSystematic(Sdf[OOAVSBND(Sdf)]),
        syst.StatSampleSystematic(Soffbeam, norm=0.1) # TODO: change after unblinding. Simulate scaling up stats by 10x.
    ]
  
    labels = [
        "Flux",
        "XSec",
        "POT Norm.",
        "Detector",
        "Dirt",
        "Beam Off",
        "Stat",
    ]

    var = "nu_E_calo"
    wgt = "glob_scale"
    cut = "selected"
    # bins = np.linspace(0, 1.5, 11)[2:]
    centers = (bins[1:] + bins[:-1]) / 2


    SCV = np.histogram(Sdf.loc[Sdf[cut], var], bins=bins, weights=Sdf.loc[Sdf[cut], wgt])[0]
    Scovs = [s.cov(var, cut, bins, SCV) for s in Ssystematics]
    Scovs.append(np.diag(SCV))


    S_detsyst_covs = [s.cov(var, cut, bins, SCV) for s in Ssystematics[labels.index("Detector")].systs]

    Schi2_detvars = [v_chi2smear(Sdf), v_chi2hi(Sdf)]
    SCHI2_DETVAR_NAMES = ["Smeared dE/dx", "Gain Hi"]

    combined = np.zeros(Scovs[0].shape)
    for c, l in zip(S_detsyst_covs, SDETVAR_NAMES[1:] + SCHI2_DETVAR_NAMES):
        _ = plt.hist(centers, bins=bins, weights=(np.sqrt(np.diag(c))/SCV), label=l, histtype="step", linewidth=2)
        combined += c
    
    for c, l in zip(Scovs[4:], labels[4:]):
        print(f"second loop {l}")
        #_ = plt.hist(centers, bins=bins, weights=(np.sqrt(np.diag(c))/SCV), label=l, histtype="step", linewidth=2)
        #combined += c

    plt.hist(centers, bins=bins, weights=(np.sqrt(np.diag(combined))/SCV), label="Total", histtype="step", color="black", linewidth=2, linestyle="--")

    add_style(plt.gca(), "Reco. Neutrino Energy [GeV]", det="SBND", ylabel="Fractional Uncertainty", 
              legend_loc="upper center", legend_ncol=2, legend_title="Uncorrelated Uncertainties")

    plt.ylim([0, 0.125])
    plt.savefig(PLOTDIR+"SBNDUncorr.png")
    plt.clf()

    for i in range(len(Sdetvars)-2):
        n, _, _ = plt.hist(Sdetvars[i].loc[Sdetvars[i][cut], var], bins=bins, weights=Sdetvars[i].loc[Sdetvars[i][cut], wgt], 
                histtype="step", linewidth=2, label=SDETVAR_NAMES[i])
 
    add_style(plt.gca(), "Reco. Neutrino Energy [GeV]", det="SBND")
    plt.savefig(PLOTDIR+"SBNDHistUncorr.png")
    plt.clf()


if __name__ == "__main__":
    main()
