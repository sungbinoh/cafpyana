import pandas as pd
import numpy as np
import h5py
from scipy.interpolate import CubicSpline

from tqdm.auto import tqdm
import pyanalib.pandas_helpers as ph
from multiprocess import Pool
from functools import partial
import syst

import gump_cuts as gc

# Dataframe names
EVT = "mcnu_%i"
WGT = "histpotdf_%i"
HDR = "trig_%i"
MC  = "stub_%i"
CRT = "hdr_%i"

EVT = "evt_%i"
WGT = "wgt_%i"
HDR = "hdr_%i"
MC  = "mcnu_%i"
CRT = "crt_%i"

xsec_syst = [
    # CCQE
    "GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape",
    'GENIEReWeight_SBN_v1_multisigma_CoulombCCQE',

    "ZExpPCAWeighter_SBNNuSyst_multisigma_MvA_ZExp_b1", 
    "ZExpPCAWeighter_SBNNuSyst_multisigma_MvA_ZExp_b2",
    "ZExpPCAWeighter_SBNNuSyst_multisigma_MvA_ZExp_b3",
    "ZExpPCAWeighter_SBNNuSyst_multisigma_MvA_ZExp_b4",
    
    'CCQETemplateReweight_SBNNuSyst_multisigma_SF_q0bin1',
    'CCQETemplateReweight_SBNNuSyst_multisigma_SF_q0bin2',
    'CCQETemplateReweight_SBNNuSyst_multisigma_SF_q0bin3',
    'CCQETemplateReweight_SBNNuSyst_multisigma_SF_q0bin4',
    # 'CCQETemplateReweight_SBNNuSyst_multisigma_SF_q0bin5',

    'CCQETemplateReweight_SBNNuSyst_multisigma_CRPA_q0bin1',
    'CCQETemplateReweight_SBNNuSyst_multisigma_CRPA_q0bin2',
    'CCQETemplateReweight_SBNNuSyst_multisigma_CRPA_q0bin3',
    'CCQETemplateReweight_SBNNuSyst_multisigma_CRPA_q0bin4',
    # 'CCQETemplateReweight_SBNNuSyst_multisigma_CRPA_q0bin5',
    
    'QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_0',
    'QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_1',
    'QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_2',
    'QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_3',
    'QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_4',
    'QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_5',

    # MEC
    'GENIEReWeight_SBN_v1_multisigma_NormCCMEC',
    'GENIEReWeight_SBN_v1_multisigma_NormNCMEC',
    "GENIEReWeight_SBN_v1_multisigma_DecayAngMEC",
    
    'MECq0q3InterpWeighting_SuSAv2ToValenica_q0binned_MECResponse_q0bin0',
    'MECq0q3InterpWeighting_SuSAv2ToValenica_q0binned_MECResponse_q0bin1',
    'MECq0q3InterpWeighting_SuSAv2ToValenica_q0binned_MECResponse_q0bin2',
    'MECq0q3InterpWeighting_SuSAv2ToValenica_q0binned_MECResponse_q0bin3',

    'MECq0q3InterpWeighting_SuSAv2ToMartini_q0binned_MECResponse_q0bin0',
    'MECq0q3InterpWeighting_SuSAv2ToMartini_q0binned_MECResponse_q0bin1',
    'MECq0q3InterpWeighting_SuSAv2ToMartini_q0binned_MECResponse_q0bin2',
    'MECq0q3InterpWeighting_SuSAv2ToMartini_q0binned_MECResponse_q0bin3',

    # RES
    "GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi",
    "GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad",
    "GENIEReWeight_SBN_v1_multisigma_MaCCRES",
    "GENIEReWeight_SBN_v1_multisigma_MaNCRES",
    "GENIEReWeight_SBN_v1_multisigma_MvCCRES",
    "GENIEReWeight_SBN_v1_multisigma_MvNCRES",
    # Non-Res

    # DIS
    # "GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse",
    'GENIEReWeight_SBN_v1_multisigma_AhtBY',
    'GENIEReWeight_SBN_v1_multisigma_BhtBY',
    'GENIEReWeight_SBN_v1_multisigma_CV1uBY',
    'GENIEReWeight_SBN_v1_multisigma_CV2uBY',

    # COH
    "GENIEReWeight_SBN_v1_multisigma_NormCCCOH",
    "GENIEReWeight_SBN_v1_multisigma_NormNCCOH",

    # FSI
    # "GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse",
    # "GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse",
    'GENIEReWeight_SBN_v1_multisigma_MFP_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrCEx_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrInel_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrAbs_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi',
    
    # 'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrG4_N',
    # 'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrINCL_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrG4LoE_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrINCLLoE_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrG4M1E_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrINCLM1E_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrG4M2E_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrINCLM2E_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrG4HiE_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrINCLHiE_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_MFPLoE_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_MFPM1E_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_MFPM2E_N',
    'GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_MFPHiE_N',

    # NCEL
    'GENIEReWeight_SBN_v1_multisigma_MaNCEL',
    'GENIEReWeight_SBN_v1_multisigma_EtaNCEL',
]

xsec_cv_rwgt = [
    "ZExpPCAWeighter_SBNNuSyst_multisigma_MvA_ZExp_b1", 
]

flux_syst = [
 'expskin_Flux',
 'horncurrent_Flux',
 'kminus_Flux',
 'kplus_Flux',
 'kzero_Flux',
 'nucleoninexsec_Flux',
 'nucleonqexsec_Flux',
 'nucleontotxsec_Flux',
 'piminus_Flux',
 'pioninexsec_Flux',
 'pionqexsec_Flux',
 'piontotxsec_Flux',
 'piplus_Flux'
]

truthvars = {
  "true_E": ("nu_E", ""),
  "true_nu_pdg": ("pdg", ""),
  "true_isnc": ("is_nc", ""),
  "genie_mode": ("genie_mode", ""),
  "true_vtx_x": ("pos_x", ""),
  "true_vtx_y": ("pos_y", ""),
  "true_vtx_z": ("pos_z", ""),
}

def scale_pot(df, pot, desired_pot):
    """Scale DataFrame by desired POT."""
    print(f"POT: {pot}\nScaling to: {desired_pot}")
    scale = desired_pot / pot
    df['glob_scale'] = scale * df.cvwgt
    return pot, scale

def load_one(fname, idf, 
    include_syst=True, nuniv=100, spline=False, xsec_univ=False, # systematic handling
    reweight_aFF=False,
    load_truth=True, load_crt=False, match_Enu=True, # load extra information
    offbeampot_SBND=False, offbeampot_ICARUS=False, # POT handling
    preselection=None, # apply preselection cut
    hdrname=HDR, evtname=EVT, wgtname=WGT, mcname=MC, crtname=CRT): # override default table names

    df =  pd.read_hdf(fname, evtname % idf)
    if preselection is not None:
        df = df[preselection(df)]

    hdr = pd.read_hdf(fname, hdrname % idf)

    match = hdr[["run", "evt"]]
    match_ind = list(match.columns)
    # if needed, include neutrino energy in matching information
    if match_Enu:
        mcdf = pd.read_hdf(fname, mcname % idf)
        match = match.merge(mcdf.nu_E.groupby(level=[0,1]).max().rename("nu_E0"), on=["__ntuple", "entry"], how="left")
        match_ind = list(match.columns)

        # Add in other meta-data to match.
        vtx = pd.DataFrame({
          "x": mcdf.pos_x,
          "y": mcdf.pos_y,
          "z": mcdf.pos_z,
        })
        any_in_AV = gc._fv_cut(vtx, "ICARUS", 0, 0, 0, 0).groupby(level=[0,1]).any().rename("AVnu")
        match = match.merge(any_in_AV, on=["__ntuple", "entry"], how="left")

    df = df.merge(match, on=["__ntuple", "entry"], how="left")

    match = match.set_index(match_ind, append=True).droplevel([0,1]).sort_index()

    # LOAD POT
    if offbeampot_SBND:
        N_GATES_ON_PER_5e12POT = 1.05104 # TODO: update
        pot = hdr.noffbeambnb.sum()*N_GATES_ON_PER_5e12POT*5e12
    elif offbeampot_ICARUS:
        trig = pd.read_hdf(fname, "trig_%i" % idf)
        N_GATES_ON_PER_5e12POT = 1.3886218026202426 # TODO: update
        pot = trig.gate_delta.sum()*(1-1/20.)*N_GATES_ON_PER_5e12POT*5e12
    else:
        pot = hdr.pot.sum()

    # LOAD TRUTH
    if load_truth:
        mcdf = pd.read_hdf(fname, mcname % idf)
        mc_tosave = {}
        for setv, load in truthvars.items():
            mc_tosave[setv] = mcdf[load]
        mcdf = pd.DataFrame(mc_tosave, mcdf.index)
        df = df.merge(mcdf, left_on=["__ntuple", "entry", "tmatch_idx"], right_index=True, how="left") 

    # LOAD CRT
    if load_crt:
        if "crthit" in df.columns: del df["crthit"]

        crtdf = pd.read_hdf(fname, crtname % idf)
        crthit = ((crt.time > -1) & (crt.time < 1.8) & (crt.plane != 50)).groupby(level=[0, 1]).any()
        crthit.name = "crthit"
        df = df.join(crthit, on=["__ntuple", "entry"])

    # LOAD AXIAL FORM FACTOR REWEIGHT
    if reweight_aFF:
        rewgt = pd.read_hdf(fname, wgtname % idf)[xsec_cv_rwgt]
        rewgt["cvwgt"] = 1.
        for w in xsec_cv_rwgt:
            rewgt["cvwgt"] = rewgt.cvwgt * rewgt[w]["cv"]
        df = df.merge(rewgt.cvwgt.rename("cvwgt"), left_on=["__ntuple", "entry", "tmatch_idx"], right_index=True, how="left")
        df.cvwgt = df.cvwgt.fillna(1.)
    else:
        df["cvwgt"] = 1.

    # EARLY RETURN IF NOT LOADING WEIGHTS
    if not include_syst:
        return df, match, pot

    # LOAD WEIGHTS
    wgt = pd.read_hdf(fname, wgtname % idf) 
    skim = {}
    for i in range(min(100, nuniv)):
        skim["flux_univ%i" % i] = np.prod([wgt[s]["univ_%i" % i] for s in flux_syst], axis=0)

    if xsec_univ:
        rng = np.random.default_rng(seed=24601) # repeatable random numbers
        rnd = np.clip(rng.normal(size=(len(xsec_syst), nuniv)), -3, 3)
        for i in range(nuniv):
            wgt_vs = []
            for j, s in enumerate(xsec_syst):
                r = rnd[j][i]
        
                if "ps1" in wgt[s]:
                    if spline:
                        w = wgt[s].fillna(1).replace([np.inf, -np.inf], 1)
                        spline_ = CubicSpline([-3, -2, -1, 0, 1, 2, 3], 
                                [w["ms3"]/w["cv"], w["ms2"]/w["cv"], w["ms1"]/w["cv"], pd.Series(1, w.index), w["ps1"]/w["cv"], w["ps2"]/w["cv"], w["ps3"]/w["cv"]])
                        s = spline_(r)
                    else:
                        s = 1 + (wgt[s]["ps1"]/wgt[s]["cv"] - 1)*r
                elif "morph" in wgt[s]:
                    s = 1 + (wgt[s]["morph"] - 1)*np.abs(np.clip(r, -1, 1))
                else:
                    assert(False)

                s = np.clip(s, 0, 10)
                wgt_vs.append(s)
            
            skim["xsec_univ%i" % i] = np.clip(np.prod(wgt_vs, axis=0), 0, 30).fillna(1.)
    else:
        for i, s in enumerate(xsec_syst):
            if "ps1" in wgt[s]:
                skim["%s_univ" % s] = np.clip(wgt[s]["ps1"]/wgt[s]["cv"], 0, 10).fillna(1.)
            elif "morph" in wgt[s]:
                skim["%s_univ" % s] = np.clip(wgt[s]["morph"], 0, 10).fillna(1.)
            else:
                assert(False)

    skim = pd.DataFrame(skim, index=wgt.index)

    mrg = df.merge(skim,
            left_on=["__ntuple", "entry", "tmatch_idx"],
            right_index=True,
            how="left") ## -- save all sllices
    mrg.loc[np.isnan(mrg[skim.columns[0]]), skim.columns] = 1

    return mrg, match, pot


def load(fname, **kwargs):
    with h5py.File(fname, "r") as f:
        ndf = len([k for k in f.keys() if k.startswith("hdr")])

    pots = 0
    dfs = []
    matches = []
    for idf in range(ndf):
        df, match, pot = load_one(fname, idf, **kwargs)
        pots += pot
        dfs.append(df)
        matches.append(match)
    df = pd.concat(dfs).reset_index(drop=True)
    match = pd.concat(matches)

    return df, match, pots
        
def loadl(flist, progress=True, njob=None, **kwargs):
    if njob is not None:
        pool = Pool(njob)
        m = pool.imap_unordered
    else:
        m = map

    # define function w/ kwargs since multiproc doesn't allow for lambdas
    doload_ = partial(load, **kwargs)

    it = m(doload_, flist)

    if progress:
        it = tqdm(it, total=len(flist))

    dfs = []
    matches = []
    pots = 0
    for df, match, pot in it:
        pots += pot
        dfs.append(df)
        matches.append(match)
    df = pd.concat(dfs).reset_index(drop=True)
    matches = pd.concat(matches)

    if njob is not None:
        pool.close()

    return df, matches, pots

def match_common_evts(mrgs, dfs, pots):
    common_ind = mrgs[0].index
    for m in mrgs[1:]:
        common_ind = common_ind.intersection(m.index)

    common_df = pd.DataFrame({"common": 1}, index=common_ind)

    outdfs = []
    outpots = []
    for m, df, p in zip(mrgs, dfs, pots):
        common_frac = common_ind.size / m.index.size
        outpots.append(common_frac*p)
        outdf = df.merge(common_df, left_on=common_ind.names, right_index=True, how="left")
        outdf["common"] = outdf["common"].fillna(0)
        outdf = outdf[outdf.common == 1]
        outdfs.append(outdf)

    return outdfs, outpots

# Systematic class helpers for what is in these files
class FluxSystematic(syst.WeightSystematic):
    def __init__(self, df, scale="glob_scale"):
        wgts = ["flux_univ%i" % i for i in range(100)]
        super().__init__(df, wgts, scale=scale)

class XSecSystematic(syst.WeightSystematic):
    def __init__(self, df, scale="glob_scale"):
        super().__init__(df, ["%s_univ" % s for s in xsec_syst], avg=False, scale=scale)

