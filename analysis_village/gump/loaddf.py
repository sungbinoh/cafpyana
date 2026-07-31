import os
import hashlib
import json

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

def tmatch(reco, mc):
    for c in mc.columns:
        if c in reco.columns:
            print(f'duplicate column found! {c}, setting right col to *_true.')
            if(isinstance(c, tuple)):
                mc.rename(columns={c:(c[0]+'_true', '')}, inplace=True)
            else:
                mc.rename(columns={c:c+'_true'}, inplace=True)

    df = ph.multicol_merge(reco.reset_index(), mc.reset_index(),
                           left_on=[("__ntuple", ""), ("entry", ""), ("tmatch_idx", "")],
                           right_on=[("__ntuple", ""), ("entry", ""), ("rec.mc.nu..index", "")],
                           how="left") # start with keeping everything...
    return df

# Dataframe names
EVT = "evt_%i"
WGT = "wgt_%i"
HDR = "hdr_%i"
MC  = "mcnu_%i"
CRT = "crt_%i"
FLASH = "flash_%i"

pot_syst = {'ms3': 0.982714, 'ms2': 0.9887274, 'ms1': 0.99474195, 'cv': 1.0, 'ps1': 1.005, 'ps2': 1.01, 'ps3': 1.015}

xsec_syst = [
    # CCQE
    "GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape",
    'GENIEReWeight_SBN_v1_multisigma_CoulombCCQE',

    # MEC
    'GENIEReWeight_SBN_v1_multisigma_NormCCMEC',
    'GENIEReWeight_SBN_v1_multisigma_NormNCMEC',
    "GENIEReWeight_SBN_v1_multisigma_DecayAngMEC",

    # RES
    "GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi",
    "GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad",
    "GENIEReWeight_SBN_v1_multisigma_MaCCRES",
    "GENIEReWeight_SBN_v1_multisigma_MaNCRES",
    "GENIEReWeight_SBN_v1_multisigma_MvCCRES",
    "GENIEReWeight_SBN_v1_multisigma_MvNCRES",
    "GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma",
    "GENIEReWeight_SBN_v1_multisigma_RDecBR1eta",

    # Non-Res
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi',
    'GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi',

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
    'GENIEReWeight_SBN_v1_multisigma_MFP_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrCEx_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrInel_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrAbs_pi',
    'GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi',

    # NCEL
    'GENIEReWeight_SBN_v1_multisigma_MaNCEL',
    'GENIEReWeight_SBN_v1_multisigma_EtaNCEL',

    # Systematics introduced by Ar23+
    "CCQETemplateReweight_SBN_v3_LFGToSF_q0bin0",
    "CCQETemplateReweight_SBN_v3_LFGToSF_q0bin1",
    "CCQETemplateReweight_SBN_v3_LFGToSF_q0bin2",
    "CCQETemplateReweight_SBN_v3_LFGToSF_q0bin3",
    "CCQETemplateReweight_SBN_v3_LFGToSF_q0bin4",

    "CCQETemplateReweight_SBN_v3_LFGToHF_q0bin0",
    "CCQETemplateReweight_SBN_v3_LFGToHF_q0bin1",
    "CCQETemplateReweight_SBN_v3_LFGToHF_q0bin2",
    "CCQETemplateReweight_SBN_v3_LFGToHF_q0bin3",
    "CCQETemplateReweight_SBN_v3_LFGToHF_q0bin4",

    "CCQETemplateReweight_SBN_v3_HFToCRPA_q0bin0",
    "CCQETemplateReweight_SBN_v3_HFToCRPA_q0bin1",
    "CCQETemplateReweight_SBN_v3_HFToCRPA_q0bin2",
    "CCQETemplateReweight_SBN_v3_HFToCRPA_q0bin3",
    "CCQETemplateReweight_SBN_v3_HFToCRPA_q0bin4",

    "QEInterference_SBN_v3_QEIntf_dial_0",
    "QEInterference_SBN_v3_QEIntf_dial_1",
    "QEInterference_SBN_v3_QEIntf_dial_2",
    "QEInterference_SBN_v3_QEIntf_dial_3",
    "QEInterference_SBN_v3_QEIntf_dial_4",
    "QEInterference_SBN_v3_QEIntf_dial_5",

    "GENIEReWeight_SBN_v3_FrG4LoE_N",
    "GENIEReWeight_SBN_v3_FrG4M1E_N",
    "GENIEReWeight_SBN_v3_FrG4M2E_N",
    "GENIEReWeight_SBN_v3_FrG4HiE_N",
    "GENIEReWeight_SBN_v3_FrINCLLoE_N",
    "GENIEReWeight_SBN_v3_FrINCLM1E_N",
    "GENIEReWeight_SBN_v3_FrINCLM2E_N",
    "GENIEReWeight_SBN_v3_FrINCLHiE_N",
    "GENIEReWeight_SBN_v3_MFPLoE_N",
    "GENIEReWeight_SBN_v3_MFPM1E_N",
    "GENIEReWeight_SBN_v3_MFPM2E_N",
    "GENIEReWeight_SBN_v3_MFPHiE_N",

    "ZExpPCAWeighter_SBN_v3_MvA_b1",
    "ZExpPCAWeighter_SBN_v3_MvA_b2",
    "ZExpPCAWeighter_SBN_v3_MvA_b3",
    "ZExpPCAWeighter_SBN_v3_MvA_b4",

    "MECq0q3InterpWeighting_SBN_v3_SuSAToVal_MECResponse_q0bin0",
    "MECq0q3InterpWeighting_SBN_v3_SuSAToVal_MECResponse_q0bin1",
    "MECq0q3InterpWeighting_SBN_v3_SuSAToVal_MECResponse_q0bin2",
    "MECq0q3InterpWeighting_SBN_v3_SuSAToVal_MECResponse_q0bin3",
    "MECq0q3InterpWeighting_SBN_v3_SuSAToMar_MECResponse_q0bin0",
    "MECq0q3InterpWeighting_SBN_v3_SuSAToMar_MECResponse_q0bin1",
    "MECq0q3InterpWeighting_SBN_v3_SuSAToMar_MECResponse_q0bin2",
    "MECq0q3InterpWeighting_SBN_v3_SuSAToMar_MECResponse_q0bin3",
]

xsec_cv_rwgt = [
    "ZExpPCAWeighter_SBN_v3_MvA_b1",
    "CCQEXSecCorr_SBN_v3_CCQEXSecCorr",
    "GENIEReWeight_SBN_v3_FrKin_PiProFix_N",
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

g4_syst = [
 'reinteractions_piminus_Geant4',
 'reinteractions_piplus_Geant4',
 'reinteractions_proton_Geant4',
]

truthvars = {
  "true_E": ("nu_E", ""),
  "true_nu_pdg": ("pdg", ""),
  "true_issig": ("is_sig", ""),
  "true_isothernumucc": ("is_other_numucc", ""),
  "true_isfv": ("is_fv", ""),
  "true_isnc": ("is_nc", ""),
  "genie_mode": ("genie_mode", ""),
  "true_vtx_x": ("pos_x", ""),
  "true_vtx_y": ("pos_y", ""),
  "true_vtx_z": ("pos_z", ""),
  "true_nmu": ("nmu", ""),
  "true_np": ("np", ""),
  "true_nn": ("nn", ""),
  "true_npi": ("npi", ""),
  "true_npi0": ("npi0", ""),
}

def scale_pot(df, pot, desired_pot):
    """Scale DataFrame by desired POT."""
    scale = desired_pot / pot
    df['glob_scale'] = scale * df.cvwgt
    return pot, scale

def _cache_key(fname, idf, **kwargs):
    """Build a deterministic hash from the input file path, split index, and all keyword args."""
    key_dict = {"fname": os.path.abspath(fname), "idf": idf, "_cache_version": _CACHE_VERSION}
    # Include the input file's identity, not just its path, so regenerating a .df in
    # place busts the cache instead of silently serving the old df *and* the old POT.
    # NB: this also means re-syncing the .df directory (which rewrites mtimes) forces
    # a full re-load even if the contents are unchanged.
    st = os.stat(fname)
    key_dict["_fsize"] = st.st_size
    key_dict["_fmtime"] = int(st.st_mtime)
    # Only include serializable kwargs (skip preselection function)
    for k, v in sorted(kwargs.items()):
        if callable(v):
            # Use the function's qualified name so different preselections bust the cache
            key_dict[k] = v.__module__ + "." + v.__qualname__
        else:
            key_dict[k] = v
    raw = json.dumps(key_dict, sort_keys=True, default=str)
    return hashlib.sha256(raw.encode()).hexdigest()[:16]

def _write_cache(cache_file, df, match, pot):
    """Write load_one output to an HDF5 cache file."""
    os.makedirs(os.path.dirname(cache_file), exist_ok=True)
    df.to_hdf(cache_file, "df", mode="w")
    match.to_hdf(cache_file, "match", mode="a")
    with h5py.File(cache_file, "a") as cf:
        cf.attrs["pot"] = pot

#@profile
def load_one(fname, idf,
    detector=None, # One of SBND, ICARUS, ICARUS Run4
    include_syst=True, nuniv=100, spline=False, xsec_univ=False, xsec_spline=False,# systematic handling
    reweight_aFF=False, pot_univ=False, flux_univ=True, sep_flux_univ=False, g4_univ=True,
    pot_spline=False, detvar_spline=False,
    load_truth=True, load_crt=False, match_Enu=True, # load extra information
    offbeampot=False, # POT handling
    preselection=None, # apply preselection cut
    cache_dir=None, # directory to cache output; None disables caching
    flashname=FLASH, hdrname=HDR, evtname=EVT, wgtname=WGT, mcname=MC, crtname=CRT, drops=None, lightmem=False): # override default table names

    assert(detector == "SBND" or detector == "ICARUS Run2" or detector == "ICARUS Run4")

    # Check cache
    if cache_dir is not None:
        cache_hash = _cache_key(fname, idf, detector=detector, include_syst=include_syst,
            nuniv=nuniv, spline=spline, xsec_univ=xsec_univ, xsec_spline=xsec_spline, reweight_aFF=reweight_aFF, pot_univ=pot_univ,
            flux_univ=flux_univ, sep_flux_univ=sep_flux_univ, g4_univ=g4_univ,
            load_truth=load_truth, load_crt=load_crt,
            match_Enu=match_Enu, offbeampot=offbeampot, preselection=preselection,
            drops=drops, lightmem=lightmem,
            flashname=flashname, hdrname=hdrname, evtname=evtname,
            wgtname=wgtname, mcname=mcname, crtname=crtname)
        cache_file = os.path.join(cache_dir, cache_hash + ".h5")
        if os.path.exists(cache_file):
            try:
                df = pd.read_hdf(cache_file, "df")
                match = pd.read_hdf(cache_file, "match")
            except Exception as err:
                print(fname, cache_file)
                raise err 
            with h5py.File(cache_file, "r") as cf:
                pot = float(cf.attrs["pot"])
            return df, match, pot

    df =  pd.read_hdf(fname, evtname % idf)
    hdr = pd.read_hdf(fname, hdrname % idf)
    ismc = hdr.ismc.iloc[0] == 1

    # set run 
    if "SBND" in fname:
        df["Run"] = 1
        Run = 1
    elif "ICARUS" in fname and "Run4" in fname:
        df["Run"] = 4
        Run = 4
    elif "ICARUS" in fname:
        df["Run"] = 2
        Run = 2
    else: assert(False)

    # apply the scaled pe flash
    if ismc: # Scale PE for MC-only
        # Best-fit MC PE scale factors from the data/MC fits in FlashMCDataComparison.ipynb
        if detector == "SBND": pe_scale = 0.642
        elif detector == "ICARUS Run2": pe_scale = 0.632
        elif detector == "ICARUS Run4": pe_scale = 0.358
    else:
        pe_scale = 1.0

    if "ICARUS" in detector:
        df["flash_maxpe"] = df["flash_maxpe_cryo0"] * pe_scale
        df.loc[df.slc_vtx_x > 0, "flash_maxpe"] = df["flash_maxpe_cryo1"] * pe_scale
    else:
        df["flash_maxpe"] = df["flash_maxpe"] * pe_scale
    df["flash_maxpe"] = df["flash_maxpe"].fillna(0.).astype(float)

    # Apply preselection
    if preselection is not None:
        df = df[preselection(df)]

    match = hdr[["run", "evt"]]
    # The columns that identify an event across samples. Everything merged onto `match`
    # after this point is metadata carried along as a column, NOT part of the key --
    # in particular AVnu, which is derived from gc._fv_cut and so depends on `detector`.
    match_ind = list(match.columns)
    # if needed, include neutrino energy in matching information
    if match_Enu:
        mcdf = pd.read_hdf(fname, mcname % idf)
        match = match.merge(mcdf.nu_E.groupby(level=[0,1]).max().rename("nu_E0"), on=["__ntuple", "entry"], how="left")
        match_ind = list(match.columns)

        # Add in other meta-data to match.
        vtx = pd.DataFrame({
          "detector": detector,
          "Run": Run,
          "x": mcdf.pos_x,
          "y": mcdf.pos_y,
          "z": mcdf.pos_z,
        })
        any_in_AV = gc._fv_cut(vtx, 0, 0, 0, 0).groupby(level=[0,1]).any().rename("AVnu")
        match = match.merge(any_in_AV, on=["__ntuple", "entry"], how="left")

    df = df.merge(match, on=["__ntuple", "entry"], how="left")

    # DROP DUPLICATED EVENTS
    # A "duplicate" is the same physical event appearing in more than one
    # (__ntuple, entry) row of the header — i.e. the same event reconstructed twice.
    # SBND MC legitimately reuses (run, evt) across distinct MC events, so when
    # match_Enu is True we include nu_E0 in the dedup key as a tie-breaker.
    # Drop ALL occurrences (keep=False), not just the extras, from both match and df.
    # Use a MultiIndex.isin mask on df rather than df.merge, so df's existing
    # MultiIndex (__ntuple, entry, rec.slc..index) is preserved.
    dedup_cols = ["run", "evt", "nu_E0"] if match_Enu else ["run", "evt"]
    dup_mask_match = match.duplicated(subset=dedup_cols, keep=False)
    n_dup_pairs = int(match.loc[dup_mask_match, dedup_cols].drop_duplicates().shape[0])
    n_dup_rows  = int(dup_mask_match.sum())
    if n_dup_rows > 0:
        bad_pairs = pd.MultiIndex.from_frame(
            match.loc[dup_mask_match, dedup_cols].drop_duplicates())
        df_pairs = pd.MultiIndex.from_arrays([df[c] for c in dedup_cols])
        df = df[~df_pairs.isin(bad_pairs)]
        match = match[~dup_mask_match]
    print(f"[{os.path.basename(fname)} idf={idf}] dedup: dropped "
          f"{n_dup_pairs} duplicated {tuple(dedup_cols)} keys ({n_dup_rows} hdr rows)")

    match = match.set_index(match_ind, append=True).droplevel([0,1]).sort_index()

    # LOAD POT
    if offbeampot:
        if detector == "SBND":
            N_GATES_ON_PER_5e12POT = 1.05104
            pot = hdr.noffbeambnb.sum()/N_GATES_ON_PER_5e12POT*5e12
        elif detector == "ICARUS Run4":
            trig = pd.read_hdf(fname, "trig_%i" % idf)
            N_GATES_ON_PER_5e12POT = 1.0631936867739828
            pot = trig.gate_delta.sum()*(1-1/20.)/N_GATES_ON_PER_5e12POT*5e12
        elif detector == "ICARUS Run2":
            trig = pd.read_hdf(fname, "trig_%i" % idf)
            N_GATES_ON_PER_5e12POT = 1.3886218026202426
            pot = trig.gate_delta.sum()*(1-1/20.)/N_GATES_ON_PER_5e12POT*5e12
    else:
        pot = hdr.pot.sum()

    # CORRECT POT FOR THE DEDUP
    # The dedup above dropped events from `match`/`df`, but `hdr` (and the ICARUS
    # trigger table) still describe every event, so the POT just computed covers
    # events that are no longer in the frame. Scale it by the surviving fraction.
    # This treats POT as uniform per event -- the same assumption `match_common_evts`
    # makes -- which is the best available: for MC the POT sits only on the
    # first_in_subrun records, so filtering `hdr` directly would charge the full
    # subrun POT to whichever record happened to be dropped.
    if n_dup_rows > 0:
        pot *= 1. - n_dup_rows / len(hdr)

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
        crthit = ((crtdf.time > -1) & (crtdf.time < 1.8) & (crtdf.plane != 50)).groupby(level=[0, 1]).any()
        crthit.name = "crthit"
        df = df.join(crthit, on=["__ntuple", "entry"])

    df["crthit"] = df.crthit.fillna(False).astype(bool) 

    # LOAD AXIAL FORM FACTOR REWEIGHT
    if reweight_aFF:
        rewgt = pd.read_hdf(fname, wgtname % idf)[xsec_cv_rwgt]
        rewgt["cvwgt"] = 1.
        for w in xsec_cv_rwgt:
            cvcol = "cv" if "cv" in rewgt[w].columns else "morph"
            rewgt["cvwgt"] = rewgt.cvwgt * rewgt[w][cvcol]
        df = df.merge(rewgt.cvwgt.rename("cvwgt"), left_on=["__ntuple", "entry", "tmatch_idx"], right_index=True, how="left")
        df.cvwgt = df.cvwgt.fillna(1.)
    else:
        df["cvwgt"] = 1.

    if drops is not None:
        df.drop(columns=drops, inplace=True, errors='ignore')

    if lightmem:
        type_map = {
            'detector': 'category',
            'Run': 'category',
            'true_isfv': 'Int8',
            'true_isothernumucc': 'Int8',
            'true_issig': 'Int8',
            'true_isnc': 'Int8'
        }
        
        valid_type_map = {col: dtype for col, dtype in type_map.items() if col in df.columns}
        
        df = df.astype(valid_type_map)
        df[df.select_dtypes(include=['float64']).columns] = df.select_dtypes(include=['float64']).astype('float32')

    # EARLY RETURN IF NOT LOADING WEIGHTS
    if not include_syst:
        if cache_dir is not None:
            _write_cache(cache_file, df, match, pot)
        return df, match, pot

    # LOAD WEIGHTS
    wgt = pd.read_hdf(fname, wgtname % idf) 
    skim = {}
    if flux_univ:
        for i in range(min(100, nuniv)):
            skim["flux_univ%i" % i] = np.prod([wgt[s]["univ_%i" % i] for s in flux_syst], axis=0)

    if g4_univ:
        for i in range(min(100, nuniv)):
            skim["g4_univ%i" % i] = np.prod([wgt[s]["univ_%i" % i] for s in g4_syst], axis=0)

    if pot_univ:
        rng = np.random.default_rng(seed=24601) # repeatable random numbers
        rnd = np.clip(rng.normal(size=nuniv), -3, 3)
        for i in range(nuniv):
            wgt_vs = []
            r = rnd[i]
        
            if "ps1" in pot_syst:
                if spline:
                    w = pot_syst
                    spline_ = CubicSpline([-3, -2, -1, 0, 1, 2, 3], 
                            [w["ms3"]/w["cv"], w["ms2"]/w["cv"], w["ms1"]/w["cv"], pd.Series(1, w.index), w["ps1"]/w["cv"], w["ps2"]/w["cv"], w["ps3"]/w["cv"]])
                    s = spline_(r)
                else:
                    s = 1 + (pot_syst["ps1"]/pot_syst["cv"] - 1)*r
            else:
                assert(False)

            wgt_vs.append(s)
            
            skim["pot_univ%i" % i] = np.prod(wgt_vs, axis=0)
    else:
        if "ps1" in pot_syst:
            skim["pot_univ"] = pot_syst["ps1"]/pot_syst["cv"]
        else:
            assert(False)

    multisim_cols = []
    multisigma_cols = []

    if sep_flux_univ:
        for j, s in enumerate(flux_syst):
            multisim_cols.append(s)
            w = wgt[s]#.fillna(1).replace([np.inf, -np.inf], 1)
            if lightmem:
                w[w.select_dtypes(include=["float64"]).columns] = w.select_dtypes(include=["float64"]).astype("float32")
            stacked_variants = np.vstack([np.nan_to_num(w["univ_%i" % i].to_numpy(), nan=1.0, posinf=1.0, neginf=1.0) for i in range(min(100, nuniv))])
            skim[s] = stacked_variants.T.tolist()
            for d in stacked_variants.T.tolist():
                if len(d) != 100:
                    print(d)

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

    if xsec_spline:
        for j, s in enumerate(xsec_syst):
            multisigma_cols.append(s)
            if "ps1" in wgt[s]:
                w = wgt[s].fillna(1).replace([np.inf, -np.inf], 1)
                stacked_variants = np.vstack([
                    np.clip((w["ms3"] / w["cv"]).to_numpy(), 0, 10),
                    np.clip((w["ms2"] / w["cv"]).to_numpy(), 0, 10),
                    np.clip((w["ms1"] / w["cv"]).to_numpy(), 0, 10),
                    np.ones(len(w)),  # Central value ratio is exactly 1.0
                    np.clip((w["ps1"] / w["cv"]).to_numpy(), 0, 10),
                    np.clip((w["ps2"] / w["cv"]).to_numpy(), 0, 10),
                    np.clip((w["ps3"] / w["cv"]).to_numpy(), 0, 10)
                ])

                # 2. Transpose to shape (n_events, 7) so each row represents an event,
                # then convert to a list of lists for uproot/awkward ingestion later
                skim[s] = stacked_variants.T.tolist()
            elif "morph" in wgt[s]:
                w = wgt[s].fillna(1).replace([np.inf, -np.inf], 1)
                if lightmem:
                    w[w.select_dtypes(include=["float64"]).columns] = w.select_dtypes(include=["float64"]).astype("float32")

                stacked_variants = np.vstack([
                    np.ones(len(w)),  # Central value ratio is exactly 1.0
                    np.clip((w["morph"]).to_numpy(), 0, 10)
                ])

                # 2. Transpose to shape (n_events, 7) so each row represents an event,
                # then convert to a list of lists for uproot/awkward ingestion later
                skim[s] = stacked_variants.T.tolist()

    else:
        for i, s in enumerate(xsec_syst):
            if "ps1" in wgt[s]:
                skim["%s_univ" % s] = np.clip(wgt[s]["ps1"]/wgt[s]["cv"], 0, 10).fillna(1.)
            elif "morph" in wgt[s]:
                skim["%s_univ" % s] = np.clip(wgt[s]["morph"], 0, 10).fillna(1.)
            elif "univ_0" in wgt[s]:
                skim["%s_univ" % s] = pd.Series(1 + np.sqrt(np.mean([(1 - wgt[s][c].clip(0, 10.))**2 for c in wgt[s].columns], axis=0)), index=wgt.index).fillna(1.)
            else:
                assert(False)

    skim = pd.DataFrame(skim, index=wgt.index)


    mrg = df.merge(skim,
            left_on=["__ntuple", "entry", "tmatch_idx"],
            right_index=True,
            how="left") ## -- save all sllices

    univ_cols = [col for col in skim.columns if "univ" in col]
    if len(multisigma_cols) > 0:
        nan_mask = mrg[multisigma_cols[0]].isna()
        for col in multisigma_cols:
            mrg.loc[nan_mask, col] = mrg.loc[nan_mask, col].apply(lambda x: [1.0] * len(mrg.loc[~nan_mask, col].iloc[0]))

    if len(multisim_cols) > 0:
        nan_mask = mrg[multisim_cols[0]].isna()
        for col in multisim_cols:
            mrg.loc[nan_mask, col] = mrg.loc[nan_mask, col].apply(lambda x: [1.0] * 100)

    if len(univ_cols) > 0:
        mrg.loc[np.isnan(mrg[univ_cols[0]]), univ_cols] = 1.0 

    if drops is not None:
        mrg.drop(columns=drops, inplace=True, errors='ignore')
    if cache_dir is not None:
        _write_cache(cache_file, mrg, match, pot)

    return mrg, match, pot


def load(fname, maxdf=None, **kwargs):
    with h5py.File(fname, "r") as f:
        ndf = len([k for k in f.keys() if k.startswith("hdr")])

    if maxdf is None:
        maxdf = ndf

    pots = 0
    dfs = []
    matches = []
    for idf in range(min(ndf, maxdf)):
        df, match, pot = load_one(fname, idf, **kwargs)
        pots += pot
        dfs.append(df)
        matches.append(match)
    df = pd.concat(dfs).reset_index(drop=True)
    match = pd.concat(matches)
    n_match_before = len(match)

    # CROSS-IDF DEDUP
    # `load_one` only sees one idf (split) at a time. The same physical event can
    # show up in more than one idf — `__ntuple` is a per-idf ordinal, not globally
    # unique, so the within-idf check can't catch this. Drop every occurrence of
    # any duplicate after the concat across idfs. Match nu_E0 in the key when it's
    # present (match_Enu=True), since SBND MC reuses (run, evt) across distinct
    # MC events and would over-drop on (run, evt) alone.
    dedup_levels = ["run", "evt"]
    if "nu_E0" in match.index.names:
        dedup_levels.append("nu_E0")
    key = pd.MultiIndex.from_arrays([
        match.index.get_level_values(name) for name in dedup_levels
    ])
    dup_mask = key.duplicated(keep=False)
    if dup_mask.any():
        bad_pairs = pd.MultiIndex.from_arrays([
            key.get_level_values(i)[dup_mask] for i in range(len(dedup_levels))
        ]).unique()
        n_dup_pairs = len(bad_pairs)
        n_dup_rows  = int(dup_mask.sum())
        match = match[~dup_mask]
        df_pairs = pd.MultiIndex.from_arrays([df[name] for name in dedup_levels])
        df = df[~df_pairs.isin(bad_pairs)]
        # As in load_one: the events are gone, so their POT must go with them.
        pots *= 1. - n_dup_rows / n_match_before
        print(f"[{os.path.basename(fname)}] cross-idf dedup: dropped "
              f"{n_dup_pairs} duplicated {tuple(dedup_levels)} keys "
              f"({n_dup_rows} match rows)")

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
    df = pd.concat(dfs, ignore_index=True)
    del dfs
    matches = pd.concat(matches)

    if njob is not None:
        pool.close()

    return df, matches, pots

def match_common_evts(mrgs, dfs, pots):
    """Restrict every sample to the events common to all of them.

    All output frames hold the *same* set of physical events, so they are all given
    the *same* POT, taken from the group's nominal (index 0 by convention -- see
    SDETVARS in the signal-box notebooks and DETVAR_FILES in mcdata_comparison.py).

    Deriving the POT per member instead (`common_frac_i * pots[i]`, what this used to
    do) is only equivalent when every sample has the same number of CAF records per
    generated POT. When it does not -- e.g. a production that lost event records but
    still books the full POT of those jobs -- the CV and the variation end up with
    different normalizations for identical events, and that shows up downstream as a
    flat normalization detector systematic that is pure bookkeeping.
    """
    common_ind = mrgs[0].index
    for m in mrgs[1:]:
        common_ind = common_ind.intersection(m.index)

    common_df = pd.DataFrame({"common": 1}, index=common_ind)

    # NB: isin(), not common_ind.size -- Index.intersection returns unique values, so
    # against a match index with repeated keys the size ratio understates the fraction
    # of rows actually kept.
    common_frac = float(mrgs[0].index.isin(common_ind).mean())
    pot_common = common_frac * pots[0]

    # The members should agree on events-per-POT: they are the same generated events
    # reconstructed differently. If they do not, one of the samples is missing event
    # records relative to its own POT bookkeeping, and the matched result is only as
    # trustworthy as that sample.
    rates = [m.index.size / p for m, p in zip(mrgs, pots) if p > 0]
    if len(rates) > 1 and max(rates) / min(rates) - 1 > 0.02:
        print("WARNING: match_common_evts members disagree on events-per-POT by %.1f%% "
              "(%s) -- check the inputs for missing events before trusting the "
              "normalization." % (100*(max(rates)/min(rates) - 1),
                                  ", ".join("%.4g" % r for r in rates)))

    outdfs = []
    for df in dfs:
        outdf = df.merge(common_df, left_on=common_ind.names, right_index=True, how="left")
        outdf["common"] = outdf["common"].fillna(0)
        outdf = outdf[outdf.common == 1]
        outdfs.append(outdf)

    return outdfs, [pot_common]*len(pots)

# Systematic class helpers for what is in these files
class FluxSystematic(syst.WeightSystematic):
    def __init__(self, df, scale="glob_scale"):
        wgts = ["flux_univ%i" % i for i in range(100)]
        super().__init__(df, wgts, scale=scale)

class G4Systematic(syst.WeightSystematic):
    def __init__(self, df, scale="glob_scale"):
        wgts = ["g4_univ%i" % i for i in range(100)]
        super().__init__(df, wgts, scale=scale)

class XSecSystematic(syst.WeightSystematic):
    def __init__(self, df, scale="glob_scale"):
        super().__init__(df, ["%s_univ" % s for s in xsec_syst], avg=False, scale=scale)

class POTSystematic(syst.WeightSystematic):
    def __init__(self, df, scale="glob_scale"):
        super().__init__(df, ["pot_univ"], avg=False, scale=scale)


