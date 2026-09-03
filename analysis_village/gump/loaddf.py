import os
import hashlib
import json
import time
import sys

import pandas as pd
import numpy as np
import h5py
from scipy.interpolate import CubicSpline

from tqdm.auto import tqdm
import pyanalib.pandas_helpers as ph
from multiprocess import Pool
from functools import partial
import syst

# rwt_map and gumple_cuts live in analysis_village/gumple/. The cwd-relative
# entry is kept for compatibility with running from analysis_village/gump; the
# file-relative one makes the import work regardless of cwd.
workspace_root = os.getcwd()
sys.path.insert(0, workspace_root + "/../gumple/")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "gumple"))
import rwt_map as rw
import gumple_cuts as gmpl

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
EVTREC = "evtrec_%i"

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
  "baseline": ("baseline", ""),
  "true_E": ("nu_E", ""),
  "true_nu_pdg": ("pdg", ""),
  "true_issig": ("is_sig", ""),
  "true_isothernumucc": ("is_other_numucc", ""),
  "true_isfv": ("is_fv", ""),
  "genie_mode": ("genie_mode", ""),
  "true_pos_x": ("pos_x", ""),# we were already getting something called true_vtx_x from recodf
  "true_pos_y": ("pos_y", ""),# this was causing *_x, *_y to get appended to columns.
  "true_pos_z": ("pos_z", ""),
  "true_nmu": ("nmu", ""),
  "true_np": ("np", ""),
  "true_nn": ("nn", ""),
  "true_npi": ("npi", ""),
  "true_npi0": ("npi0", ""),
}


# ---------------------------------------------------------------------------
# GENIE event record (evtrec) -> pre-FSI truth kinematics
#
# The evtrec table (makedf.makedf.make_genie_evtrec_df) is the raw GHEP particle
# stack: one row per GENIE particle, index (__ntuple, entry, pindex), momenta and
# energies in GeV, vtx_* in metres.
# ---------------------------------------------------------------------------

# genie::EGHepStatus codes used below
_GHEP_INITIAL     = 0   # kIStInitialState        -- the probe and the target nucleus
_GHEP_NUCLEON_TGT = 11  # kIStNucleonTarget       -- struck nucleon, or a 2p2h cluster
_GHEP_PREFSI      = 14  # kIStHadronInTheNucleus  -- the hadrons that enter FSI

_NEUTRON_MASS = 0.939565
_PROTON_MASS = 0.938272

# INTRANUKE fate of a hadron that went through FSI, from the GHEP rescatter code.
# Surfaced per species as genie_prefsi_<species>_fsi, alongside the GHEP status code as
# genie_prefsi_<species>_status. Only status-14 particles carry a fate, so leptons and
# photons are always FSI_NONE. FSI_NOINT means the cascade ran and left the hadron
# alone: "FSI was applied" is `fsi >= FSI_NOINT`, "FSI changed the particle" is
# `fsi > FSI_NOINT`. NB there is no separate absorption fate -- a pion absorbed on the
# nucleus is FSI_INELAS with no pion among its daughters.
FSI_NONE   = -1  # no FSI applied: leptons, photons, and the hadrons INTRANUKE does not
                 # transport (eta, Lambda, rho, K0, Sigma+, ...)
FSI_NOINT  = 1   # cascade ran, hadron did not interact
FSI_CEX    = 2   # charge exchange, e.g. pi+ n -> pi0 p
FSI_ELAS   = 3   # elastic: identity kept, one nucleon knocked out
FSI_INELAS = 4   # inelastic, including pion absorption
FSI_PIPROD = 7   # pion production

# Only the codes above occur in the GUMP productions; anything else decodes to "codeN".
GENIE_FSI_NAMES = {FSI_NONE: "none", FSI_NOINT: "noint", FSI_CEX: "cex",
                   FSI_ELAS: "elas", FSI_INELAS: "inelas", FSI_PIPROD: "piprod"}

def genie_fsi_name(code):
    """Name of an FSI fate code, for labelling cuts and plots."""
    if code is None or (isinstance(code, float) and np.isnan(code)):
        return "nan"
    return GENIE_FSI_NAMES.get(int(code), "code%d" % int(code))

# pre-FSI hadron species: output moniker -> pdg selection. Monikers follow the
# make_mcdf convention (mu/p/p2/cpi/e), extended with n/n2 for the neutron, which has
# no make_mcdf counterpart. The photon is NOT here -- see below. NB the neutron IS
# transported by INTRANUKE, so unlike the photon it carries a real fate: its status is
# 14 and its fsi is FSI_NOINT/CEX/ELAS/INELAS/..., never FSI_NONE.
_PREFSI_PDG = {
    "p":   lambda pdg: pdg == 2212,
    "n":   lambda pdg: pdg == 2112,
    "cpi": lambda pdg: np.abs(pdg) == 211,
    "pi0": lambda pdg: pdg == 111,
}

# Photons are a special case. They never carry status 14 -- measured on the -13
# files, 0 of 42717 ICARUS Run2 events have one -- because a photon does not
# rescatter and so never enters INTRANUKE; every photon in the record is status 1.
# Their pre-FSI analogue is a photon born from a primary-vertex state, i.e. one
# whose mother is a decayed resonance (status 3, e.g. the RDecBR1gamma radiative
# Delta decay) or a pre-fragmentation hadronic state (status 12), or a pre-FSI
# hadron itself. Photons whose mother is the target nucleus (status 0) are nuclear
# de-excitation of the residual nucleus -- emitted after the interaction, not part
# of the primary hadronic system -- and are excluded here. That cut matters: 98%
# of the photons in the record (14388 of 14689) are de-excitation photons.
_PREFSI_GAMMA_MOTHER = [3, 12, 14]

# species carrying a momentum, in output order. "lep" is the primary lepton and
# "p2"/"n2" the sub-leading pre-FSI proton and neutron; all are handled specially below.
_GENIE_SPECIES = ["lep", "p", "p2", "n", "n2", "cpi", "g", "pi0"]

_GENIE_SCALARS = ["genie_Enu", "genie_q0", "genie_q3", "genie_W",
                  "genie_pmiss", "genie_emiss"]

# per species: the rotated momentum, the GHEP status code, and the FSI fate code
_GENIE_PER_SPECIES = ["px", "py", "pz", "status", "fsi"]

GENIE_COLS = _GENIE_SCALARS + ["genie_prefsi_%s_%s" % (s, q)
                               for s in _GENIE_SPECIES for q in _GENIE_PER_SPECIES]

def _p3(d):
    """(N,3) momentum array from a frame with px/py/pz columns."""
    return np.c_[d.px.to_numpy(float), d.py.to_numpy(float), d.pz.to_numpy(float)]

def _evtrec_link(mcdf):
    """Index of each mcnu row's entry in the GENIE event record.

    Prefer the stored link (rec.mc.nu.genie_evtrec_idx) where the production kept
    it -- maple does, gump's make_gump_nudf does not. Otherwise reconstruct it: the
    GenieEvtRecTree entry number is a running counter over the neutrinos of the
    input file, so it is the position of the mcnu row within its __ntuple.

    That reconstruction is verified against the interaction vertex by
    _evtrec_kinematics; on the -13 files it resolves every evtrec entry to exactly
    one mcnu row, with the vertex and Enu agreeing exactly. NB the *recTree* entry
    is NOT the evtrec entry -- assuming it is picks the wrong neutrino for 14% of
    ICARUS and 90% of SBND records, because a record holding two neutrinos advances
    the GENIE counter by two while advancing the recTree entry by one.
    """
    flat = [c[0] if isinstance(c, tuple) else c for c in mcdf.columns]
    if "genie_evtrec_idx" in flat:
        return mcdf[mcdf.columns[flat.index("genie_evtrec_idx")]]
    srt = mcdf.sort_index()
    return pd.Series(srt.groupby(level=0).cumcount().to_numpy(),
                     index=srt.index).reindex(mcdf.index)

def _evtrec_kinematics(er, mcdf):
    """Pre-FSI GENIE truth kinematics, indexed like `mcdf` so it can be joined onto
    the slice frame through tmatch_idx.

    All momenta are given in an event-by-event frame built from the record itself:

        z_hat = p_nu / |p_nu|                (the initial-state neutrino direction)
        y_hat = the outgoing lepton's momentum transverse to z_hat, normalised
        x_hat = y_hat x z_hat

    so the neutrino is (0, 0, Enu) and the primary lepton is (0, +pT, pL) -- the
    lepton px is identically zero and its py is positive by construction, and x is
    the out-of-plane direction.

    Returns (frame, stats) where stats carries the diagnostics load_one prints.
    """
    nan = pd.DataFrame(np.nan, index=mcdf.index, columns=GENIE_COLS)
    stats = {"n_mcnu": len(mcdf), "n_resolved": 0, "vtx_ok": np.nan, "n_degen": 0}
    if er is None or not len(er) or not len(mcdf):
        return nan, stats

    # --- per-GENIE-event pieces, all indexed by (__ntuple, evtrec entry) ---------
    # The probe sits at pindex 0 and its (single) daughter is the primary lepton.
    probe = er[er.index.get_level_values("pindex") == 0].droplevel("pindex")
    probe = probe[probe.status == _GHEP_INITIAL]

    lidx = pd.MultiIndex.from_arrays(
        [probe.index.get_level_values(0), probe.index.get_level_values(1),
         probe.fdaughter.to_numpy()], names=er.index.names)
    lep = er.reindex(lidx)
    lep.index = probe.index

    tgt = er[er.status == _GHEP_NUCLEON_TGT].groupby(level=[0, 1]).first().reindex(probe.index)

    # --- rotation basis ---------------------------------------------------------
    pnu, plep = _p3(probe), _p3(lep)
    with np.errstate(invalid="ignore", divide="ignore"):
        zh = pnu / np.linalg.norm(pnu, axis=1, keepdims=True)
        pt = plep - np.sum(plep*zh, axis=1, keepdims=True)*zh
        yh = pt / np.linalg.norm(pt, axis=1, keepdims=True)

    # A lepton collinear with the neutrino defines no transverse direction: |pt| is then
    # pure round-off and pt/|pt| is noise, giving a non-orthonormal frame and inflated
    # momenta. The azimuth is arbitrary for such an event, so fix it by convention --
    # detector y, made transverse to the neutrino -- which keeps the frame orthonormal,
    # so magnitudes and relative angles stay right and the lepton lands at (0, ~0, |p|).
    # Safe because the beam is far from detector y (|z_hat.y_det| < 0.01). The threshold
    # is uncritical: affected events sit at |pt|/|p| <= 4e-16, the next one up at 7e-4.
    # `~( > )` also catches the exactly-collinear rows, where yh is already NaN.
    degen = ~(np.linalg.norm(pt, axis=1) > 1e-9*np.linalg.norm(plep, axis=1))
    if degen.any():
        yfb = np.zeros_like(zh)
        yfb[:, 1] = 1.
        yfb -= np.sum(yfb*zh, axis=1, keepdims=True)*zh
        yfb /= np.linalg.norm(yfb, axis=1, keepdims=True)
        yh = np.where(degen[:, None], yfb, yh)
    xh = np.cross(yh, zh)

    def rot(d):
        """Project a species' momentum onto (x_hat, y_hat, z_hat)."""
        v = _p3(d.reindex(probe.index))
        return np.c_[np.sum(v*xh, 1), np.sum(v*yh, 1), np.sum(v*zh, 1)]

    # --- scalars ----------------------------------------------------------------
    Enu = probe.E.to_numpy(float)
    q0 = Enu - lep.E.to_numpy(float)
    q3 = np.linalg.norm(pnu - plep, axis=1)
    Q2 = q3**2 - q0**2

    # W in GENIE's own convention: an on-shell nucleon at rest, so this is directly
    # comparable to rec.mc.nu.w. The target pdg is a di-nucleon cluster (2000000200
    # /201/202/300) for 2p2h, where the neutron/proton average is the best available
    # nucleon mass.
    tpdg = tgt.pdg.to_numpy(float)
    is_nucleon = np.isin(tpdg, [2212, 2112])
    M = np.where(tpdg == 2212, _PROTON_MASS,
                 np.where(tpdg == 2112, _NEUTRON_MASS, 0.5*(_PROTON_MASS + _NEUTRON_MASS)))
    W2 = M**2 + 2*M*q0 - Q2
    W = np.sqrt(np.where(W2 > 0, W2, np.nan))

    # Fermi momentum and removal energy of the struck nucleon -- the quantities the
    # LFG/SF/HF CCQE template dials move. Undefined for a 2p2h cluster target, so
    # those events are left NaN rather than quietly mixing in a cluster momentum.
    ptgt = np.linalg.norm(_p3(tgt), axis=1)
    pmiss = np.where(is_nucleon, ptgt, np.nan)
    emiss = np.where(is_nucleon, M - tgt.E.to_numpy(float), np.nan)

    out = pd.DataFrame({"genie_Enu": Enu, "genie_q0": q0, "genie_q3": q3,
                        "genie_W": W, "genie_pmiss": pmiss, "genie_emiss": emiss},
                       index=probe.index)

    # --- pre-FSI species, leading by |p| ----------------------------------------
    pre = er[er.status == _GHEP_PREFSI]
    pre = pre.assign(_pmag=np.linalg.norm(_p3(pre), axis=1)).sort_values("_pmag", ascending=False)

    parts = {"lep": lep}
    for name, sel in _PREFSI_PDG.items():
        s = pre[sel(pre.pdg.to_numpy())]
        parts[name] = s.groupby(level=[0, 1]).head(1).droplevel("pindex")
    # sub-leading proton: second row of the momentum-ordered proton list
    sp = pre[pre.pdg == 2212]
    sp = sp.groupby(level=[0, 1]).head(2)
    parts["p2"] = sp[sp.groupby(level=[0, 1]).cumcount() == 1].droplevel("pindex")
    # sub-leading neutron: the same, on the momentum-ordered neutron list
    sn = pre[pre.pdg == 2112]
    sn = sn.groupby(level=[0, 1]).head(2)
    parts["n2"] = sn[sn.groupby(level=[0, 1]).cumcount() == 1].droplevel("pindex")

    # photons: status-1, but only those from a primary-vertex parent (see above)
    gam = er[er.pdg == 22]
    if len(gam):
        gmom = er.reindex(pd.MultiIndex.from_arrays(
            [gam.index.get_level_values(0), gam.index.get_level_values(1),
             gam.fmother.to_numpy()], names=er.index.names))
        gam = gam[np.isin(gmom.status.to_numpy(float), _PREFSI_GAMMA_MOTHER)]
        gam = gam.assign(_pmag=np.linalg.norm(_p3(gam), axis=1)).sort_values("_pmag", ascending=False)
    parts["g"] = gam.groupby(level=[0, 1]).head(1).droplevel("pindex")

    for name in _GENIE_SPECIES:
        d = parts[name].reindex(probe.index)
        v = rot(d)
        for i, c in enumerate("xyz"):
            out["genie_prefsi_%s_p%s" % (name, c)] = v[:, i]
        # status is 14 for the pre-FSI hadrons and 1 for the lepton and photon, so it
        # says which species saw FSI at all without knowing how each was selected.
        out["genie_prefsi_%s_status" % name] = d.status.to_numpy(float)
        out["genie_prefsi_%s_fsi" % name] = d.rescat.to_numpy(float)

    # events whose azimuth came from the detector-y fallback rather than the lepton
    stats["n_degen"] = int(degen.sum())

    # --- map onto the mcnu index via the evtrec link ----------------------------
    link = _evtrec_link(mcdf).to_numpy(float)
    link = np.where(np.isfinite(link), link, -1).astype(np.int64)
    gidx = pd.MultiIndex.from_arrays([mcdf.index.get_level_values(0), link],
                                     names=out.index.names)
    res = out.reindex(gidx)
    res.index = mcdf.index

    # --- self-check: the link is reconstructed, so verify it against the vertex --
    # evtrec vtx_* is in metres, mcnu pos_* in cm. Both name the same interaction
    # point, so on a correct link they agree to round-off on every resolved row.
    vtx = probe[["vtx_x", "vtx_y", "vtx_z"]].reindex(gidx)
    got = np.isfinite(vtx.vtx_x.to_numpy(float))
    stats["n_resolved"] = int(got.sum())
    if got.any():
        pos = np.c_[mcdf.pos_x.to_numpy(float), mcdf.pos_y.to_numpy(float),
                    mcdf.pos_z.to_numpy(float)]
        d = np.abs(vtx.to_numpy(float)*100. - pos).max(axis=1)
        stats["vtx_ok"] = float((d[got] < 1e-3).mean())

    return res, stats

detvar_rwt_files = [
  'SBND_WMXThetaXW.txt',
  'SBND_WMYZ.txt',
  'ICARUSRun2_WMYZ.txt',
  'ICARUSRun4_WMYZ.txt',
  'SBND_DENT.txt',
  ['SBND_0xSCE.txt', 'SBND_2xSCE.txt'],
  'ICARUSRun2_SCE.txt',
  'ICARUSRun4_SCE.txt',
  'SBND_SmeareddEdx.txt',
  'ICARUSRun2_SmeareddEdx.txt',
  'ICARUSRun2_WMXThetaXW.txt',
  'ICARUSRun4_SmeareddEdx.txt',
  'ICARUSRun4_WMXThetaXW.txt',
  'SBND_GainHi.txt',
  'ICARUSRun2_GainHi.txt',
  'ICARUSRun4_GainHi.txt',
  'SBND_EMBAlpha.txt',
  'ICARUSRun2_EMBAlpha.txt',
  'ICARUSRun4_EMBAlpha.txt',
  'SBND_EMBBeta.txt',
  'ICARUSRun2_EMBBeta.txt',
  'ICARUSRun4_EMBBeta.txt',
  'SBND_EMBR.txt',
  'ICARUSRun2_EMBR.txt',
  'ICARUSRun4_EMBR.txt',
  ['SBND_TrigEffMin.txt', 'SBND_TrigEffPls.txt'],
  ['ICARUSRun2_TrigEffMin.txt', 'ICARUSRun2_TrigEffPls.txt'],
  ['ICARUSRun4_TrigEffMin.txt', 'ICARUSRun4_TrigEffPls.txt'],
  'SBND_BIND.txt',
  'ICARUSRun2_BIND.txt',
  'ICARUSRun4_BIND.txt',
  'ICARUSRun2_Z=0_TRKSPLT.txt',
  'ICARUSRun4_Z=0_TRKSPLT.txt',
  'ICARUSRun2_EastCathode_TRKSPLT.txt',
  'ICARUSRun4_EastCathode_TRKSPLT.txt',
  'ICARUSRun2_WestCathode_TRKSPLT.txt',
  'ICARUSRun4_WestCathode_TRKSPLT.txt',
]

detvar_rwt_lbls = [
  'WireMod_SBND_multisigma_WMXThetaXW',
  'WireMod_SBND_multisigma_WMYZ',
  'WireMod_ICARUSRun2_multisigma_WMYZ',
  'WireMod_ICARUSRun4_multisigma_WMYZ',
  'DENT_SBND_multisigma_DENT',
  'SCE_SBND_multisigma_SCE',
  'SCE_ICARUSRun2_multisigma_SCE',
  'SCE_ICARUSRun4_multisigma_SCE',
  'SBND_PID_Smear',
  'ICARUSRun2_PID_Smear',
  'WireMod_ICARUSRun2_multisigma_WMXThetaXW',
  'ICARUSRun4_PID_Smear',
  'WireMod_ICARUSRun4_multisigma_WMXThetaXW',
  'SBND_PID_Gain',
  'ICARUSRun2_PID_Gain',
  'ICARUSRun4_PID_Gain',
  'SBND_PID_Alpha',
  'ICARUSRun2_PID_Alpha',
  'ICARUSRun4_PID_Alpha',
  'SBND_PID_Beta',
  'ICARUSRun2_PID_Beta',
  'ICARUSRun4_PID_Beta',
  'SBND_PID_R',
  'ICARUSRun2_PID_R',
  'ICARUSRun4_PID_R',
  'SBND_TrigEff',
  'ICARUSRun2_TrigEff',
  'ICARUSRun4_TrigEff',
  'BIND',
  'BIND',
  'BIND',
  'ICARUSRun2_Z=0_TRKSPLT',
  'ICARUSRun4_Z=0_TRKSPLT',
  'ICARUSRun2_EastCathode_TRKSPLT',
  'ICARUSRun4_EastCathode_TRKSPLT',
  'ICARUSRun2_WestCathode_TRKSPLT',
  'ICARUSRun4_WestCathode_TRKSPLT',
]

std_drops = ['is_clear_cosmic', 'crlongtrkdiry', 'p_len', 'has_stub',
             'true_pcand_pdg', 'true_p_dir_x', 'true_p_dir_y', 'true_p_dir_z', 
             'true_pcand_dir_x', 'true_pcand_dir_y', 'true_pcand_dir_z', 
             'true_pcand_end_x', 'true_pcand_end_y', 'true_pcand_end_z',
             'true_mucand_pdg', 'true_mucand_dir_x', 'true_mucand_dir_y', 
             'true_mucand_dir_z', 'true_mucand_end_x', 'true_mucand_end_y', 
             'true_mucand_end_z', 'stub_l0_5cm_dedx','stub_l0_5cm_charge',
             'stub_l1cm_dedx','stub_l1cm_charge','stub_l2cm_dedx',
             'stub_l2cm_charge','stub_l3cm_dedx','stub_l3cm_charge',
             'stub_l4cm_dedx','stub_l4cm_charge','prot_chi2smear5_of_prot_cand', 
             'prot_chi2smear5_of_mu_cand', 'mu_chi2smear5_of_mu_cand', 
             'mu_chi2smear5_of_prot_cand', 'tmatch_pur', 'tmatch_eff', 
             'true_baseline', 'true_nu_pdg_x', 'true_nu_pdg_y',
             'true_nmu_27MeV', 'true_np_20MeV', 'true_np_50MeV', 
             'true_npi_30MeV', 'flash_sumpe', 'true_mucand_p', 
             'true_pcand_p', 'p_true_p', 'true_mu_end_x', 
             'true_p_end_x', 'true_mu_end_y', 'true_p_end_y', 'true_mu_end_z', 
             'true_p_end_z', 'true_nu_E', 'p_true_pdg', 'mu_true_pdg', 
             'mu_chi22lo_of_mu_cand', 'mu_chi22hi_of_mu_cand', 
             'prot_chi22lo_of_mu_cand', 'prot_chi22hi_of_mu_cand',
             'mu_chi22lo_of_prot_cand', 'mu_chi22hi_of_prot_cand', 
             'prot_chi22lo_of_prot_cand', 'prot_chi22hi_of_prot_cand', 
             'true_mu_p', 'true_p_p', 'pot_univ']

def get_std_drops():
    return std_drops

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

# Track-splitting planes for the split_tracks argument: region -> (dim, coord).
# Same planes as TrackSplittingCorrection.py / the signal-box track-split
# systematic: the cathodes sit at the center of each ICARUS drift volume.
SPLIT_REGIONS = {
    "Z=0":          ("z", 0.0),
    "East Cathode": ("x", 0.5*(gmpl.ICARUSRun4FVCuts["C0"]["x"]["min"] + gmpl.ICARUSRun4FVCuts["C0"]["x"]["max"])),
    "West Cathode": ("x", 0.5*(gmpl.ICARUSRun4FVCuts["C1"]["x"]["min"] + gmpl.ICARUSRun4FVCuts["C1"]["x"]["max"])),
}

# Fraction of the plane-crossing muons that split_tracks actually splits, per
# region: the Run2+Run4 combined values measured by TrackSplittingCorrection.py
# (TrackSplittingCorrection-2026-08-01.md), as used by the signal-box track-split
# systematic.
SPLIT_FRAC = {"Z=0": 0.0960, "East Cathode": 0.0447, "West Cathode": 0.0398}

# Run periods each plane applies to. The east cathode lies outside the Run 2 muon
# fiducial volume (the East-East TPC was off in Run 2), so it has Run 4 crossers
# only -- same table as TrackSplittingCorrection.split_points().
SPLIT_RUNS = {"Z=0": [2, 4], "East Cathode": [4], "West Cathode": [2, 4]}

# Nominal binding-energy shift [GeV] for shift_binding_E, as in the signal-box
# BE systematic
BE_SHIFT = 0.025

# Fraction of events the binding-energy shift is applied to, as in the signal-box
# BE systematic (eres_ar23_ar25.BE_FRACTION, mcdata_comparison)
BE_FRACTION = 0.5

# Columns syst.split_tracks recomputes on the plane-crossing rows. mu_T is no
# longer stored in the GUMPLE dfs (recompute as mu_E - MUON_MASS if needed).
_SPLIT_COLS = ["mu_end_x", "mu_end_y", "mu_end_z", "mu_len",
               "nu_E_calo", "del_p", "del_Tp", "del_phi", "mu_E"]

def _weight_col(df):
    """The column the mixture weight is folded into.

    _apply_variations runs inside load_one, before scale_pot creates glob_scale;
    cvwgt is always there by then and scale_pot computes glob_scale = scale*cvwgt,
    so folding into cvwgt propagates. Prefer glob_scale when it already exists so
    this is also correct on an already-scaled frame.
    """
    return "glob_scale" if "glob_scale" in df.columns else "cvwgt"

def _mix(df, varied, rows, frac):
    """Deterministic f-weighted mixture of a variation applied to `rows`.

    The affected rows enter twice -- unvaried at (1-f)*w and varied at f*w --
    and every other row is left alone, so histogramming the result with the
    weight column reproduces (1-f)*nominal + f*varied exactly, with no added MC
    statistical noise. Same convention as syst.TrackSplittingSystematic.
    (syst.shift_binding_energy instead shifts a deterministic fraction of the
    rows in place, keeping the row count.)

    `rows` is positional (index labels are not guaranteed unique), and `varied`
    holds exactly those rows, already carrying their unscaled weights.
    """
    w = _weight_col(df)
    df = df.copy()
    wi = df.columns.get_loc(w)
    df.iloc[rows, wi] = df.iloc[rows, wi].to_numpy()*(1 - frac)
    varied = varied.copy()
    varied[w] = varied[w].to_numpy()*frac
    return pd.concat([df, varied])

def _apply_variations(df, shift_binding_E, split_tracks,
                      shift_fraction=None, split_fraction=None):
    """Apply the requested variations to a loaded df.

    split_tracks names a plane in SPLIT_REGIONS; muons crossing it are truncated
    there (syst.split_tracks), restricted to the runs the plane applies to
    (SPLIT_RUNS). shift_binding_E recomputes the reco kinematics under
    BE -> BE + BE_SHIFT (syst.shift_binding_energy; needs genie_mode, i.e.
    load_truth on MC).

    Each variation is applied to only a fraction of the events -- SPLIT_FRAC[region]
    and BE_FRACTION, the values measured for / used by the signal-box systematics --
    as the deterministic f-weighted mixture described in _mix. Pass
    split_fraction / shift_fraction to override; 1.0 varies every affected event
    in place, as before the fractions existed.

    With a split fraction below 1 the returned frame has MORE ROWS than the
    input and duplicate index labels (the split adds a copy of the crossing
    rows via _mix), though the total weight is unchanged. The BE shift keeps
    the row count: syst.shift_binding_energy recomputes the shifted kinematics
    in place on a deterministic, evenly-interleaved fraction of the rows and
    leaves the rest at their stored CV values. Cut columns must be evaluated
    downstream of this (loaddf computes none itself): the split moves
    mu_end_*, which the FV cuts read.

    Run after the cache read/write on purpose: the cache always holds the
    unvaried df, so one cached load serves every variant and adding a variation
    does not bust the existing cache.
    """
    if split_tracks is not None:
        dim, coord = SPLIT_REGIONS[split_tracks]
        f = SPLIT_FRAC[split_tracks] if split_fraction is None else split_fraction
        s, crosses = syst.split_tracks(df, dim, coord, runs=SPLIT_RUNS[split_tracks])
        # positional assignment: index labels are not guaranteed unique
        rows = np.flatnonzero(crosses)
        if len(rows) == 0:
            pass # no crossers -> nothing to vary (and nothing to concat)
        elif f >= 1.0:
            for c in _SPLIT_COLS:
                df.iloc[rows, df.columns.get_loc(c)] = s[c].to_numpy()
        else:
            df = _mix(df, s, rows, f)
    if shift_binding_E:
        f = BE_FRACTION if shift_fraction is None else shift_fraction
        # hand off to shift_binding_energy so this stays the same universe as
        # the signal-box systematic: it shifts a deterministic fraction f of
        # the rows in place (unshifted rows keep their stored CV columns; row
        # count and weights unchanged). Costs a whole-frame copy, which is why
        # the caller should slim the frame first if memory is tight.
        df = syst.shift_binding_energy(df, BE_SHIFT, fraction=f, scale=_weight_col(df))
    return df

_DETECTOR_RUN = {"SBND": 1, "ICARUS Run2": 2, "ICARUS Run4": 4}

def load_one(fname, idf,
    detector=None, # One of SBND, ICARUS, ICARUS Run4
    include_syst=True, nuniv=100, spline=False, xsec_univ=False, xsec_spline=False,# systematic handling
    reweight_aFF=False, pot_univ=False, flux_univ=True, sep_flux_univ=False, g4_univ=True, sep_g4_univ=False,
    pot_spline=False, detvar_spline=False, spline_dir="rwt_outputs",
    load_truth=True, load_crt=False, load_evtrec=False, match_Enu=True, # load extra information
    offbeampot=False, # POT handling
    preselection=None, # apply preselection cut
    shift_binding_E=False, split_tracks=None, # variations applied to the output df (see _apply_variations)
    shift_fraction=None, split_fraction=None, # fraction of events each variation is applied to (None -> BE_FRACTION / SPLIT_FRAC)
    cache_dir=None, # directory to cache output; None disables caching
    flashname=FLASH, hdrname=HDR, evtname=EVT, wgtname=WGT, mcname=MC, crtname=CRT, evtrecname=EVTREC, drops=None, lightmem=False): # override default table names
    assert(detector == "SBND" or detector == "ICARUS Run2" or detector == "ICARUS Run4")
    # Check cache
    if cache_dir is not None:
        cache_hash = _cache_key(fname, idf, detector=detector, include_syst=include_syst,
            nuniv=nuniv, spline=spline, xsec_univ=xsec_univ, xsec_spline=xsec_spline, reweight_aFF=reweight_aFF, pot_univ=pot_univ,
            flux_univ=flux_univ, sep_flux_univ=sep_flux_univ, g4_univ=g4_univ,
            load_truth=load_truth, load_crt=load_crt, load_evtrec=load_evtrec,
            match_Enu=match_Enu, offbeampot=offbeampot, preselection=preselection,
            drops=drops, lightmem=lightmem,
            flashname=flashname, hdrname=hdrname, evtname=evtname,
            wgtname=wgtname, mcname=mcname, crtname=crtname, evtrecname=evtrecname)
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
            return _apply_variations(df, shift_binding_E, split_tracks, shift_fraction, split_fraction), match, pot

    df =  pd.read_hdf(fname, evtname % idf)
    hdr = pd.read_hdf(fname, hdrname % idf)
    ismc = hdr.ismc.iloc[0] == 1

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
    # in particular AVnu, which is derived from gmpl._fv_cut and so depends on `detector`.
    match_ind = list(match.columns)
    # if needed, include neutrino energy in matching information
    if match_Enu:
        mcdf = pd.read_hdf(fname, mcname % idf)
        mcdf["detector"] = detector
        if "Run2" in detector.replace(" ", ""):
            mcdf["Run"] = 2
        elif "Run4" in detector.replace(" ", ""):
            mcdf["Run"] = 4
        elif "SBND" in detector.replace(" ", ""):
            mcdf["Run"] = 1
        else:
            print("Run ambiguous!")

        match = match.merge(mcdf.nu_E.groupby(level=[0,1]).max().rename("nu_E0"), on=["__ntuple", "entry"], how="left")
        match_ind = list(match.columns)

        # Add in other meta-data to match.
        # Newer makedf productions stamp detector/Run onto the mcnu frame
        # (makedf.py:769-770); the sbn-rewgted-13/-14 files predate that, so fall
        # back to the detector this load was told to use. gmpl._fv_cut only consults
        # Run when the label is the ambiguous "ICARUS", which the fallback never is.
        vtx = pd.DataFrame({
          "detector": mcdf.detector if "detector" in mcdf.columns else detector,
          "Run": mcdf.Run if "Run" in mcdf.columns else _DETECTOR_RUN[detector],
          "x": mcdf.pos_x,
          "y": mcdf.pos_y,
          "z": mcdf.pos_z,
        })
        any_in_AV = gmpl._fv_cut(vtx, 0, 0, 0, 0).groupby(level=[0,1]).any().rename("AVnu")
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

    # LOAD GENIE EVENT RECORD
    # Pre-FSI truth kinematics from the raw GHEP stack, joined on tmatch_idx exactly
    # like the truth block above: tmatch_idx names the mcnu row, and _evtrec_kinematics
    # returns its columns on the mcnu index. Slices with no truth match (cosmics) fall
    # out of the left join as NaN, as they already do for truthvars.
    if load_evtrec:
        with h5py.File(fname, "r") as f:
            has_evtrec = (evtrecname % idf) in f
        if has_evtrec:
            er = pd.read_hdf(fname, evtrecname % idf)
            mcdf = pd.read_hdf(fname, mcname % idf)
            gdf, gstats = _evtrec_kinematics(er, mcdf)
            del er
            # The evtrec link is reconstructed, not stored (see _evtrec_link), so say
            # out loud how much of the sample it reached and shout if the vertex
            # cross-check fails -- a broken link would silently attach another
            # neutrino's kinematics rather than raise.
            frac = gstats["n_resolved"] / max(gstats["n_mcnu"], 1)
            print(f"[{os.path.basename(fname)} idf={idf}] evtrec: "
                  f"{gstats['n_resolved']}/{gstats['n_mcnu']} neutrinos resolved "
                  f"({100*frac:.1f}%), vertex check {gstats['vtx_ok']:.5f}, "
                  f"{gstats['n_degen']} collinear-lepton frame(s)")
            if gstats["n_resolved"] > 0 and not (gstats["vtx_ok"] > 0.999):
                print(f"WARNING: {os.path.basename(fname)} idf={idf}: the GENIE event "
                      f"record link does not reproduce the mcnu vertex "
                      f"({gstats['vtx_ok']:.5f} agree) -- genie_* columns are NOT "
                      f"trustworthy for this file.")
        else:
            # No evtrec in this file (data, detvar, dirt, ...). Emit the columns as NaN
            # so a mixed file list still concatenates to one schema.
            gdf = pd.DataFrame(np.nan, index=pd.read_hdf(fname, mcname % idf).index,
                               columns=GENIE_COLS)
        df = df.merge(gdf, left_on=["__ntuple", "entry", "tmatch_idx"], right_index=True, how="left")

    # LOAD CRT
    if load_crt:
        if "crthit" in df.columns: del df["crthit"]

        crtdf = pd.read_hdf(fname, crtname % idf)
        crthit = ((crtdf.time > -1) & (crtdf.time < 1.8) & (crtdf.plane != 50)).groupby(level=[0, 1]).any()
        crthit.name = "crthit"
        df = df.join(crthit, on=["__ntuple", "entry"])

    df["crthit"] = df.crthit.fillna(False).astype(bool) 

    # LOAD WEIGHTS
    if include_syst:
        wgt = pd.read_hdf(fname, wgtname % idf) 

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
        df["total_pot"] = pot
        return _apply_variations(df, shift_binding_E, split_tracks, shift_fraction, split_fraction), match, pot

    # LOAD WEIGHTS
    wgt = pd.read_hdf(fname, wgtname % idf) 
    skim = {}

    if flux_univ:
        num_to_process = min(100, nuniv)
        
        # Pre-cache the system lookups to avoid doing it inside the inner loops
        system_data = [wgt[s] for s in flux_syst]
        
        new_columns_dict = {}
        for i in range(num_to_process):
            univ_key = "univ_%i" % i
            # np.prod over the pre-cached systems list
            new_columns_dict["flux_univ%i" % i] = np.prod([sys[univ_key] for sys in system_data], axis=0)
            
        # --- FIX HERE: Merging two dictionaries ---
        skim.update(new_columns_dict)

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

    if pot_spline:
        for d in ["SBND", "ICARUS Run2", "ICARUS Run4"]:
            col_str = f"multisigma_{d.replace(' ', '')}_POT"
            multisigma_cols.append(col_str)
            if detector == d:
                skim[f"{col_str}"] = [list(pot_syst.values()) for _ in range(len(wgt))]
            else:
                skim[f"{col_str}"] = [[1.0]*7 for _ in range(len(wgt))]

    if sep_flux_univ:
        for j, s in enumerate(flux_syst):
            w = wgt[s]
            col_str = s if 'multisim' in s else 'multisim_' + s
            n_univs = min(100, nuniv)
            univ_cols = [f"univ_{i}" for i in range(n_univs)]
            stacked_variants = w[univ_cols].to_numpy(dtype=np.float32).T
            np.nan_to_num(stacked_variants, copy=False, nan=1.0, posinf=1.0, neginf=1.0)
            skim[col_str] = stacked_variants.T.tolist()
            multisim_cols.append(col_str)

            #if not 'multisim' in s:
            #    col_str = 'multisim_'+s
            #else:
            #    col_str = s
            #multisim_cols.append(col_str)
            #w = wgt[s]
            #if lightmem:
            #    float64_cols = w.select_dtypes(include=["float64"]).columns
            #    w.loc[:, float64_cols] = w.loc[:, float64_cols].astype("float32")
            #stacked_variants = np.vstack([np.nan_to_num(w["univ_%i" % i].to_numpy(), nan=1.0, posinf=1.0, neginf=1.0) for i in range(min(100, nuniv))])
            #skim[col_str] = stacked_variants.T.tolist()

    if sep_g4_univ:
        for j, s in enumerate(g4_syst):
            if not 'multisim' in s:
                col_str = 'multisim_'+s
            else:
                col_str = s
            multisim_cols.append(col_str)
            w = wgt[s]
            if lightmem:
                float64_cols = w.select_dtypes(include=["float64"]).columns
                if len(float64_cols) > 0:
                            w.loc[:, float64_cols] = w.loc[:, float64_cols].astype("float32")
            stacked_variants = np.vstack([np.nan_to_num(w["univ_%i" % i].to_numpy(), nan=1.0, posinf=1.0, neginf=1.0) for i in range(min(100, nuniv))])
            skim[col_str] = stacked_variants.T.tolist()

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
                if not 'multisigma' in s:
                    col_str = 'multisigma_'+s
                else:
                    col_str = s
                skim[col_str] = stacked_variants.T.tolist()
                multisigma_cols.append(col_str)
            elif "morph" in wgt[s]:
                w = wgt[s].fillna(1).replace([np.inf, -np.inf], 1)
                if lightmem:
                    float64_cols = w.select_dtypes(include=["float64"]).columns
                    w.loc[:, float64_cols] = w.loc[:, float64_cols].astype("float32")

                stacked_variants = np.vstack([
                    np.ones(len(w)),  # Central value ratio is exactly 1.0
                    np.clip((w["morph"]).to_numpy(), 0, 10)
                ])
                if not 'multisigma' in s:
                    col_str = 'multisigma_'+s
                else:
                    col_str = s
                skim[col_str] = stacked_variants.T.tolist()
                multisigma_cols.append(col_str)
            elif "multisim" in s:
                w = wgt[s]
                n_univs = min(100, nuniv)
                univ_cols = [f"univ_{i}" for i in range(n_univs)]
                stacked_variants = w[univ_cols].to_numpy(dtype=np.float32).T
                np.nan_to_num(stacked_variants, copy=False, nan=1.0, posinf=1.0, neginf=1.0)
                skim[s] = stacked_variants.T.tolist()
                col_str = s if 'multisim' in s else 'multisim_' + s
                multisim_cols.append(s)

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

    if detvar_spline:
        for s, f in zip(detvar_rwt_lbls, detvar_rwt_files):
            if isinstance(f, (str, bytes)):
                fs = [spline_dir + '/' + f]
            else:
                fs = [spline_dir + '/' + fi for fi in f]
            
            allowed_substrings = ["ICARUSRun4", "ICARUSRun2", "SBND"]

            if not all(any(sub in s for sub in allowed_substrings) for s in fs):
                # Find the specific offender to make the error message helpful
                invalid_string = next(s for s in fs if not any(sub in s for sub in allowed_substrings))
                raise ValueError(f"Validation failed: '{invalid_string}' is invalid. Check that your reweight files are all for the same detector.")

            if not 'multisigma' in s:
                col_str = 'multisigma_' + s
            else:
                col_str = s

            # allow for f
            if detector.replace(' ', '') in fs[0]:
                s_df = rw.apply_map(mrg, fs, s)
                mrg[col_str] = s_df
            elif not col_str in mrg.columns:
            # Don't overwrite existing columns with defaults (eg BIND)
                mrg[col_str] = [[1.0]*(len(fs)+1) for _ in range(len(mrg))]
            multisigma_cols.append(col_str)

    univ_cols = [col for col in skim.columns if "univ" in col]
    if len(multisigma_cols) > 0:
        nan_mask = mrg[multisigma_cols[0]].isna()
        n_missing = nan_mask.sum()
        for col in multisigma_cols:
            valid_rows = mrg.loc[~nan_mask, col]
            if len(valid_rows) > 0:
                col_len = len(valid_rows.iloc[0])
            else:
                col_len = 7  # Fallback to standard 7-knot default if the whole block is NaN
            
            # 2. Vectorized assignment: Create the block of lists all at once
            default_val = [1.0] * col_len
            mrg.loc[nan_mask, col] = pd.Series([default_val] * n_missing, index=mrg.index[nan_mask])

    if len(multisim_cols) > 0:

        for col in multisim_cols:
            nan_mask = mrg[col].isna()
            n_missing = nan_mask.sum()
            valid_rows = mrg.loc[~nan_mask, col]
            col_len = 100#len(mrg[col].iloc[0]) 

            # 2. Vectorized assignment: Create the block of lists all at once
            default_val = [1.0] * col_len
            mrg.loc[nan_mask, col] = pd.Series([default_val] * n_missing, index=mrg.index[nan_mask])

    if len(univ_cols) > 0:
        mrg.loc[np.isnan(mrg[univ_cols[0]]), univ_cols] = 1.0 

    if drops is not None:
        mrg.drop(columns=drops, inplace=True, errors='ignore')
    if cache_dir is not None:
        _write_cache(cache_file, mrg, match, pot)

    mrg["total_pot"] = pot
    return _apply_variations(mrg, shift_binding_E, split_tracks, shift_fraction, split_fraction), match, pot


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

    df["total_pot"] = pots
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

    df["total_pot"] = pots
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
