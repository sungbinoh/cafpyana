"""Systematic weights per true neutrino (rec.mc.nu.wgt), as gump's
make_gump_nuwgtdf / make_gump_nuslimwgtdf (analysis_village/gump/makedf.py,
makedf/makedf.make_mcnudf): flux (makedf.bnbsyst), GENIE (makedf.geniesyst
with GENIE_SYSTEMATICS) and Geant4 reinteractions (makedf.g4syst).

Index (entry, rec.mc.nu..index) as the mcnu table (make_truth.make_sp_truthdf), so the
weights join to the candidates on (file, __ntuple, entry, tmatch_idx) like the
truth. Columns are (systematic, variation): multisim (syst, "univ_<i>"),
+/- sigma unisims (syst, "ps1"/"ms1"/..., "cv"), morph unisims (syst,
"morph"); slim=True instead keeps the products over each group, ("Flux",
"univ_<i>") and ("GENIE", "univ_<i>"), plus the Geant4 multisims.

Only neutrinos that can enter the analysis are kept, to bound the size: a true
vertex in the active volume, or matched (T_match) to a Wire-Cell candidate.
Neutrinos without a row get weight 1 in every universe.

The imports avoid makedf.makedf (it needs CVMFS for chi2pid).
"""

import numpy as np
import pandas as pd

from makedf import bnbsyst, geniesyst, g4syst, getsyst
from makedf import wirecell as wc

from analysis_village.SinglePhotonWC import cuts
from analysis_village.SinglePhotonWC.make_truth import _load

# gump_genie_systematics (analysis_village/gump/makedf.py), plus the radiative and
# eta resonance-decay branching ratios (in later gump versions) that move the NC
# Delta -> N gamma signal
GENIE_SYSTEMATICS = [
    # CCQE
    "GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape",
    "GENIEReWeight_SBN_v1_multisigma_RPA_CCQE",
    "GENIEReWeight_SBN_v1_multisigma_CoulombCCQE",
    "GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE",
    "GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE",
    "GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE",
    "GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE",
    # MEC
    "GENIEReWeight_SBN_v1_multisigma_NormCCMEC",
    "GENIEReWeight_SBN_v1_multisigma_NormNCMEC",
    "GENIEReWeight_SBN_v1_multisigma_DecayAngMEC",
    # RES
    "GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi",
    "GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad",
    "GENIEReWeight_SBN_v1_multisigma_MaCCRES",
    "GENIEReWeight_SBN_v1_multisigma_MaNCRES",
    "GENIEReWeight_SBN_v1_multisigma_MvCCRES",
    "GENIEReWeight_SBN_v1_multisigma_MvNCRES",
    "GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma",   # Delta -> N gamma branching ratio
    "GENIEReWeight_SBN_v1_multisigma_RDecBR1eta",
    # DIS
    "GENIEReWeight_SBN_v1_multisigma_AhtBY",
    "GENIEReWeight_SBN_v1_multisigma_BhtBY",
    "GENIEReWeight_SBN_v1_multisigma_CV1uBY",
    "GENIEReWeight_SBN_v1_multisigma_CV2uBY",
    # COH
    "GENIEReWeight_SBN_v1_multisigma_NormCCCOH",
    "GENIEReWeight_SBN_v1_multisigma_NormNCCOH",
    # FSI
    "GENIEReWeight_SBN_v1_multisigma_MFP_pi",
    "GENIEReWeight_SBN_v1_multisigma_FrCEx_pi",
    "GENIEReWeight_SBN_v1_multisigma_FrInel_pi",
    "GENIEReWeight_SBN_v1_multisigma_FrAbs_pi",
    "GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi",
    "GENIEReWeight_SBN_v1_multisigma_MFP_N",
    "GENIEReWeight_SBN_v1_multisigma_FrCEx_N",
    "GENIEReWeight_SBN_v1_multisigma_FrInel_N",
    "GENIEReWeight_SBN_v1_multisigma_FrAbs_N",
    "GENIEReWeight_SBN_v1_multisigma_FrPiProd_N",
    # NCEL
    "GENIEReWeight_SBN_v1_multisigma_MaNCEL",
    "GENIEReWeight_SBN_v1_multisigma_EtaNCEL",
]

MULTISIM_NUNIV = 250   # universes per flux / Geant4 multisim, and slim GENIE throws


def _kept_neutrinos(f, N):
    """Mask over the true neutrinos N: true vertex in the AV, or the truth
    match of a Wire-Cell candidate."""
    in_av = cuts.dist_to_active_boundary(N.position_x, N.position_y, N.position_z) > 0
    keep = pd.Series(in_av, index=N.index)
    if wc.has_wc(f) and "wc.T_match.level" in f[wc.WC_TREE].keys():
        tm = wc.make_wc_tmatchdf(f)
        matched = set(zip(tm.index.get_level_values("entry"), tm.tmatch_idx[tm.tmatch_idx >= 0]))
        keep |= pd.Series([(e, n) in matched for e, n in zip(N.index.get_level_values(0), N.index.get_level_values(1))],
                          index=N.index)
    return keep.to_numpy()


def make_sp_wgtdf(f, slim=False, wgt_types=("bnb", "genie", "g4")):
    if "globalTree" not in f:
        return pd.DataFrame()
    N = _load(f, "rec.mc.nu.", ["position.x", "position.y", "position.z"])
    if not len(N):
        return pd.DataFrame()
    # all neutrinos go in (getsyst's slim products assume the full list), then the kept rows are selected
    nuind = pd.Series(N.index.get_level_values(1), index=N.index)
    parts = []
    if "bnb" in wgt_types:
        parts.append(bnbsyst.bnbsyst(f, nuind, multisim_nuniv=MULTISIM_NUNIV, slim=slim))
    if "genie" in wgt_types:
        parts.append(geniesyst.geniesyst(f, nuind, multisim_nuniv=MULTISIM_NUNIV, slim=slim, systematics=GENIE_SYSTEMATICS))
    if "g4" in wgt_types:
        # g4syst.g4syst has no universe argument (getsyst's default, 100)
        parts.append(getsyst.getsyst(f, g4syst.g4_systematics, nuind, multisim_nuniv=MULTISIM_NUNIV))
    w = pd.concat(parts, axis=1).astype(np.float32)
    return w[_kept_neutrinos(f, N)]   # index (entry, rec.mc.nu..index), as the mcnu table


def make_sp_slimwgtdf(f):
    return make_sp_wgtdf(f, slim=True)
