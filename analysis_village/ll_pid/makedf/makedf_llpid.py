from makedf.makedf import *
from makedf import llpid, chi2pid
from pyanalib.pandas_helpers import *
from makedf.util import *
import warnings
warnings.filterwarnings('ignore')

def make_llpiddf_calo_update(f):
    return make_llpiddf(f, updatecalo="cv")

def make_llpiddf(f, updatecalo=None, variation="cv", rrcorr=True, cutClearCosmic=False, tracks=None):
    hitdfs = [make_trkhitdf(f, plane) for plane in range(0, 3)]

    if tracks is not None:
        hitdfs = [df[df.index.droplevel(-1).isin(tracks)].copy() for df in hitdfs]

    if cutClearCosmic:
        is_clear_cosmic = loadbranches(f["recTree"], ["rec.slc.is_clear_cosmic"]).rec.slc.is_clear_cosmic
        keep = is_clear_cosmic.index[is_clear_cosmic == 0]
        hitdfs = [df[df.index.droplevel([-2, -1]).isin(keep)].copy() for df in hitdfs]

    dedxname = "dedx"
    if updatecalo is not None:
        det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det
        det = "SBND" if (1 == det.unique()) else "ICARUS"
        ismc = make_mchdrdf(f).ismc.iloc[0]
        for plane, df in enumerate(hitdfs):
            df["dedx_redo"] = chi2pid.dedx(df, gain=det, calibrate=det, plane=plane, isMC=ismc, new_calo_params=chi2pid.CALO_VARIATIONS[updatecalo])
        dedxname = "dedx_redo"

    lldf = llpid.likelihood(hitdfs, dedxname=dedxname, rrcorr=rrcorr, params=llpid.LLPID_VARIATIONS[variation])
    return lldf

def make_pandora_llpid_df_calo_update(f, **trkArgs):
    return make_pandora_llpid_df(f, trkScoreCut=False, trkDistCut=50., cutClearCosmic=True, requireFiducial=False, updatecalo=True, **trkArgs)

def make_pandora_llpid_df(f, updatecalo=False, variation="cv", rrcorr=True, **pandoraArgs):
    # make_pandora_df plus the llpid columns under pfp.trk.llpid, e.g. df.pfp.trk.llpid.ell_mu
    # updatecalo: redo dedx for the chi2pid and the llpid alike
    pandoradf = make_pandora_df(f, updatecalo=updatecalo, **pandoraArgs)

    # only fit the tracks that survived make_pandora_df's cuts
    lldf = make_llpiddf(f, updatecalo="cv" if updatecalo else None, variation=variation, rrcorr=rrcorr, tracks=pandoradf.index)

    # pfps with no usable hits get no llpid row: NaN likelihoods, 0 hits, False flags
    flags = lldf.columns[lldf.dtypes == bool]
    nhits = [c for c in lldf.columns if c.startswith("nhits_")]
    lldf = lldf.reindex(pandoradf.index)
    lldf[flags] = lldf[flags].fillna(False).astype(bool)
    lldf[nhits] = lldf[nhits].fillna(0.)
    lldf["version"] = llpid.llpid_version   

    lldf.columns = pd.MultiIndex.from_tuples([("pfp", "trk", "llpid", c) for c in lldf.columns])
    pandoradf = multicol_merge(pandoradf, lldf, left_index=True, right_index=True, how="left", validate="one_to_one")
    return pandoradf


