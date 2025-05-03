import pandas as pd
import numpy as np
import sys

import util
import kinematics


# CUT VALUES
MUSEL_MUSCORE_TH = 25
MUSEL_PSCORE_TH = 100
MUSEL_LEN_TH = 50
PSEL_MUSCORE_TH = 0
PSEL_PSCORE_TH = 90
STUB_DQDXS = [5.5e5, 4e5, 3e5, 0]
DELP_TH = 0.25

# HELPER FUNCTIONS
def cut_all(cuts):
    ret = cuts[0]
    for c in cuts[1:]:
        sum(c)
        ret = ret & c
    return ret

def InFV(df, DETECTOR):
    return util.InFV(df, 50, det=DETECTOR)

def InBeam(t):
    return (t > 0.) & (t < 1.800)

def FVCut(df, DETECTOR):
    return InFV(df.slc.vertex, DETECTOR)

def PIDCut(df):
    # muon cut on muon candidates
    
    # TODO: use scores of all 3 planes
    # muon_chi2 = (Avg(df, "muon", drop_0=True) < MUSEL_MUSCORE_TH) & (Avg(df, "proton", drop_0=True) > MUSEL_PSCORE_TH)
    
    # TODO: used BDT scores
    # len_cut = (masterdf.len.squeeze() > MUSEL_LEN_TH)
    # dazzle_muon = (masterdf.dazzle.muonScore > 0.6)
    # muon_cut = (muon_chi2) & (len_cut | dazzle_muon)
    
    mu_score_cut = (df.mu.pfp.trk.chi2pid.I2.chi2_muon < MUSEL_MUSCORE_TH) & (df.mu.pfp.trk.chi2pid.I2.chi2_proton > MUSEL_PSCORE_TH)
    mu_len_cut = (df.mu.pfp.trk.len > MUSEL_LEN_TH)
    mu_cut = (mu_score_cut) & (mu_len_cut)
    
    # proton cut on proton candidates
    p_score_cut = (df.p.pfp.trk.chi2pid.I2.chi2_muon > PSEL_MUSCORE_TH) & (df.p.pfp.trk.chi2pid.I2.chi2_proton < PSEL_PSCORE_TH) 
    p_cut = p_score_cut
    
    # select slices with mu+p
    slc_mu_cut = mu_cut.groupby(level=[0,1,2]).any()
    slc_p_cut = p_cut.groupby(level=[0,1,2]).any()
    
    return slc_mu_cut & slc_p_cut

def StubCut(df):
    cut_list = []
    for i,l in enumerate(["0_5", "1", "2", "3"]):
        this_cut = np.invert((df.stub["l"+l+"cm"].Q/df.stub["l"+l+"cm"].length) > STUB_DQDXS[i])
        cut_list.append(this_cut)
    
    return cut_all(cut_list)

# transverse momentum cut
def dPtCut(df):
    return df.del_p < DELP_TH
    
def is_weightcol(c):
    return c.startswith("GENIEReWeight") or c.endswith("Flux")

# Reconstructed variables we want to save
def recodf(df):
    dPt = df.del_p.rename("dPt")
    caloE = kinematics.neutrino_energy(df.mu.pfp.trk.P.p_muon, df.mu.pfp.trk.dir, df.p.pfp.trk.P.p_proton, df.p.pfp.trk.dir).rename("caloE")
    # Lookup truth match
    # TODO: make this better, include cut on purity
    tmatch = df['rec.mc.nu..index'].fillna(-1).astype(int).rename("tmatch")

    return pd.concat([dPt, caloE, tmatch], axis=1) # .droplevel(0, 0)

def tmatchdf(df, mcdf, DETECTOR):
    # filter the columns from the mcdf we want
    # map old name to new name
    tosave = {
      "E": "trueE",
      "Q2": "trueQ2",
      "pdg": "truepdg",
      "iscc": "iscc",
      "genie_mode": "genie_mode",
      "np_50MeV": "np_50MeV",
      "np_20MeV": "np_20MeV",
    }

    def savecol(c):
        return c[0] in list(tosave.keys()) or is_weightcol(c[0])

    mcdf_cols = [c for c in mcdf.columns if savecol(c)]
    mcdf = mcdf[mcdf_cols]

    # Get the column depth matching on both sides
    mcdf.columns = pd.MultiIndex.from_tuples([(tosave[c[0]], c[1]) if c[0] in tosave else (c[0], c[1]) for c in mcdf.columns])
    df.columns = pd.MultiIndex.from_tuples([(c, "") for c in df.columns])

    tmatch_df = pd.merge(df, mcdf, how="left", left_on=["__ntuple", "entry", "tmatch"], right_index=True) 

    # Fill nan's for not-matching columns
    for c in tosave.values():
        tmatch_df[(c, "")] = tmatch_df[(c, "")].fillna(0.)

    # Systematic weights are all 1 for no truth match
    tmatch_df.fillna(1, inplace=True)

    # TODO: get actual baseline
    if DETECTOR == "SBND":
        tmatch_df[("baseline", "")] = 110.
        print("Using SBND Baseline of 110")
    elif DETECTOR == "ICARUS":
        tmatch_df[("baseline", "")] = 600.
        print("Using ICARUS Baseline of 600")
    else:
        print("No detector configured!")
        sys.exit()

    return tmatch_df

def main(f):
    if "icarus" in f:
        DETECTOR = "ICARUS"
    elif "sbnd" in f:
        DETECTOR = "SBND"
    df = pd.read_hdf(f, key="evt")
    print("len:", len(df))
    mcdf = pd.read_hdf(f, key="mcnu")
    hdrdf = pd.read_hdf(f, key="hdr")
    POT = sum(hdrdf.pot.tolist())
    # TODO: lookup in hdr df
    # about 1M neutrinos / 1e20 POT
    # POT = 1e20*(mcdf.index.size / 1e6)

    # For now, don't include dPt cut. Instead, we will bin it while fitting
    cuts = [FVCut(df, DETECTOR), PIDCut(df), StubCut(df), dPtCut(df)]

    # Apply selection
    seldf = df[cut_all(cuts)]

    # Get reconstructed variables we want
    ret = recodf(seldf)
    # Get reconstructed variables we want
    # Lookup truth info
    ret = tmatchdf(ret, mcdf, DETECTOR)
    # save
    savefile = f.split(".")[0] + ".gump.df"
    
    # Separately save the weights and the variables
    weightcol = [c for c in ret.columns if is_weightcol(c[0])]
    varcol = [c for c in ret.columns if not is_weightcol(c[0])]

    # Prune unneeded multi-indexes and save
    ret[varcol].droplevel(1,1).reset_index(drop=True).to_hdf(savefile, key="var")
    ret[weightcol].reset_index(drop=True).to_hdf(savefile, key="wgt")

    # Save POT
    pd.DataFrame([POT]).to_hdf(savefile, key="pot")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("USAGE: python gump_selection.py input.df")
        sys.exit(1)
    main(sys.argv[1])
