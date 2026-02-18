from pyanalib.pandas_helpers import *
import pyanalib.calo_helpers as caloh
from .branches import *
from .util import *
from .calo import *
from .chi2pid import SBND_CALO_PARAMS
from . import numisyst, g4syst, geniesyst, bnbsyst

PDG = {
    "muon": [13, "muon", 0.105,],
    "proton": [2212, "proton", 0.938272,],
    "neutron": [2112, "neutron", 0.9395654,],
    "pizero": [111, "pizero", 0.1349768],
    "pipm": [211, "piplus", 0.13957039],
    "argon": [1000180400, "argon", (18*0.938272 + 22*0.9395654)],
    "gamma": [22, "gamma", 0 ],
    "lambda": [3122, "lambda", 1.115683],
    "kaon_p": [321, "kaon_p",  0.493677],
    "sigma_p": [3222, "sigma_p", 1.18936],
    "kaon_0": [311, "kaon_0", 0.497648],
    "sigma_0": [3212, "sigma_0", 1.19246],
    "lambda_p_c": [4122, "lambda_p_c", 2.28646],
    "sigma_pp_c": [4222, "sigma_pp_c", 2.45397],
    "electron": [11, "electron", 0.510998950],
    "sigma_p_c": [4212, "sigma_p_c", 2.4529],
}

## == For additional column in mcdf with primary particle multiplicities
## ==== "<column name>": ["<particle name>", <KE cut in GeV>]
## ==== <particle name> is used to collect PID and mass from the "PDG" dictionary
TRUE_KE_THRESHOLDS = {"nmu_27MeV": ["muon", 0.027],
                      "np_20MeV": ["proton", 0.02],
                      "nmu_25MeV": ["muon", 0.025],
                      "nmu_100MeV": ["muon", 0.1],
                      "np_50MeV": ["proton", 0.05],
                      "nn_0MeV": ["neutron", 0.0],
                      "ne_25MeV": ["electron", 0.025],
                      "npi_25MeV": ["pipm", 0.025],
                      "npi0_0MeV": ["pizero", 0.0],
                      "nk_25MeV": ["kaon_p", 0.025],
                      }

def make_hdrdf(f):
    hdr = loadbranches(f["recTree"], hdrbranches).rec.hdr
    return hdr

def make_mchdrdf(f):
    hdr = loadbranches(f["recTree"], mchdrbranches).rec.hdr
    return hdr

def make_potdf_bnb(f):
    pot = loadbranches(f["recTree"], bnbpotbranches).rec.hdr.bnbinfo
    return pot

def make_potdf_numi(f):
    pot = loadbranches(f["recTree"], numipotbranches).rec.hdr.numiinfo
    return pot

def make_mcnuwgtdf(f):
    return make_mcnudf(f, include_weights=True, multisim_nuniv=1000)

def make_mcnuwgtdf_slim(f):
    return make_mcnudf(f, include_weights=True, multisim_nuniv=1000, slim=True)

def make_mcnudf(f, include_weights=False, multisim_nuniv=250, wgt_types=["bnb","genie"], slim=False, **kwargs):
    # ----- sbnd or icarus? -----
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det
    if (1 == det.unique()):
        det = "SBND"
    else:
        det = "ICARUS"

    mcdf = make_mcdf(f,**kwargs)
    mcdf["ind"] = mcdf.index.get_level_values(1)
    if include_weights:
        if len(wgt_types) == 0:
            print("include_weights is set to True, pass at least one type of wgt to save")

        else:
            if det == "ICARUS":
                wgtdf = pd.concat([numisyst.numisyst(mcdf.pdg, mcdf.E), geniesyst.geniesyst(f, mcdf.ind), g4syst.g4syst(f, mcdf.ind)], axis=1)
            elif det == "SBND":
                df_list = []
                if "bnb" in wgt_types:
                    bnbwgtdf = bnbsyst.bnbsyst(f, mcdf.ind, multisim_nuniv=multisim_nuniv, slim=slim)
                    df_list.append(bnbwgtdf)
                if "genie" in wgt_types:
                    geniewgtdf = geniesyst.geniesyst(f, mcdf.ind, slim=slim)
                    df_list.append(geniewgtdf)
                if "g4" in wgt_types:
                    g4wgtdf = g4syst.g4syst(f, mcdf.ind, slim=slim)
                    df_list.append(g4wgtdf)
                wgtdf = pd.concat(df_list, axis=1)
            mcdf = multicol_concat(mcdf, wgtdf)
    return mcdf

def make_mchdf(f, include_weights=False):
    mcdf = loadbranches(f["recTree"], mchbranches).rec.mc.prtl
    if include_weights:
        wgtdf = numisyst.numisyst(14, mcdf.E) # TODO: what PDG?
        mcdf = pd.concat([mcdf, wgtdf], axis=1)
    return mcdf

def make_crtspdf(f):
    crtspdf = loadbranches(f["recTree"], crtspbranches).rec
    return crtspdf

def make_opflashdf(f):
    opflashdf = loadbranches(f["recTree"], opflashbranches).rec.opflashes
    return opflashdf

def make_histpotdf(f):
    # get the value from the TotalPOT histogram
    pot = f['TotalPOT'].values()
    histpot = pd.DataFrame(data={'TotalPOT':pot})
    histpot.index.name = 'entry'
    return histpot

def make_histgenevtdf(f):
    # get the value from the TotalPOT histogram
    genevt = f['TotalGenEvents'].values()
    histgenevt = pd.DataFrame(data={'TotalGenEvents':genevt})
    histgenevt.index.name = 'entry'
    return histgenevt

def make_trkdf(f, scoreCut=False, requiret0=False, requireCosmic=False, mcs=False, trackScoreCut=0.5, updaterecomb=False):
    trkdf = loadbranches(f["recTree"], trkbranches)
    if scoreCut:
        trkdf = trkdf.rec.slc.reco[trkdf.rec.slc.reco.pfp.trackScore > trackScoreCut]
    else:
        trkdf = trkdf.rec.slc.reco

    if requiret0:
        trkdf = trkdf[~np.isnan(trkdf.pfp.t0)]

    if requireCosmic:
        trkdf = trkdf[trkdf.pfp.parent == -1]

    if mcs:
        mcsdf = loadbranches(f["recTree"], [trkmcsbranches[0]]).rec.slc.reco.pfp.trk.mcsP
        mcsdf_angle = loadbranches(f["recTree"], [trkmcsbranches[1]]).rec.slc.reco.pfp.trk.mcsP
        mcsdf_angle.index.set_names(mcsdf.index.names, inplace=True)

        mcsdf = mcsdf.merge(mcsdf_angle, how="left", left_index=True, right_index=True)
        mcsgroup = list(range(mcsdf.index.nlevels-1))
        cumlen = mcsdf.seg_length.groupby(level=mcsgroup).cumsum()*14 # convert rad length to cm
        maxlen = (cumlen*(mcsdf.seg_scatter_angles >= 0)).groupby(level=mcsgroup).max()
        trkdf[("pfp", "trk", "mcsP", "len", "", "")] = maxlen
    
    # Just for muons and just on the collection plane
    if updaterecomb:
        calo_params = SBND_CALO_PARAMS.copy()
        calo_params_recomb = SBND_CALO_PARAMS.copy()
        hdrdf = make_hdrdf(f)
        ismc = hdrdf.ismc.iloc[0]
        idx = 0 if ismc else 1
        # Map uncertainty dict keys to actual SBND_CALO_PARAMS keys
        unc_dict = {
            "alpha_emb": 0.008,
            "beta_90": 0.008,  # Fixed: was "beta_emb"
            "R_emb": 0.02,  # Fixed: was "R_90"
            "c_cal_frac": 0.02, #Assign 2% uncertainty on c_cal_frac
        }
        sign_tag = {-1: "m", 0: "0", 1: "p"}  # Define sign_tag mapping
        recomb_keys = list(unc_dict.keys())
        for recomb_key in recomb_keys:
            for plane in range(2, 3): #only collection
                hitdf = make_trkhitdf(f, plane)
                trk_keys = trkdf.index.unique()
                hit_keys3 = hitdf.index.droplevel(-1)  # -1 since hitdf has an additional index
                mask_match = hit_keys3.isin(trk_keys)
                hitdf = hitdf[mask_match]
                this_etau = calo_params_recomb["etau"][idx]

                for sig in [-1, 0, 1]:
                    #Only keep alpha_emb00, no other null variations
                    if sig == 0 and recomb_key != "alpha_emb":
                        continue
                    # Select MC (idx=0) or data (idx=1) value, then add uncertainty
                    if recomb_key == "c_cal_frac":
                        base_value = calo_params[recomb_key][plane]
                        perturbed_value = base_value + sig * unc_dict[recomb_key]
                    else:
                        base_value = calo_params[recomb_key][idx]
                        perturbed_value = base_value + sig * unc_dict[recomb_key]

                    # Update the parameter list (keep both MC and data, but use perturbed value for current idx)
                    calo_params_recomb[recomb_key] = calo_params[recomb_key].copy()
                    if recomb_key == "c_cal_frac":
                        calo_params_recomb[recomb_key][plane] = perturbed_value
                    else:
                        calo_params_recomb[recomb_key][idx] = perturbed_value
                    
                    recomb_tag = "_{}{}{}".format(recomb_key, sign_tag[sig], np.abs(sig))

                    # Pass scalar values to new_dedx (not lists)
                    new_dedx = caloh.new_dedx(hitdf, calo_params_recomb["c_cal_frac"][plane], plane, 
                                             calo_params_recomb["alpha_emb"][idx],
                                             calo_params_recomb["beta_90"][idx],
                                             calo_params_recomb["R_emb"][idx],
                                             this_etau, ismc)
                    hitdf[('dedx', '')] = new_dedx

                    for par in ['muon', 'proton']:
                        this_chi2_new = hitdf.groupby(level=['entry', 'rec.slc..index', 'rec.slc.reco.pfp..index']).apply(lambda group: caloh.calculate_chi2_for_entry(group, par))
                        this_chi2_new = this_chi2_new.apply(pd.Series).rename(columns={0: "chi2", 1: "ndof"})
                        this_chi2_col = ('pfp', 'trk', 'chi2pid', 'I' + str(plane), 'chi2_' + par + recomb_tag, '')
                        this_ndof_col = ('pfp', 'trk', 'chi2pid', 'I' + str(plane), 'ndof_' + par + recomb_tag, '')
                        try:
                            trkdf[this_chi2_col] = this_chi2_new.chi2
                            trkdf[this_ndof_col] = this_chi2_new.ndof
                            trkdf[this_chi2_col] = trkdf[this_chi2_col].fillna(0.)
                            trkdf[this_ndof_col] = trkdf[this_ndof_col].fillna(0)
                        except:
                            print("no selected tracks")


    trkdf[("pfp", "tindex", "", "", "", "")] = trkdf.index.get_level_values(2)

    return trkdf

def make_trkhitdf(f,plane):
    df = loadbranches(f["recTree"], trkhitbranches)[("rec","slc","reco","pfp","trk","calo",f"I{plane}","points")]

    # Firsthit and Lasthit info
    ihit = df.index.get_level_values(-1)
    df["firsthit"] = ihit == 0

    lasthit = df.groupby(level=list(range(df.index.nlevels-1))).tail(1).copy()
    lasthit["lasthit"] = True
    df["lasthit"] = lasthit.lasthit
    df.lasthit = df.lasthit.fillna(False)

    return df

def make_slcdf(f, **kwargs):
    slcdf = loadbranches(f["recTree"], slcbranches)
    slcdf = slcdf.rec

    slc_mcdf = make_mcdf(f, slc_mcbranches, slc_mcprimbranches, **kwargs)
    slc_mcdf.columns = pd.MultiIndex.from_tuples([tuple(["slc", "truth"] + list(c)) for c in slc_mcdf.columns])
    slcdf = multicol_merge(slcdf, slc_mcdf, left_index=True, right_index=True, how="left", validate="one_to_one")

    return slcdf

def make_mcdf(f, branches=mcbranches, primbranches=mcprimbranches, include_mu=True, include_p=True, include_pi=True):
    # load the df
    mcdf = loadbranches(f["recTree"], branches)
    while mcdf.columns.nlevels > 2:
        mcdf.columns = mcdf.columns.droplevel(0)
    # Add in primary particle info
    mcprimdf = loadbranches(f["recTree"], primbranches)
    while mcprimdf.columns.nlevels > 2:
        mcprimdf.columns = mcprimdf.columns.droplevel(0)

    mcprimdf.index = mcprimdf.index.rename(mcdf.index.names[:2] + mcprimdf.index.names[2:])

    max_proton_KE = mcprimdf[np.abs(mcprimdf.pdg)==PDG["proton"][0]].genE.groupby(level=[0,1]).max() - PDG["proton"][2]
    mcdf = multicol_add(mcdf, max_proton_KE.rename("max_proton_ke"), default=0.)

    mcdf.max_proton_ke = mcdf.max_proton_ke.fillna(0.)

    # particle counts
    mcdf = multicol_add(mcdf, (np.abs(mcprimdf.pdg)==2112).groupby(level=[0,1]).sum().rename("nn"))
    mcdf = multicol_add(mcdf, (np.abs(mcprimdf.pdg)==2212).groupby(level=[0,1]).sum().rename("np"))
    mcdf = multicol_add(mcdf, (np.abs(mcprimdf.pdg)==13).groupby(level=[0,1]).sum().rename("nmu"))
    mcdf = multicol_add(mcdf, (np.abs(mcprimdf.pdg)==11).groupby(level=[0,1]).sum().rename("ne"))
    mcdf = multicol_add(mcdf, (np.abs(mcprimdf.pdg)==211).groupby(level=[0,1]).sum().rename("npi"))
    mcdf = multicol_add(mcdf, (np.abs(mcprimdf.pdg)==111).groupby(level=[0,1]).sum().rename("npi0"))
    mcdf = multicol_add(mcdf, (np.abs(mcprimdf.pdg)==22).groupby(level=[0,1]).sum().rename("ng"))
    mcdf = multicol_add(mcdf, (np.abs(mcprimdf.pdg)==321).groupby(level=[0,1]).sum().rename("nk"))
    mcdf = multicol_add(mcdf, (np.abs(mcprimdf.pdg)==310).groupby(level=[0,1]).sum().rename("nk0"))
    mcdf = multicol_add(mcdf, (np.abs(mcprimdf.pdg)==3112).groupby(level=[0,1]).sum().rename("nsm"))
    mcdf = multicol_add(mcdf, (np.abs(mcprimdf.pdg)==3222).groupby(level=[0,1]).sum().rename("nsp"))

    # particle counts w/ threshold
    for identifier, (particle, threshold) in TRUE_KE_THRESHOLDS.items():
        this_KE = mcprimdf[np.abs(mcprimdf.pdg)==PDG[particle][0]].genE - PDG[particle][2]
        mcdf = multicol_add(mcdf, ((np.abs(mcprimdf.pdg)==PDG[particle][0]) & (this_KE > threshold)).groupby(level=[0,1]).sum().rename(identifier))
 
    if include_mu:
        # lepton info
        mudf = mcprimdf[np.abs(mcprimdf.pdg)==13].sort_values(mcprimdf.index.names[:2] + [("genE", "")]).groupby(level=[0,1]).last()
        mudf.columns = pd.MultiIndex.from_tuples([tuple(["mu"] + list(c)) for c in mudf.columns])
        mcdf = multicol_merge(mcdf, mudf, left_index=True, right_index=True, how="left", validate="one_to_one")
        # primary track variables
        mcdf.loc[:, ('mu','totp','')] = np.sqrt(mcdf.mu.genp.x**2 + mcdf.mu.genp.y**2 + mcdf.mu.genp.z**2)
        # opening angles
        mcdf.loc[:, ('mu','dir','x')] = mcdf.mu.genp.x/mcdf.mu.totp
        mcdf.loc[:, ('mu','dir','y')] = mcdf.mu.genp.y/mcdf.mu.totp
        mcdf.loc[:, ('mu','dir','z')] = mcdf.mu.genp.z/mcdf.mu.totp
        # Is contained
        mcdf.loc[:, ('mu','is_contained','')] = InFV(mcdf.mu.start, 0, det="SBND AV") & InFV(mcdf.mu.end, 0, det="SBND AV")

    if include_pi:
        cpidf = mcprimdf[np.abs(mcprimdf.pdg)==211].sort_values(mcprimdf.index.names[:2] + [("genE", "")]).groupby(level=[0,1]).last()
        cpidf.columns = pd.MultiIndex.from_tuples([tuple(["cpi"] + list(c)) for c in cpidf.columns])
        mcdf = multicol_merge(mcdf, cpidf, left_index=True, right_index=True, how="left", validate="one_to_one")

    if include_p:
        pdf = mcprimdf[mcprimdf.pdg==2212].sort_values(mcprimdf.index.names[:2] + [("genE", "")]).groupby(level=[0,1]).last()
        pdf.columns = pd.MultiIndex.from_tuples([tuple(["p"] + list(c)) for c in pdf.columns])

        mcdf = multicol_merge(mcdf, pdf, left_index=True, right_index=True, how="left", validate="one_to_one")

        mcdf.loc[:, ('p','totp','')] = np.sqrt(mcdf.p.genp.x**2 + mcdf.p.genp.y**2 + mcdf.p.genp.z**2)

        mcdf.loc[:, ('p','dir','x')] = mcdf.p.genp.x/mcdf.p.totp
        mcdf.loc[:, ('p','dir','y')] = mcdf.p.genp.y/mcdf.p.totp
        mcdf.loc[:, ('p','dir','z')] = mcdf.p.genp.z/mcdf.p.totp

    return mcdf

def make_mcprimdf(f):
    mcprimdf = loadbranches(f["recTree"], mcprimbranches)
    return mcprimdf

def make_pandora_df(f, trkScoreCut=False, trkDistCut=10., cutClearCosmic=False, requireFiducial=False, **trkArgs):
    # load
    trkdf = make_trkdf(f, trkScoreCut, **trkArgs)
    slcdf = make_slcdf(f)

    # merge in tracks
    slcdf = multicol_merge(slcdf, trkdf, left_index=True, right_index=True, how="right", validate="one_to_many")

    # distance from vertex to track start
    slcdf = multicol_add(slcdf, dmagdf(slcdf.slc.vertex, slcdf.pfp.trk.start).rename(("pfp", "dist_to_vertex")))

    if trkDistCut > 0:
        slcdf = slcdf[slcdf.pfp.dist_to_vertex < trkDistCut]
    if cutClearCosmic:
        slcdf = slcdf[slcdf.slc.is_clear_cosmic==0]
    # require fiducial verex
    if requireFiducial:
        slcdf = slcdf[InFV(slcdf.slc.vertex, 50)]

    return slcdf

def make_spine_df(f, trkDistCut=-1, requireFiducial=True, **trkArgs):
    # load
    partdf = make_spinepartdf(f, **trkArgs)
    partdf.columns = pd.MultiIndex.from_tuples([tuple(["particle"] + list(c)) for c in partdf.columns])
    eslcdf = make_spineslcdf(f)

    # merge in tracks
    eslcdf = multicol_merge(eslcdf, partdf, left_index=True, right_index=True, how="right", validate="one_to_many")
    eslcdf = multicol_add(eslcdf, dmagdf(eslcdf.vertex, eslcdf.particle.start_point).rename("dist_to_vertex"))

    if trkDistCut > 0:
        eslcdf = eslcdf[eslcdf.dist_to_vertex < trkDistCut]
    # require fiducial verex
    if requireFiducial:
        eslcdf = eslcdf[InFV(eslcdf.vertex, 50)]

    return eslcdf

def make_stubs(f, det="ICARUS"):
    stubdf = loadbranches(f["recTree"], stubbranches)
    stubdf = stubdf.rec.slc.reco.stub

    stubpdf = loadbranches(f["recTree"], stubplanebranches)
    stubpdf = stubpdf.rec.slc.reco.stub.planes

    stubdf["nplane"] = stubpdf.groupby(level=[0,1,2]).size()
    stubdf["plane"] = stubpdf.p.groupby(level=[0,1,2]).first()

    stubhitdf = loadbranches(f["recTree"], stubhitbranches)
    stubhitdf = stubhitdf.rec.slc.reco.stub.planes.hits

    stubhitdf = stubhitdf.join(stubpdf)
    stubhitdf = stubhitdf.join(stubdf.efield_vtx)
    stubhitdf = stubhitdf.join(stubdf.efield_end)

    hdrdf = make_mchdrdf(f)
    def dEdx2dQdx(dEdx): # MC parameters
        return recombination_sbnd(dEdx, np.pi/2) if det == "SBND" else recombination_icarus(dEdx, np.pi/2)

    MIP_dqdx = dEdx2dQdx(1.7) 

    stub_end_charge = stubhitdf.charge[stubhitdf.wire == stubhitdf.hit_w].groupby(level=[0,1,2,3]).first().groupby(level=[0,1,2]).first()
    stub_end_charge.name = ("endp_charge", "", "")

    stub_pitch = stubpdf.pitch.groupby(level=[0,1,2]).first()
    stub_pitch.name = ("pitch", "", "")

    stubdir_is_pos = (stubhitdf.hit_w - stubhitdf.vtx_w) > 0.
    when_sum = ((stubhitdf.wire > stubhitdf.vtx_w) == stubdir_is_pos) & (((stubhitdf.wire < stubhitdf.hit_w) == stubdir_is_pos) | (stubhitdf.wire == stubhitdf.hit_w)) 
    stubcharge = (stubhitdf.charge[when_sum]).groupby(level=[0,1,2,3]).sum().groupby(level=[0,1,2]).first()
    stubcharge.name = ("charge", "", "")

    stubinccharge = (stubhitdf.charge).groupby(level=[0,1,2,3]).sum().groupby(level=[0,1,2]).first()
    stubinccharge.name = ("inc_charge", "", "")

    hit_before_start = ((stubhitdf.wire < stubhitdf.vtx_w) == stubdir_is_pos)
    stub_inc_sub_charge = (stubhitdf.charge - MIP_dqdx*stubhitdf.ontrack*(~hit_before_start)*stubhitdf.trkpitch).groupby(level=[0,1,2,3]).sum().groupby(level=[0,1,2]).first()
    stub_inc_sub_charge.name = ("inc_sub_charge", "", "")

    stubdf = stubdf.join(stubcharge)
    stubdf = stubdf.join(stubinccharge)
    stubdf = stubdf.join(stub_inc_sub_charge)
    stubdf = stubdf.join(stub_end_charge)
    stubdf = stubdf.join(stub_pitch)
    stubdf["length"] = magdf(stubdf.vtx - stubdf.end)
    stubdf["Q"] = stubdf.inc_sub_charge
    stubdf["truth_pdg"] = stubdf.truth.p.pdg
    stubdf["truth_interaction_id"] = stubdf.truth.p.interaction_id 
    stubdf["truth_gen_E"] = stubdf.truth.p.genE 

    # TODO: convert charge to energy
    stubdf["ke"] = np.nan # Q2KE(stubdf.Q)
    # TODO: also do calorimetric variations
    stubdf["ke_callo"] = np.nan # Q2KE_mc_callo(stubdf.Q)
    stubdf["ke_calhi"] = np.nan # Q2KE_mc_calhi(stubdf.Q)

    stubdf.ke = stubdf.ke.fillna(0)
    stubdf.Q = stubdf.Q.fillna(0)

    stubdf["dedx"] = stubdf.ke / stubdf.length

    stubdf["dedx_callo"] = stubdf.ke_callo / stubdf.length
    stubdf["dedx_calhi"] = stubdf.ke_calhi / stubdf.length

    dqdx = stubdf.inc_sub_charge / stubdf.length
    length = stubdf.length
    hasstub = (length < 4.) & \
        (((length > 0.) & (dqdx > 5.5e5)) |\
        ((length > 0.5) & (dqdx > 3.5e5)) |\
        ((length > 1) & (dqdx > 3e5)) |\
        ((length > 2) & (dqdx > 2e5)))

    stubdf["dqdx"] = dqdx 
    stubdf['pass_proton_stub'] = hasstub
    return stubdf

    ## It seems there is a bug. First stub in each length is included for a slice...
    # only take collection plane
    #stubdf = stubdf[stubdf.plane == 2]

    #stub_length_bins = [0, 0.5, 1, 2, 3, 4]
    #stub_length_name = ["l0_5cm", "l1cm", "l2cm", "l3cm", "l4cm"]
    #tosave = ["dedx", "dedx_callo", "dedx_calhi", "Q", "length", "charge", "inc_charge"]

    #df_tosave = []
    #for blo, bhi, name in zip(stub_length_bins[:-1], stub_length_bins[1:], stub_length_name):
    #    stub_tosave = stubdf.dedx[(stubdf.length > blo) & (stubdf.length < bhi)].groupby(level=[0,1]).idxmax()
    #    for col in tosave:
    #        s = stubdf.loc[stub_tosave, col]
    #        s.name = ("stub", name, col, "", "", "")
    #        s.index = s.index.droplevel(-1)
    #        df_tosave.append(s)

    #return pd.concat(df_tosave, axis=1)

def make_spineslcdf(f):
    eslcdf = loadbranches(f["recTree"], eslcbranches)
    eslcdf = eslcdf.rec.dlp

    etintdf = loadbranches(f["recTree"], etruthintbranches)
    etintdf = etintdf.rec.dlp_true
    
    # match to the truth info
    mcdf = make_mcdf(f)
    # mc is truth
    mcdf.columns = pd.MultiIndex.from_tuples([tuple(["truth"] + list(c)) for c in mcdf.columns])

    # Do matching
    # 
    # First get the ML true particle IDs matched to each reco particle
    eslc_matchdf = loadbranches(f["recTree"], eslcmatchedbranches)
    eslc_match_overlap_df = loadbranches(f["recTree"], eslcmatchovrlpbranches)
    eslc_match_overlap_df.index.names = eslc_matchdf.index.names

    eslc_matchdf = multicol_merge(eslc_matchdf, eslc_match_overlap_df, left_index=True, right_index=True, how="left", validate="one_to_one")
    eslc_matchdf = eslc_matchdf.rec.dlp

    # Then use bestmatch.match to get the nu ids in etintdf
    eslc_matchdf_wids = pd.merge(eslc_matchdf, etintdf, left_on=["entry", "match"], right_on=["entry", "id"], how="left")
    eslc_matchdf_wids.index = eslc_matchdf.index

    # Now use nu_ids to get the true interaction information
    eslc_matchdf_trueints = multicol_merge(eslc_matchdf_wids, mcdf, left_on=["entry", "nu_id"], right_index=True, how="left")
    eslc_matchdf_trueints.index = eslc_matchdf_wids.index

    # delete unnecesary matching branches
    del eslc_matchdf_trueints[("match", "")]
    del eslc_matchdf_trueints[("nu_id", "")]
    del eslc_matchdf_trueints[("id", "")]

    # first match is best match
    bestmatch = eslc_matchdf_trueints.groupby(level=list(range(eslc_matchdf_trueints.index.nlevels-1))).first()

    # add extra levels to eslcdf columns
    eslcdf.columns = pd.MultiIndex.from_tuples([tuple(list(c) + [""]*2) for c in eslcdf.columns])

    eslcdf_withmc = multicol_merge(eslcdf, bestmatch, left_index=True, right_index=True, how="left")

    # Fix position names (I0, I1, I2) -> (x, y, z)
    def mappos(s):
        if s == "I0": return "x"
        if s == "I1": return "y"
        if s == "I2": return "z"
        return s
    def fixpos(c):
        if c[0] not in ["end_point", "start_point", "start_dir", "vertex", "momentum"]: return c
        return tuple([c[0]] + [mappos(c[1])] + list(c[2:]))

    eslcdf_withmc.columns = pd.MultiIndex.from_tuples([fixpos(c) for c in eslcdf_withmc.columns])

    return eslcdf_withmc

def make_spinepartdf(f):
    epartdf = loadbranches(f["recTree"], eparticlebranches)
    epartdf = epartdf.rec.dlp.particles

    tpartdf = loadbranches(f["recTree"], trueparticlebranches)
    tpartdf = tpartdf.rec.true_particles
    # cut out EMShowerDaughters
    # tpartdf = tpartdf[(tpartdf.parent == 0)]

    etpartdf = loadbranches(f["recTree"], etrueparticlebranches)
    etpartdf = etpartdf.rec.dlp_true.particles
    etpartdf.columns = [s for s in etpartdf.columns]
    
    # Do matching
    # 
    # First get the ML true particle IDs matched to each reco particle
    epart_matchdf = loadbranches(f["recTree"], eparticlematchedbranches)
    epart_match_overlap_df = loadbranches(f["recTree"], eparticlematchovrlpbranches)
    epart_match_overlap_df.index.names = epart_matchdf.index.names
    epart_matchdf = multicol_merge(epart_matchdf, epart_match_overlap_df, left_index=True, right_index=True, how="left", validate="one_to_one")
    epart_matchdf = epart_matchdf.rec.dlp.particles
    # get the best match (highest match_overlap), assume it's sorted
    bestmatch = epart_matchdf.groupby(level=list(range(epart_matchdf.index.nlevels-1))).first()
    bestmatch.columns = [s for s in bestmatch.columns]

    # Then use betmatch.match to get the G4 track IDs in etpartdf
    bestmatch_wids = pd.merge(bestmatch, etpartdf, left_on=["entry", "match"], right_on=["entry", "id"], how="left")
    bestmatch_wids.index = bestmatch.index

    # Now use the G4 track IDs to get the true particle information
    bestmatch_trueparticles = multicol_merge(bestmatch_wids, tpartdf, left_on=["entry", "track_id"], right_on=["entry", ("G4ID", "")], how="left")
    bestmatch_trueparticles.index = bestmatch_wids.index

    # delete unnecesary matching branches
    del bestmatch_trueparticles[("match", "")]
    del bestmatch_trueparticles[("track_id", "")]
    del bestmatch_trueparticles[("id", "")]

    # add extra level to epartdf columns
    epartdf.columns = pd.MultiIndex.from_tuples([tuple(list(c) + [""]) for c in epartdf.columns])

    # put everything in epartdf
    for c in bestmatch_trueparticles.columns:
        epartdf[tuple(["truth"] + list(c))] = bestmatch_trueparticles[c]

    # Fix position names (I0, I1, I2) -> (x, y, z)
    def mappos(s):
        if s == "I0": return "x"
        if s == "I1": return "y"
        if s == "I2": return "z"
        return s
    def fixpos(c):
        if c[0] not in ["end_point", "start_point", "start_dir", "vertex"]: return c
        return tuple([c[0]] + [mappos(c[1])] + list(c[2:]))

    epartdf.columns = pd.MultiIndex.from_tuples([fixpos(c) for c in epartdf.columns])

    return epartdf
