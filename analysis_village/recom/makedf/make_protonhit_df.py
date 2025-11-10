from makedf.makedf import *
from pyanalib.pandas_helpers import *
from makedf.util import *

pd.options.mode.chained_assignment = None

###############################################
#### slice and pfp level selections
###############################################

def InFV_nohiyz(data):
    xmin = 10.
    xmax = 190.
    zmin = 10.
    zmax = 450.
    ymax_highz = 100.
    pass_xz = (np.abs(data.x) > xmin) & (np.abs(data.x) < xmax) & (data.z > zmin) & (data.z < zmax)
    pass_y = ((data.z < 250) & (np.abs(data.y) < 190.)) | ((data.z > 250) & (data.y > -190.) & (data.y < ymax_highz))
    return pass_xz & pass_y

def InFV_nohiyz_trk(data):
    xmax = 190.
    zmin = 10.
    zmax = 450.
    ymax_highz = 100.
    pass_xz = (np.abs(data.x) < xmax) & (data.z > zmin) & (data.z < zmax)
    pass_y = ((data.z < 250) & (np.abs(data.y) < 190.)) | ((data.z > 250) & (data.y > -190.) & (data.y < ymax_highz))
    return pass_xz & pass_y

def pass_slc_with_n_pfps(df, n = 2):
    group_levels = ['entry', 'rec.slc..index']

    # Count how many pfps per slice
    pfp_counts = df.groupby(level=group_levels).size()

    # Get only slices with at least 2 2pfps
    valid_slices = pfp_counts[pfp_counts == n].index

    # Apply the mask to original
    filtered_df = df.loc[df.index.droplevel('rec.slc.reco.pfp..index').isin(valid_slices)]

    return filtered_df

def select_topology(df, nprong=2, max_len_cut=25, min_len_cut=25, contained=True):

    npfps = df.groupby(level=[0,1]).count().slc.producer

    ## exact number of pfps
    twopfps_idx = npfps[npfps == nprong].index
    twopfps = df.reset_index(level=[2]).loc[twopfps_idx]
    twopfps = twopfps.reset_index().set_index(["entry", "rec.slc..index", "rec.slc.reco.pfp..index"])

    ## exact number of pfps passing trk score cut in contained cut
    twopfps["is_track"] = (twopfps.pfp.trackScore > 0.5)
    twopfps["is_contained_track"] = (twopfps.pfp.trackScore > 0.5) & (InFV_nohiyz_trk(twopfps.pfp.trk.start) & InFV_nohiyz_trk(twopfps.pfp.trk.end))
    if contained:
        ntrks = twopfps.groupby(level=[0,1]).is_contained_track.sum()
    else:
        ntrks = twopfps.groupby(level=[0,1]).is_track.sum()

    ## exact number of pfpfs passing len cut
    nprong_slc_idx = ntrks[(ntrks == nprong)].index
    nprong_df = df.reset_index().set_index(["entry", "rec.slc..index"]).loc[nprong_slc_idx]
    nprong_df = nprong_df.reset_index().set_index(["entry", "rec.slc..index", "rec.slc.reco.pfp..index"])

    nprong_df['is_len'] = nprong_df.pfp.trk.len > 10. ## default len cut value based on 1-prong selection
    if nprong == 2:
        nprong_df['is_len'] = nprong_df.pfp.trk.len > 25.
    ntrks_len = nprong_df.groupby(level=[0,1]).is_len.sum()

    nprong_len_slc_idx = ntrks_len[(ntrks_len == nprong)].index
    nprong_len_df = df.reset_index().set_index(["entry", "rec.slc..index"]).loc[nprong_len_slc_idx]
    nprong_len_df = nprong_len_df.reset_index().set_index(["entry", "rec.slc..index", "rec.slc.reco.pfp..index"])

    longer_df = nprong_len_df.sort_values(by=("pfp", "trk","len","", "", ""), ascending=False).groupby(level=[0,1]).nth(0)
    shorter_df = nprong_len_df.sort_values(by=("pfp", "trk", "len", "", "", ""), ascending=True).groupby(level=[0,1]).nth(0)
    return longer_df, shorter_df

def select_calorimetry(df, muscore_cut_val=20, pscore_cut_val=80, pid=2212):
    if pid == 2212:
        muscore_cut = (df.pfp.trk.chi2pid.I2.chi2_muon > muscore_cut_val)
        pscore_cut = (df.pfp.trk.chi2pid.I2.chi2_proton < pscore_cut_val)

    elif pid == 13:
        muscore_cut = (df.pfp.trk.chi2pid.I2.chi2_muon < muscore_cut_val)
        pscore_cut = (df.pfp.trk.chi2pid.I2.chi2_proton > pscore_cut_val)

    score_selected = df[muscore_cut & pscore_cut]

    return score_selected

def check_flipped(hitdf, plane=2):
    trk_hitdfs = []
    hitdf["flipped"] = False

    trks = hitdf.reset_index(level=[3]).index.unique()
    for pidx in range(len(trks)):
        this_trk_hits = hitdf.reset_index(level=[3]).loc[trks[pidx]]

        if len(this_trk_hits) == 0:
            continue

        first_half = this_trk_hits.iloc[:len(this_trk_hits)//2]
        last_half = this_trk_hits.iloc[len(this_trk_hits)//2:]
        first_avg = np.mean(first_half[("dqdx")])
        last_avg = np.mean(last_half[("dqdx")])
        if first_avg > last_avg:
            max_rr = this_trk_hits[("rr")].max()
            this_trk_hits[("rr")] = max_rr - this_trk_hits[("rr")]
            this_trk_hits["flipped"] = True

        trk_hitdfs.append(this_trk_hits)

    hitdf = pd.concat(trk_hitdfs)
    hitdf = hitdf.reset_index().set_index(["entry", "rec.slc..index", "rec.slc.reco.pfp..index", "rec.slc.reco.pfp.trk.calo.{}.points..index".format(plane)])
    return hitdf

def thetazx_cut(df):
    dirx = df.pfp.trk.dir.x
    dirz = df.pfp.trk.dir.z
    coszx = np.abs(dirx / np.hypot(dirx, dirz))
    coszx_pass = coszx < 0.75
    
    return df[coszx_pass]

###############################################
#### Hit level selections
###############################################
def drop_endhits(hitdf, plane=2, ndrop=3):
    trk_hitdfs = []

    trks = hitdf.reset_index(level=[3]).index.unique()
    for pidx in range(len(trks)):  # removed tqdm
        this_trk_hits = hitdf.reset_index(level=[3]).loc[trks[pidx]]
        this_trk_hits = this_trk_hits.sort_values(by=("rr"), ascending=False)
        this_trk_hits = this_trk_hits.iloc[ndrop:]
        trk_hitdfs.append(this_trk_hits)

    hitdf = pd.concat(trk_hitdfs)
    return hitdf

def check_badorder(df, plane=2):
    trk_hitdfs = []
    df["badorder"] = False

    trks = df.index.unique()
    for pidx in range(len(trks)):
        this_trk_hits = df.loc[trks[pidx]] #.reset_index()

        idx_rr = this_trk_hits.sort_values(by=("rr"), ascending=False).index
        idx_wire = this_trk_hits.sort_values(by=("wire"), ascending=False).index
        if (idx_rr[0] > idx_rr[-1]):
            idx_rr = idx_rr[::-1]
        if (idx_wire[0] > idx_wire[-1]):
            idx_wire = idx_wire[::-1]

        if ((idx_rr != idx_wire).any()):
            this_trk_hits["badorder"] = True

        trk_hitdfs.append(this_trk_hits)

    df = pd.concat(trk_hitdfs)
    return df

def check_ntpcs(df, plane=2):
    trk_hitdfs = []
    df["ntpcs"] = False

    trks = df.index.unique()
    for pidx in range(len(trks)):
        this_trk_hits = df.loc[trks[pidx]] #.reset_index()

        if len(this_trk_hits) == 0:
            continue

        tpcs = this_trk_hits.tpc.unique()

        if (len(tpcs) == 2):
            this_trk_hits["ntpcs"] = 2
        elif (len(tpcs) == 1):
            this_trk_hits["ntpcs"] = 1
        else:
            raise ValueError(f"Number of TPCs is {len(tpcs)}")

        trk_hitdfs.append(this_trk_hits)

    df = pd.concat(trk_hitdfs)
    return df

def check_wireskip(df, plane=2):
    trk_hitdfs = []
    df["wireskip"] = 0

    trks = df.index.unique()
    for pidx in range(len(trks)):
        this_trk_hits = df.loc[trks[pidx]]
        if len(this_trk_hits) == 0:
            continue

        this_trk_hits = this_trk_hits.reset_index()
        nwires = this_trk_hits.wire.max() - this_trk_hits.wire.min()
        nidxs = this_trk_hits.index.max() - this_trk_hits.index.min()
        if (nwires != nidxs):
            this_trk_hits["wireskip"] = np.abs(nwires - nidxs)

        trk_hitdfs.append(this_trk_hits)

    df = pd.concat(trk_hitdfs)
    return df

###############################################
#### main function
###############################################

def make_protonhit_df(f):

    hdrdf = make_hdrdf(f)
    run = hdrdf.run
    subrun = hdrdf.subrun
    evt = hdrdf.evt

    pandoradf = make_pandora_df(f)
    pandoradf = pandoradf[InFV_nohiyz(pandoradf.slc.vertex)]
    pandoradf = pandoradf[(pandoradf.pfp.trk.len > 4.) & (pandoradf.pfp.dist_to_vertex < 6.)]

    trkhitdf = make_trkhitdf_plane2(f)
    trktruehitdf = make_trktruehitdf_plane2(f)
    trkhitdf = trkhitdf.join(trktruehitdf)
    trkhitdf = trkhitdf.join(subrun)
    trkhitdf = trkhitdf.join(evt)
    
    #print(pandoradf)
    #print(pandoradf.pfp.trk.columns)

    # 1 prong 
    longer_data, shorter_data = select_topology(pandoradf, nprong=1, max_len_cut=25, min_len_cut=25, contained=True)
    proton_candidate = select_calorimetry(shorter_data, muscore_cut_val=20, pscore_cut_val=80)
    proton_candidate = thetazx_cut(proton_candidate)

    proton_candidate_ture_pdg = proton_candidate.pfp.trk.truth.p.pdg
    
    proton_candidate_chi2_muon = proton_candidate.pfp.trk.chi2pid.I2.chi2_muon
    proton_candidate_chi2_proton = proton_candidate.pfp.trk.chi2pid.I2.chi2_proton
    
    proton_hits = trkhitdf.reset_index(level=[3]).loc[proton_candidate.index]
    proton_hits = proton_hits.reset_index().set_index(["entry", "rec.slc..index", "rec.slc.reco.pfp..index", "rec.slc.reco.pfp.trk.calo.{}.points..index".format(2)])
    proton_hits['chi2_muon'] = proton_candidate_chi2_muon
    proton_hits['chi2_proton'] = proton_candidate_chi2_proton
    proton_hits['true_pdg'] = proton_candidate_ture_pdg
    
    if len(proton_hits) == 0:
        return None ### FIXME

    proton_hits_1 = check_flipped(proton_hits, plane=2)
    proton_hits_1["selection"] = 1

    # 2 prong
    longer_data, shorter_data = select_topology(pandoradf, nprong=2, max_len_cut=25, min_len_cut=25, contained=True)
    proton_candidate = select_calorimetry(shorter_data, muscore_cut_val=20, pscore_cut_val=80)
    proton_candidate = thetazx_cut(proton_candidate)

    proton_candidate_ture_pdg = proton_candidate.pfp.trk.truth.p.pdg
    
    proton_candidate_chi2_muon = proton_candidate.pfp.trk.chi2pid.I2.chi2_muon
    proton_candidate_chi2_proton = proton_candidate.pfp.trk.chi2pid.I2.chi2_proton
    
    proton_hits = trkhitdf.reset_index(level=[3]).loc[proton_candidate.index]
    proton_hits = proton_hits.reset_index().set_index(["entry", "rec.slc..index", "rec.slc.reco.pfp..index", "rec.slc.reco.pfp.trk.calo.{}.points..index".format(2)])
    proton_hits['chi2_muon'] = proton_candidate_chi2_muon
    proton_hits['chi2_proton'] = proton_candidate_chi2_proton
    proton_hits['true_pdg'] = proton_candidate_ture_pdg

    if len(proton_hits) == 0:
        return None ### FIXME

    proton_hits_2 = check_flipped(proton_hits, plane=2)
    proton_hits_2["selection"] = 2

    proton_hits = pd.concat([proton_hits_1, proton_hits_2], ignore_index=False)

    # apply hit-level selections
    proton_hits = proton_hits[proton_hits.flipped == False]
    proton_hits = drop_endhits(proton_hits, ndrop=3, plane=2)
    proton_hits = proton_hits[proton_hits.rr > 0.6]

    proton_hits = check_badorder(proton_hits)
    proton_hits = check_ntpcs(proton_hits)
    proton_hits = check_wireskip(proton_hits)

    proton_hits = proton_hits[(proton_hits.pitch < 1) & (proton_hits.badorder == False)]

    #print(proton_hits)

    # output dfs
    proton_hits_selected = proton_hits.loc[:, [
        ("selection"),
        ("run"),
        ("subrun"),
        ("evt"),
        ("tpc"),
        ("plane"),
        ("rr"),
        ("dqdx"),
        ("integral"),
        ("sumadc"),
        ("h_e"),
        ("h_nelec"),
        ("x"),
        ("y"),
        ("z"),
        ("t"),
        ("pitch"),
        ("phi"),
        ("efield"),
        ("chi2_muon"),
        ("chi2_proton"),
        ("true_pdg"),
    ]].copy()

    proton_hits_selected.columns = [
        "selection",
        "run",
        "subrun",
        "evt",
        "tpc",
        "plane",
        "rr",
        "dqdx",
        "integral",
        "sumadc",
        "h_e",
        "h_nelec",
        "x",
        "y",
        "z",
        "t",
        "pitch",
        "phi",
        "efield",
        "trk_chi2_mu",
        "trk_chi2_p",
        "true_pdg",
    ]
    proton_hits_selected = proton_hits_selected.reset_index()
    #print(proton_hits_selected)

    return proton_hits_selected
