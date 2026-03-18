from makedf.makedf import *
from pyanalib.pandas_helpers import *
from makedf.util import *
import warnings
warnings.filterwarnings('ignore')

def make_mcnudf_hnl(f,**args):
    mcdf = make_mcnudf(f,**args)
    # drop mcdf columns not relevant for this analysis
    if 'mu'  in list(zip(*list(mcdf.columns)))[0]:  mcdf = mcdf.drop('mu', axis=1,level=0)
    if 'p'   in list(zip(*list(mcdf.columns)))[0]:  mcdf = mcdf.drop('p',  axis=1,level=0)
    if 'cpi' in list(zip(*list(mcdf.columns)))[0]:  mcdf = mcdf.drop('cpi',axis=1,level=0)
    return mcdf

def make_hnldf_mc_wgt(f):
    df = make_hnldf_mc(f,include_weights=True)
    return df

def make_hnldf_mc(f, include_weights=False,multisim_nuniv=100,slim=True):
    
    slcdf = make_hnldf(f)
    mcdf = make_mcnudf_hnl(f,include_weights=include_weights,multisim_nuniv=multisim_nuniv,slim=slim)
    mcdf.columns = pd.MultiIndex.from_tuples([tuple(["slc", "truth"] + list(c)) for c in mcdf.columns])
    df = multicol_merge(slcdf.reset_index(), 
                        mcdf.reset_index(),
                        left_on=[('entry', '', '', '', '', ''), 
                                ('slc', 'tmatch', 'idx', '', '', '')], 
                        right_on=[('entry', '', '', '', '', ''), 
                                ('rec.mc.nu..index', '', '')], 
                        how="left")
    df = df.set_index(slcdf.index.names, verify_integrity=True)
    return df

def make_hnldf_data(f):
    slcdf = make_hnldf(f)
    # drop truth cols for data
    slcdf = slcdf.drop('tmatch', axis=1,level=1) # slc level
    slcdf = slcdf.drop('truth',  axis=1,level=2) # pfp level
   
    framedf = make_framedf(f)
    timingdf = make_timingdf(f)
    ftdf = multicol_merge(framedf, timingdf, left_index=True, right_index=True, how="left", validate="one_to_one")

    df = multicol_merge(slcdf.reset_index(), 
                        ftdf.reset_index(),
                        left_on=[('entry', '', '', '', '', '')],
                        right_on=[('entry', '', '', '', '', '')], 
                        how="left")
    df = df.set_index(slcdf.index.names, verify_integrity=True)
    return df
    
def p_to_energy(p, mass):
    return np.sqrt(p**2 + mass**2)

def make_hnldf(f):
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det
    DETECTOR = "SBND"
    
    pfpdf = make_pfpdf(f)
    slcdf = loadbranches(f["recTree"], slcbranches + barycenterFMbranches + correctedflashbranches)
    slcdf = slcdf.rec
    
    pfpdf = pfpdf.drop('pfochar',axis=1,level=1)

    #SHOWER
    #count number of showers per event before selecting primary and secondary shower
    nshwdf = pfpdf[(pfpdf.pfp.trackScore < 0.5) & (pfpdf.pfp.shw.maxplane_energy > 0)].groupby(level=[0,1]).size().to_frame('n_shws') #count number of showers
    nshwdf.columns = pd.MultiIndex.from_tuples([tuple(['slc'] + list(nshwdf.columns))]) #change nshws column to have multiindex with same structure as shwdf
    slcdf = multicol_merge(slcdf, nshwdf,left_index=True,right_index=True,how="left",validate="one_to_one")
    slcdf['slc','n_shws'] = slcdf['slc','n_shws'].fillna(0)

    del nshwdf

    ## primary shw candidate is shw pfp with highest energy, valid energy, and score < 0.5
    shwdf = pfpdf[(pfpdf.pfp.trackScore < 0.5) & (pfpdf.pfp.shw.maxplane_energy > 0)].sort_values(pfpdf.pfp.index.names[:-1] + [('pfp','shw','maxplane_energy','','','')]).groupby(level=[0,1]).nth(-1)
    # drop all columns that are from trk attributes
    shwdf = shwdf.drop('trk',axis=1,level=1)
    shwdf.columns = shwdf.columns.set_levels(['primshw'],level=0)
    slcdf = multicol_merge(slcdf, shwdf.droplevel(-1),left_index=True,right_index=True,how="left",validate="one_to_one")

    del shwdf

    ## secondary shower is shw pfp with second highest energy, valid energy, and score < 0.5 
    shwsecdf = pfpdf[(pfpdf.pfp.trackScore < 0.5) & (pfpdf.pfp.shw.maxplane_energy > 0)].sort_values(pfpdf.pfp.index.names[:-1] + [('pfp','shw','maxplane_energy','','','')]).groupby(level=[0,1]).nth(-2)
    shwsecdf = shwsecdf.drop('trk',axis=1,level=1)
    shwsecdf.columns = shwsecdf.columns.set_levels(['secshw'],level=0)
    slcdf = multicol_merge(slcdf, shwsecdf.droplevel(-1),left_index=True,right_index=True,how="left",validate="one_to_one")

    del shwsecdf

    #TRACKS
    #GENERAL TRACKS: track score > 0.5 & track length > 0 
    trkdf = pfpdf[(pfpdf.pfp.trackScore > 0.5) & (pfpdf.pfp.trk.len > 0)]
    ntrkdf = trkdf.groupby(level=[0,1]).size().to_frame('n_trks') #count number of tracks
    ntrkdf.columns = pd.MultiIndex.from_tuples([tuple(['slc'] + list(ntrkdf.columns))]) #change ntrks column to have multiindex with same structure as trkdf
    slcdf = multicol_merge(slcdf, ntrkdf,left_index=True,right_index=True,how="left",validate="one_to_one")
    slcdf['slc','n_trks'] = slcdf['slc','n_trks'].fillna(0)

    del ntrkdf

    #save only longest track
    trkdf = trkdf.sort_values(pfpdf.pfp.index.names[:-1] + [('pfp','trk','len','','','')]).groupby(level=[0,1]).nth(-1)
    # drop all columns that are from shw attributes
    trkdf = trkdf.drop('shw',axis=1,level=1)
    trkdf.columns = trkdf.columns.set_levels(['primtrk'],level=0)
    slcdf = multicol_merge(slcdf, trkdf.droplevel(-1),left_index=True,right_index=True,how="left",validate="one_to_one")

    del trkdf

    #MUON-LIKE TRACKS:
    mudf = pfpdf[(pfpdf.pfp.trackScore > 0.5) & (pfpdf.pfp.trk.len > 10) & (pfpdf.pfp.trk.chi2pid.I2.chi2_muon < 20) & (pfpdf.pfp.trk.chi2pid.I2.chi2_proton > 85)]
    nmudf = mudf.groupby(level=[0,1]).size().to_frame('n_mu_trks') #count number of tracks
    nmudf.columns = pd.MultiIndex.from_tuples([tuple(['slc'] + list(nmudf.columns))])
    slcdf = multicol_merge(slcdf, nmudf,left_index=True,right_index=True,how="left",validate="one_to_one")
    slcdf['slc','n_mu_trks'] = slcdf['slc','n_mu_trks'].fillna(0)

    del nmudf

    #save only primary muon with the highest energy
    mu_mass = 0.1056583745 # GeV/c^2
    mup_column = ('pfp','trk','rangeP','p_muon','','')
    energy_column = ('pfp','trk','energy','','','')
    mudf[energy_column] = mudf.apply(lambda row: p_to_energy(row[mup_column], mu_mass), axis=1)
    mudf = mudf.sort_values(pfpdf.pfp.index.names[:-1] + [energy_column]).groupby(level=[0,1]).nth(-1)
    # drop all columns that are from shw attributes
    mudf = mudf.drop('shw',axis=1,level=1)
    mudf.columns = mudf.columns.set_levels(['primmu'],level=0)
    slcdf = multicol_merge(slcdf, mudf.droplevel(-1),left_index=True,right_index=True,how="left",validate="one_to_one")

    del mudf

    #PROTON-LIKE TRACKS:
    prodf = pfpdf[(pfpdf.pfp.trackScore > 0.5) & (pfpdf.pfp.trk.len > 0) & (pfpdf.pfp.trk.chi2pid.I2.chi2_muon > 20) & (pfpdf.pfp.trk.chi2pid.I2.chi2_proton < 85)]
    nprodf = prodf.groupby(level=[0,1]).size().to_frame('n_proton_trks') #count number of tracks
    nprodf.columns = pd.MultiIndex.from_tuples([tuple(['slc'] + list(nprodf.columns))])
    slcdf = multicol_merge(slcdf, nprodf,left_index=True,right_index=True,how="left",validate="one_to_one")
    slcdf['slc','n_proton_trks'] = slcdf['slc','n_proton_trks'].fillna(0)

    del nprodf

    #save only primary proton with the highest energy
    pro_mass = 0.9382720813 # GeV/c^2
    prop_column = ('pfp','trk','rangeP','p_proton','','')
    prodf[energy_column] = prodf.apply(lambda row: p_to_energy(row[prop_column], pro_mass), axis=1)
    prodf = prodf.sort_values(pfpdf.pfp.index.names[:-1] + [energy_column]).groupby(level=[0,1]).nth(-1)
    # drop all columns that are from shw attributes
    prodf = prodf.drop('shw',axis=1,level=1)
    prodf.columns = prodf.columns.set_levels(['primproton'],level=0)
    slcdf = multicol_merge(slcdf, prodf.droplevel(-1),left_index=True,right_index=True,how="left",validate="one_to_one")

    del prodf

    # pre-selection cuts
    slcdf = slcdf[slcdf.slc.is_clear_cosmic==0]
    slcdf = slcdf[slcdf.slc.nu_score > 0.5]
    slcdf = slcdf[InFV(df=slcdf.slc.vertex, inzback=0, det=DETECTOR)]    
    
    # recalc dEdx for the primary shower
    # do after pre-selection to speed up
    for plane in range(3):
        trkhitdf = make_trkhitdf(f,plane)
        slchitdf = multicol_merge(slcdf.reset_index(), 
                                   trkhitdf.reset_index(),
                                   left_on=[('entry', '', '', '', '', ''), 
                                           ('rec.slc..index', '', '', '', '', ''),
                                           ('primshw','tindex', '', '', '', '')], 
                                   right_on=[('entry', '', '', '', '', ''), 
                                           ('rec.slc..index', '', '', '', '', ''),
                                           ('rec.slc.reco.pfp..index', '', '', '', '', ''),],
                                   how="left")
        slchitdf = slchitdf.set_index(trkhitdf.index.names,verify_integrity=True)
        slchitdf = multicol_add(slchitdf,dmagdf(slchitdf.primshw.shw.start,slchitdf).rename("sp_to_start"))
        
        # require that the spacepoints are between 0.5 cm and 5 cm of the shower start
        # require that the spacepoints are within the AV
        slchitdf = slchitdf[(slchitdf.sp_to_start < 5) & (slchitdf.sp_to_start > 0.5)]
        slchitdf = slchitdf[InAV(slchitdf)]

        slchitdf['dedx_reco'] = chi2pid.dedx(slchitdf,gain="SBND",calibrate="SBND",plane=plane)
        this_dedx_col = ('primshw','shw','plane',f'I{plane}','dEdx_new')
        this_hits_col = ('primshw','shw','plane',f'I{plane}','nHits_dEdx')

        slcdf = multicol_add(slcdf,slchitdf[('dedx_reco', '', '', '', '', '')].groupby(slchitdf.index.names[:-2]).median().rename(this_dedx_col), default=-999)
        slcdf = multicol_add(slcdf,slchitdf[('dedx_reco', '', '', '', '', '')].groupby(slchitdf.index.names[:-2]).count().rename(this_hits_col),default=-999)
    
    slcdf['primshw','shw','maxplane_dEdx_new','','',''] = np.nan
    slcdf['primshw','shw','maxplane_idx','','',''] = slcdf.loc(axis=1)['primshw','shw','plane',:,"nHits_dEdx"].idxmax(axis=1).apply(lambda x: x[3])

    conditions = [slcdf['primshw','shw','maxplane_idx','','','']=="I2",slcdf['primshw','shw','maxplane_idx','','','']=="I1",slcdf['primshw','shw','maxplane_idx','','','']=="I0"]
    choices = [slcdf['primshw','shw','plane','I2','dEdx_new',''],slcdf['primshw','shw','plane','I1','dEdx_new',''],slcdf['primshw','shw','plane','I0','dEdx_new','']]
    slcdf['primshw','shw','maxplane_dEdx_new','','',''] = np.select(conditions,choices,default=np.nan)
    
    return slcdf 
