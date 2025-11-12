from makedf.makedf import *
from pyanalib.pandas_helpers import *
from makedf.util import *

def make_mcnudf_nuecc(f,**args):
    mcdf = make_mcnudf(f,**args)
    # drop mcdf columns not relevant for this analysis
    if 'mu'  in list(zip(*list(mcdf.columns)))[0]:  mcdf = mcdf.drop('mu', axis=1,level=0)
    if 'p'   in list(zip(*list(mcdf.columns)))[0]:  mcdf = mcdf.drop('p',  axis=1,level=0)
    if 'cpi' in list(zip(*list(mcdf.columns)))[0]:  mcdf = mcdf.drop('cpi',axis=1,level=0)
    return mcdf

def make_nueccdf_mc_wgt(f):
    df = make_nueccdf_mc(f,include_weights=True)
    return df

def make_nueccdf_mc(f, include_weights=False,multisim_nuniv=100,slim=True):
    
    slcdf = make_nueccdf(f)
    mcdf = make_mcnudf_nuecc(f,include_weights=include_weights,multisim_nuniv=multisim_nuniv,slim=slim)
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

def make_nueccdf_data(f):
    slcdf = make_nueccdf(f)
    # drop truth cols for data
    slcdf = slcdf.drop('tmatch', axis=1,level=1) # slc level
    slcdf = slcdf.drop('truth',  axis=1,level=2) # pfp level
    
    ## keep the only relevant column (for now)
    framedf = make_framedf(f)[['frameApplyAtCaf']]
    
    df = multicol_merge(slcdf.reset_index(), 
                        framedf.reset_index(),
                        left_on=[('entry', '', '', '', '', '')],
                        right_on=[('entry', '', '', '', '', '')], 
                        how="left")
    df = df.set_index(slcdf.index.names, verify_integrity=True)
    return df

def make_nueccdf(f):
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det
    if (1 == det.unique()):
        DETECTOR = "SBND"
    else:
        DETECTOR = "ICARUS"

    assert DETECTOR == "SBND"
    
    pfpdf = make_pfpdf(f)
    slcdf = loadbranches(f["recTree"], slcbranches)
    slcdf = slcdf.rec
    
    pfpdf = pfpdf.drop('pfochar',axis=1,level=1)
    ## primary shw candidate is shw pfp with highest energy, valid energy, and score < 0.5
    shwdf = pfpdf[(pfpdf.pfp.trackScore < 0.5) & (pfpdf.pfp.shw.maxplane_energy > 0)].sort_values(pfpdf.pfp.index.names[:-1] + [('pfp','shw','maxplane_energy','','','')]).groupby(level=[0,1]).nth(-1)
    # drop all columns that are from trk attributes
    shwdf = shwdf.drop('trk',axis=1,level=1)
    shwdf.columns = shwdf.columns.set_levels(['primshw'],level=0)
    slcdf = multicol_merge(slcdf, shwdf.droplevel(-1),left_index=True,right_index=True,how="right",validate="one_to_one")

    ## primary trk is track pfp with the longest length
    trkdf = pfpdf[(pfpdf.pfp.trackScore > 0.5) & (pfpdf.pfp.trk.len > 0)].sort_values(pfpdf.pfp.index.names[:-1] + [('pfp','trk','len','','','')]).groupby(level=[0,1]).nth(-1)
    # drop all columns that are from shw attributes
    trkdf = trkdf.drop('shw',axis=1,level=1)
    trkdf.columns = trkdf.columns.set_levels(['primtrk'],level=0)
    slcdf = multicol_merge(slcdf, trkdf.droplevel(-1),left_index=True,right_index=True,how="left",validate="one_to_one")

    ## secondary shower is shw pfp with second highest energy, valid energy, and score < 0.5 
    shwsecdf = pfpdf[(pfpdf.pfp.trackScore < 0.5) & (pfpdf.pfp.shw.maxplane_energy > 0)].sort_values(pfpdf.pfp.index.names[:-1] + [('pfp','shw','maxplane_energy','','','')]).groupby(level=[0,1]).nth(-2)
    shwsecdf = shwsecdf.drop('trk',axis=1,level=1)
    shwsecdf.columns = shwsecdf.columns.set_levels(['secshw'],level=0)
    slcdf = multicol_merge(slcdf, shwsecdf.droplevel(-1),left_index=True,right_index=True,how="left",validate="one_to_one")
    
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
