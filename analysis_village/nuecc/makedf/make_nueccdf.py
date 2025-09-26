from makedf.makedf import *
from pyanalib.pandas_helpers import *
from makedf.util import *

def make_nueccdf(f, 
                       include_weights=True, 
                       multisim_nuniv=100,
                       wgt_types=["bnb","genie"],
                       slim=True,
                       cutClearCosmics=True,
                       ):
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det
    if (1 == det.unique()):
        DETECTOR = "SBND"
    else:
        DETECTOR = "ICARUS"

    assert DETECTOR == "SBND"
    
    mcdf = make_mcnudf(f)
    pfpdf = make_trkdf(f)
    slcdf = make_slcdf(f)
    
    # set shower energy as the one with the plane that has the most number of hits (maxplane)
    pfpdf['pfp','shw','maxplane','','',''] = pfpdf.loc(axis=1)['pfp','shw','plane',:,"nHits"].idxmax(axis=1).apply(lambda x: x[3])
    pfpdf['pfp','shw','maxplane_energy','','',''] = np.nan
    conditions = [pfpdf['pfp','shw','maxplane','','','']=="I2",pfpdf['pfp','shw','maxplane','','','']=="I1",pfpdf['pfp','shw','maxplane','','','']=="I0"]
    choices = [pfpdf['pfp','shw','plane','I2','energy',''],pfpdf['pfp','shw','plane','I1','energy',''],pfpdf['pfp','shw','plane','I0','energy','']]
    pfpdf['pfp','shw','maxplane_energy','','',''] = np.select(conditions,choices,default=np.nan)

    ## primary shw candidate is shw pfp with highest energy, valid energy, and score < 0.5
    shwdf = pfpdf[(pfpdf.pfp.trackScore < 0.5) & (pfpdf.pfp.shw.maxplane_energy > 0)].sort_values(pfpdf.pfp.index.names[:-1] + [('pfp','shw','maxplane_energy','','','')]).groupby(level=[0,1]).nth(-1)
    # drop all columns that are from trk attributes
    shwdf = shwdf.drop('trk',axis=1,level=1)
    shwdf.columns = shwdf.columns.set_levels(['primshw'],level=0)
    slcdf = multicol_merge(slcdf, shwdf,left_index=True,right_index=True,how="right",validate="one_to_one")

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
    slcdf = multicol_merge(slcdf, shwsecdf,left_index=True,right_index=True,how="left",validate="one_to_one")
    
    mcdf.columns = pd.MultiIndex.from_tuples([tuple(list(c) +["", ""]) for c in mcdf.columns])     # match # of column levels

    # pre-selection cuts
    slcdf = slcdf[slcdf.slc.is_clear_cosmic==0]
    slcdf = slcdf[slcdf.slc.nu_score > 0.5]
    slcdf = slcdf[InFV(df=slcdf.slc.vertex, inzback=0, det=DETECTOR)]
    
    df = multicol_merge(slcdf.reset_index(), 
                        mcdf.reset_index(),
                        left_on=[('entry', '', '', '', '', ''), 
                                ('slc', 'tmatch', 'idx', '', '', '')], 
                        right_on=[('entry', '', '', '', '', ''), 
                                ('rec.mc.nu..index', '', '')], 
                        how="left")

    df = df.set_index(slcdf.index.names, verify_integrity=True)
    
    return df
