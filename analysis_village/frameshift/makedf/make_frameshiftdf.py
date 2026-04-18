from makedf.makedf import *
from pyanalib.pandas_helpers import *
from makedf.util import *
import warnings
warnings.filterwarnings('ignore')

def make_frameshiftdf_slc_FV(f, FVcut = True):
    df = make_frameshiftdf_slc(f, FVcut)
    return df
    
def make_frameshiftdf_slc(f, FVcut = False):

    DETECTOR = "SBND"
    
    slcdf = loadbranches(f["recTree"], slcbranches + barycenterFMbranches + correctedflashbranches)
    slcdf = slcdf.rec
    # drop truth cols for data
    slcdf = slcdf.drop('tmatch', axis=1,level=1) # slc level

    # pre-selection cuts
    slcdf = slcdf[slcdf.slc.is_clear_cosmic==0]
    slcdf = slcdf[slcdf.slc.nu_score > 0.5]
    if FVcut:
        slcdf = slcdf[InFV(df=slcdf.slc.vertex, inzback=0, det=DETECTOR)]    

    #save frameshift/timing info columns
    framedf = make_framedf(f)
    #framedf = make_framedf_new(f)
    timingdf = make_timingdf(f)
    ftdf = multicol_merge(framedf, timingdf, left_index=True, right_index=True, how="left", validate="one_to_one")
    
    df = multicol_merge(slcdf.reset_index(), 
                        ftdf.reset_index(),
                        left_on=[('entry')],
                        right_on=[('entry')], 
                        how="left")
    df = df.set_index(slcdf.index.names, verify_integrity=True)
    return df

def make_frameshiftdf_opflash(f):

    DETECTOR = "SBND"
    
    opdf = loadbranches(f["recTree"], opflashbranches)
    opdf = opdf[opdf.rec.opflashes.totalpe > 1000]    

    #save frameshift/timing info columns
    #framedf = make_framedf(f)
    framedf = make_framedf_new(f)
    timingdf = make_timingdf(f)
    ftdf = multicol_merge(framedf, timingdf, left_index=True, right_index=True, how="left", validate="one_to_one")
    
    df = multicol_merge(opdf.reset_index(), 
                        ftdf.reset_index(),
                        left_on=[('entry')],
                        right_on=[('entry')], 
                        how="left")
    df = df.set_index(opdf.index.names, verify_integrity=True)
    
    return df

def make_frameshiftdf_crtsp(f):

    DETECTOR = "SBND"

    crtbranches = [
                    'rec.crt_spacepoints.time'
                  ,'rec.crt_spacepoints.time_err'
                  ,'rec.crt_spacepoints.pe'
                  ,'rec.crt_spacepoints.complete'
                  ,'rec.crt_spacepoints.tagger'
                  ]
    crtdf = loadbranches(f["recTree"], crtbranches)
    crtdf = crtdf[crtdf.rec.crt_spacepoints.complete == 1]   

    #save frameshift/timing info columns
    framedf = make_framedf(f)
    #framedf = make_framedf_new(f)
    timingdf = make_timingdf(f)
    ftdf = multicol_merge(framedf, timingdf, left_index=True, right_index=True, how="left", validate="one_to_one")
    
    df = multicol_merge(crtdf.reset_index(), 
                        ftdf.reset_index(),
                        left_on=[('entry')],
                        right_on=[('entry')], 
                        how="left")
    df = df.set_index(crtdf.index.names, verify_integrity=True)
    
    return df
    
def make_frameshiftdf_crttrk(f):

    DETECTOR = "SBND"

    crtbranches = ['rec.sbnd_crt_tracks.time']
    crtdf = loadbranches(f["recTree"], crtbranches)

    #save frameshift/timing info columns
    #framedf = make_framedf(f)
    framedf = make_framedf_new(f)    
    timingdf = make_timingdf(f)
    ftdf = multicol_merge(framedf, timingdf, left_index=True, right_index=True, how="left", validate="one_to_one")
    
    df = multicol_merge(crtdf.reset_index(), 
                        ftdf.reset_index(),
                        left_on=[('entry')],
                        right_on=[('entry')], 
                        how="left")
    df = df.set_index(crtdf.index.names, verify_integrity=True)
    
    return df

def make_frameshiftdf_crtspmatch(f):

    DETECTOR = "SBND"
    
    pfpdf = make_pfpdf(f)
    pfpdf = pfpdf.drop('pfochar',axis=1,level=1)
    pfpdf = pfpdf.drop('shw',axis=1,level=1)

    trk_score = ('pfp', 'trackScore', '', '', '', '')
    pfpdf = pfpdf[pfpdf[trk_score] > 0.0] #save valid pfp

    crt_branches = [('pfp', 'trk', 'crtspacepoint', 'score','','')
                      , ('pfp', 'trk', 'crtspacepoint', 'spacepoint', 'time','')]
    crtdf = pfpdf[crt_branches]  
    crtdf = crtdf.dropna()
    
    #save frameshift/timing info columns
    #framedf = make_framedf(f)
    framedf = make_framedf_new(f)
    timingdf = make_timingdf(f)
    ftdf = multicol_merge(framedf, timingdf, left_index=True, right_index=True, how="left", validate="one_to_one")
    
    df = multicol_merge(crtdf.reset_index(), 
                        ftdf.reset_index(),
                        left_on=[('entry')],
                        right_on=[('entry')], 
                        how="left")
    df = df.set_index(crtdf.index.names, verify_integrity=True)

    return df


def make_frameshiftdf_crttrkmatch(f):

    DETECTOR = "SBND"
    
    pfpdf = make_pfpdf(f)
    pfpdf = pfpdf.drop('pfochar',axis=1,level=1)
    pfpdf = pfpdf.drop('shw',axis=1,level=1)

    trk_score = ('pfp', 'trackScore', '', '', '', '')
    pfpdf = pfpdf[pfpdf[trk_score] > 0.0] #save valid pfp

    crt_branches = [('pfp', 'trk', 'crtsbndtrack', 'score','','')
                      , ('pfp', 'trk', 'crtsbndtrack', 'track', 'time','')]
    crtdf = pfpdf[crt_branches]  
    crtdf = crtdf.dropna()
    
    #save frameshift/timing info columns
    #framedf = make_framedf(f)
    framedf = make_framedf_new(f)   
    timingdf = make_timingdf(f)
    ftdf = multicol_merge(framedf, timingdf, left_index=True, right_index=True, how="left", validate="one_to_one")
    
    df = multicol_merge(crtdf.reset_index(), 
                        ftdf.reset_index(),
                        left_on=[('entry')],
                        right_on=[('entry')], 
                        how="left")
    df = df.set_index(crtdf.index.names, verify_integrity=True)
    
    return df