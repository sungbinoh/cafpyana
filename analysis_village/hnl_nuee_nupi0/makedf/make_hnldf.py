from makedf.makedf import *
from pyanalib.pandas_helpers import *
from makedf.util import *
import warnings
warnings.filterwarnings('ignore')

#----------------------------------------------------------------------------------#
# Functions to read MC neutrinos
#----------------------------------------------------------------------------------#

def make_mcnudf_hnl(f,**args):
    mcdf = make_mcnudf(f,**args)
    # drop mcdf columns not relevant for this analysis
    if 'mu'  in list(zip(*list(mcdf.columns)))[0]:  mcdf = mcdf.drop('mu', axis=1,level=0)
    if 'p'   in list(zip(*list(mcdf.columns)))[0]:  mcdf = mcdf.drop('p',  axis=1,level=0)
    if 'cpi' in list(zip(*list(mcdf.columns)))[0]:  mcdf = mcdf.drop('cpi',axis=1,level=0)
    return mcdf

def make_hnldf_mcnu_wgt(f):
    df = make_hnldf_mcnu(f,include_weights=True)
    return df

def make_hnldf_mcnu_nopreselect_savepfp(f):
    df = make_hnldf_mcnu(f,applyPreselection=False, savePfp=True)
    return df

def make_hnldf_mcnu(f, include_weights=False,multisim_nuniv=100,slim=True,applyPreselection=True, savePfp=False):
    slcdf = make_hnldf(f, applyPreselection=applyPreselection, savePfp=savePfp)
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

#----------------------------------------------------------------------------------#
# Functions to read MC MeVPrtl HNLs
#----------------------------------------------------------------------------------#

def make_hnldf_mevprtl_wgt(f):
    df = make_hnldf_mevprtl(f,include_weights=True)
    return df

def make_hnldf_mevprtl_nopreselect_savepfp(f):
    df = make_hnldf_mevprtl(f,applyPreselection=False, savePfp=True)
    return df

def make_hnldf_mevprtl(f, include_weights=False, multisim_nuniv=100,slim=True,applyPreselection=True, savePfp=False):
    slcdf = make_hnldf(f, applyPreselection=applyPreselection, savePfp=savePfp)
    mcdf = make_mevprtldf(f,include_weights=include_weights,multisim_nuniv=multisim_nuniv,slim=slim)
    mcdf.columns = pd.MultiIndex.from_tuples([tuple(["slc", "prtl"] + list(c)) for c in mcdf.columns])
    df = multicol_merge(slcdf.reset_index(), 
                        mcdf.reset_index(),
                        left_on=[('entry', '', '', '', '', ''), 
                                ('slc', 'tmatch', 'idx', '', '', '')], 
                        right_on=[('entry', '', '', '', '', ''), 
                                ('rec.mc.prtl..index', '', '')], 
                        how="left")
    df = df.set_index(slcdf.index.names, verify_integrity=True)
    return df

#----------------------------------------------------------------------------------#
# Functions to read data
#----------------------------------------------------------------------------------#

def make_hnldf_data(f, applyPreselection=True, savePfp=False):
    slcdf = make_hnldf(f, applyPreselection=applyPreselection, savePfp=savePfp)
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

def make_hnldf_data_nopreselect_savepfp(f):
    df = make_hnldf_data(f, applyPreselection=False, savePfp=True)
    return df
    
#----------------------------------------------------------------------------------#
# Functions to read reconstruction (both MC/data)
#----------------------------------------------------------------------------------#
    
def p_to_energy(p, mass):
    return np.sqrt(p**2 + mass**2)

def make_hnldf(f, applyPreselection=False, savePfp=False, trackScore=0.51):
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det
    DETECTOR = "SBND"
    
    pfpdf = make_pfpdf(f)
    slcdf = loadbranches(f["recTree"], slcbranches + barycenterFMbranches + correctedflashbranches)
    slcdf = slcdf.rec
    
    pfpdf = pfpdf.drop('pfochar',axis=1,level=1)

    #Save all PFPs in a slice
    if savePfp:
        # Keep slcdf as the left/base table and expand rows with pfp-level info.
        slcdf = multicol_merge(slcdf,
                           pfpdf.reset_index(level='rec.slc.reco.pfp..index'),
                           left_index=True,
                           right_index=True,
                           how="left",
                           validate="one_to_many")

        pfp_idx_col = next(
            c for c in slcdf.columns
            if c == 'rec.slc.reco.pfp..index' or (isinstance(c, tuple) and len(c) > 0 and c[0] == 'rec.slc.reco.pfp..index')
        )
        slcdf = slcdf.set_index(pfp_idx_col, append=True)
    
    #Save only 2 primary showers in a slice
    else:
        #Make sure to count only valid PFPs so apply a PFP sanity filter
        #1. valid trackScore > 0
        pfpdf = pfpdf[(pfpdf.pfp.trackScore > 0)]
        #2. valid shower Energy
        pfpdf = pfpdf[(pfpdf.pfp.shw.bestplane_energy > 0)]
        #3. valid start x/y/z
        pfpdf = pfpdf[(pfpdf.pfp.shw.start.x != -999) & (pfpdf.pfp.shw.start.y != -999) & (pfpdf.pfp.shw.start.z != -999)]
        #4. valid track length > 0
        pfpdf = pfpdf[(pfpdf.pfp.trk.len > 0)]

        #SHOWER
        #count number of showers per event before selecting primary and secondary shower
        nshwdf = pfpdf[(pfpdf.pfp.trackScore < trackScore)].groupby(level=[0,1]).size().to_frame('n_shws') #count number of showers
        nshwdf.columns = pd.MultiIndex.from_tuples([tuple(['slc'] + list(nshwdf.columns))]) #change nshws column to have multiindex with same structure as shwdf
        slcdf = multicol_merge(slcdf, nshwdf,left_index=True,right_index=True,how="left",validate="one_to_one")
        slcdf['slc','n_shws'] = slcdf['slc','n_shws'].fillna(0)
        del nshwdf

        # primary shw candidate is shw pfp with highest energy, valid energy, and score < trackScore
        shwdf = pfpdf[(pfpdf.pfp.trackScore < trackScore)].sort_values(pfpdf.pfp.index.names[:-1] + [('pfp','shw','bestplane_energy','','','')]).groupby(level=[0,1]).nth(-1)
        # drop all columns that are from trk attributes
        shwdf = shwdf.drop('trk',axis=1,level=1)
        shwdf.columns = shwdf.columns.set_levels(['primshw'],level=0)
        slcdf = multicol_merge(slcdf, shwdf.droplevel(-1),left_index=True,right_index=True,how="left",validate="one_to_one")
        del shwdf

        # secondary shower is shw pfp with second highest energy, valid energy, and score < 0.5 
        shwsecdf = pfpdf[(pfpdf.pfp.trackScore < trackScore)].sort_values(pfpdf.pfp.index.names[:-1] + [('pfp','shw','bestplane_energy','','','')]).groupby(level=[0,1]).nth(-2)
        shwsecdf = shwsecdf.drop('trk',axis=1,level=1)
        shwsecdf.columns = shwsecdf.columns.set_levels(['secshw'],level=0)
        slcdf = multicol_merge(slcdf, shwsecdf.droplevel(-1),left_index=True,right_index=True,how="left",validate="one_to_one")
        del shwsecdf

        #TRACKS
        #GENERAL TRACKS: track score > 0.5 & track length > 0 
        trkdf = pfpdf[(pfpdf.pfp.trackScore > trackScore)]
        ntrkdf = trkdf.groupby(level=[0,1]).size().to_frame('n_trks') #count number of tracks
        ntrkdf.columns = pd.MultiIndex.from_tuples([tuple(['slc'] + list(ntrkdf.columns))]) #change ntrks column to have multiindex with same structure as trkdf
        slcdf = multicol_merge(slcdf, ntrkdf,left_index=True,right_index=True,how="left",validate="one_to_one")
        slcdf['slc','n_trks'] = slcdf['slc','n_trks'].fillna(0)
        del ntrkdf

        ##save only longest track
        #trkdf = trkdf.sort_values(pfpdf.pfp.index.names[:-1] + [('pfp','trk','len','','','')]).groupby(level=[0,1]).nth(-1)
        ## drop all columns that are from shw attributes
        #trkdf = trkdf.drop('shw',axis=1,level=1)
        #trkdf.columns = trkdf.columns.set_levels(['primtrk'],level=0)
        #slcdf = multicol_merge(slcdf, trkdf.droplevel(-1),left_index=True,right_index=True,how="left",validate="one_to_one")
        #del trkdf

        ##MUON-LIKE TRACKS:
        #mudf = pfpdf[(pfpdf.pfp.trackScore > 0.5) & (pfpdf.pfp.trk.len > 10) & (pfpdf.pfp.trk.chi2pid.I2.chi2_muon < 20) & (pfpdf.pfp.trk.chi2pid.I2.chi2_proton > 85)]
        #nmudf = mudf.groupby(level=[0,1]).size().to_frame('n_mu_trks') #count number of tracks
        #nmudf.columns = pd.MultiIndex.from_tuples([tuple(['slc'] + list(nmudf.columns))])
        #slcdf = multicol_merge(slcdf, nmudf,left_index=True,right_index=True,how="left",validate="one_to_one")
        #slcdf['slc','n_mu_trks'] = slcdf['slc','n_mu_trks'].fillna(0)
        #del nmudf

        #save only primary muon with the highest energy
        #mu_mass = 0.1056583745 # GeV/c^2
        #mup_column = ('pfp','trk','rangeP','p_muon','','')
        #energy_column = ('pfp','trk','energy','','','')
        #mudf[energy_column] = mudf.apply(lambda row: p_to_energy(row[mup_column], mu_mass), axis=1)
        #mudf = mudf.sort_values(pfpdf.pfp.index.names[:-1] + [energy_column]).groupby(level=[0,1]).nth(-1)
        # drop all columns that are from shw attributes
        #mudf = mudf.drop('shw',axis=1,level=1)
        #mudf.columns = mudf.columns.set_levels(['primmu'],level=0)
        #slcdf = multicol_merge(slcdf, mudf.droplevel(-1),left_index=True,right_index=True,how="left",validate="one_to_one")
        #del mudf

        #PROTON-LIKE TRACKS:
        #prodf = pfpdf[(pfpdf.pfp.trackScore > 0.5) & (pfpdf.pfp.trk.len > 0) & (pfpdf.pfp.trk.chi2pid.I2.chi2_muon > 20) & (pfpdf.pfp.trk.chi2pid.I2.chi2_proton < 85)]
        #nprodf = prodf.groupby(level=[0,1]).size().to_frame('n_proton_trks') #count number of tracks
        #nprodf.columns = pd.MultiIndex.from_tuples([tuple(['slc'] + list(nprodf.columns))])
        #slcdf = multicol_merge(slcdf, nprodf,left_index=True,right_index=True,how="left",validate="one_to_one")
        #slcdf['slc','n_proton_trks'] = slcdf['slc','n_proton_trks'].fillna(0)
        #del nprodf

        #save only primary proton with the highest energy
        #pro_mass = 0.9382720813 # GeV/c^2
        #prop_column = ('pfp','trk','rangeP','p_proton','','')
        #prodf[energy_column] = prodf.apply(lambda row: p_to_energy(row[prop_column], pro_mass), axis=1)
        #prodf = prodf.sort_values(pfpdf.pfp.index.names[:-1] + [energy_column]).groupby(level=[0,1]).nth(-1)
        # drop all columns that are from shw attributes
        #prodf = prodf.drop('shw',axis=1,level=1)
        #prodf.columns = prodf.columns.set_levels(['primproton'],level=0)
        #slcdf = multicol_merge(slcdf, prodf.droplevel(-1),left_index=True,right_index=True,how="left",validate="one_to_one")

        #del prodf

    # pre-selection cuts
    if applyPreselection:
        slcdf = slcdf[slcdf.slc.is_clear_cosmic==0]
        slcdf = slcdf[slcdf.slc.nu_score > 0.5]
        slcdf = slcdf[InFV(df=slcdf.slc.vertex, inzback=0, det=DETECTOR)]
    
    return slcdf 
