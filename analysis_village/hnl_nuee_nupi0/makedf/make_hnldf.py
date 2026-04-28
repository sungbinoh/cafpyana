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

def make_slcdf(f):

    slcdf = loadbranches(f["recTree"], presel_branches)

    
    return slcdf 

def ensure_lexsorted(frame, axis):
    """Ensure DataFrame axes are fully lexsorted when using MultiIndex.
    
    This avoids pandas PerformanceWarning about indexing past lexsort depth.
    
    Parameters
    ----------
    frame : pandas.DataFrame
        DataFrame to check and sort if needed.
    axis : int
        Axis to check (0 for index, 1 for columns).
    
    Returns
    -------
    pandas.DataFrame
        DataFrame with sorted index/columns if MultiIndex, otherwise unchanged.
    """
    # axis: 0 -> index, 1 -> columns
    idx = frame.index if axis == 0 else frame.columns
    if isinstance(idx, pd.MultiIndex) and getattr(idx, "lexsort_depth", 0) < idx.nlevels:
        # sort by all levels (returns a new frame)
        return frame.sort_index(axis=axis)
    return frame

def add_variables(df, beam_x: float = -74.0, beam_y: float = 0.0):
    """Add derived kinematic columns

    Columns added:

    Per shower (primshw / secshw if present):
      (shw, 'shw', 'angle_z', '', '', '')          -- angle w.r.t. beam axis z [deg]

    Slice-level:
      ('slc', 'vertex', 'transverse_distance_beam_2', '', '', '')
                                                    -- d_T^2 = (vtx_x-beam_x)^2 + (vtx_y-beam_y)^2 [cm^2]


    Two-shower only (when secshw columns are present):
      ('slc', 'm_alt', '', '', '', '')              -- transverse mass of the shower pair [MeV]
                                                       m_alt = sqrt(2*ET1*ET2*(1-cos_theta))*1000

    Parameters
    ----------
    df : pd.DataFrame
        Slice-level DataFrame produced by topology.make_topo_df.
    beam_x, beam_y : float
        Beam centre x and y position [cm]. Default: (-74, 0).

    Returns
    -------
    pd.DataFrame
        Same DataFrame with new columns appended.
    """
    df = ensure_lexsorted(df, axis=1)

    def _angle_z(shw):
        dx = df[(shw, 'shw', 'dir', 'x', '', '')].values
        dy = df[(shw, 'shw', 'dir', 'y', '', '')].values
        dz = df[(shw, 'shw', 'dir', 'z', '', '')].values
        n  = np.sqrt(dx**2 + dy**2 + dz**2)
        with np.errstate(invalid='ignore'):
            return np.degrees(np.arccos(np.clip(dz / np.where(n > 0, n, np.nan), -1, 1)))

    for shw in ('primshw', 'secshw'):
        if (shw, 'shw', 'dir', 'z', '', '') in df.columns:
            df[(shw, 'shw', 'angle_z', '', '', '')] = _angle_z(shw)

    vtx_x_col = ('slc', 'vertex', 'x', '', '', '')
    vtx_y_col = ('slc', 'vertex', 'y', '', '', '')
    vtx_z_col = ('slc', 'vertex', 'z', '', '', '')
    if vtx_x_col in df.columns and vtx_y_col in df.columns:
        df[('slc', 'vertex', 'transverse_distance_beam_2', '', '', '')] = (
            (df[vtx_x_col].values - beam_x) ** 2 +
            (df[vtx_y_col].values - beam_y) ** 2
        )

    has_prim = ('primshw', 'shw', 'bestplane_energy', '', '', '') in df.columns
    has_sec  = ('secshw',  'shw', 'bestplane_energy', '', '', '') in df.columns

    if has_prim and has_sec:
        E1 = df[('primshw', 'shw', 'bestplane_energy', '', '', '')].values
        E2 = df[('secshw',  'shw', 'bestplane_energy', '', '', '')].values

        def _unit(shw):
            dx = df[(shw, 'shw', 'dir', 'x', '', '')].values
            dy = df[(shw, 'shw', 'dir', 'y', '', '')].values
            dz = df[(shw, 'shw', 'dir', 'z', '', '')].values
            n  = np.sqrt(dx**2 + dy**2 + dz**2)
            n  = np.where(n > 0, n, np.nan)
            return dx/n, dy/n, dz/n

        ux1, uy1, uz1 = _unit('primshw')
        ux2, uy2, uz2 = _unit('secshw')

        with np.errstate(invalid='ignore'):
            ET1 = E1 * np.sqrt(np.clip(1 - uz1**2, 0, None))
            ET2 = E2 * np.sqrt(np.clip(1 - uz2**2, 0, None))
            cos_theta = np.clip(ux1*ux2 + uy1*uy2 + uz1*uz2, -1, 1)
            df[('slc', 'm_alt', '', '', '', '')] = (
                np.sqrt(2 * ET1 * ET2 * (1 - cos_theta)) * 1000  # GeV -> MeV
            )

    return df

def make_hnldf(f, applyPreselection=False, savePfp=False, trackScore=0.51):
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det
    DETECTOR = "SBND"

    slcdf = loadbranches(f["recTree"], slcbranches + barycenterFMbranches + correctedflashbranches)
    slcdf = slcdf.rec

    # ── Preselection on slices (reduces rows before any merge) ───────────────
    if applyPreselection:
        slcdf = slcdf[slcdf.slc.is_clear_cosmic == 0]
        slcdf = slcdf[slcdf.slc.nu_score > 0.5]
        slcdf = slcdf[InFV(df=slcdf.slc.vertex, inzback=0, det=DETECTOR)]

    # ── Load PFPs only after slcdf is reduced ────────────────────────────────
    pfpdf = make_pfpdf(f)
    pfpdf = pfpdf.drop('pfochar', axis=1, level=1)

    # ── Filter pfpdf to surviving slices ─────────────────────────
    pfpdf = pfpdf[pfpdf.index.droplevel(-1).isin(slcdf.index)]

    if savePfp:
        slcdf = multicol_merge(slcdf,
                               pfpdf.reset_index(level='rec.slc.reco.pfp..index'),
                               left_index=True, right_index=True,
                               how="left", validate="one_to_many")
        pfp_idx_col = next(
            c for c in slcdf.columns
            if c == 'rec.slc.reco.pfp..index' or
               (isinstance(c, tuple) and len(c) > 0 and c[0] == 'rec.slc.reco.pfp..index')
        )
        slcdf = slcdf.set_index(pfp_idx_col, append=True)

    else:
        # ── PFP sanity filter ─────────────────────────────────────────────────
        pfpdf = pfpdf[
            (pfpdf.pfp.trackScore > 0) &
            (pfpdf.pfp.shw.bestplane_energy > 0) &
            (pfpdf.pfp.shw.start.x != -999) &
            (pfpdf.pfp.shw.start.y != -999) &
            (pfpdf.pfp.shw.start.z != -999) &
            (pfpdf.pfp.trk.len > 0) &
            (pfpdf.pfp.shw.bestplane_dEdx > 0)
            
        ]

        shw_mask = pfpdf.pfp.trackScore < trackScore
        trk_mask = pfpdf.pfp.trackScore > trackScore

        # ── Count tracks first ────────────────────────────
        ntrkdf = pfpdf[trk_mask].groupby(level=[0, 1]).size().to_frame('n_trks')
        ntrkdf.columns = pd.MultiIndex.from_tuples([('slc', 'n_trks')])
        slcdf = multicol_merge(slcdf, ntrkdf, left_index=True, right_index=True,
                               how="left", validate="one_to_one")
        slcdf['slc', 'n_trks'] = slcdf['slc', 'n_trks'].fillna(0)
        del ntrkdf

        # ── Topology cut: 0 tracks ────────────────────────────────────────────
        slcdf = slcdf[slcdf['slc', 'n_trks'] == 0]
        pfpdf = pfpdf[pfpdf.index.droplevel(-1).isin(slcdf.index)]
        shw_mask = pfpdf.pfp.trackScore < trackScore

        # ── Count showers ─────────────────────────────────────────────────────
        nshwdf = pfpdf[shw_mask].groupby(level=[0, 1]).size().to_frame('n_shws')
        nshwdf.columns = pd.MultiIndex.from_tuples([('slc', 'n_shws')])
        slcdf = multicol_merge(slcdf, nshwdf, left_index=True, right_index=True,
                               how="left", validate="one_to_one")
        slcdf['slc', 'n_shws'] = slcdf['slc', 'n_shws'].fillna(0)
        del nshwdf

        # ── Topology cut: 1 or 2 showers ─────────────────────────────────────
        slcdf = slcdf[slcdf['slc', 'n_shws'].isin([1, 2])]
        pfpdf = pfpdf[pfpdf.index.droplevel(-1).isin(slcdf.index)]
        shw_mask = pfpdf.pfp.trackScore < trackScore

        # ── Primary shower (highest energy) ───────────────────────────────────
        sort_cols = pfpdf.pfp.index.names[:-1] + [('pfp', 'shw', 'bestplane_energy', '', '', '')]
        shwdf = (pfpdf[shw_mask]
                 .sort_values(sort_cols)
                 .groupby(level=[0, 1])
                 .nth(-1)
                 .drop('trk', axis=1, level=1))
        shwdf.columns = shwdf.columns.set_levels(['primshw'], level=0)
        slcdf = multicol_merge(slcdf, shwdf.droplevel(-1), left_index=True,
                               right_index=True, how="left", validate="one_to_one")
        del shwdf

        # ── Secondary shower (2nd highest energy) ─────────────────────────────
        shwsecdf = (pfpdf[shw_mask]
                    .sort_values(sort_cols)
                    .groupby(level=[0, 1])
                    .nth(-2)
                    .drop('trk', axis=1, level=1))
        shwsecdf.columns = shwsecdf.columns.set_levels(['secshw'], level=0)
        slcdf = multicol_merge(slcdf, shwsecdf.droplevel(-1), left_index=True,
                               right_index=True, how="left", validate="one_to_one")
        del shwsecdf
        
    # ── Add derived variables ─────────────────────────────────────────────────
    slcdf = add_variables(slcdf)
    
    return slcdf