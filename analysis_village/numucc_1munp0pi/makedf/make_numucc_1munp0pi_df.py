from makedf.makedf import *
from pyanalib.pandas_helpers import *
from makedf.util import *

pd.options.mode.chained_assignment = None

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

def InTrueContain(data):
    xmax = 190.
    zmin = 10.
    zmax = 490.
    ymax = 190.
    return (np.abs(data.x) < xmax) & (np.abs(data.y) < ymax) & (data.z > zmin) & (data.z < zmax)

def is_true_numucc(mcnudf):
    iscc = mcnudf.iscc
    isnumu = (np.abs(mcnudf.pdg) == 14)
    return iscc & isnumu

def is_signal(mcnudf, mcnuprimdf, mcprimvisEbranches, mcprimdaughtersdf, tpartdf):
    is_numucc = is_true_numucc(mcnudf)
    isfv = InFV_nohiyz(mcnudf.position)

    ## muon len cut
    muon_len_cut = 50. ## cm
    first_prim_df = (mcnuprimdf.groupby(level=[0,1], sort=False).head(1).droplevel(2)) ## by the is_numucc this will be muon
    is_muon_len_pass = first_prim_df.length > muon_len_cut

    ## energy depo for p, pipm, and gamma
    p_visE_cut = 0.05 ## GeV
    pipm_visE_cut = 0.025 ## GeV
    gamma_visE_cut = 0.025 ## GeV
    total_visE = mcnuprimdf.I0.I2.visE + mcnuprimdf.I1.I2.visE
    mcnuprimdf[("tot_visE", "", "")] = total_visE
    
    n_p = ((mcnuprimdf.pdg==2212) & (mcnuprimdf.tot_visE > p_visE_cut)).groupby(level=[0,1]).sum()
    n_pipm = ((np.abs(mcnuprimdf.pdg)==211) & (mcnuprimdf.tot_visE > pipm_visE_cut)).groupby(level=[0,1]).sum()
    n_gamma = ((np.abs(mcnuprimdf.pdg)==22) & (mcnuprimdf.tot_visE > gamma_visE_cut)).groupby(level=[0,1]).sum()

    is_pass_n_p = n_p > 0
    is_pass_n_pipm = n_pipm == 0
    is_pass_n_gamma = n_gamma == 0

    ## energy depo for pizero
    pizero = mcnuprimdf[mcnuprimdf.pdg==111]
    pizero_daughters = mcprimdaughtersdf[mcprimdaughtersdf.index.droplevel(-1).isin(pizero.index)]

    pizero_daughter_tparts = tpartdf[
        tpartdf.reset_index().set_index(['entry', 'G4ID']).index.isin(
            pd.MultiIndex.from_frame(pizero_daughters.reset_index()[['entry', 'daughters']])
        )
    ]

    total_visE_pizero_daughter_tparts = pizero_daughter_tparts.plane.I0.I2.visE + pizero_daughter_tparts.plane.I1.I2.visE
    pizero_daughter_tparts[("tot_visE", "", "", "")] = total_visE_pizero_daughter_tparts

    #### match tot_visE to pizero_daughters
    mapper = pizero_daughter_tparts.reset_index().set_index(['entry', 'G4ID'])['tot_visE']
    pizero_daughters['tot_visE'] = mapper.reindex(
        pd.MultiIndex.from_frame(pizero_daughters.reset_index()[['entry', 'daughters']])
    ).to_numpy()

    n_pizero_gamma = ((pizero_daughters.tot_visE > gamma_visE_cut)).groupby(level=[0,1]).sum() ## add all photons for a neutrino interaction
    mcnudf['n_pizero_gamma'] = n_pizero_gamma
    mcnudf['n_pizero_gamma'] = mcnudf['n_pizero_gamma'].fillna(0)
    n_pizero_gamma = mcnudf.n_pizero_gamma

    is_pass_n_pizero_gamma = n_pizero_gamma == 0

    ## containment cut
    #### prim muon
    is_muon_contained = InTrueContain(first_prim_df.start) & InTrueContain(first_prim_df.end)
    #### I don't apply containment to secondaries of the prim muon since they are actually secondaries of the neutrino interaction == primes
    
    ## I will apply containment requirement for particles passing the visE criteria
    #### primary protons
    n_uncontained_p = ((mcnuprimdf.pdg==2212) & (mcnuprimdf.tot_visE > p_visE_cut) & ((~InTrueContain(mcnuprimdf.start)) | (~InTrueContain(mcnuprimdf.end)))).groupby(level=[0,1]).sum()
    is_pass_n_uncontained_p = n_uncontained_p == 0

    #### secondaries of the primatry protons
    primp = mcnuprimdf[(mcnuprimdf.pdg==2212) & (mcnuprimdf.tot_visE > p_visE_cut)]
    primp_daughters = mcprimdaughtersdf[mcprimdaughtersdf.index.droplevel(-1).isin(primp.index)]
    
    primp_daughter_tparts = tpartdf[
        tpartdf.reset_index().set_index(['entry', 'G4ID']).index.isin(
            pd.MultiIndex.from_frame(primp_daughters.reset_index()[['entry', 'daughters']])
        )
    ]
    primp_daughter_tparts['contain'] = InTrueContain(primp_daughter_tparts.start) & InTrueContain(primp_daughter_tparts.end)
    primp_daughter_tparts['charged'] = (primp_daughter_tparts.pdg == 2212) | (np.abs(primp_daughter_tparts.pdg) == 211) | (np.abs(primp_daughter_tparts.pdg) == 321)## consider only proton, charged pion and charged kaon secondaries
    mapper = primp_daughter_tparts.reset_index().set_index(['entry', 'G4ID'])['contain']
    primp_daughters['contain'] = mapper.reindex(
        pd.MultiIndex.from_frame(primp_daughters.reset_index()[['entry', 'daughters']])
    ).to_numpy()
    mapper = primp_daughter_tparts.reset_index().set_index(['entry', 'G4ID'])['charged']
    primp_daughters['charged'] = mapper.reindex(
        pd.MultiIndex.from_frame(primp_daughters.reset_index()[['entry', 'daughters']])
    ).to_numpy()
 
    n_bad_primp_daughters = ((primp_daughters.charged) & (~primp_daughters.contain)).groupby(level=[0,1]).sum()
    mcnudf['n_bad_primp_daughters'] = n_bad_primp_daughters
    mcnudf['n_bad_primp_daughters'] = mcnudf['n_bad_primp_daughters'].fillna(0)
    n_bad_primp_daughters = mcnudf.n_bad_primp_daughters
    is_pass_n_bad_primp_daughters = n_bad_primp_daughters == 0

    ## final decision
    return is_numucc & isfv & is_muon_len_pass & is_pass_n_p & is_pass_n_pipm & is_pass_n_gamma & is_pass_n_pizero_gamma & is_muon_contained & is_pass_n_uncontained_p & is_pass_n_bad_primp_daughters
    
def make_numucc_1munp0pi_mcnu_df(f):
    hdrdf = make_hdrdf(f)
    tpartdf = loadbranches(f["recTree"], trueparticlebranches).rec.true_particles
    #tpartdf = tpartdf[tpartdf.startE > 0.]

    mcnudf = make_mcnudf(f)
    
    mcnuprimdf = make_mcprimdf(f).rec.mc.nu.prim
    mcnuprimdf.columns = pd.MultiIndex.from_tuples(
        [tup + ("",) for tup in mcnuprimdf.columns],
        names=(list(mcnuprimdf.columns.names or []) + [""])
    )
    mcprimvisEbranches = make_mcprimvisEdf(f).rec.mc.nu.prim.plane
    mcnuprimdf = pd.concat([mcnuprimdf, mcprimvisEbranches], axis=1)
    
    mcprimdaughtersdf = make_mcprimdaughtersdf(f).rec.mc.nu.prim
    
    mcnudf['is_signal'] = is_signal(mcnudf, mcnuprimdf, mcprimvisEbranches, mcprimdaughtersdf, tpartdf)

    return mcnudf

    
