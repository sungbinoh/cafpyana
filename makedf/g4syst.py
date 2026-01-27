from . import getsyst

g4_systematics = [
     'reinteractions_kminus_Geant4',
    'reinteractions_kplus_Geant4',
    'reinteractions_neutron_Geant4',
    'reinteractions_piminus_Geant4',
    'reinteractions_piplus_Geant4',
    'reinteractions_proton_Geant4'
]


def g4syst(f, nuind, multisim_nuniv=250, slim=False):
    g4wgtdf = getsyst.getsyst(f, g4_systematics, nuind, multisim_nuniv=multisim_nuniv, slim=slim, slimname="G4")

    if slim:  # keep only the multiplied "g4.univ_" columns
        g4_cols = [c for c in g4wgtdf.columns if c[0] == "G4"]
        g4wgtdf = g4wgtdf[g4_cols]
        
    return g4wgtdf

