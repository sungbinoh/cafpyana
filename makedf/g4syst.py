from . import getsyst

g4_systematics = [
    "reinteractions_kplus_Geant4",
    "reinteractions_kminus_Geant4",
    "reinteractions_neutron_Geant4",
    "reinteractions_piminus_Geant4",
    "reinteractions_piplus_Geant4",
    "reinteractions_proton_Geant4"
]

def g4syst(f, nuind, multisim_nuniv=250,slim=False):
    return getsyst.getsyst(f, g4_systematics, nuind, multisim_nuniv=multisim_nuniv, slim=slim)

