from . import getsyst

g4_systematics = [
    "reinteractions_piminus_Geant4",
    "reinteractions_piplus_Geant4",
    "reinteractions_proton_Geant4"
]

def g4syst(f, nuind, slim=False):
    return getsyst.getsyst(f, g4_systematics, nuind, slim=slim)

