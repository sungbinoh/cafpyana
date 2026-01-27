from . import getsyst

g4_systematics = [
     'reinteractions_kminus_Geant4',
    'reinteractions_kplus_Geant4',
    'reinteractions_neutron_Geant4',
    'reinteractions_piminus_Geant4',
    'reinteractions_piplus_Geant4',
    'reinteractions_proton_Geant4'
]

def g4syst(f, nuind):
    return getsyst.getsyst(f, g4_systematics, nuind)

