from . import getsyst
import pandas as pd

# Regen systematic variations
regen_systematics = [
    'expskin_Flux',
    'horncurrent_Flux',
    'kminus_Flux',
    'kplus_Flux',
    'kzero_Flux',
    'nucleoninexsec_Flux',
    'nucleonqexsec_Flux',
    'nucleontotxsec_Flux',
    'piminus_Flux',
    'pioninexsec_Flux',
    'pionqexsec_Flux',
    'piontotxsec_Flux',
    'piplus_Flux'
]

def bnbsyst(f, nuind, multisim_nuniv=250, slim=False):
    bnbwgtdf = getsyst.getsyst(f, regen_systematics, nuind, multisim_nuniv=multisim_nuniv, slim=slim, slimname="Flux")

    if slim:  # keep only the multiplied "Flux.univ_" columns
        flux_cols = [c for c in bnbwgtdf.columns if c[0] == "Flux"]
        bnbwgtdf = bnbwgtdf[flux_cols]
        
    return bnbwgtdf

