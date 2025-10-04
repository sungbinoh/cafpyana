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
    bnbwgtdf = getsyst.getsyst(f, regen_systematics, nuind, multisim_nuniv=multisim_nuniv)

    # multiply all knobs and save to "Flux.univ_"
    # dummy df to hold the product of all syst knobs -- iterative inserting causes PerformanceWarning
    flux_cols = pd.MultiIndex.from_product(
        [["Flux"], [f"univ_{i}" for i in range(multisim_nuniv)]],
    )
    flux_wgt = pd.DataFrame(
        1.0,
        index=bnbwgtdf.index,
        columns=flux_cols,
    )
    for syst in regen_systematics:
        flux_wgt *= bnbwgtdf[syst].to_numpy()
    bnbwgtdf = pd.concat([bnbwgtdf, flux_wgt], axis=1)

    if slim:  # keep only the multiplied "Flux.univ_" columns
        bnbwgtdf = bnbwgtdf[flux_cols]
        
    return bnbwgtdf

