"""
Kaon analysis data frame maker
"""

import pandas as pd

import makedf.makedf as makedf
import pyanalib.pandas_helpers as pd_helpers
import makedf.util as util


KPDG = {
    'kplus': makedf.PDG["kaon_p"][0],
    'kzero': makedf.PDG["kaon_0"][0]
}
KMASS = {
    'kplus': makedf.PDG["kaon_p"][2],
    'kzero': makedf.PDG["kaon_0"][2]
}

# TODO pick a better number
TRUE_KE_CUT = 0.


def make_kaon_mcdf(f: pd.DataFrame) -> pd.DataFrame:
    mcdf = makedf.make_mcdf(f)
    mcprimdf = makedf.make_mcprimdf(f)
    while mcprimdf.columns.nlevels > 2:
        mcprimdf.columns = mcprimdf.columns.droplevel(0)
    mcprimdf.index = mcprimdf.index.rename(mcdf.index.names[:2] + mcprimdf.index.names[2:])

    # add kaon info: Number above threshold & primary info
    for kname in ('kplus', 'kzero'):
        ke = mcprimdf[mcprimdf.pdg==KPDG[kname]].genE - KMASS[kname] 
        mcdf = pd_helpers.multicol_add(mcdf, ((mcprimdf.pdg==KPDG[kname]) \
                                              & (ke > TRUE_KE_CUT)).groupby(level=[0,1]).sum().rename(f'n{kname}'))
        kdf = mcprimdf[mcprimdf.pdg==KPDG[kname]].sort_values(mcprimdf.index.names[:2] + [("genE", "")]).groupby(level=[0,1]).last()
        kdf.columns = pd.MultiIndex.from_tuples([tuple([kname] + list(c)) for c in kdf.columns])
        mcdf = pd_helpers.multicol_merge(mcdf, kdf, left_index=True, right_index=True, how="left", validate="one_to_one")

    return mcdf
