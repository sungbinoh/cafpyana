import pandas as pd
import os
import sys

workspace_root = os.getcwd()
sys.path.insert(0, workspace_root + "/../../")

import pyanalib.pandas_helpers as ph
from pyanalib.split_df_helpers import *

def debug(dfs, keys=['nu_E']):

    pd.options.display.max_rows =1000 
    pd.set_option("display.max_rows", 1000) 
    pd.options.display.min_rows =500 
    pd.set_option("display.min_rows", 500)

    set_0 = (dfs[0].drop_duplicates().sort_values(keys))
    set_1 = (dfs[1].drop_duplicates().sort_values(keys))

    print(set_0.iloc[30:50])
    print(set_1.iloc[30:50])
    print(common_ids.sort_values(keys))

def main():
    """Main analysis pipeline."""

    Nom = load_dfs("/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/SBND_SpringMC_Nom.df", ['mcnu'])['mcnu']
    Var = load_dfs("/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/SBND_SpringMC_0xSCE.df", ['mcnu'])['mcnu']
    debug([Nom, Var])
    
