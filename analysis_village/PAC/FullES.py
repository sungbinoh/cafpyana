import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, sys
# Add the head direcoty to sys.path
workspace_root = os.getcwd()  
sys.path.insert(0, workspace_root + "/../../")

# import this repo's classes
from pyanalib.split_df_helpers import *
import pyanalib.pandas_helpers as ph

def main():
    """Main analysis pipeline."""
    f = "test.df"

    dfs = load_dfs(f, ['mcnu'])
    mcdf = dfs['mcnu']

    for c in mcdf.columns:
        print(mcdf[c])

if __name__ == "__main__":
    main()
