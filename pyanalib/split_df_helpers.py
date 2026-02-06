"""
Helper functions for working with the HDF5 df files.  

Each dataset (identified by a *key*) is stored as a set of smaller dfs (one per split),
so that very large samples can be handled in chunks.  
"""

import pandas as pd
from os import path

# utils for the case where events are split into keys within a single file
def get_n_split(file):
    this_split_df = pd.read_hdf(file, key="split")
    this_n_split = this_split_df.n_split.iloc[0]
    return this_n_split

def print_keys(file):
    with pd.HDFStore(file, mode='r') as store:
        keys = store.keys()       # list of all keys in the file
        print("Keys:", keys)

def load_dfs(file, keys2load, n_max_concat=100):
    out_df_dict = {}
    this_n_keys = get_n_split(file)
    n_concat = min(n_max_concat, this_n_keys)
    for key in keys2load:
        dfs = []  # collect all splits for this key
        for i in range(n_concat):
            this_df = pd.read_hdf(file, key=f"{key}_{i}")
            dfs.append(this_df)
        out_df_dict[key] = pd.concat(dfs, ignore_index=False)

    return out_df_dict


# utils for the case where events are split into different files (which also has split keys like above)
def load_and_concat_mc_dfs(
    file_dir,
    chunk_tags=None,
    df_tag="",
    keys2load=['hdr', 'evt'],
    n_max_concat=3,
    sub_dir="MC",
    sample_dir="BNB_cosmics",
):
    """
    Loops over chunk_tags, loads mc_hdr_df and mc_evt_df, then concats them.
    Keeps the first level of multiindex value unique (__ntuple) by bumping it up by 
    the previous dfs lengths' summed.
    """

    df_lists = {k:[] for k in keys2load}
    ntuple_offset = 0

    for tag in chunk_tags:
        mc_file = path.join(file_dir, sub_dir, sample_dir, tag+df_tag+".df")
        mc_n_split = get_n_split(mc_file)
        print(f"Reading file with tag {tag}, mc_n_split: {mc_n_split}")
        mc_dfs = load_dfs(mc_file, keys2load, n_max_concat=n_max_concat)

        # Make __ntuple unique by offsetting first level index (if exists)
        for df_key in keys2load:
            df = mc_dfs[df_key]
            if isinstance(df.index, pd.MultiIndex):
                levels = list(df.index.levels)
                names = df.index.names
                # __ntuple should be at level 0
                if "__ntuple" in names:
                    idx_loc = names.index("__ntuple")
                else:
                    idx_loc = 0
                # Add offset
                new_tuples = []
                for tup in df.index:
                    tup = list(tup)
                    tup[idx_loc] = tup[idx_loc] + ntuple_offset
                    new_tuples.append(tuple(tup))
                df.index = pd.MultiIndex.from_tuples(new_tuples, names=names)
            else:
                # If single index, and it's __ntuple
                if df.index.name == "__ntuple":
                    df.index = df.index + ntuple_offset

            df_lists[df_key].append(df)

        # bump index for next file
        if isinstance(df.index, pd.MultiIndex):
            ntuple_vals = df.index.get_level_values(0)
        else:
            ntuple_vals = df.index
        ntuple_offset += ntuple_vals.max()+1

    concat_dfs = {k:pd.concat(df_lists[k], axis=0, sort=False) for k in df_lists.keys()}
    return concat_dfs
