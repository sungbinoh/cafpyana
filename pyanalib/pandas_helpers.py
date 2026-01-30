import pandas as pd
import numpy as np
import awkward as ak

def broadcast(v, df):
    for vi, ii in zip(v.index.names, df.index.names):
        if vi != ii:
            raise ValueError("Value index (%s) does not match index (%s)." % (str(vi), str(ii)))
    if len(v.index.names) > len(df.index.names):
        raise ValueError("Value index too long.")
    if len(v.index.names) == len(df.index.names):
        return v

    rpt = df.groupby(level=list(range(v.index.nlevels))).size()
    has_value = v.index.intersection(rpt.index)
    v_rpt = np.repeat(v.loc[has_value].values, rpt)

    return pd.Series(v_rpt, df.index).rename(v.name) 

def multicol_concat(lhs, rhs):
    # Fix the columns
    lhs_col = lhs.columns
    rhs_col = rhs.columns

    nlevel = max(lhs_col.nlevels, rhs_col.nlevels)

    def pad(c):
       return tuple(list(c) + [""]*(nlevel - len(c))) 

    lhs.columns = pd.MultiIndex.from_tuples([pad(c) for c in lhs_col])
    rhs.columns = pd.MultiIndex.from_tuples([pad(c) for c in rhs_col])

    return pd.concat([lhs, rhs], axis=1)

def multicol_add(df, s, default=None, **panda_kwargs):
    # if both the series and the df is one level, we can do a simple join()
    if isinstance(s.name, str) and df.columns.nlevels == 1:
        return df.join(s, **panda_kwargs)

    if isinstance(s.name, str):
        s.name = (s.name,)

    nlevel = max(df.columns.nlevels, len(s.name))
    def pad(col, c=""):
       return tuple(list(col) + [c]*(nlevel - len(col))) 

    if df.columns.nlevels < nlevel:
        df.columns = pd.MultiIndex.from_tuples([pad(c) for c in df.columns])

    setname = None
    if len(s.name) < nlevel:
        setname = pad(s.name)
        s.name = pad(s.name, ".")

    # Work around bugs in pandas

    # Reindex s to match index of df
    s = s.reindex(df.index)
    # fill default
    if default is not None:
        s.fillna(default, inplace=True)

    ret = df.join(s,  **panda_kwargs)
    # Another pandas bug work around -- we can't pad with ""'s. So pad with "."'s and map "." -> ""
    if setname is not None:
        ret.columns = pd.MultiIndex.from_tuples([c if i != len(ret.columns) - 1 else setname for i,c in enumerate(ret.columns)])

    return ret

def multicol_merge(lhs, rhs, **panda_kwargs):
    # Fix the columns
    lhs_col = lhs.columns
    rhs_col = rhs.columns

    nlevel = max(lhs_col.nlevels, rhs_col.nlevels)

    def pad(c):
       nc = 1 if isinstance(c, str) else len(c)
       c0 = [c] if isinstance(c, str) else list(c)
       return tuple(c0 + [""]*(nlevel - nc)) 

    lhs.columns = pd.MultiIndex.from_tuples([pad(c) for c in lhs_col])
    rhs.columns = pd.MultiIndex.from_tuples([pad(c) for c in rhs_col])

    return lhs.merge(rhs, **panda_kwargs)

def _detect_vectors(branch, tree_keys):
    ret = []
    hierarchy = branch.split(".")
    for i in range(len(hierarchy)):
        subbranch = ".".join(hierarchy[:i+1])
        lenbranch = subbranch + "..length"
        if lenbranch in tree_keys:
            ret.append(subbranch)
    return ret

def idarray(ids, lens):
    return np.repeat(ids.values, lens.values)

def loadbranches(tree, branches, **uprargs):
    vectors = []
    keys = list(tree.keys())
    for i,branch in enumerate(branches):
        this_vectors = _detect_vectors(branch, keys)
        if i == 0:
            vectors = this_vectors
        elif len(this_vectors) == 0: # This case is ok since it will automatically broadcast
            pass
        # All the branches must have the same vector structure for this to work
        elif vectors != this_vectors:
            raise ValueError("Branches %s and %s have different vector structures in the CAF." % (branches[0], branch))

    lengths = [ak.to_dataframe(tree.arrays([v+"..length"], library="ak", **uprargs), how="inner") for v in vectors]
    data = ak.to_dataframe(tree.arrays(branches, library="ak", **uprargs), how=None)

    # If there's no vectors, we can just return the top guy
    if len(lengths) == 0:
        df = data[0]
    else:
        tomerge = lengths + data
        # Otherwise, iteratively merge the branches
        df = tomerge[0]
        df.index.name = "entry"

        # handle the rest
        for i in range(1, len(tomerge)):
            thismerge = tomerge[i]
            v_ind = i - 1

            # Build the information in the right-hand table needed to do the join
            # The "upidx" will be matched to the index vector-by-vector
            for i in range(v_ind):
                thismerge[vectors[v_ind] + "..upidx" + str(i)] = idarray(df[vectors[i]+ "..index"], df[vectors[v_ind] + "..length"])

            # Inner join! Throw away rows in the right-hand with no match in the left-hand
            df = pd.merge(df, thismerge, how="inner",
                         left_on = ["entry"] + [v+"..index" for v in vectors[:v_ind]],
                         right_on = ["entry"] + [vectors[v_ind] + "..upidx" + str(i) for i in range(v_ind)],
                         validate="one_to_many")

            # Make sure no rows in the right-hand were dropped
            assert(df.shape[0] == thismerge.shape[0])

            # postprocess: build the index
            df[vectors[v_ind] + "..index"] = df.groupby(["entry"] + [v+"..index" for v in vectors[:v_ind]]).cumcount()

        # Set the index
        df.set_index([v+"..index" for v in vectors], append=True, verify_integrity=True, inplace=True)

        # Drop all the metadata info we don't need anymore
        df = df[branches]

    # Setup branch names so df reflects structure of CAF file
    bsplit = [b.split(".") for b in branches]
    # Replace any reserved names
    def unreserve(s):
        if s == "index":
            return "idx"
        if s[0].isdigit(): # make the name a legal field 
            return "I" + s
        return s

    bsplit = [[unreserve(s) for s in b] for b in bsplit]

    depth = max([len(b) for b in bsplit])

    def pad(b):
        return tuple(b + [""]*(depth - len(b)))

    df.columns = pd.MultiIndex.from_tuples([pad(b) for b in bsplit])

    return df

def pad_column_name(col, pad_ref): # TODO: merge with the above "pad" function
    if isinstance(pad_ref, pd.DataFrame):
        nlevels = pad_ref.columns.nlevels
    elif isinstance(pad_ref, int):
        nlevels = pad_ref
    else:
        raise ValueError("pad_ref must be a pandas DataFrame or an integer")
    ndummies = nlevels - len(col)
    extended_name = col + ("",) * ndummies
    return extended_name



