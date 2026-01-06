import uproot
import numpy as np
import pandas as pd
import awkward as ak

def get_env(f):
    if "env" not in f or "envtree" not in f["env"]:
        return pd.DataFrame({"caf_name": [np.nan]})

    env_tree = f["env"]["envtree"]

    # Be defensive in case branches are missing
    if "key" not in env_tree or "value" not in env_tree:
        return pd.DataFrame({"caf_name": [np.nan]})

    keys = env_tree["key"].array(library="np")
    vals = env_tree["value"].array(library="np")

    mask = (keys == "output")
    output_vals = vals[mask]
    first_output = output_vals[0] if getattr(output_vals, "size", 0) else ""

    value = np.nan if first_output is None else first_output
    return pd.DataFrame({"caf_name": [value]})
