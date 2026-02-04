import numpy as np
import pandas as pd
from tqdm import tqdm

def get_MCstat_unc(evt_df, hdr_df, n_universes=100):
    """
    Create a unique seed based on event metadata
    Using a hash function that's deterministic
    """

    meta_seeds = []
    for i in tqdm(range(len(evt_df))):
        this_hdr_df = hdr_df.loc[evt_df.reset_index(level=[2]).index[i]]
        runno = this_hdr_df.run
        subrunno = this_hdr_df.subrun
        evtno = this_hdr_df.evt
        slcid = evt_df.loc[evt_df.index[i]].slc.self
        unique_seed = hash(f"run_{runno}_subrun_{subrunno}_evt_{evtno}_slcid_{slcid}") % (2**32)  # Ensure it's a 32-bit integer
        if unique_seed in meta_seeds:
            print("duplicate seed found", unique_seed)
            break
        meta_seeds.append(unique_seed)
    # make sure the seeds are unique!
    assert len(meta_seeds) == len(set(meta_seeds))

    # generate universes
    MCstat_univ_events = np.zeros((n_universes, len(evt_df)))
    poisson_mean = 1.0
    # get Poisson weights and save to "MCstat.univ_"
    # TODO: generalize column level padding
    mcstat_univ_cols = pd.MultiIndex.from_product(
        [["MCstat"], [f"univ_{i}" for i in range(n_universes)], [""], [""], [""], [""], [""]],
    )
    mcstat_univ_wgt = pd.DataFrame(
        1.0,
        index=evt_df.index,
        columns=mcstat_univ_cols,
    )
    for uidx in range(n_universes):
        universe_seed = hash(f"universe_{uidx}") % (2**32)
        
        poisson_weights = []
        for meta_seed in meta_seeds:
            # Combine universe seed with event seed for unique randomness -- per event, per universe
            combined_seed = (universe_seed + meta_seed) % (2**32)
            np.random.seed(combined_seed)
            
            poisson_val = np.random.poisson(poisson_mean)
            poisson_weights.append(poisson_val)
        
        mcstat_univ_wgt[("MCstat", "univ_{}".format(uidx), "", "", "", "", "")] = np.array(poisson_weights)
        MCstat_univ_events[uidx, :] = np.array(poisson_weights)

    evt_df = evt_df.join(mcstat_univ_wgt)

    return evt_df, MCstat_univ_events
