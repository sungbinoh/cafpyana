import pandas as pd
import numpy as np


def mcstat(evt_df: pd.DataFrame,
                   hdr_df: pd.DataFrame,
                   n_universes: int = 100):
    """
    Generate Poisson-fluctuated statistical universes for MC, to be treated as multisim uncertainties

    - A deterministic seed for each event is constructed from its metadata (run, sub-run, event number and slice ID)

    - For every requested universe an additional seed that depends on the universe index is added, 
      to ensure that each universe is statistically independent

    - In every universe a Poisson random number with mean mu = 1 is drawn for each event

    - Generated universes are appended to evt_df under the branch name ("MCstat", "univ_{universe_index}")

    Inputs:
    - evt_df : pandas df. Must contain a ``slc.self`` column
    - hdr_df : pandas df. Must contain run, subrun and evt columns
    - n_universes : int, default 100. Number of independent statistical universes to generate

    Returns:
    - evt_df with the newly added weight columns 
    - np array with shape (n_universes, len(evt_df)) holding the generated Poisson weights
    """

    # build a deterministic seed per event from its metadata
    meta_seeds = []
    for i in range(len(evt_df)):
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

    # extra safety check -- seeds should be unique!
    assert len(meta_seeds) == len(set(meta_seeds))

    # generate universes
    n_universes = 100
    MCstat_univ_events = np.zeros((n_universes, len(evt_df)))
    poisson_mean = 1.0

    # get Poisson weights and save to "MCstat.univ_"
    # dummy df to hold the weights -- iterative inserting causes PerformanceWarning
    mcstat_univ_cols = pd.MultiIndex.from_product(
        [["MCstat"], [f"univ_{i}" for i in range(n_universes)], [""], [""], [""], [""], [""], [""]],
    )
    mcstat_univ_wgt = pd.DataFrame(
        1.0,
        index=evt_df.index,
        columns=mcstat_univ_cols,
    )

    for uidx in range(n_universes):
        universe_seed = hash(f"universe_{uidx}") % (2**32)
        
        poisson_weights = []
        for sidx, meta_seed in enumerate(meta_seeds):
            # Combine universe seed with event seed for unique randomness -- per event, per universe
            combined_seed = (universe_seed + meta_seed) % (2**32)
            np.random.seed(combined_seed)
            
            poisson_val = np.random.poisson(poisson_mean)
            poisson_weights.append(poisson_val)
        
        mcstat_univ_wgt[("MCstat", "univ_{}".format(uidx), "", "", "", "", "", "")] = np.array(poisson_weights)
        MCstat_univ_events[uidx, :] = np.array(poisson_weights)

    evt_df = evt_df.join(mcstat_univ_wgt)
    return evt_df, MCstat_univ_events