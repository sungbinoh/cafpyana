import uproot
import numpy as np
import pandas as pd
import awkward as ak
import hashlib

CHUNK_SIZE = "10 MB"

def _report_entry_start(report):
    """Get global entry start from an uproot iterate report object."""
    for attr in ("global_entry_start", "entry_start", "start"):
        if hasattr(report, attr):
            return int(getattr(report, attr))
    if isinstance(report, dict):
        for key in ("global_entry_start", "entry_start", "start"):
            if key in report:
                return int(report[key])
    return 0


def _read_wgt_metadata(f):
    wgt_names = list(f["globalTree"]["global/wgts/wgts.name"].arrays(library="np")["wgts.name"][0])
    wgt_types = f["globalTree"]["global/wgts/wgts.type"].arrays(library="np")["wgts.type"][0]
    wgt_nuniv = f["globalTree"]["global/wgts/wgts.nuniv"].arrays(library="np")["wgts.nuniv"][0]
    return wgt_names, wgt_types, wgt_nuniv


def _build_isyst_map(wgt_nuniv):
    isyst_map = pd.Series(np.repeat(list(range(len(wgt_nuniv))), wgt_nuniv), name="isyst")
    isyst_map.index.name = "iwgt"
    return isyst_map


def _iter_filtered_wgts(rec_tree, nuniv, isyst_map, needed_syst_indices, step_size):
    """Yield filtered weight DataFrames chunk-by-chunk to cap peak memory."""
    for batch, report in rec_tree.iterate(
        filter_name=["rec.mc.nu.wgt.univ"],
        library="ak",
        step_size=step_size,
        report=True,
    ):
        entry_offset = _report_entry_start(report)
        wgts = ak.to_dataframe(batch["rec.mc.nu.wgt.univ"], how=None)[0]
        wgts["inu"] = wgts.index.get_level_values(1) // nuniv
        wgts["iwgt"] = wgts.index.get_level_values(1) % nuniv
        wgts = wgts.reset_index()
        wgts["entry"] = wgts["entry"] + entry_offset
        wgts = wgts.set_index(["entry", "inu", "iwgt"]).drop(columns="subentry")
        wgts.columns = ["wgt"]
        wgts = wgts.join(isyst_map, on="iwgt")
        wgts = wgts[wgts.isyst.isin(needed_syst_indices)]
        if len(wgts):
            yield wgts


def getsyst(f, systematics, nuind, multisim_nuniv=100, slim=False, slimname="slim"):
    if "globalTree" not in f:
        return pd.DataFrame(index=nuind.index)

    nuidx = pd.MultiIndex.from_arrays([nuind.index.get_level_values(0), nuind])
    wgt_names, wgt_types, wgt_nuniv = _read_wgt_metadata(f)
    isyst_map = _build_isyst_map(wgt_nuniv)
    nuniv = wgt_nuniv.sum()

    needed_syst_indices = {wgt_names.index(s) for s in systematics if s in wgt_names}
    if len(needed_syst_indices) != len(systematics):
        missing = [s for s in systematics if s not in wgt_names]
        raise Exception("Missing systematic weights in file: %s" % ", ".join(missing))

    if slim:
        cols = pd.MultiIndex.from_product([[slimname], [f"univ_{i}" for i in range(multisim_nuniv)]])
        systs_slim = pd.DataFrame(1.0, index=nuidx, columns=cols)
    else:
        systs_slim = None

    step_size = CHUNK_SIZE
    syst_chunks = {}

    for wgts in _iter_filtered_wgts(f["recTree"], nuniv, isyst_map, needed_syst_indices, step_size):
        for s in systematics:
            isyst = wgt_names.index(s)
            this_systs = []
            filtered = wgts[wgts.isyst == isyst]
            if filtered.empty:
                continue

            # Get weight type
            # +/- 1,2,3 sigma
            if wgt_types[isyst] == 3 and wgt_nuniv[isyst] == 1:  # morph unisim
                s_morph = filtered.wgt.groupby(level=[0, 1]).first()
                s_morph.name = (s, "morph")

                if slim:
                    for i in range(multisim_nuniv):
                        seed_input = s + str(i)
                        seed_int = int(hashlib.md5(seed_input.encode('utf-8')).hexdigest(), 16) % (2**32)
                        np.random.seed(seed_int)
                        wgt = 1 + (s_morph - 1) * 2 * np.abs(np.random.normal(0, 1))
                        wgt = np.maximum(wgt, 0)
                        systs_slim[(slimname, f"univ_{i}")] = systs_slim[(slimname, f"univ_{i}")].values * wgt
                else:
                    this_systs.append(s_morph)

            elif wgt_types[isyst] == 3 and wgt_nuniv[isyst] > 1:  # +/- sigma unisim
                nwgt = filtered.wgt.groupby(level=[0, 1]).size().values[0]
                nsigma = nwgt // 2

                pivoted = (
                    filtered.wgt
                    .reset_index(level=2)
                    .pivot_table(values="wgt", index=["entry", "inu"], columns="iwgt")
                )

                for isigma in range(nsigma):
                    s_ps = pivoted.iloc[:, 2 * isigma].rename((s, "ps%i" % (isigma + 1)))
                    s_ms = pivoted.iloc[:, 2 * isigma + 1].rename((s, "ms%i" % (isigma + 1)))

                    if slim and isigma == 0:  # use ps1
                        for i in range(multisim_nuniv):
                            seed_input = s + str(i)
                            seed_int = int(hashlib.md5(seed_input.encode('utf-8')).hexdigest(), 16) % (2**32)
                            np.random.seed(seed_int)
                            wgt = 1 + (s_ps - 1) * np.random.normal(0, 1)
                            wgt = np.maximum(wgt, 0)
                            col = (slimname, f"univ_{i}")
                            systs_slim.loc[wgt.index, col] = systs_slim.loc[wgt.index, col] * wgt
                    else:
                        this_systs.append(s_ps)
                        this_systs.append(s_ms)

                # check if we also saved the 0-sigma weight. This is conventionally put last
                if nwgt % 2 != 0:
                    s_cv = pivoted.iloc[:, nwgt - 1].rename((s, "cv"))
                    this_systs.append(s_cv)
                # otherwise, assume it's one
                else:
                    this_systs.append(pd.Series(1, index=pivoted.index, name=(s, "cv")))

            elif wgt_types[isyst] == 0:  # multisim
                this_wgts = filtered.wgt.groupby(level=[0, 1]).head(multisim_nuniv)
                this_wgts = this_wgts.reset_index(level=2)
                this_wgts = this_wgts.pivot_table(values="wgt", index=["entry", "inu"], columns="iwgt")
                this_wgts.columns = pd.MultiIndex.from_tuples(
                    [(s, "univ_%i" % i) for i in range(len(this_wgts.columns))]
                )

                if slim:
                    for i in range(multisim_nuniv):
                        col = (slimname, f"univ_{i}")
                        src = (s, f"univ_{i}")
                        if src in this_wgts.columns:
                            systs_slim.loc[this_wgts.index, col] = (
                                systs_slim.loc[this_wgts.index, col].values * this_wgts[src].values
                            )

                for c in this_wgts.columns:
                    this_systs.append(this_wgts[c])
            else:
                raise Exception("Cannot decode systematic uncertainty: %s" % s)

            for syst in this_systs:
                syst_chunks.setdefault(syst.name, []).append(syst)

    if slim:
        s_idx = systs_slim.index.get_indexer(nuidx)
        systs_slim.loc[s_idx < 0, :] = 1.
        systs_slim.index = nuind.index
        return systs_slim

    if not syst_chunks:
        return pd.DataFrame(index=nuind.index)

def print_syst_all(f):
    if "globalTree" not in f:
        print("no globalTree")

    globalTree = f["globalTree"]
    wgt_names = [n for n in f["globalTree"]['global/wgts/wgts.name'].arrays(library="np")['wgts.name'][0]]
    wgt_types = f["globalTree"]['global/wgts/wgts.type'].arrays(library="np")['wgts.type'][0]
    wgt_nuniv = f["globalTree"]['global/wgts/wgts.nuniv'].arrays(library="np")['wgts.nuniv'][0]
    for i in range(len(wgt_names)):
        print(f"Index: {i:<3} | Name: {wgt_names[i]:<30} | Type: {wgt_types[i]:<5} | Univ: {wgt_nuniv[i]}")

def get_all_syst_df(f, multisim_nuniv=1000):
    ## adding this function to open globaltree only once for each flat.caf
    if "globalTree" not in f:
        print("no globalTree")
        #return pd.DataFrame(index=nuind.index)

    globalTree = f["globalTree"]
    wgt_names = [n for n in globalTree['global/wgts/wgts.name'].arrays(library="np")['wgts.name'][0]]
    wgt_types = globalTree['global/wgts/wgts.type'].arrays(library="np")['wgts.type'][0]
    wgt_nuniv = globalTree['global/wgts/wgts.nuniv'].arrays(library="np")['wgts.nuniv'][0]

    nuniv = wgt_nuniv.sum()

    # 1. Load DataFrame and explicitly extract the first column as a pure Series
    wgts_df = ak.to_dataframe(f["recTree"]['rec.mc.nu.wgt.univ'].arrays(library="ak"), how=None)[0]
    wgts_series = wgts_df.iloc[:, 0].copy()

    # 2. Calculate indices
    idx_entry = wgts_series.index.get_level_values(0)
    idx_subentry = wgts_series.index.get_level_values(1)

    inu = idx_subentry // nuniv
    iwgt = idx_subentry % nuniv

    # 3. Apply MultiIndex and unstack ONCE.
    # Since wgts_series is a 1D Series, unstacking creates a DataFrame with flat integer columns (0, 1, 2...)
    wgts_series.index = pd.MultiIndex.from_arrays(
        [idx_entry, inu, iwgt],
        names=["entry", "inu", "iwgt"]
    )
    df_wgts = wgts_series.unstack("iwgt")

    # 4. Dictionary construction
    systs_dict = {}
    start_col = 0

    for isyst in range(len(wgt_names)):
        s = wgt_names[isyst]
        w_type = wgt_types[isyst]
        w_nuniv = wgt_nuniv[isyst]

        if w_type == 3 and w_nuniv == 1: # morph unisim
            systs_dict[(s, "morph")] = df_wgts[start_col]

        elif w_type == 3 and w_nuniv > 1: # +/- sigma unisim
            nsigma = w_nuniv // 2
            for isigma in range(nsigma):
                systs_dict[(s, f"ps{isigma+1}")] = df_wgts[start_col + 2*isigma]
                systs_dict[(s, f"ms{isigma+1}")] = df_wgts[start_col + 2*isigma + 1]

            # check if we also saved the 0-sigma weight (conventionally put last)
            if w_nuniv % 2 != 0:
                systs_dict[(s, "cv")] = df_wgts[start_col + w_nuniv - 1]
            else:
                systs_dict[(s, "cv")] = pd.Series(1.0, index=df_wgts.index)

        elif w_type == 0: # multisim
            limit = min(w_nuniv, multisim_nuniv)
            for i in range(limit):
                systs_dict[(s, f"univ_{i}")] = df_wgts[start_col + i]

        else:
            raise Exception(f"Cannot decode systematic uncertainty: {s}")

        # Advance the pointer
        start_col += w_nuniv

    # 5. Build final DataFrame
    systs = pd.DataFrame(systs_dict)

    return systs

def filter_systs_nuind(f, systs, nuind):
    ## use outputs of the "get_wgts_df_and_list" above as input for this function
    if "globalTree" not in f:
        return pd.DataFrame(index=nuind.index)

    nuidx = pd.MultiIndex.from_arrays([nuind.index.get_level_values(0), nuind])

    systs_match = systs.reindex(nuidx, fill_value=1.0)
    systs_match.index = nuind.index

    return systs_match
