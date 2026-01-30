import uproot
import numpy as np
import pandas as pd
import awkward as ak

def getsyst(f, systematics, nuind, multisim_nuniv=100, slim=False, slimname="slim"):
    if "globalTree" not in f:
        return pd.DataFrame(index=nuind.index)

    nuidx = pd.MultiIndex.from_arrays([nuind.index.get_level_values(0), nuind])

    if slim:
        # one column to save them all
        cols = pd.MultiIndex.from_product(
            [[slimname], [f"univ_{i}" for i in range(multisim_nuniv)]],
        )
        systs_slim = pd.DataFrame(
            1.0,
            index=nuidx,
            columns=cols,
        )

    globalTree = f["globalTree"]
    wgt_names = [n for n in f["globalTree"]['global/wgts/wgts.name'].arrays(library="np")['wgts.name'][0]]
    wgt_types = f["globalTree"]['global/wgts/wgts.type'].arrays(library="np")['wgts.type'][0]
    wgt_nuniv = f["globalTree"]['global/wgts/wgts.nuniv'].arrays(library="np")['wgts.nuniv'][0]

    isyst = pd.Series(np.repeat(list(range(len(wgt_nuniv))), wgt_nuniv), name="isyst")
    isyst.index.name = "iwgt"
    nuniv = wgt_nuniv.sum()

    wgts = ak.to_dataframe(f["recTree"]['rec.mc.nu.wgt.univ'].arrays(library="ak"), how=None)[0]
    wgts["inu"] = wgts.index.get_level_values(1) // nuniv
    wgts["iwgt"] = wgts.index.get_level_values(1) % nuniv
    wgts = wgts.reset_index().set_index(["entry", "inu", "iwgt"]).drop(columns="subentry")
    wgts.columns = ["wgt"]
    wgts = wgts.join(isyst)

    systs = []
    for s in systematics:
        isyst = wgt_names.index(s)
        this_systs = []

        # Get weight type
        # +/- 1,2,3 sigma
        if wgt_types[isyst] == 3 and wgt_nuniv[isyst] == 1: # morph unisim
            s_morph = wgts[wgts.isyst == isyst].wgt.groupby(level=[0,1]).first()
            s_morph.name = (s, "morph")

            if slim:
                for i in range(multisim_nuniv):
                    np.random.seed(hash(s+str(i)) % (2**32))
                    wgt = 1 + (s_morph - 1) * 2 * np.abs(np.random.normal(0, 1)) # std -> unc.
                    systs_slim[(slimname, f"univ_{i}")] = systs_slim[(slimname, f"univ_{i}")].values * wgt

            else:
                this_systs.append(s_morph)

        elif wgt_types[isyst] == 3 and wgt_nuniv[isyst] > 1: # +/- sigma unisim
            nwgt = wgts[wgts.isyst == isyst].wgt.groupby(level=[0,1]).size().values[0]
            nsigma = nwgt // 2
            for isigma in range(nsigma):
                s_ps = wgts[wgts.isyst == isyst].wgt.groupby(level=[0,1]).nth(2*isigma)
                s_ps.name = (s, "ps%i" % (isigma+1))
                s_ms = wgts[wgts.isyst == isyst].wgt.groupby(level=[0,1]).nth(2*isigma+1)
                s_ms.name = (s, "ms%i" % (isigma+1))

                if slim and isigma == 0: # use ps1
                    for i in range(multisim_nuniv):
                        np.random.seed(hash(s+str(i)) % (2**32))
                        wgt = 1 + (s_ps - 1) * np.random.normal(0, 1)
                        wgt = wgt.reset_index(level=2, drop=True)  # Drop the 'iwgt' level to match systs_slim index
                        systs_slim[(slimname, f"univ_{i}")] = systs_slim[(slimname, f"univ_{i}")].values * wgt
    
                else:
                    this_systs.append(s_ps.droplevel(2))
                    this_systs.append(s_ms.droplevel(2))

            # check if we also saved the 0-sigma weight. This is conventionally put last
            if nwgt % 2 != 0:
                s_cv = wgts[wgts.isyst == isyst].wgt.groupby(level=[0,1]).nth(nwgt-1)
                s_cv.name = (s, "cv")
                this_systs.append(s_cv.droplevel(2))
            # otherwise, assume it's one
            else:
                this_systs.append(pd.Series(1, index=this_systs[-1].index, name=(s, "cv")))

        elif wgt_types[isyst] == 0: # multisim
            this_wgts = wgts[wgts.isyst == isyst].wgt.groupby(level=[0,1]).head(multisim_nuniv) # limit to 250 universes
            this_wgts = this_wgts.reset_index(level=2)
            this_wgts = this_wgts.pivot_table(values="wgt", index=["entry", "inu"], columns="iwgt")
            this_wgts.columns = pd.MultiIndex.from_tuples([(s, "univ_%i"% i) for i in range(len(this_wgts.columns))])

            if slim:
                for i in range(multisim_nuniv):
                    systs_slim[(slimname, f"univ_{i}")] = systs_slim[(slimname, f"univ_{i}")].values * this_wgts[(s, f"univ_{i}")]

            for c in this_wgts.columns:
                this_systs.append(this_wgts[c])

        else:
            raise Exception("Cannot decode systematic uncertainty: %s" % s)

        for syst in this_systs:
            systs.append(syst)

    if slim:
        s_idx = systs_slim.index.get_indexer(nuidx)
        systs_slim.loc[s_idx < 0, :] = 1.
        systs_slim.index = nuind.index
        return systs_slim

    else:
        systs = pd.DataFrame(systs).T
        s_idx = systs.index.get_indexer(nuidx)
        systs_match = systs.iloc[s_idx]
        systs_match.loc[s_idx < 0, :] = 1.
        systs_match.index = nuind.index
        return systs_match

