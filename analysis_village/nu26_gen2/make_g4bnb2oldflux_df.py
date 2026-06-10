from makedf.makedf import *
from pyanalib.pandas_helpers import *
import uproot

INPUT_G4BNB = "/exp/sbnd/app/users/sungbino/framework_dev/g4bnb2old/nuMom_immParent_merged_200.root"
INPUT_OLD = "/exp/sbnd/app/users/sungbino/framework_dev/g4bnb2old/nuMom_immParent_oldFiles_200_merged.root"

nu_parent_pids = [211, -211, 321, -321, 130, 13, -13]
nu_types = { 211: ["numu", "nue"],
             -211: ["numubar", "nuebar"],
             321: ["numu", "nue"],
             -321: ["numubar", "nuebar"],
             130: ["numu", "nue", "numubar", "nuebar"],
             13: ["numu", "nuebar"],
             -13: ["numubar", "nue"]
            }
nu_pids = {"numu": 14, "numubar": -14, "nue": 12, "nuebar": -12}

px_edges = None
py_edges = None
pz_edges = None

def InAV(data):
    xmin = -200.
    xmax = 200.
    ymin = -200.
    ymax = 200.
    zmin = 0.
    zmax = 500.
    pass_xyz = (np.abs(data.x) > xmin) & (np.abs(data.x) < xmax) & (data.y > ymin) & (data.y < ymax) & (data.z > zmin) & (data.z < zmax)
    return pass_xyz

def call_flux_hists(f_root):
    maps = []
    px_edges = py_edges = pz_edges = None

    with uproot.open(f_root) as f_in:
        for nu_par in nu_parent_pids:
            for nu_type in nu_types[nu_par]:
                nu_pid = nu_pids[nu_type]
                hist_prefix = f"hSBND_{nu_par}_{nu_type}"

                # Fetch histograms
                h_3mom = f_in[hist_prefix + "_3mom"]

                # Extract bin edges once
                if px_edges is None:
                    px_edges = h_3mom.axis(0).edges()
                    py_edges = h_3mom.axis(1).edges()
                    pz_edges = h_3mom.axis(2).edges()

                # Get the 3D numpy array
                vals = h_3mom.values()

                # OPTIMIZATION 1: Find only the indices where the value is NOT zero
                xbins, ybins, zbins = np.nonzero(vals)

                # Extract only the non-zero scales
                non_zero_scales = vals[xbins, ybins, zbins]

                # OPTIMIZATION 2: Create DataFrame directly with downcasted datatypes
                vals_df = pd.DataFrame({
                    "px_index": xbins.astype(np.uint16),
                    "py_index": ybins.astype(np.uint16),
                    "pz_index": zbins.astype(np.uint16),
                    "scale": non_zero_scales.astype(np.float32),
                    "parent_pdg": np.int16(nu_par),
                    "pdg": np.int16(nu_pid)
                })

                maps.append(vals_df)

    # Concatenate the lists
    out_df = pd.concat(maps, ignore_index=True)

    #print(out_df)
    #print(out_df.info())
    return out_df, px_edges, py_edges, pz_edges

g4bnb_df, g4bnb_x_bins, g4bnb_y_bins, g4bnb_z_bins = call_flux_hists(INPUT_G4BNB)
oldflux_df, oldflux_x_bins, oldflux_y_bins, oldflux_z_bins = call_flux_hists(INPUT_OLD)

g4bnb_df.columns = pd.MultiIndex.from_tuples([(col, '') for col in g4bnb_df.columns])
oldflux_df.columns = pd.MultiIndex.from_tuples([(col, '') for col in oldflux_df.columns])

#print(g4bnb_dfs)

def make_g4bnb2oldflux_df(f):
    this_mcdf = loadbranches(f["recTree"], mcbranches).rec.mc.nu
    is_av = InAV(this_mcdf.position)
    this_mcdf = this_mcdf[is_av]

    selected_cols = [
        ('E', ''),
        ('momentum', 'x'),
        ('momentum', 'y'),
        ('momentum', 'z'),
        ('pdg', ''),
        ('parent_pdg', ''),
        ('genie_mode', ''),
    ]
    this_mcdf = this_mcdf[selected_cols]

    # Cap px and py between -5 and +5, and pz beween 0 and 5: template histograms are only available in this region
    this_mcdf[('momentum', 'x')] = this_mcdf[('momentum', 'x')].clip(lower=-5.0, upper=5.0)
    this_mcdf[('momentum', 'y')] = this_mcdf[('momentum', 'y')].clip(lower=-5.0, upper=5.0)
    this_mcdf[('momentum', 'z')] = this_mcdf[('momentum', 'z')].clip(lower=0.0, upper=5.0)

    this_mcdf['px_index'] = pd.cut(this_mcdf[('momentum', 'x')], bins=g4bnb_x_bins, labels=False)
    this_mcdf['py_index'] = pd.cut(this_mcdf[('momentum', 'y')], bins=g4bnb_y_bins, labels=False)
    this_mcdf['pz_index'] = pd.cut(this_mcdf[('momentum', 'z')], bins=g4bnb_z_bins, labels=False)

    merge_cols = [('px_index', ''), ('py_index', ''), ('pz_index', ''), ('parent_pdg', ''), ('pdg', '')]

    # add g4bnb entry
    this_mcdf = (this_mcdf.reset_index().merge(g4bnb_df, on=merge_cols, how='left').set_index(['entry', 'rec.mc.nu..index']))
    this_mcdf = this_mcdf.rename(columns={'scale': 'g4bnb_entry'}, level=0)
    this_mcdf.g4bnb_entry = this_mcdf.g4bnb_entry.fillna(-1)

    #print(this_mcdf)

    # add oldflux entry
    this_mcdf = (this_mcdf.reset_index().merge(oldflux_df, on=merge_cols, how='left').set_index(['entry', 'rec.mc.nu..index']))
    this_mcdf = this_mcdf.rename(columns={'scale': 'oldflux_entry'}, level=0)
    this_mcdf.oldflux_entry = this_mcdf.oldflux_entry.fillna(-1)


    # Calculate the individual 1D weights
    this_mcdf[('old_over_g4bnb', '')] = this_mcdf[('oldflux_entry', '')] / this_mcdf[('g4bnb_entry', '')]

    this_mcdf[('event_weight_old_over_g4bnb_per_entry', '')] = (
        this_mcdf[('old_over_g4bnb', '')].groupby(level=0).transform('prod')
    )

    hdrdf = make_hdrdf(f)

    # 1. Grab only the columns you want from the header dataframe
    header_info = hdrdf[['run', 'subrun', 'evt']].copy()

    # 2. Convert its flat columns into MultiIndex tuples to match this_mcdf
    # 'run' -> ('run', ''), 'subrun' -> ('subrun', ''), etc.
    header_info.columns = pd.MultiIndex.from_tuples([(c, '') for c in header_info.columns])

    # 3. Join them together!
    # Because both dataframes share 'entry' in their row index,
    # Pandas will automatically duplicate the run/subrun/evt for neutrinos in the same entry.
    this_mcdf = this_mcdf.join(header_info, how='left')

    #print(this_mcdf)

    return this_mcdf
