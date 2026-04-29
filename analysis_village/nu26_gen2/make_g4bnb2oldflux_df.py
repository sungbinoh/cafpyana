from makedf.makedf import *
from pyanalib.pandas_helpers import *
import uproot

INPUT_G4BNB = "/exp/sbnd/app/users/sungbino/framework_dev/g4bnb2old/nuMom_immParent_G4BNB.root"
INPUT_OLD = "/exp/sbnd/app/users/sungbino/framework_dev/g4bnb2old/nuMom_immParent_oldFiles.root"

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

def call_flux_hists(f_root):
    maps_px, maps_py, maps_pz = [], [], []
    px_edges = py_edges = pz_edges = None

    with uproot.open(f_root) as f_in:
        for nu_par in nu_parent_pids:
            for nu_type in nu_types[nu_par]:
                nu_pid = nu_pids[nu_type]
                hist_prefix = f"hSBND_{nu_par}_{nu_type}"
                
                # Fetch histograms
                h_px = f_in[hist_prefix + "_px"]
                h_py = f_in[hist_prefix + "_py"]
                h_pz = f_in[hist_prefix + "_pz"]

                # Extract bin edges once
                if px_edges is None:
                    px_edges = h_px.axis(0).edges()
                    py_edges = h_py.axis(0).edges()
                    pz_edges = h_pz.axis(0).edges()

                # Process each component into a temporary list for concatenation
                for h, container in zip([h_px, h_py, h_pz], [maps_px, maps_py, maps_pz]):
                    vals = h.values()
                    
                    # Create the DataFrame for this specific parent/neutrino combination
                    df = pd.DataFrame({
                        "bin_idx": range(len(vals)),
                        "entry": vals.ravel(),
                        "nu_parent_pid": nu_par,
                        "nu_pid": nu_pid,
                    })
                    
                    # Reorder to match your preferred style
                    new_order = ["nu_parent_pid", "nu_pid", "bin_idx", "entry"]
                    container.append(df[new_order])

    # Combine into final DataFrames
    out_px = pd.concat(maps_px, ignore_index=True)
    out_py = pd.concat(maps_py, ignore_index=True)
    out_pz = pd.concat(maps_pz, ignore_index=True)

    return {"px": out_px, "py": out_py, "pz": out_pz}, px_edges, py_edges, pz_edges
                
g4bnb_dfs, g4bnb_x_bins, g4bnb_y_bins, g4bnb_z_bins = call_flux_hists(INPUT_G4BNB)
oldflux_dfs, oldflux_x_bins, oldflux_y_bins, oldflux_z_bins = call_flux_hists(INPUT_OLD)

#print(g4bnb_dfs)
                  
def make_g4bnb2oldflux_df(f):
    this_mcdf = loadbranches(f["recTree"], mcbranches).rec.mc.nu

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

    for comp in ['px', 'py', 'pz']:
        idx_col_tuple = (f'{comp}_index', '')
        
        # ---------------------------------------------------------
        # 1. Prepare G4BNB dataframe for joining
        # ---------------------------------------------------------
        g4_merge_df = g4bnb_dfs[comp].rename(columns={
            'nu_parent_pid': 'parent_pdg',
            'nu_pid': 'pdg',
            'bin_idx': f'{comp}_index'
        })
        
        # Set the columns we want to match on as the INDEX of the flux dataframe
        g4_merge_df = g4_merge_df.set_index(['parent_pdg', 'pdg', f'{comp}_index'])
        
        # Make the actual value column a MultiIndex so it stacks nicely
        g4_merge_df.columns = pd.MultiIndex.from_tuples([('g4bnb_entry', comp)])
        
        # .join() automatically preserves your left row index!
        this_mcdf = this_mcdf.join(
            g4_merge_df, 
            on=[('parent_pdg', ''), ('pdg', ''), idx_col_tuple]
        )
        
        # ---------------------------------------------------------
        # 2. Prepare OLD flux dataframe for joining
        # ---------------------------------------------------------
        old_merge_df = oldflux_dfs[comp].rename(columns={
            'nu_parent_pid': 'parent_pdg',
            'nu_pid': 'pdg',
            'bin_idx': f'{comp}_index'
        })
        
        # Set the index to use as a lookup table
        old_merge_df = old_merge_df.set_index(['parent_pdg', 'pdg', f'{comp}_index'])
        
        # Make the column a MultiIndex
        old_merge_df.columns = pd.MultiIndex.from_tuples([('oldflux_entry', comp)])
        
        # Join!
        this_mcdf = this_mcdf.join(
            old_merge_df, 
            on=[('parent_pdg', ''), ('pdg', ''), idx_col_tuple]
        )

    # Calculate the individual 1D weights
    this_mcdf[('px_old_over_g4bnb', '')] = this_mcdf[('oldflux_entry', 'px')] / this_mcdf[('g4bnb_entry', 'px')]
    this_mcdf[('py_old_over_g4bnb', '')] = this_mcdf[('oldflux_entry', 'py')] / this_mcdf[('g4bnb_entry', 'py')]
    this_mcdf[('pz_old_over_g4bnb', '')] = this_mcdf[('oldflux_entry', 'pz')] / this_mcdf[('g4bnb_entry', 'pz')]

    # Calculate the total weight
    this_mcdf[('total_old_over_g4bnb', '')] = (
        this_mcdf[('px_old_over_g4bnb', '')] * this_mcdf[('py_old_over_g4bnb', '')] * this_mcdf[('pz_old_over_g4bnb', '')]
    )

    this_mcdf[('event_weight_old_over_g4bnb_per_entry', '')] = (
        this_mcdf[('total_old_over_g4bnb', '')].groupby(level=0).transform('prod')
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
