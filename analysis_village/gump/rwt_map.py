import pandas as pd
import os
import sys
from cycler import cycler
import argparse
from functools import reduce

workspace_root = os.getcwd()
sys.path.insert(0, workspace_root + "/../../")

import pyanalib.pandas_helpers as ph
import warnings
from pyanalib.split_df_helpers import *
from makedf.util import *
from analysis_village.gump.gump_cuts import *
import analysis_village.gump.PID as PID 

def clean_pot(df):
    # 1. Identify rows that are duplicates of a (run, subrun) pair.
    # keep='first' marks the very first occurrence as False (not a duplicate)
    # and all subsequent occurrences as True.
    id_cols = [col for col in df.columns if col[0] in ['run', 'subrun', '__ntuple']]

    is_duplicate = df.duplicated(subset=id_cols, keep='first')
    
    # 2. Set 'pot' to 0.0 for every row that is a duplicate
    df.loc[is_duplicate, 'pot'] = 0.0
    
    print(f"POT: {df['pot'].astype('float64').sum()}")
    return df

class FileHistogramFunction:
    def __init__(self, filename):
        with open(filename, 'r') as f:
            line1 = f.readline().strip('# ').split(',')[:-1]
            line2 = f.readline().strip('# ').split(',')[:-1]

            # Extract x metadata
            self.x_edges = np.array([float(l) for l in line1])
            self.y_edges = np.array([float(l) for l in line2]) 

        # 2. Load the actual data grid (skipping the header lines)
        self.grid = np.loadtxt(filename, delimiter=",")

    def __call__(self, x_arr, y_arr):
        # x_arr and y_arr are now numpy arrays (e.g., df.nu_E_calo.values)
        
        # Use digitize to find bin indices for all points at once
        ix = np.digitize(x_arr, self.x_edges) - 1
        iy = np.digitize(y_arr, self.y_edges) - 1
        
        # Handle out-of-bounds (set to a default or clip)
        mask = (ix >= 0) & (ix < self.grid.shape[0]) & \
               (iy >= 0) & (iy < self.grid.shape[1])
        
        # Pre-fill result with 1.0 (your default)
        result = np.ones_like(x_arr, dtype=float)
        
        # Apply grid values where mask is True
        # We use indexing with arrays here
        result[mask] = self.grid[ix[mask], iy[mask]]
        
        return np.nan_to_num(result, nan=1.0)

def save_histogram(filename, hist_values, x_edges, y_edges):
    # Extract metadata to store in the header
    nx, ny = hist_values.shape
    
    # Create a header string
    header=""
    for x in x_edges:
        header += f"{x},"
    header +="\n"
    for y in y_edges:
        header += f"{y},"

    # Save the 2D grid
    print(f"Saving: {filename}")
    np.savetxt(filename, hist_values, header=header, delimiter=",")

def make_hists_from_files(files, output):
    b_E = np.array([0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.25, 1.5, 3.0])
    b_p = np.array([0.0, 0.2, 0.4, 0.6])

    dataframes = match_across_detvar(files)

    hist2ds = []
    for df in dataframes:
        p = df['pot'].astype('float64').sum()
        df = all_cuts(df, df.detector.iloc[0])
        counts, _, _ = np.histogram2d(df.nu_E_calo, df.del_p, bins=[b_E, b_p])
        hist2ds.append(counts/p)

    hist2ds = np.array(hist2ds)
    hist2ds_rat = hist2ds[1]/hist2ds[0]

    save_histogram(output, hist2ds_rat, b_E, b_p)

def make_hists_from_var(var_dataframes, cv_dataframe, output):
    b_E = np.array([0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.25, 1.5, 3.0])
    b_p = np.array([0.0, 0.2, 0.4, 0.6])

    hist2ds = []
    for df in var_dataframes:
        counts, _, _ = np.histogram2d(df.nu_E_calo, df.del_p, bins=[b_E, b_p])
        hist2ds.append(counts)

    hist2ds = np.array(hist2ds)

    upper = np.max(hist2ds, axis=0)
    lower = np.min(hist2ds, axis=0)

    cv_hist2d, _, _  = np.histogram2d(cv_dataframe.nu_E_calo, cv_dataframe.del_p, bins=[b_E, b_p])

    upper_hist2ds_rat = upper/cv_hist2d
    save_histogram('rwt_outputs/pls_'+output, upper_hist2ds_rat, b_E, b_p)

    lower_hist2ds_rat = lower/cv_hist2d
    save_histogram('rwt_outputs/min_'+output, lower_hist2ds_rat, b_E, b_p)

def match_reco_to_mc(recomatchdf, mcdf):
    recomatchdf = recomatchdf[['run', 'subrun', 'evt', 'tmatch_idx', 'pot']].drop_duplicates() # get rid of most cols, to avoid *_x, *_y when merging later
                                                                                        # also drop duplicates, because we just watch to attach things to mc
                                                                                        # i.e. if run 1183, sr 85, evt 1 has tmatch_idx 0.0 for multiple slices
                                                                                        # don't care how many slices, just care about matching hdr info the mc truth
    
    recomatchdf = recomatchdf[recomatchdf.tmatch_idx != -999] # get rid of everything that doesn't match to mc
    
    recomatchdf.columns = pd.MultiIndex.from_tuples([(col, '') for col in recomatchdf.columns])

    mcmatchdf = ph.multicol_merge(recomatchdf.reset_index(), mcdf.reset_index(),
                           left_on=[("__ntuple", ""), ("entry", ""), ("tmatch_idx", "")],
                           right_on=[("__ntuple", ""), ("entry", ""), ("rec.mc.nu..index", "")],
                           how="right").set_index(['__ntuple', 'entry', "rec.mc.nu..index"])
    return mcmatchdf

def match_across_detvar(files):

    keep_cols = ['slc_vtx_x', 'slc_vtx_y', 'slc_vtx_z', 'nu_E_calo', 'detector', 'tmatch_idx', 'del_p', 'nu_score',
                 'mu_chi2_of_mu_cand', 'mu_chi2_of_prot_cand', 'prot_chi2_of_mu_cand', 'prot_chi2_of_prot_cand', 'mu_len',
                 'other_trk_length', 'other_shw_length', 'mu_end_x', 'mu_end_y', 'mu_end_z', 'p_end_x', 'p_end_y', 'p_end_z', 'crthit', 'nu_E_true']
    dataframes = []
    unfilt_mcdataframes = []
    unfilt_recodataframes = []
    
    for f in files:
        dfs = load_dfs(f, ['evt', 'mcnu', 'hdr']) 
        mcdf = dfs['mcnu']
        recodf = dfs['evt']
        recodf = recodf[keep_cols]
        hdrdf = dfs['hdr']

        if 'del_p' in mcdf.columns:
            mcdf.rename(columns={'del_p':'del_p_true'}, inplace=True)

        DETECTOR = recodf.detector.iloc[0]
        
        # Keeps all hdrdf rows, can't be slices with out header
        merged = recodf.reset_index().merge(
            hdrdf.reset_index(),
            on=['__ntuple', 'entry'],
            how='right'  
        )

        # when we have a header without a slice, give it an index to preserve header info
        merged['rec.slc..index'] = merged['rec.slc..index'].fillna(0) 

        recomatchdf = merged.set_index(['__ntuple', 'entry']) # use this for matching header to truth info
        recodf = merged.set_index(['__ntuple', 'entry', 'rec.slc..index']) # use this for matching common truth to reco

        unfilt_mcdataframes.append(match_reco_to_mc(recomatchdf, mcdf))
        unfilt_recodataframes.append(recodf)

    # now match up all of our mc dfs with hdr info added
    filt_mcdataframes = filter_n_common_events(unfilt_mcdataframes, keys=["run", "subrun", "evt", "nu_E"])    

    dataframes = [apply_mc_filt_to_reco(filt_mc, reco) for filt_mc, reco in zip(filt_mcdataframes, unfilt_recodataframes)]
    return dataframes

def apply_mc_filt_to_reco(filt_mc, reco):

    valid_df = filt_mc[['run', 'subrun', 'pot']].dropna().drop_duplicates()
    #valid_df.columns = pd.MultiIndex.from_tuples([(col, '') for col in valid_df.columns])

    new_filt_mc = filt_mc.drop(columns=[("run", ""), ("subrun", ""), ("pot", ""), ("evt", ""), ("nu_E", "")])
    df = ph.multicol_merge(reco.reset_index(), new_filt_mc.reset_index(),
                           left_on=[("__ntuple", ""), ("entry", ""), ("tmatch_idx", "")],
                           right_on=[("__ntuple", ""), ("entry", ""), ("rec.mc.nu..index", "")],
                           how="left") # start with keeping everything...

    df = df.set_index(['__ntuple', 'entry', 'rec.slc..index'])
    
    # drop cases where run, subrun not in in matched filter.
    df = df.reset_index().merge(valid_df, on=['run', 'subrun', 'pot'], how='inner').set_index(df.index.names)       
    #df = df.reset_index().merge(valid_df, on=[('run', ''), ('subrun', ''), ('pot', '')], how='inner').set_index(df.index.names)       

    # Identify slices with truth
    valid_tmatch = df['tmatch_idx'] != -999
    
    # Tag events that have at least one valid truth match
    event_has_truth = valid_tmatch.groupby(['__ntuple', 'entry']).transform('any')

    # Tag events where EVERY slice has a NaN PDG
    event_all_pdg_nan = df['pdg'].isna().groupby(['__ntuple', 'entry']).transform('all')

    # Create a mask for rows belonging to events that meet both criteria
    to_drop = event_has_truth & event_all_pdg_nan

    # copy pot everywhere so we don't lose it to mask
    df['pot'] = df.groupby([pd.Grouper(level='__ntuple'), 'run', 'subrun'])['pot'].transform('max')
    
    # apply mask 
    df_filtered = df[~to_drop]
   
    # clean pot after mask 
    df = clean_pot(df_filtered.reset_index())

    return df.set_index(['__ntuple', 'entry', 'rec.slc..index'])

def apply_double_map(df, min_map_file, pls_map_file, col_name):
    min_func = FileHistogramFunction(min_map_file)
    pls_func = FileHistogramFunction(pls_map_file)
    data = {
            (col_name, 'ps1') : [pls_func(E, del_p) for E, del_p in zip(df.nu_E_calo, df.del_p)],
            (col_name, 'ms1') : [min_func(E, del_p) for E, del_p in zip(df.nu_E_calo, df.del_p)]
            }

    new_df = pd.DataFrame(data)
    return new_df

def apply_map(df, map_file, col_name):
    func = FileHistogramFunction(map_file)
    weights = func(df.nu_E_calo.values, df.del_p.values)
    return pd.DataFrame({(col_name, 'ps1'): weights}, index=df.index)

def filter_n_common_events(dfs, keys=['run', 'subrun', 'evt', 'nu_E_true']):
    if not dfs:
        return []

    # we use an additional col to count the number of duplicate columns
    # matching with this also helps catch instances of multiple matching truth slices
    # in same event.
    
    processed_dfs = []
    for df in dfs:
        temp_df = df.copy()
        temp_df['occ_id'] = temp_df.groupby(keys, dropna=False).cumcount()
        processed_dfs.append(temp_df)
    
    matching_keys = keys + ['occ_id']

    id_sets = (df[matching_keys].drop_duplicates() for df in processed_dfs)
    common_ids = reduce(lambda left, right: pd.merge(left, right, on=matching_keys), id_sets)
    
    filtered_dfs = []

    for df in processed_dfs:
        index_names = df.index.names
        
        f_df = df.reset_index() \
                 .merge(common_ids, on=matching_keys) \
                 .set_index(index_names) \
                 .drop(columns=['occ_id']) # Clean up the helper column
        
        f_df.sort_values(keys, inplace=True)
        
        filtered_dfs.append(f_df)

    return filtered_dfs

def make_plots(dataframes, norms=[], file_titles=["0xSCE", "CV", "2xSCE"], global_title="GUMP Selection SCE Detector Variations", output_tag='SCE'):
    if len(norms) == 0:
        norms = np.ones(len(dataframes))

    b_E = np.array([0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.25, 1.5, 3.0])
    b_p = np.array([0.0, 0.2, 0.4, 0.6])

    bins = [b_E, b_p]
    reco_vars = ['nu_E_calo', 'del_p']

    # 1D dists

    colors = ['red', 'blue', 'green']
    for var, b in zip(reco_vars, bins):
        plt.figure(figsize=(10, 6))
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        for df, t, norm, c in zip(dataframes, file_titles, norms, colors):
            n, _, _ = plt.hist(df[var], bins=b, label=t, weights=(1/norm)*np.ones(len(df.nu_E_calo)), histtype='step', linewidth=1, color=c)
            nerr, _ = np.histogram(df[var], bins=b, weights=(1/(norm)**2)*np.ones(len(df.nu_E_calo)))
            mcstaterr = np.sqrt(nerr)
            plt.bar(b[:-1], height=2*mcstaterr, bottom=n-mcstaterr, width=(b[1:]-b[:-1]), align='edge', color=c, alpha=0.5)[0]

        plt.xlabel(r'Reconstructed Energy $E_{calo}$ [GeV]')
        plt.ylabel('Events (Normalized)')
        plt.title(global_title, fontsize=20)
        plt.legend()
        plt.savefig('rwt_outputs/'+var+'_'+output_tag+'.png', dpi=300)
        plt.clf()

    # 2D dist
    max_counts = []
    for df, norm in zip(dataframes, norms):
        counts, _, _ = np.histogram2d(df.nu_E_calo, df.del_p, weights=(1/norm)*np.ones(len(df.nu_E_calo)), bins=[b_E, b_p])
        max_counts.append(counts.max())
    global_vmax = max(max_counts)

    fig, axes = plt.subplots(nrows=1, ncols=len(dataframes), figsize=(18, 5), constrained_layout=True)
    hist2ds = []
    for df, a, t, norm in zip(dataframes, axes, file_titles, norms):
        im = a.hist2d(df.nu_E_calo, df.del_p, bins=[b_E, b_p], vmin=0.0, vmax=global_vmax, weights = (1/norm)*np.ones((len(df.nu_E_calo))))
        hist2ds.append(np.transpose(im[0]))

        a.set_xlabel(r'Reconstructed Energy $E_{calo}$ [GeV]')
        a.set_ylabel(r'$\delta p$ [GeV/c]')
        a.set_title(t)

    hist2ds = np.array(hist2ds)

    cbar = fig.colorbar(im[3], ax=axes.ravel().tolist(), label='Events (Normalized)', format='%.0e')
    fig.suptitle(global_title, fontsize=40)

    fig.savefig('rwt_outputs/2d'+output_tag+'.png', dpi=300)
    plt.close(fig) 

def plot_2d_hist_from_file(filename, plot_title, output_tag):
    x_edges = []
    y_edges = []
    data_rows = []

    with open(filename, 'r') as f:
        lines = f.readlines()
        x_edges = [float(x) for x in lines[0].strip('# ').split(',') if x.strip()]
        y_edges = [float(y) for y in lines[1].strip('# ').split(',') if y.strip()]
        
        for line in lines[2:]:
            if line.strip():
                row = [float(val) for val in line.strip().split(',') if val.strip()]
                data_rows.append(row)

    z_values = np.array(data_rows)

    plt.figure(figsize=(10, 6))
    X, Y = np.meshgrid(x_edges, y_edges)
    mesh = plt.pcolormesh(x_edges, y_edges, z_values.T, cmap='seismic', linewidth=0.1, vmin=0.8, vmax=1.2)
    
    plt.colorbar(mesh, label='Value')
    plt.title(plot_title)
    plt.xlabel(r'Reconstructed Energy $E_{calo}$ [GeV]')
    plt.ylabel(r'$\delta p$ [GeV/c]')
    
    mesh.get_cmap().set_bad(color='gray')
    plt.savefig('rwt_outputs/2d_ratio_'+output_tag+'.png', dpi=300)
    plt.clf() 

def main():
    """Main analysis pipeline."""

    make_hists_from_files(["/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/SBND_SpringMC_Nom.df", "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/SBND_SpringMC_WMYZ.df"], "rwt_outputs/YZ.txt")
    make_hists_from_files(["/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/SBND_SpringMC_Nom.df", "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/SBND_SpringMC_WMXThetaXW.df"], "rwt_outputs/XThetaXW.txt")

    make_hists_from_files(["/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/SBND_SpringMC_Nom.df", "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/SBND_SpringMC_0xSCE.df"], "rwt_outputs/min_SCE.txt")
    make_hists_from_files(["/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/SBND_SpringMC_Nom.df", "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/SBND_SpringMC_2xSCE.df"], "rwt_outputs/pls_SCE.txt")

    dfs = []
    for i in range(20):
        dfs.append(load_dfs(f"/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/no-flashes/SBND_SpringMC_rewgt_{i}.df", ['evt'])['evt'])
    recodf = pd.concat(dfs)
    recodf_precut = recodf.copy()
    recodf = all_cuts(recodf, recodf.detector.iloc[0])
    make_hists_from_var(PID.get_smear_vars(recodf_precut), recodf.copy(), 'SBND_smear.txt')
    make_hists_from_var(PID.get_gain_vars(recodf_precut), recodf.copy(), 'SBND_gain.txt')

    recodf = load_dfs(f"/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/ICARUS_SpringMCOverlay_rewgt.df", ['evt'])['evt']
    recodf_precut = recodf.copy()
    recodf = all_cuts(recodf, "ICARUS")
    make_hists_from_var(PID.get_smear_vars(recodf_precut), recodf.copy(), 'ICARUS_smear.txt')
    make_hists_from_var(PID.get_gain_vars(recodf_precut), recodf.copy(), 'ICARUS_gain.txt')

    #just a random dataframe to test map with
    apply_double_map(recodf.copy(), 'rwt_outputs/min_SCE.txt', 'rwt_outputs/pls_SCE.txt', 'CAFPYANA_SBN_v1_multisigma_SCE')
    apply_map(recodf, 'rwt_outputs/YZ.txt', 'CAFPYANA_SBN_v1_multisigma_WMYZ')
    apply_map(recodf, 'rwt_outputs/XThetaXW.txt', 'CAFPYANA_SBN_v1_multisigma_WMXThetaXW')


def plot():
    # SCE
    prefix = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/"
    dataframes = match_across_detvar([prefix+"SBND_SpringMC_0xSCE.df", prefix+"SBND_SpringMC_Nom.df", prefix+"SBND_SpringMC_2xSCE.df"])
    pots = []

    for i in range(len(dataframes)):
        pots.append(dataframes[i]['pot'].astype('float64').sum())
        dataframes[i] = all_cuts(dataframes[i], dataframes[i].detector.iloc[0])

    make_plots(dataframes, norms=pots, file_titles=["0xSCE", "CV", "2xSCE"], global_title="SCE Detector Variations", output_tag='SCE')

    plot_2d_hist_from_file("rwt_outputs/min_SCE.txt", "0xSCE Reweight", 'min_SCE')
    plot_2d_hist_from_file("rwt_outputs/pls_SCE.txt", "2xSCE Reweight", 'pls_SCE')

    # WireMod
    prefix = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/"
    dataframes = match_across_detvar([prefix+"SBND_SpringMC_WMXThetaXW.df", prefix+"SBND_SpringMC_WMYZ.df", prefix+"SBND_SpringMC_Nom.df",])
    pots = []
    for i in range(len(dataframes)):
        pots.append(dataframes[i]['pot'].astype('float64').sum())
        dataframes[i] = all_cuts(dataframes[i], dataframes[i].detector.iloc[0])

    make_plots(dataframes, norms=pots, file_titles=['XThetaXW', 'YZ', 'Nominal'], global_title="WireMod Detector Variations", output_tag="WM")

    plot_2d_hist_from_file("rwt_outputs/YZ.txt", "WireMod YZ Variation", 'WMYZ')
    plot_2d_hist_from_file("rwt_outputs/XThetaXW.txt", "WireMod XThetaXW", 'WMXThetaXW')

    # SBND, smear 
    plot_2d_hist_from_file("rwt_outputs/min_SBND_smear.txt", "SBND Smear Minimum", 'min_SBND_smear')
    plot_2d_hist_from_file("rwt_outputs/pls_SBND_smear.txt", "SBND Smear Maximum", 'pls_SBND_smear')

    # SBND, gain
    plot_2d_hist_from_file("rwt_outputs/min_SBND_gain.txt", "SBND Gain Minimum", 'min_SBND_gain')
    plot_2d_hist_from_file("rwt_outputs/pls_SBND_gain.txt", "SBND Gain Maximum", 'pls_SBND_gain')

    # ICARUS, smear
    plot_2d_hist_from_file("rwt_outputs/min_ICARUS_smear.txt", "ICARUS Smear Minimum", 'min_ICARUS_smear')
    plot_2d_hist_from_file("rwt_outputs/pls_ICARUS_smear.txt", "ICARUS Smear Maximum", 'pls_ICARUS_smear')

    # ICARUS, gain
    plot_2d_hist_from_file("rwt_outputs/min_ICARUS_gain.txt", "ICARUS Gain Minimum", 'min_ICARUS_gain')
    plot_2d_hist_from_file("rwt_outputs/pls_ICARUS_gain.txt", "ICARUS Gain Maximum", 'pls_ICARUS_gain')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog = 'rwt_map',
                                     description = 'Reweight map generator script.')
    parser.add_argument('-m','--main', action='store_true')
    parser.add_argument('-p','--plot', action='store_true')

    args = parser.parse_args()

    if args.main:
        main()
    if args.plot:
        plot()
