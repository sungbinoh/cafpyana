import pandas as pd
import os
import sys
from cycler import cycler
import argparse

workspace_root = os.getcwd()
sys.path.insert(0, workspace_root + "/../../")

import pyanalib.pandas_helpers as ph
import warnings
from pyanalib.split_df_helpers import *
from makedf.util import *
from analysis_village.gump.gump_cuts import *
import analysis_village.gump.PID as PID 

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

    def __call__(self, x, y):
        # Check boundaries
        if not (min(self.x_edges) <= x <= max(self.x_edges) and min(self.y_edges) <= y <= max(self.y_edges)):
            return 1.0
             
        # Convert coordinate to index
        centered_x = self.x_edges - x
        centered_y = self.y_edges - y
        ix = len(centered_x[centered_x < 0]) - 1
        iy = len(centered_y[centered_y < 0]) - 1

        # Return value (Transpose if necessary depending on your save orientation)
        return np.nan_to_num(self.grid[ix, iy], nan=1.0)

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

    dataframes, pot = match_across_detvar(files)

    hist2ds = []
    for df, p in zip(dataframes, pot):
        #df = all_cuts(df, df.detector.iloc[0])
        counts, _, _ = np.histogram2d(df.nu_E_calo, df.del_p, bins=[b_E, b_p])
        hist2ds.append(counts*(1e20/p))

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

def match_across_detvar(files):
    dataframes = []
    unfilt_dataframes = []
    pot = []
    keep_cols = ['slc_vtx_x', 'slc_vtx_y', 'slc_vtx_z', 'nu_E_calo', 'detector', 'tmatch_idx', 'del_p']
    for f in files:
        dfs = load_dfs(f, ['evt', 'mcnu', 'hdr'])
        mcdf = dfs['mcnu']
        recodf = dfs['evt']
        recodf = recodf[keep_cols]
        hdrdf = dfs['hdr']
        pot.append(sum(hdrdf.pot.fillna(0.0)))
        print(f"pot:{pot[-1]}")
        if 'del_p' in mcdf.columns:
            mcdf.rename(columns={'del_p':'del_p_true'}, inplace=True)

        DETECTOR = recodf.detector.iloc[0]
        merged = recodf.reset_index().merge(
            hdrdf.reset_index(),
            on=['__ntuple', 'entry'],
            how='right'  # Keeps all hdrdf rows
        )
        merged['rec.slc..index'] = merged['rec.slc..index'].fillna(-1).astype('int64')
        merged.loc[merged['rec.slc..index'] != 0, 'pot'] = 0.0
       
        recodf = merged.set_index(['__ntuple', 'entry', 'rec.slc..index'])

        matchdf = recodf.copy()
        matchdf.columns = pd.MultiIndex.from_tuples([(col, '') for col in matchdf.columns])

        matchdf = ph.multicol_merge(matchdf.reset_index(), mcdf.reset_index(),
                               left_on=[("__ntuple", ""), ("entry", ""), ("tmatch_idx", "")],
                               right_on=[("__ntuple", ""), ("entry", ""), ("rec.mc.nu..index", "")],
                               #how="outer").set_index(['rec.mc.nu..index', "rec.slc..index"])
                               how="outer").set_index(['__ntuple', 'entry', "rec.slc..index"])
    
        pd.options.display.max_rows =100 
        pd.set_option("display.max_rows", 100) 
        pd.options.display.min_rows =50 
        pd.set_option("display.min_rows", 50)
        unfilt_dataframes.append(matchdf)
    
    dataframes, frac = filter_n_common_events(unfilt_dataframes, keys=["run", "subrun", "evt", "nu_E"])

    return dataframes, [frac[i]*pot[i] for i in range(len(pot))]

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
    data = {
            (col_name, 'ps1') : [func(E, del_p) for E, del_p in zip(df.nu_E_calo, df.del_p)],
            }

    new_df = pd.DataFrame(data)
    return new_df

import pandas as pd
from functools import reduce

def filter_n_common_events(dfs, keys=['run', 'subrun', 'evt', 'nu_E_true']):
    if not dfs:
        return []

    id_sets = (df[keys].drop_duplicates() for df in dfs)
    print("0:")
    print(dfs[0].sort_values(keys))
    print("1:")
    print(dfs[1].sort_values(keys))
    print("merge:")
    print(pd.merge(dfs[0], dfs[1], on=keys).sort_values(keys))
    common_ids = reduce(lambda left, right: pd.merge(left, right, on=keys), id_sets)
    n_common = len(common_ids)
    filtered_dfs = []

    allindices = [set(df[keys].to_records(index=False).tolist()) for df in dfs]

    frac = []
    for df, indices in zip(dfs, allindices):
        index_names = df.index.names
        n_original = len(df[keys].drop_duplicates())
        frac.append(n_common/n_original) 
        print(f"n_common:{n_common}, n_original: {n_original}, frac: {frac[-1]}")
        f_df = df.reset_index() \
                 .merge(common_ids, on=keys) \
                 .set_index(index_names)
        f_df.sort_values(keys, inplace=True)
        filtered_dfs.append(f_df)

    return filtered_dfs,frac 

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

    # make_hists_from_files(["/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/SBND_SpringMC_Nom.df", "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/SBND_SpringMC_WMYZ.df"], "rwt_outputs/YZ.txt")
    # make_hists_from_files(["/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/SBND_SpringMC_Nom.df", "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/SBND_SpringMC_WMXThetaXW.df"], "rwt_outputs/XThetaXW.txt")

    #make_hists_from_files(["/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/SBND_SpringMC_0xSCE.df", "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/SBND_SpringMC_0xSCE.df"], "rwt_outputs/dummy.txt")
    make_hists_from_files(["/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/SBND_SpringMC_Nom.df", "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/SBND_SpringMC_0xSCE.df"], "rwt_outputs/min_SCE.txt")
    #make_hists_from_files(["/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/SBND_SpringMC_Nom.df", "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/SBND_SpringMC_2xSCE.df"], "rwt_outputs/pls_SCE.txt")

    # dfs = []
    # for i in range(20):
    #     dfs.append(load_dfs(f"/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/SBND_SpringMC_rewgt_{i}.df", ['evt'])['evt'])
    # recodf = pd.concat(dfs)
    # recodf_precut = recodf.copy()
    # recodf = all_cuts(recodf, recodf.detector.iloc[0])
    # make_hists_from_var(PID.get_smear_vars(recodf_precut), recodf.copy(), 'SBND_smear.txt')
    # make_hists_from_var(PID.get_gain_vars(recodf_precut), recodf.copy(), 'SBND_gain.txt')

    # recodf = load_dfs(f"/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/ICARUS_SpringMCOverlay_rewgt.df", ['evt'])['evt']
    # recodf_precut = recodf.copy()
    # recodf = all_cuts(recodf, "ICARUS")
    # make_hists_from_var(PID.get_smear_vars(recodf_precut), recodf.copy(), 'ICARUS_smear.txt')
    # make_hists_from_var(PID.get_gain_vars(recodf_precut), recodf.copy(), 'ICARUS_gain.txt')

    # #just a random dataframe to test map with
    # apply_double_map(recodf.copy(), 'min_SCE.txt', 'pls_SCE.txt', 'CAFPYANA_SBN_v1_multisigma_SCE')
    # apply_map(recodf, 'YZ.txt', 'CAFPYANA_SBN_v1_multisigma_WMYZ')
    # apply_map(recodf, 'XThetaXW.txt', 'CAFPYANA_SBN_v1_multisigma_WMXThetaXW')


def plot():
    # SCE
    prefix = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/"
    dataframes, pot = match_across_detvar([prefix+"SBND_SpringMC_0xSCE.df", prefix+"SBND_SpringMC_0xSCE.df", ])
    #dataframes, pot = match_across_detvar([prefix+"SBND_SpringMC_0xSCE.df", prefix+"SBND_SpringMC_Nom.df", prefix+"SBND_SpringMC_2xSCE.df"])
    #for i in range(len(dataframes)):
    #    dataframes[i] = all_cuts(dataframes[i], dataframes[i].detector.iloc[0])
    #make_plots(dataframes, norms=pot, file_titles=["0xSCE", "CV", "2xSCE"], global_title="SCE Detector Variations", output_tag='SCE')

    # plot_2d_hist_from_file("rwt_outputs/min_SCE.txt", "0xSCE Reweight", 'min_SCE')
    # plot_2d_hist_from_file("rwt_outputs/pls_SCE.txt", "2xSCE Reweight", 'pls_SCE')

    # # WireMod
    # prefix = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/"
    # dataframes, pot = match_across_detvar([prefix+"SBND_SpringMC_WMXThetaXW.df", prefix+"SBND_SpringMC_WMYZ.df", prefix+"SBND_SpringMC_Nom.df",])

    # make_plots(dataframes, norms=pot, file_titles=['XThetaXW', 'YZ', 'Nominal'], global_title="WireMod Detector Variations", output_tag="WM")

    # plot_2d_hist_from_file("rwt_outputs/YZ.txt", "WireMod YZ Variation", 'WMYZ')
    # plot_2d_hist_from_file("rwt_outputs/XThetaXW.txt", "WireMod XThetaXW", 'WMXThetaXW')

    # # SBND, smear 
    # plot_2d_hist_from_file("rwt_outputs/min_SBND_smear.txt", "SBND Smear Minimum", 'min_SBND_smear')
    # plot_2d_hist_from_file("rwt_outputs/pls_SBND_smear.txt", "SBND Smear Maximum", 'pls_SBND_smear')

    # # SBND, gain
    # plot_2d_hist_from_file("rwt_outputs/min_SBND_gain.txt", "SBND Gain Minimum", 'min_SBND_gain')
    # plot_2d_hist_from_file("rwt_outputs/pls_SBND_gain.txt", "SBND Gain Maximum", 'pls_SBND_gain')

    # # ICARUS, smear
    # plot_2d_hist_from_file("rwt_outputs/min_ICARUS_smear.txt", "ICARUS Smear Minimum", 'min_ICARUS_smear')
    # plot_2d_hist_from_file("rwt_outputs/pls_ICARUS_smear.txt", "ICARUS Smear Maximum", 'pls_ICARUS_smear')

    # # ICARUS, gain
    # plot_2d_hist_from_file("rwt_outputs/min_ICARUS_gain.txt", "ICARUS Gain Minimum", 'min_ICARUS_gain')
    # plot_2d_hist_from_file("rwt_outputs/pls_ICARUS_gain.txt", "ICARUS Gain Maximum", 'pls_ICARUS_gain')

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
