import pandas as pd
import os
import sys
from cycler import cycler
from plot_tools import *
# Add the head direcoty to sys.path
workspace_root = os.getcwd()
sys.path.insert(0, workspace_root + "/../../")
import pyanalib.pandas_helpers as ph
from pyanalib.split_df_helpers import *
from makedf.util import *

def all_cuts(recodf, DETECTOR):
    ## Apply cuts
    slc_vtx = pd.DataFrame({'x':recodf.slc_vtx_x,
                            'y':recodf.slc_vtx_y,
                            'z':recodf.slc_vtx_z})

    recodf = recodf[fv_cut(slc_vtx, DETECTOR)]

    ### NuScore cut
    recodf = recodf[cosmic_cut(recodf)]

    ### Two prong cut
    recodf = recodf[twoprong_cut(recodf)]

    ### containment cut
    recodf = recodf[mufv_cut(recodf, DETECTOR)]
    recodf = recodf[pfv_cut(recodf, DETECTOR)]

    ### PID cut
    recodf = recodf[pid_cut(recodf.mu_chi2_of_mu_cand, recodf.mu_chi2_of_prot_cand,
                            recodf.prot_chi2_of_mu_cand, recodf.prot_chi2_of_prot_cand,
                            recodf.mu_len)]

    ### crthitveto cut
    if DETECTOR == "ICARUS":
        recodf = recodf[crthitveto_cut(recodf)]

    return recodf

# def plot_sce(files):
#     
#     b_p = np.array([0.0, 0.2, 0.4, 0.6])
# 
#     fig, axes = plt.subplots(nrows=1, ncols=3)
# 
#     for f, a in zip(files, axes):
#         prefix = "/exp/sbnd/data/users/nrowe/GUMP/det_syst/"
# 
#         nsplits = get_n_split(prefix+f)
# 
#         for i in range(nsplits):
#             recodf = pd.read_hdf(prefix+f, key='evt_'+str(i))
# 
#             ## Figure out which detector this is
#             DETECTOR = recodf.detector.iloc[0]
# 
#             recodf = all_cuts(recodf, DETECTOR)
# 
#             if i == 0: filedf = recodf.copy()
#             else: filedf = pd.concat([filedf, recodf])
#             del recodf
# 
#         a.hist2d(filedf.nu_E_calo, filedf.del_p, bins=[b_E, b_p])
# 
#     fig.savefig('2dSCE.png')
#     fig.clf()

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
            return 0.0
             
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
    np.savetxt(filename, hist_values, header=header, delimiter=",")

def make_sce_hists(files):
    b_E = np.array([0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.25, 1.5, 2.0, 3.0])
    b_p = np.array([0.0, 0.2, 0.4, 0.6])
    
    prefix = "/exp/sbnd/data/users/nrowe/GUMP/det_syst/"
    dataframes = []

    # First pass: Load data to find global min/max for the colorbar
    for f in files:
        nsplits = get_n_split(prefix + f)
        temp_list = []
        for i in range(nsplits):
            recodf = pd.read_hdf(prefix + f, key='evt_' + str(i))
            DETECTOR = recodf.detector.iloc[0]
            recodf = all_cuts(recodf, DETECTOR)
            temp_list.append(recodf)
        dataframes.append(pd.concat(temp_list))

    hist2ds = []
    for df in dataframes:
        im = plt.hist2d(df.nu_E_calo, df.del_p, bins=[b_E, b_p])
        hist2ds.append(im[0])
        # hist2ds.append(np.transpose(im[0]))
    plt.clf()

    hist2ds = np.array(hist2ds)
    hist2ds_rat = [hist2ds[0]/hist2ds[1], hist2ds[2]/hist2ds[1]]

    save_histogram('min_SCE.txt', hist2ds_rat[0], b_E, b_p)
    save_histogram('pls_SCE.txt', hist2ds_rat[1], b_E, b_p)

def apply_sce_map(df, min_map_file, pls_map_file):
    min_func = FileHistogramFunction(min_map_file)
    pls_func = FileHistogramFunction(pls_map_file)
    data = {
            ('CAFPYANA_SBN_v1_multisigma_SCE', 'ps1') : [pls_func(E, del_p) for E, del_p in zip(df.nu_E_calo, df.del_p)],
            ('CAFPYANA_SBN_v1_multisigma_SCE', 'ms1') : [min_func(E, del_p) for E, del_p in zip(df.nu_E_calo, df.del_p)]
            }

    new_df = pd.DataFrame(data)

    w = np.ones_like(df.nu_E_calo)
    w[df.del_p > 0.6] = 0.0

    plt.hist(df.nu_E_calo, bins = min_func.x_edges, weights=w, label='cv', histtype='step')
    plt.hist(df.nu_E_calo, bins = min_func.x_edges, weights=new_df['CAFPYANA_SBN_v1_multisigma_SCE']['ms1'], label='ms1-rwgt', histtype='step')
    plt.hist(df.nu_E_calo, bins = min_func.x_edges, weights=new_df['CAFPYANA_SBN_v1_multisigma_SCE']['ps1'], label='ps1-rwgt', histtype='step')
    plt.legend()
    plt.savefig('sce_debug.png')
    return new_df

def main():
    """Main analysis pipeline."""

    make_sce_hists(["SCE_0.df", "SCE_1.df", "SCE_2.df"])

    prefix = "/exp/sbnd/data/users/nrowe/GUMP/det_syst/"

    cv_f = "/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-3-fix-zexp-again-again/SBND_SpringMC_rewgt_1.df"
    nsplits = get_n_split(cv_f)
    temp_list = []
    for i in range(nsplits):
        recodf = pd.read_hdf(cv_f, key='evt_' + str(i))
        DETECTOR = recodf.detector.iloc[0]
        recodf = all_cuts(recodf, DETECTOR)
        temp_list.append(recodf)

    recodf = pd.concat(temp_list)
    apply_sce_map(recodf, 'min_SCE.txt', 'pls_SCE.txt')

if __name__ == "__main__":
    main()
