import pandas as pd
import os
import sys
from cycler import cycler

#from analysis_village.gump.plot_tools import *
# Add the head direcoty to sys.path
workspace_root = os.getcwd()
sys.path.insert(0, workspace_root + "/../../")
import pyanalib.pandas_helpers as ph
from pyanalib.split_df_helpers import *
from makedf.util import *
from analysis_village.gump.gump_cuts import *

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
    return new_df

def main():
    """Main analysis pipeline."""

    make_sce_hists(["SCE_0.df", "SCE_1.df", "SCE_2.df"])

    prefix = "/exp/sbnd/data/users/nrowe/GUMP/det_syst/"

    # just a random dataframe to test map with
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
