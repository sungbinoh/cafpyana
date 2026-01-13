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
#     b_E = np.array([0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.25, 1.5, 2.0, 3.0])
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

def plot_sce(files):
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
            recodf = recodf[del_p_cut(recodf)]
            temp_list.append(recodf)
        dataframes.append(pd.concat(temp_list))

    # 2. Determine global max for shared colorbar scale
    # We calculate the histograms once to find the highest bin value
    max_counts = []
    for df in dataframes:
        counts, _, _ = np.histogram2d(df.nu_E_calo, df.del_p, bins=[b_E, b_p])
        max_counts.append(counts.max())
    global_vmax = max(max_counts)

    titles = ["0xSCE", "CV", "2xSCE"]
    for df, t in zip(dataframes, titles):
        plt.hist(df.nu_E_calo, bins=b_E, label=t, histtype='step')

    plt.xlabel(r'Reconstructed Energy $E_{calo}$ [GeV]')
    plt.ylabel('Events')
    plt.title("GUMP Selection SCE Detector Variations", fontsize=20)
    plt.legend()
    plt.savefig('E_SCE.png', dpi=300)
    plt.clf() # Better than fig.clf() for memory management

    for df, t in zip(dataframes, titles):
        plt.hist(df.del_p, bins=b_p, label=t, histtype='step')

    plt.xlabel(r'$\delta p$ [GeV/c]')
    plt.ylabel('Events')
    plt.title("GUMP Selection SCE Detector Variations", fontsize=20)
    plt.legend()
    plt.savefig('p_SCE.png', dpi=300)
    plt.clf() # Better than fig.clf() for memory management

    # 1. Change dimensions: figsize controls the width and height in inches
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(18, 5), constrained_layout=True)

    hist2ds = []
    # Second pass: Plotting
    for df, a, t in zip(dataframes, axes, titles):
        # 3. Apply shared scale using vmin and vmax
        im = a.hist2d(df.nu_E_calo, df.del_p, bins=[b_E, b_p], vmin=0, vmax=global_vmax)
        hist2ds.append(np.transpose(im[0]))

        # 4. Add individual x and y axis labels
        a.set_xlabel(r'Reconstructed Energy $E_{calo}$ [GeV]')
        a.set_ylabel(r'$\delta p$ [GeV/c]')
        a.set_title(t)
    hist2ds = np.array(hist2ds)

    # 5. Add a colorbar with the same scale
    # 'im[3]' refers to the QuadMesh object returned by hist2d
    cbar = fig.colorbar(im[3], ax=axes.ravel().tolist(), label='Events')
    fig.suptitle("GUMP Selection SCE Detector Variations", fontsize=40)

    fig.savefig('2dSCE.png', dpi=300)
    plt.close(fig) # Better than fig.clf() for memory management

    hist2ds_del = [hist2ds[1] - hist2ds[0], hist2ds[2] - hist2ds[1]]

    diff1 = hist2ds[1] - hist2ds[0]
    diff2 = hist2ds[2] - hist2ds[1]
    diffs = [diff1, diff2]

    # 4. Determine Symmetric Scale for Diverging Colorbar
    # To keep 0 as white, vmax and vmin must be equal and opposite
    abs_max = max(np.abs(diff1).max(), np.abs(diff2).max())

    # 5. Plotting using pcolormesh
    titles = ["Difference (2 - 1)", "Difference (3 - 2)"]
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5), constrained_layout=True)
    for i, (ax, data) in enumerate(zip(axes, diffs)):
        # pcolormesh needs the transpose of the counts array if using [x, y]
        im = ax.pcolormesh(b_E, b_p, data, cmap='seismic', vmin=-abs_max, vmax=abs_max)
        
        ax.set_title(titles[i])
        ax.set_xlabel(r'$E_{calo}$ [GeV]')
        ax.set_ylabel(r'$\Delta p$ [GeV/c]')

    # 6. Global formatting
    fig.suptitle('Residual Differences between SCE Iterations', fontsize=16)
    cbar = fig.colorbar(im, ax=axes, label='$\Delta$ Events')


    #titles = ['CV-0xSCE','2xSCE-CV']
    #for hd, a, t in zip(hist2ds_del, axes, titles):
    #    a.imshow(hd, aspect='auto')
    #    a.set_xlabel(r'Reconstructed Energy $E_{calo}$ [GeV]')
    #    a.set_ylabel(r'$\delta p$ [GeV/c]')
    #    a.set_title(t)

    #cbar = fig.colorbar(im[3], ax=axes.ravel().tolist(), label=r'$\Delta$ Events')
    fig.suptitle("GUMP Selection SCE Detector Variations", fontsize=40)
    fig.savefig('2ddeltaSCE.png', dpi=300)
    plt.close(fig)

def main():
    """Main analysis pipeline."""

    # plot_sce(["SCE_0.df"])
    plot_sce(["SCE_0.df", "SCE_1.df", "SCE_2.df"])

if __name__ == "__main__":
    main()


