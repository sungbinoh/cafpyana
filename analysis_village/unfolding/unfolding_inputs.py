import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl


cmap = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)

## A function is to make text in 2D plots more visible
def get_text_color(value):
    rgba = cmap(norm(value))
    # Compute luminance (perceived brightness)
    luminance = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
    return "black" if luminance > 0.5 else "white"


def bin_range_labels(edges):
    return [f"{edges[i]:.2f}–{edges[i+1]:.2f}" for i in range(len(edges)-1)]


## get smearing matrix
def get_smear_matrix(true_var_signal_sel_df, var_signal_df, bins_2d, weights=None, 
                     plot=True, var_labels=None,
                     save_fig=False, save_fig_name=None):
    if weights is not None:
        reco_vs_true, bins, _ = np.histogram2d(true_var_signal_sel_df, var_signal_df, bins=bins_2d, weights=weights)
    else:
        reco_vs_true, bins, _ = np.histogram2d(true_var_signal_sel_df, var_signal_df, bins=bins_2d)
    Response = reco_vs_true.T

    if plot:
        print(reco_vs_true)
        unif_bin = np.linspace(0., float(len(bins) - 1), len(bins))
        extent = [unif_bin[0], unif_bin[-1], unif_bin[0], unif_bin[-1]]
        fig, ax = plt.subplots(figsize=(10, 10))
        plt.imshow(Response, extent=extent, origin="lower", cmap="viridis")
        plt.colorbar(label="Entries", shrink=0.7)

        x_edges = np.array(bins)
        y_edges = np.array(bins)
        x_tick_positions = (unif_bin[:-1] + unif_bin[1:]) / 2
        y_tick_positions = (unif_bin[:-1] + unif_bin[1:]) / 2

        x_labels = bin_range_labels(x_edges)
        y_labels = bin_range_labels(y_edges)

        plt.xticks(x_tick_positions, x_labels, rotation=45, ha="right")
        plt.yticks(y_tick_positions, y_labels)

        if var_labels is not None:
            plt.xlabel(var_labels[2], fontsize=20)
            plt.ylabel(var_labels[1], fontsize=20)

        for i in range(Response.shape[0]):      # rows (y)
            for j in range(Response.shape[1]):  # columns (x)
                value = Response[i, j]
                if not np.isnan(value):  # skip NaNs
                    plt.text(
                        j + 0.5, i + 0.5,
                        f"{value:.1f}",
                        ha="center", va="center",
                        color=get_text_color(value),
                        fontsize=10
                    )

        plt.title("Smearing", fontsize=20)
        plt.tight_layout()

        if save_fig and save_fig_name is not None:
            plt.savefig("{}.pdf".format(save_fig_name), bbox_inches="tight")
        plt.show()

    return reco_vs_true


# efficiency
def get_eff(reco_vs_true, true_signal):
    reco_vs_true.T.sum(axis=0)
    eff = reco_vs_true.T.sum(axis=0) / true_signal  # efficiency per truth bin
    return eff


def get_response_matrix(reco_vs_true, eff, bins, 
                     plot=True, var_labels=None,
                     save_fig=False, save_fig_name=None):
    denom = reco_vs_true.T.sum(axis=0)
    num = reco_vs_true.T
    Response = np.divide(
        num * eff, denom,
        out=np.zeros_like(num, dtype=float),  # fill with 0 where invalid
        where=denom != 0
    )

    if plot:
        fig, ax = plt.subplots(figsize=(10, 10))
        unif_bin = np.linspace(0., float(len(bins) - 1), len(bins))
        extent = [unif_bin[0], unif_bin[-1], unif_bin[0], unif_bin[-1]]
        plt.imshow(Response, extent=extent, origin="lower", cmap="viridis")
        plt.colorbar(label="Response", shrink=0.7)

        x_edges = np.array(bins)
        y_edges = np.array(bins)
        x_tick_positions = (unif_bin[:-1] + unif_bin[1:]) / 2
        y_tick_positions = (unif_bin[:-1] + unif_bin[1:]) / 2

        x_labels = bin_range_labels(x_edges)
        y_labels = bin_range_labels(y_edges)

        plt.xticks(x_tick_positions, x_labels, rotation=45, ha="right")
        plt.yticks(y_tick_positions, y_labels)

        if var_labels is not None:
            plt.xlabel(var_labels[2], fontsize=20)
            plt.ylabel(var_labels[1], fontsize=20)

        for i in range(Response.shape[0]):      # rows (y)
            for j in range(Response.shape[1]):  # columns (x)
                value = Response[i, j]
                if not np.isnan(value):  # skip NaNs
                    plt.text(
                        j + 0.5, i + 0.5,
                        f"{value:.3f}",
                        ha="center", va="center",
                        color=get_text_color(value),
                        fontsize=10
                    )

        plt.title("Response", fontsize=20)
        plt.tight_layout()

        if save_fig and save_fig_name is not None:
            plt.savefig("{}.pdf".format(save_fig_name), bbox_inches="tight")
        plt.show()
    return Response

