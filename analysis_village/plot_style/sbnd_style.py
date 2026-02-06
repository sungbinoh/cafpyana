import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap

# =========================
# Okabe-Ito palette
# =========================
OKABE_ITO_COLORS = {
    "orange": (0.90, 0.60, 0.0),
    "sky_blue": (0.35, 0.70, 0.90),
    "blue_green": (0.0, 0.60, 0.50),
    "yellow": (0.95, 0.90, 0.25),
    "blue": (0.0, 0.45, 0.70),
    "vermilion": (0.80, 0.30, 0.0),
    "red_purple": (0.80, 0.60, 0.70),
}

COLOR_CYCLES = {
    "okabe_ito": [
        "black",
        OKABE_ITO_COLORS["vermilion"],
        OKABE_ITO_COLORS["sky_blue"],
        OKABE_ITO_COLORS["orange"],
        OKABE_ITO_COLORS["blue_green"],
        OKABE_ITO_COLORS["red_purple"],
        OKABE_ITO_COLORS["blue"],
        OKABE_ITO_COLORS["yellow"]
    ],
    "sbnd_logo": [
        OKABE_ITO_COLORS["vermilion"],
        OKABE_ITO_COLORS["blue_green"],
        OKABE_ITO_COLORS["blue"],
        OKABE_ITO_COLORS["sky_blue"]
    ]
}

# =========================
# Core style functions
# =========================
def apply_color_cycle(cycle_name="okabe_ito"):
    colors = COLOR_CYCLES.get(cycle_name, COLOR_CYCLES["okabe_ito"])
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors)

def apply_cvd_palette():
    plt.set_cmap("cividis")

def apply_sea_palette():
    colors = [(1, 1, 1), (0.0, 0.45, 0.70)]
    if "SeaPalette" not in plt.colormaps():
        cmap = LinearSegmentedColormap.from_list("SeaPalette", colors)
        plt.register_cmap("SeaPalette", cmap)
    plt.set_cmap("SeaPalette")


def apply_symmetric_palette():
    colors = [
        (0.0, 0.45, 0.70),
        (1.0, 1.0, 1.0),
        (0.8, 0.3, 0.0)
    ]
    if "SymmetricPalette" not in plt.colormaps():
        cmap = LinearSegmentedColormap.from_list("SymmetricPalette", colors)
        plt.register_cmap(name="SymmetricPalette", cmap=cmap)
    plt.set_cmap("SymmetricPalette")

_color_counter = {c: 0 for c in COLOR_CYCLES.keys()}

def next_color(cycle="okabe_ito", start=None):
    """Return next color in a given SBND color cycle."""
    colors = COLOR_CYCLES[cycle]
    if start is not None:
        _color_counter[cycle] = start % len(colors)
    idx = _color_counter[cycle]
    _color_counter[cycle] = (idx + 1) % len(colors)
    return colors[idx]

# =========================
# SBND text helpers
# =========================

def sbnd_watermark():
    return r"$\bf{SBND}$"

def sbnd_wip(ax, x, y, fontsize=30):
    ax.text(x, y, f"{sbnd_watermark()} Work in Progress", horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, fontsize=fontsize,  **{'fontname':'Helvetica'})

def sbnd_preliminary(ax, x, y, ha='left', va='top', fontsize=30):
    ax.text(x, y, f"{sbnd_watermark()} Preliminary", horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, fontsize=fontsize,  **{'fontname':'Helvetica'})

def sbnd_data(ax, x, y, ha='left', va='top', fontsize=30):
    ax.text(x, y, f"{sbnd_watermark()} Data", horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, fontsize=fontsize,  **{'fontname':'Helvetica'})

def sbnd_official(ax, x, y, fontsize=30):
    ax.text(x, y, f"{sbnd_watermark()}", horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, fontsize=fontsize,  **{'fontname':'Helvetica'})

# =========================
# Layout helpers
# =========================
def split_canvas(fig, ysplit=0.23):
    gs = fig.add_gridspec(2, 1, height_ratios=[1 - ysplit, ysplit], hspace=0.)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    return ax1, ax2

def center_titles(ax):
    ax.title.set_ha('center')
    ax.set_title(ax.get_title(), ha='center')
    ax.xaxis.label.set_ha('center')
    ax.yaxis.label.set_va('center')
    ax.yaxis.labelpad = 22

# =========================
# DEFAULT STYLE (auto-apply)
# =========================
def _apply_default_style():
    apply_color_cycle("okabe_ito")

    plt.rcParams.update({
        "figure.figsize": (8, 6),
        "errorbar.capsize": 3,
        "savefig.dpi": 300,
        "figure.dpi": 100,
        "font.family": "sans-serif",
        "font.size": 14,
        "axes.labelsize": 25,
        "axes.titlesize": 25,
        "axes.linewidth": 1,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "xtick.direction": "in",
        "xtick.major.size": 6,
        "ytick.major.size": 6,
        "xtick.minor.size": 3,
        "ytick.minor.size": 3,
        "xtick.top": False,
        "ytick.right": False,
        "lines.linewidth": 2.0,
        "lines.markersize": 6,
        "legend.frameon": False,
        "legend.fontsize": 20,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "figure.subplot.left": 0.15,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.15,
        "figure.subplot.top": 0.92

    })

_apply_default_style()
