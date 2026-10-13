"""
rwt_map.py -- detector-systematic outputs for GUMP / GUMPLE.

For each detector variation this writes either or both of:
  * a reweight map (default): the 2D (nu_E_calo, del_p) bin-by-bin ratio
    (variation / CV) of selected, POT-scaled event counts, saved as a text grid
    with the bin edges in the header;
  * sbruce trees (--sbruce-trees): the selected events of the variation, and of
    the CV it is compared to, converted with the sbruce tools.

Variations that come from a dedicated sample (WireMod, SCE, DENT, ...) are
matched to the CV event-by-event first. Matching is done per *group*: by default
each variation is its own group, so it is matched against the CV on its own and
gets its own matched CV. This mirrors load_detvar() in
SignalBoxSystematics-GUMPLE.ipynb. Put several variations in one group only if
they should share a single jointly-matched CV.

Variations derived from the CV itself (binding energy, track splitting, chi2 /
dE/dx, trigger) need no matching and are compared to the full, unmatched CV.
"""

import argparse
import glob
import importlib
import os
import sys
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

import matplotlib
matplotlib.use("Agg")  # files only; avoids slow/hanging GUI backends over X forwarding
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

# Paths are resolved relative to this file (analysis_village/gumple/), so the
# script behaves the same regardless of the working directory.
_HERE = os.path.dirname(os.path.abspath(__file__))
for _path in (os.path.join(_HERE, "..", ".."),      # repo root
              os.path.join(_HERE, "..", "gump"),    # loaddf, syst
              _HERE):                               # gumple_cuts
    sys.path.insert(0, os.path.normpath(_path))

import loaddf
import syst
import gumple_cuts as gmpl


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Output locations are relative to the directory the script is run from, so
# each selection's run directory gets its own outputs.
DEFAULT_OUTDIR = "rwt_outputs"
DEFAULT_TREEOUTDIR = "rwt_outputs"
DEFAULT_PLOTDIR = "rwt_output_plots"
DEFAULT_DFDIR = "/exp/sbnd/data/users/gputnam/GUMPLE/sbn-rewgted-22/"

# [nu_E_calo edges, del_p edges].
#   "1D" : single effective del_p bin; reco-E bins chosen for the
#          statistics of the 1p vs N>1p channels
#   "2D" : production binning
BINNINGS = {
    "1D": [np.array([0.3, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.10, 1.25, 1.5]),
           np.array([-1000.0, 0.0, 1000.0])],
    "2D": [np.array([0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.25, 1.5]),
           np.array([0.0, 0.2, 0.4, 0.6])],
}

# Variations built directly from the (unmatched) CV dataframe: (name, function).
CV_DERIVED_VARIATIONS = [
    ("Smeared dE/dx", syst.v_chi2smear),
    ("Biased dE/dx",  syst.v_chi2dedxbias),
    ("Gain Low",       syst.v_chi2lo),
    ("Gain Hi",       syst.v_chi2hi),
    ("EMB Alpha p",   syst.v_chi2alphap),
    ("EMB Beta p",    syst.v_chi2betap),
    ("EMB R p",       syst.v_chi2Rp),
    ("EMB Alpha m",   syst.v_chi2alpham),
    ("EMB Beta m",    syst.v_chi2betam),
    ("EMB R m",       syst.v_chi2Rm),
    ("TrigEffPls",    lambda df: syst.v_flashscale(df, 1)),
    ("TrigEffMin",    lambda df: syst.v_flashscale(df, -1)),
]

# --- sbruce tree settings ---
SBRUCE_SUBDIR = "sbruce_trees"   # created inside the output directory
DEFAULT_WEIGHT_COL = "cvwgt"
# Tree name for the unmatched CV, i.e. the nominal for BIND, TRKSPLT and the
# CV-derived variations. Each matched group's CV is written as <names>Nominal.
UNMATCHED_CV_TREE_NAME = "BINDNominal"


def _icarus_config(run, goal_pot):
    prefix = f"ICARUSRun{run}_Spring_Overlay_"
    return {
        "goal_pot": goal_pot,
        "cv": f"ICARUSRun{run}_SpringMCOverlay_rewgt_*.df",
        "trksplt_regions": ["Z=0", "East Cathode", "West Cathode"],
        "ratio_scales": {},
        "match_sets": [
            {
                "cv": None,  # None -> use the main CV above
                "groups": [
                    {"WMXThetaXW": [prefix + "WMXThXW.df"]},
                    {"WMYZ":       [prefix + "WMYZ.df"]},
                    {"SCE":        [prefix + "SCE.df"]},
                ],
            },
        ],
    }


# Per-detector settings. File entries are relative to --dfdir and may be globs.
#   goal_pot        : POT every sample is scaled to
#   cv              : main CV sample
#   trksplt_regions : track-splitting regions (ICARUS only)
#   ratio_scales    : per-variation multiplier on (var - cv) in the maps, e.g. 2x smearing
#   match_sets      : groups of dedicated-sample variations, see module docstring.
#                     Optional keys: "light" (lightmem load, default True) and
#                     "binning" (override the --binning choice).
DETECTOR_CONFIG = {
    "SBND": {
        "goal_pot": 1e20,
        "cv": "SBNDMCCV_*.df",
        "trksplt_regions": [],
        # SBND dE/dx smearing systematic: 13% -> 2x13%. Applied to the map only;
        # the sbruce tree holds the 1x-smeared events and the fitter applies the 2x.
        "ratio_scales": {"Smeared dE/dx": 2.0},
        "match_sets": [
            {
                "cv": None,
                "groups": [
                    {"WMXThetaXW": ["SBND_SpringMC_WMXThetaXW.df"]},
                    {"WMYZ":       ["SBND_SpringMC_WMYZ.df"]},
                ],
            },
            {
                # The SBND SCE / DENT samples have their own nominal sample.
                # 2xSCE / 0xSCE are the +/-1 sigma SCE variations; each is
                # matched against Nom separately (one group each).
                # Previously loaded with loaddf.load (no lightmem) and always
                # binned in 2D regardless of --binning; both kept here.
                "cv": ["SBND_SpringMC_Nom.df"],
                "light": False,
                "binning": "2D",
                "groups": [
                    {"2xSCE": ["SBND_SpringMC_2xSCE.df"]},
                    {"0xSCE": ["SBND_SpringMC_0xSCE.df"]},
                    {"DENT":  ["SBND_SpringMC_DENT.df"]},
                ],
            },
        ],
    },
    "ICARUS Run2": _icarus_config(2, 2e20),
    "ICARUS Run4": _icarus_config(4, 3e20),
}


# ---------------------------------------------------------------------------
# Map file I/O
# ---------------------------------------------------------------------------

def _parse_edges(line):
    return np.array([float(v) for v in line.strip().lstrip("#").split(",") if v.strip()])


def read_map_file(filename):
    """Return (x_edges, y_edges, grid) from a file written by save_histogram."""
    with open(filename, "r") as f:
        x_edges = _parse_edges(f.readline())
        y_edges = _parse_edges(f.readline())
    grid = np.loadtxt(filename, delimiter=",", ndmin=2)  # '#' header lines are skipped
    return x_edges, y_edges, grid


def save_histogram(filename, hist_values, x_edges, y_edges):
    header = ",".join(str(x) for x in x_edges) + ",\n" + ",".join(str(y) for y in y_edges) + ","
    print(f"Saving: {filename}")
    np.savetxt(filename, hist_values, header=header, delimiter=",")


class FileHistogramFunction:
    """Look up per-event weights from a saved map. Out-of-range or NaN -> 1."""

    def __init__(self, filename):
        self.x_edges, self.y_edges, self.grid = read_map_file(filename)

    def __call__(self, x_arr, y_arr):
        ix = np.digitize(x_arr, self.x_edges) - 1
        iy = np.digitize(y_arr, self.y_edges) - 1
        mask = (ix >= 0) & (ix < self.grid.shape[0]) & (iy >= 0) & (iy < self.grid.shape[1])

        result = np.ones(len(x_arr), dtype=float)
        result[mask] = self.grid[ix[mask], iy[mask]]
        return np.nan_to_num(result, nan=1.0)


def apply_map(df, map_file, col_name):
    """Column of per-event weight lists: [1 (nominal), w_map1, w_map2, ...]."""
    map_files = [map_file] if isinstance(map_file, (str, bytes)) else map_file
    weights = [np.ones(len(df))]
    for mf in map_files:
        weights.append(FileHistogramFunction(mf)(df.nu_E_calo.values, df.del_p.values))
    return pd.DataFrame({col_name: np.column_stack(weights).tolist()}, index=df.index)


# ---------------------------------------------------------------------------
# Debug plots
# ---------------------------------------------------------------------------

def plot_2d_hist_from_file(filename, plot_title, output_tag, plotdir=DEFAULT_PLOTDIR):
    x_edges, y_edges, grid = read_map_file(filename)
    z = np.ma.masked_invalid(grid.T)

    # Colour scale symmetric about 1 so white always means "no change".
    dev = float(np.abs(z - 1).max()) if z.count() else 0.0
    dev = dev if dev > 0 else 0.01

    cmap = plt.get_cmap("seismic").copy()
    cmap.set_bad(color="gray")  # empty CV bins

    fig, ax = plt.subplots(figsize=(10, 6))
    mesh = ax.pcolormesh(x_edges, y_edges, z, cmap=cmap, vmin=1 - dev, vmax=1 + dev)
    fig.colorbar(mesh, ax=ax, label="Ratio (Variation / CV)")
    ax.set_title(plot_title)
    ax.set_xlabel(r"Reconstructed Energy $E_{calo}$ [GeV]")
    ax.set_ylabel(r"$\delta p$ [GeV/c]")

    os.makedirs(plotdir, exist_ok=True)
    fig.savefig(os.path.join(plotdir, f"2d_ratio_{output_tag}.png"), dpi=300)
    plt.close(fig)


def plot_del_p_slices(filename, plot_title, output_tag, plotdir=DEFAULT_PLOTDIR):
    """Ratio vs nu_E_calo, one curve per del_p bin of the map."""
    x_edges, y_edges, grid = read_map_file(filename)
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])

    fig, ax = plt.subplots(figsize=(10, 6))
    for iy in range(len(y_edges) - 1):
        label = rf"$\delta p \in [{y_edges[iy]:.1f}, {y_edges[iy+1]:.1f}]$ GeV/c"
        ax.plot(x_centers, grid[:, iy], marker="o", markersize=4, label=label)

    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1)
    ax.set_xlabel(r"Reconstructed Energy $E_{calo}$ [GeV]")
    ax.set_ylabel("Ratio (Variation / CV)")
    ax.set_title(plot_title)
    ax.legend()
    ax.grid(True, alpha=0.3)

    os.makedirs(plotdir, exist_ok=True)
    fig.savefig(os.path.join(plotdir, f"slice_{output_tag}.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_map(filename, plotdir=DEFAULT_PLOTDIR):
    """Make the 2D map and del_p-slice debug plots for one saved map file."""
    tag = os.path.splitext(os.path.basename(filename))[0]
    title = tag.replace("_", " ")
    plot_2d_hist_from_file(filename, title, tag, plotdir)
    plot_del_p_slices(filename, title + r" -- $\delta p$ slices", tag, plotdir)


# ---------------------------------------------------------------------------
# sbruce trees
# ---------------------------------------------------------------------------

class SbruceTreeWriter:
    """
    Write the selected events of a sample as an sbruce tree: a flat TTree via
    export_dataframe_to_uproot, converted by run_makesbruce_macro. The flat file
    is removed if the conversion succeeds and kept for inspection if it fails.

    tree_vars=None (the default) writes every column of the dataframe; otherwise
    only tree_vars plus weight_col are written.
    """

    def __init__(self, treeoutdir, tree_vars=None, weight_col=DEFAULT_WEIGHT_COL):
        try:  # only needed when trees are requested, so imported here
            from sbruce import export_dataframe_to_uproot, run_makesbruce_macro
        except ImportError as e:
            raise ImportError(f"--sbruce-trees needs the sbruce module on the Python path: {e}")
        self._export = export_dataframe_to_uproot
        self._convert = run_makesbruce_macro
        self.weight_col = weight_col
        self.columns = None if tree_vars is None else list(dict.fromkeys(list(tree_vars) + [weight_col]))
        self.outdir = os.path.join(treeoutdir, SBRUCE_SUBDIR)
        os.makedirs(self.outdir, exist_ok=True)

    def __call__(self, selected_df, tag):
        columns = self.columns if self.columns is not None else list(selected_df.columns)

        # A missing weight column only warns; the tree is still written without it.
        if self.weight_col not in selected_df.columns:
            print(f"  [!] {tag}: weight column '{self.weight_col}' not in dataframe, "
                  "writing the tree without it")
            columns = [c for c in columns if c != self.weight_col]

        # Explicitly requested (-v) variables must all be present.
        missing = [c for c in columns if c not in selected_df.columns]
        if missing:
            print(f"  [!] Skipping sbruce tree for {tag}: missing columns {missing}")
            return

        flat_file = os.path.join(self.outdir, f"{tag}_flat.root")
        sbruce_file = os.path.join(self.outdir, f"{tag}_sbruce.root")
        print(f"Writing sbruce tree: {sbruce_file}")
        self._export(selected_df[columns].copy(), flat_file, tree_name="SelectedEvents")

        if self._convert(flat_file, sbruce_file) == 0:
            os.remove(flat_file)
        else:
            print(f"  [!] sbruce conversion failed for {tag}, keeping {flat_file} for inspection")


# ---------------------------------------------------------------------------
# Loading / histogramming helpers
# ---------------------------------------------------------------------------

def _copy(x):
    return x.copy() if hasattr(x, "copy") else x


def _set_total_pot(df, pot):
    """After matching, make loaddf's total_pot column (if present) the sample's matched POT."""
    if "total_pot" not in df.columns:
        return
    if np.ndim(pot) != 0:
        raise TypeError(f"Expected a single POT value after matching, got {type(pot).__name__} "
                        f"with shape {np.shape(pot)}; total_pot can't be updated")
    df["total_pot"] = float(pot)


def _clean_name(name):
    return name.replace("/", "").replace(" ", "")


def _resolve_files(df_dir, patterns):
    """Expand file names / globs relative to df_dir. Fails loudly if nothing matches."""
    if isinstance(patterns, str):
        patterns = [patterns]
    files = []
    for p in patterns:
        matched = sorted(glob.glob(os.path.join(df_dir, p)))
        if not matched:
            raise FileNotFoundError(f"No files match {os.path.join(df_dir, p)}")
        files += matched
    return files


def load_sample(files, detector, light=True, **kwargs):
    slim_drops = [
    # Misc / matching
    "charge_center_z", "tmatch_idx", "tmatch_eff", "slice_index", "crthi_ismc",

    # True muon
    "true_mu_end_x", "true_mu_end_y", "true_mu_end_z",

    # True proton
    "true_p_end_x", "true_p_end_y", "true_p_end_z",

    # True second proton
    "true_p2_p", "true_p2_dir_x", "true_p2_dir_y", "true_p2_dir_z",

    # True charged pion
    "true_cpi_dir_x", "true_cpi_dir_y", "true_cpi_pdg",

    # True gamma
    "true_g_p", "true_g_dir_x", "true_g_dir_y", "true_g_dir_z",
    "true_g_end_x", "true_g_end_y", "true_g_end_z",

    # True pi0
    "true_pi0_p", "true_pi0_dir_x", "true_pi0_dir_y", "true_pi0_dir_z",
    "true_pi0_end_x", "true_pi0_end_y", "true_pi0_end_z",

    # Spill info
    "spill_TOR860", "spill_TOR875", "spill_LM875A", "spill_LM875B",
    "spill_LM875C", "spill_THCURR", "spill_FOM",

    # Early cuts and PFP counts
    "cut_contained", "cut_cathode",
    "n_pfp_no_calo", "n_shower", "n_other", "has_muon",
    "cut_np", "cut_0shwother",

    # Momentum sum
    "psum_p", "psum_ke", "psum_E",
    "psum_dir_x", "psum_dir_y", "psum_dir_z",

    # Muon candidate (reco)
    "mu_dist_start", "mu_prim_pfp", "mu_contained10",
    "mu_true_p", "mu_true_pdg",
    "true_mucand_p",  # NOTE: likely a typo of "true_mucand_p" below

    # Muon candidate (truth)
    "true_mucand_p", "true_mucand_dir_x", "true_mucand_dir_y", "true_mucand_dir_z",
    "true_mucand_end_x", "true_mucand_end_y", "true_mucand_end_z",

    # Proton candidate
    "p_dist_to_vertex", "true_pcand_pdg",
    "true_pcand_dir_x", "true_pcand_dir_y", "true_pcand_dir_z",
    "true_pcand_end_x", "true_pcand_end_y", "true_pcand_end_z",

    # Selection cuts
    "cut_presel", "cut_cosmic", "cut_flash", "cut_trk",
    "cut_muon", "cut_protons", "cut_far_shw",

    # Selection flags / weights
    "gump_sel", "maple_sel", "common", "glob_scale",
]
    """Load a sample with the standard preselection. Returns (df, match, pot)."""
    common = dict(preselection=gmpl.slcfv_cut, include_syst=False, detector=detector)
    if light:
        return loaddf.loadl(files, lightmem=True, drops=slim_drops, **common, **kwargs)
    if len(files) != 1:
        raise ValueError(f"Non-light loads take exactly one file, got {files}")
    return loaddf.load(files[0], **common, **kwargs)


def hist2d(selected_df, bins):
    """POT-weighted (nu_E_calo, del_p) histogram of already-selected events."""
    return np.histogram2d(selected_df.nu_E_calo.to_numpy(), selected_df.del_p.to_numpy(),
                          bins=bins, weights=selected_df.glob_scale.to_numpy())[0]


def ratio_hist(var_hist, cv_hist, scale=1.0):
    """(cv + scale*(var - cv)) / cv. Bins with an empty CV are NaN (-> weight 1)."""
    numerator = cv_hist + scale * (var_hist - cv_hist)
    return np.divide(numerator, cv_hist,
                     out=np.full(cv_hist.shape, np.nan), where=cv_hist > 0)


# ---------------------------------------------------------------------------
# Processing
# ---------------------------------------------------------------------------

class _Outputs:
    """Sends each CV / variation comparison to the enabled outputs (maps and/or trees)."""

    def __init__(self, detector, selection, bins, outdir, treeoutdir, plotdir, write_maps, tree_writer):
        self.tag = detector.replace(" ", "")
        self.selection = selection
        self.bins = bins
        self.outdir = outdir
        self.treeoutdir = treeoutdir
        self.plotdir = plotdir
        self.write_maps = write_maps
        self.tree_writer = tree_writer

    def _name(self, name):
        return f"{self.tag}_{_clean_name(name)}"

    def compare(self, cv_df, variations, cv_tree_name, bins=None):
        """
        cv_df: POT-scaled CV. variations: iterable of (name, POT-scaled df, map scale),
        consumed one at a time so a generator can build each variation lazily.
        """
        bins = self.bins if bins is None else bins

        cv_sel = cv_df.loc[self.selection(cv_df)]
        cv_hist = hist2d(cv_sel, bins) if self.write_maps else None
        if self.tree_writer is not None:
            self.tree_writer(cv_sel, self._name(cv_tree_name))
        del cv_sel

        for name, var_df, scale in variations:
            var_sel = var_df.loc[self.selection(var_df)]
            if self.write_maps:
                path = os.path.join(self.outdir, f"{self._name(name)}.txt")
                save_histogram(path, ratio_hist(hist2d(var_sel, bins), cv_hist, scale), bins[0], bins[1])
                if self.plotdir is not None:
                    plot_map(path, self.plotdir)
            if self.tree_writer is not None:
                self.tree_writer(var_sel, self._name(name))
            del var_df, var_sel


def _run_match_set(match_set, main_cv, detector, df_dir, goal_pot, out, max_workers=None):
    """Process one match set, matching each group to its own copy of the CV."""
    light = match_set.get("light", True)
    bins = BINNINGS[match_set["binning"]] if match_set.get("binning") else None

    if match_set["cv"] is None:
        cv_df, cv_match, cv_pot = main_cv
    else:
        cv_df, cv_match, cv_pot = load_sample(_resolve_files(df_dir, match_set["cv"]), detector, light=light)

    def _load_group(group):
        names = list(group)
        loaded = [load_sample(_resolve_files(df_dir, files), detector, light=light)
                  for files in group.values()]
        return names, loaded

    # Pipeline: prefetch the next group's files in the background while processing
    # the current one. Only one HDF5 file is ever open at a time, avoiding
    # thread-safety issues, but I/O is hidden behind CPU work.
    groups = list(match_set["groups"])
    with ThreadPoolExecutor(max_workers=1) as pool:
        future = pool.submit(_load_group, groups[0])
        for i, group in enumerate(tqdm(groups, desc=f"{detector} matched variations")):
            names, loaded = future.result()
            if i + 1 < len(groups):
                future = pool.submit(_load_group, groups[i + 1])

            # Copies so the shared CV stays unmatched and unscaled for later groups
            dfs = [cv_df.copy()] + [l[0] for l in loaded]
            matches = [_copy(cv_match)] + [l[1] for l in loaded]
            pots = [_copy(cv_pot)] + [l[2] for l in loaded]
            del loaded

            print("Number of dfs being handed in: ", len(dfs))
            for indD, D in enumerate(dfs):
                print(f"Length of DF{indD} handed in: {len(D)}")

            dfs, pots = loaddf.match_common_evts(matches, dfs, pots)
            print("Number of dfs being coming out: ", len(dfs))
            for indD, D in enumerate(dfs):
                print(f"Length of DF{indD} coming out: {len(D)}")

            for d, p in zip(dfs, pots):
                _set_total_pot(d, p)  # loaddf's value is the pre-matching POT
                loaddf.scale_pot(d, p, goal_pot)

            cv_tree_name = "_".join(_clean_name(n) for n in names) + "Nominal"
            out.compare(dfs[0], [(n, d, 1.0) for n, d in zip(names, dfs[1:])], cv_tree_name, bins)
            del dfs


def remake_detvar_maps(detector, df_dir, selection=gmpl.all_gump_cuts, binning="2D", outdir=DEFAULT_OUTDIR, treeoutdir=DEFAULT_TREEOUTDIR,
                       plotdir=None, write_maps=True, tree_writer=None, max_workers=None):
    """
    Process all detector variations for one detector.
      write_maps  : write reweight maps to outdir
      plotdir     : if given (and maps are written), also make debug plots there
      tree_writer : an SbruceTreeWriter to also write sbruce trees, or None
    """
    cfg = DETECTOR_CONFIG[detector]
    goal_pot = cfg["goal_pot"]
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(treeoutdir, exist_ok=True)
    out = _Outputs(detector, selection, BINNINGS[binning], outdir, treeoutdir, plotdir, write_maps, tree_writer)

    products = (["maps"] if write_maps else []) + (["sbruce trees"] if tree_writer is not None else [])
    print(f"=== {detector}: {' + '.join(products)}, binning {binning}, writing to {os.path.abspath(outdir)} ===")
    print(out.bins)

    cv_files = _resolve_files(df_dir, cfg["cv"])
    main_cv = load_sample(cv_files, detector)  # (df, match, pot), unscaled

    # 1. Dedicated-sample variations, each matched to its own CV copy
    for match_set in cfg["match_sets"]:
        _run_match_set(match_set, main_cv, detector, df_dir, goal_pot, out, max_workers=max_workers)

    # 2. Everything else is compared to the full, unmatched CV
    cv_df, _, cv_pot = main_cv
    del main_cv
    loaddf.scale_pot(cv_df, cv_pot, goal_pot)

    def unmatched_variations():
        def _load_bind():
            df, _, pot = load_sample(cv_files, detector, shift_binding_E=True)
            loaddf.scale_pot(df, pot, goal_pot)
            return "BIND", df, 1.0
    
        yield _load_bind()

    def unmatched_variations():
        # Pipeline: each file is prefetched in a single background thread while
        # the previous result is being yielded and processed. Only one HDF5 file
        # is ever open at a time, avoiding thread-safety issues.
        def _load_bind():
            df, _, pot = load_sample(cv_files, detector, shift_binding_E=True)
            print(df.del_Tp)
            loaddf.scale_pot(df, pot, goal_pot)
            return "BIND", df, 1.0

        def _load_trksplt(region):
            df, _, pot = load_sample(cv_files, detector, split_tracks=region)
            loaddf.scale_pot(df, pot, goal_pot)
            return f"{region}_TRKSPLT", df, 1.0

        disk_tasks = [_load_bind] + [
            (lambda r=region: _load_trksplt(r)) for region in cfg["trksplt_regions"]
        ]
        with ThreadPoolExecutor(max_workers=1) as pool:
            future = pool.submit(disk_tasks[0])
            for i, task in enumerate(disk_tasks):
                name, df, scale = future.result()
                if i + 1 < len(disk_tasks):
                    future = pool.submit(disk_tasks[i + 1])
                yield name, df, scale

        for name, make_variation in CV_DERIVED_VARIATIONS:
            yield name, make_variation(cv_df), cfg["ratio_scales"].get(name, 1.0)

    out.compare(cv_df, unmatched_variations(), UNMATCHED_CV_TREE_NAME)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def resolve_function(func_string):
    """
    Resolve strings like 'gmpl.all_gump_cuts', 'all_gump_cuts', or
    'my_cuts_module.custom_cut' into callable Python functions.
    """
    if "." not in func_string:
        if hasattr(gmpl, func_string):
            return getattr(gmpl, func_string)
        raise argparse.ArgumentTypeError(
            f"Function name must be 'module.function' (e.g. 'gmpl.all_gump_cuts'), got '{func_string}'")

    module_name, func_name = func_string.rsplit(".", 1)
    if module_name == "gmpl":
        module_name = "analysis_village.gumple.gumple_cuts"

    try:
        return getattr(importlib.import_module(module_name), func_name)
    except (ImportError, AttributeError) as e:
        raise argparse.ArgumentTypeError(f"Could not import '{func_string}': {e}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Build detector-systematic reweight maps and/or sbruce trees.")
    parser.add_argument("-s", "--selection", type=resolve_function, default=gmpl.all_gump_cuts,
                        help="Selection function (e.g. 'gmpl.all_gump_cuts' or 'gmpl.all_maplemp_cuts')")
    parser.add_argument("-o", "--outdir", type=str, default=DEFAULT_OUTDIR,
                        help="Output directory (relative to the current directory) for maps.")
    parser.add_argument("--treeoutdir", type=str, default=DEFAULT_OUTDIR,
                        help="Output directory (relative to the current directory) for trees. "
                             f"sbruce trees go in its '{SBRUCE_SUBDIR}' subdirectory")
    parser.add_argument("-d", "--dfdir", type=str, default=DEFAULT_DFDIR,
                        help="Input directory containing the dataframe (.df) files")
    parser.add_argument("-b", "--binning", type=str, default="2D", choices=list(BINNINGS),
                        help="Reweight map binning")
    parser.add_argument("-D", "--detectors", nargs="+", default=list(DETECTOR_CONFIG),
                        choices=list(DETECTOR_CONFIG),
                        help="Detectors to process (quote names with spaces, e.g. 'ICARUS Run2')")
    parser.add_argument("--no-maps", dest="maps", action="store_false",
                        help="Don't write reweight maps (e.g. to make only sbruce trees)")
    parser.add_argument("-p", "--plot", action="store_true",
                        help="Also make debug plots (2D map + del_p slices) for every map")
    parser.add_argument("--plotdir", type=str, default=DEFAULT_PLOTDIR,
                        help="Directory for debug plots, used with --plot (relative to the current directory)")
    parser.add_argument("-t", "--sbruce-trees", action="store_true",
                        help="Also write sbruce trees of the selected events for every variation and its CV")
    parser.add_argument("-v", "--tree-vars", type=str, nargs="+", default=None,
                        help="Variables in each sbruce tree, in addition to --weight-col "
                             "(default: every column of the dataframe)")
    parser.add_argument("-w", "--weight-col", type=str, default=DEFAULT_WEIGHT_COL,
                        help="Event weight column saved in the sbruce trees")
    parser.add_argument("-j", "--jobs", type=int, default=None,
                        help="Number of parallel worker threads (default: one per CPU)")
    args = parser.parse_args()

    if not args.maps and not args.sbruce_trees:
        parser.error("--no-maps without --sbruce-trees would produce no output")
    if args.plot and not args.maps:
        parser.error("--plot makes plots of the maps, so it can't be combined with --no-maps")

    tree_writer = (SbruceTreeWriter(args.treeoutdir, args.tree_vars, args.weight_col)
                   if args.sbruce_trees else None)

    with ProcessPoolExecutor(max_workers=len(args.detectors)) as pool:
        futures = [pool.submit(remake_detvar_maps, det, args.dfdir, selection=args.selection,
                               outdir=args.outdir, treeoutdir=args.treeoutdir,
                               binning=args.binning, plotdir=args.plotdir if args.plot else None,
                               write_maps=args.maps, tree_writer=tree_writer,
                               max_workers=args.jobs) for det in args.detectors]
        for f in as_completed(futures):
            f.result()
