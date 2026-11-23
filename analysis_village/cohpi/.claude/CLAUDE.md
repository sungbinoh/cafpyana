# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Context

This directory is the SBND charged-current coherent pion (CCCOH, ν_μ Ar → μ π Ar) analysis inside the `cafpyana` repo (repo root is `../..`). cafpyana converts SBN flat CAF ROOT files into pandas DataFrames saved as HDF5 (`.df`); everything here is either a DataFrame maker plugged into that framework or a Jupyter notebook that consumes its output. There is no build, lint, or test suite.

## Running

All commands run from the repo root (`../..`) after `source setup.sh` (sets up spack ROOT/xrootd/hdf5, the Python venv in `envs/`, `CAFPYANA_DIR`, and grid scratch dirs on sbnd/icarus/jupyter hosts).

```bash
# Local run on a few files
python run_df_maker.py -c analysis_village/cohpi/configs/cohpi_slcdf.py -o <out_prefix> -i a.root,b.root   # or -l filelist.txt
# Grid submission (writes under $CAFPYANA_GRID_OUT_DIR)
python run_df_maker.py -c analysis_village/cohpi/configs/cohpi_all_weight_signal.py -o <out_prefix> -l filelist.txt -ngrid <N>
```
Useful flags: `-nfile N` (only first N files, for quick tests), `-ncpu`, `-split <GB>` (HDF5 split size).

Notebooks do `sys.path.append('../../../')` (unfolding/) or `'../../../../'` (plotting/*) to reach the repo root, so run them from their own directory. They use `%autoreload 2`, so edits to `makedf`/`pyanalib`/`analysis_village/*` modules are picked up live. Input `.df` paths are hardcoded in the notebooks (often under `/data/sungbino/...`).

## Architecture

**Config → maker pattern.** Each file in `configs/` defines `DFS` (list of functions taking an uproot file `f`) and `NAMES` (HDF5 keys, same order). `run_df_maker.py` calls each function per input file and stores results under those keys. Configs map to makers:

| config | maker module | purpose |
|---|---|---|
| `cohpi_slcdf.py` | `make_cohpidf_slc_ana.py` | per-slice selection variables + truth-matched `mcnu` (event-selection plots) |
| `cohpi_all_weight_signal.py` | `make_cohpidf_all_weight_signal.py` | same selection plus **all** systematic weights (GENIE/Flux/G4 multisims) attached to selected slices and to true signal (`mcnu`); used for syst/xsec studies |
| `cohpi_wienersvd.py` | `make_cohpidf_wienersvd_dfs.py` | inputs for Wiener-SVD unfolding |
| `cohpi_finaldf.py` / `_data.py` | `make_cohpidf_v2.py` | older final-df maker (MC / data) |
| `cohpi_ttree.py` | `makettree/cohpi` | TTree output via `run_ttree_maker.py` |

Makers build on framework pieces in the repo root: `makedf/makedf.py` (`make_mcdf`, `make_mchdrdf`, `make_hdrdf`, `make_potdf_bnb`, `make_pandora_df[_calo_update]`), `makedf/getsyst.py` (`get_all_syst_df`, `filter_systs_nuind`), and `pyanalib/pandas_helpers.py` (`multicol_merge`, `multicol_concat`, …). DataFrames use pandas MultiIndex columns padded with `''` (e.g. `('pfp','trk','chi2pid','I2','mu_over_p','')`); new columns must match the existing depth.

**Selection flow (in `make_cohpidf_slc`)**: FV cut on slice vertex (`InFV_nohiyz`, excludes high-y at z>250) → not clear cosmic → exactly 2 track candidates (dist_to_vertex<6, len>4, trackScore>0.5), 0 proton candidates (chi2_muon[_new]>25), 0 showers → direction/containment and kinematic vars (`reco_t` = |t| from μ+π range momenta, opening angle, etc.). Slice-level outputs are reduced with `make_slc_var` (groupby `entry, rec.slc..index`). Slices are truth-matched to `mcnu` via `slc.tmatch.idx` ↔ `rec.mc.nu..index`. `is_recalo = True` switches to recalibrated calorimetry (`chi2_muon_new`).

**Truth categories (`nuint_categ`, in `make_cohpi_nudf`)**: -1 out-of-FV, 1 signal, 0 NC, 2 non-signal CCCOH, 3 CCQE, 4 2p2h, 5 RES, 6 DIS, 8 other. Plot labels/colors for these live in `unfolding/cohpi_topologies.py` (`mode_list`, `mode_labels`, plus a `-3` data-driven non-beam category). Binning and reco/truth column names for unfolded variables (P_μ, cosθ_μ, cosθ_π) are in `unfolding/var_config_cohpi.py`, using `analysis_village/unfolding/variable_configs.VariableConfig`.

**`all_weight_signal` global-state coupling**: `make_cohpi_signal_nudf` (listed first in `DFS`) computes all syst weights once and stores them in the module global `all_systs_weights_df`, which `make_cohpidf_slc` then reuses. The order in `DFS` is therefore load-bearing.

## Gotchas

- **Selection lives in `makedf/cohpi_selection.py`**: FV, truth `Signal` (GENIE-mode based: `iscc & genie_mode == 3`, 1 μ >27 MeV, 1 π± >30 MeV, in FV), pfp candidate definitions, dir-z flags, preselection, and the final cuts (`FINAL_CUTS`, `SPHERE_CUTS`/`APPLY_SPHERE_CUT`, `final_selection_masks` for cut-flows). `make_cohpidf_slc_ana.py` (preselection only) and `make_cohpidf_all_weight_signal.py` (full selection) import it; change cuts there, not in the makers. `make_cohpidf_v2.py` and `make_cohpidf_wienersvd_dfs.py` still carry their own copies with the old 1μ1π0p0n0π0 signal definition. The notebooks still hard-code cut values.
- Kinematic helpers (`reco_t`, `Avg`, opening angle, …) and `make_cohpi_nudf` are still copy-pasted in each maker.
- **Calo variations (all_weight_signal, MC only)**: `calo_variation_chi2` recomputes plane-2 chi2 for every `chi2pid.CALO_VARIATIONS` key (same recipe as `make_pandora_df`'s CV recalo, hits loaded once, no hit FV cut, so `cv` == `*_new`) and the maker stores `n_proton_<var>`, `{long,short}_chi2_{mu,pro}_<var>`. On MC the chi2-dependent cuts (`n_proton == 0`, chi2_proton) are **not** applied in the maker; apply offline: `final_selection_mask(evt, True, apply_chi2=False) & calo_var_masks(evt, var)` (nominal: `final_selection_mask(evt, True)`). Data still gets the full selection. `makedf/make_trkdf(updatecalo=...)` is not used (its InFV call is broken) and shared `makedf/` code is intentionally left untouched.
- Notebooks hold the downstream analysis (`plotting/evt_sel`: selection plots, efficiencies; `plotting/syst`: merge grid output → covariances/event rates; `unfolding/run_*_syst.ipynb`: per-syst-source unfolding, chosen by `syst_name`/`output_str` in the second cell). Sample names encode production (e.g. `SBND2026A_gen2 ... v10_14_02_*`).
- Branch policy from the repo README: Gen1 (v10_06) → `dev_sbnd_gen1`, Gen2 (v10_14) → `dev_sbnd_gen2`; `main` tracks the latest SBND samples.
