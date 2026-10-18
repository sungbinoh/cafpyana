# NCDelta selection: local optimization workspace

This directory is a standalone copy of the SBND NC Δ→Nγ (single photon)
selection from cafpyana (`analysis_village/NCDelta`). It needs only numpy,
pandas, pytables and matplotlib (`pip install -r requirements.txt`).

The goal is to improve the selection's signal/background separation by
changing the cuts in `cuts.py`.

## Files

| File | Role |
|---|---|
| `cuts.py` | All selection definitions: topology, pre-selection, cut-based stage, energy window, truth categories. **The file being optimized.** |
| `dfio.py` | Reads `.df` files (HDF5, keys `slc_<i>`, `mcnu_<i>`, `hdr_<i>`). |
| `evaluate.py` | Efficiency, purity, S/√B, S/√(S+B) per channel. The objective. |
| `plot_composition.py` | Stacked truth-composition plots, including N−1 plots of every cut. |
| `fetch_dataframes.sh` | Copies the dataframes from the Fermilab GPVM. |
| `data/` | The dataframes (not in git). |

## Running

```bash
python evaluate.py -i data/<dir>/ --pot 1e20                      # current cuts
python evaluate.py -i data/<dir>/ --pot 1e20 --config cuts.json   # with overrides
python evaluate.py -i data/<dir>/ --pot 1e20 --split even         # half the files
python plot_composition.py -i data/<dir>/ --pot 1e20 -o plots/
```

`--config` is a JSON file of `CUTBASED` overrides per channel, e.g.
`{"1g1p": {"trk_proton_like": false, "nu_score_min": 0.4}}`. `null` disables
a cut. For tuning loops, call `evaluate.evaluate(slc, nu, scale, overrides)`
directly from Python instead of the CLI: loading the data is the slow part.

## Physics

- **Signal:** NC Δ(1232) → Nγ with the true vertex in the active volume
  (`cuts.is_signal`, category 0 of `cuts.truth_category`).
- **Channels:** `1g1p` is one shower and one track (the proton); `1g0p` is
  one shower and no tracks. Each slice is in at most one channel.
- **Dominant background:** NC π⁰ where only one photon is reconstructed.
  The other photon escapes, is too soft, merges with the first, or converts in a dead region.
- **Units:** energies in MeV, lengths in cm.
- One row per slice. `shw_*` columns describe the leading primary shower and
  `trk_*` the leading primary track; NaN means none. `n_shw`, `n_trk` count
  primary showers and tracks; `n_nonprim_shw` counts non-primary showers.

## Rules

1. **Never cut on truth.** Columns starting with `true_`, plus `shw_true_*`,
   `trk_true_*` and `tmatch_*`, are labels for evaluation only. Cuts may use
   only reconstructed columns.
2. **Pre-selection is already applied** to these dataframes (unless they were
   made with the `nopresel` config). Pre-selection cuts can be tightened here,
   but loosening them has no effect and gives misleading results.
3. **Don't change the truth definitions** (`is_signal`, `is_ncdelta_rad`,
   `truth_category`, `CATEGORY_LABELS`) to improve the numbers.
4. **Guard against overfitting.** The signal sample is small. Tune on
   `--split even` and report results on `--split odd`, or vice versa, and
   always show the raw `mc_signal` / `mc_bkg` counts next to the scaled
   numbers. A gain smaller than the MC statistical uncertainty
   (roughly 1/√mc_signal) is not a gain.
5. **Keep the public interface of `cuts.py`.** The cafpyana dataframe maker,
   which is not in this directory, imports `TRK_SCORE_CUT`, `topology_cut`,
   `preselection` and `dist_to_active_boundary`. The plotting script relies on
   `stage_cuts`, `nminus1_masks`, `PRESEL_*`, `CUTBASED_*` and `ENERGY_WINDOW`.
   New cuts go into `cutbased_cuts` with a config key in `CUTBASED_*`. Add a
   matching entry in `plot_composition.cut_var_specs` so the cut gets an N−1 plot.
6. The fiducial volume in `cuts.py` is a copy of the one in cafpyana's
   `mcs_pid/cuts.py`; don't change it here.
7. Report the figure of merit you optimized, the cut values, and the
   efficiency and purity of each channel, before and after.

`cuts.py` is copied back into cafpyana unchanged, so it must stay importable
with only numpy and pandas.
