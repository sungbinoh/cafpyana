# CAFANA-MAPLE vs cafpyana-maple: differences and validation

## Validation summary (2026-08-11)

Both frameworks were run on the same 94 ICARUS Run4 ReCAF2026 overlay flat
CAFs (`/pnfs/sbn/scratch/users/twester/Run4_ReCAF2026/flatcaf/000071/000710/`,
~4,700 events; list in `validation_files.list`):

- CAFANA-MAPLE: `cafana-maple/run_cafana_maple.sh` (SL7 container, sbnana
  v10_01_04, the unmodified selection/variable headers) → `SB_GUMP_run4val.root`
- cafpyana-maple: `run_df_maker.py -c analysis_village/maple/configs/maple_evt_cafanapid.py`
  → `maple_run4val_cafanapid.df`
- Comparison: `compare_cafana_cafpyana.py` (matches on Run/Subrun/Evt/Slice_index)

Result — with the CAFANA-compat PID mode, the two frameworks are **identical**:

| Quantity | CAFANA | cafpyana | mismatches |
|---|---|---|---|
| selected slices (`selectedReco`) | 40 | 40 | 0 |
| truth 1muNp interactions (`selectedNu`) | 58 | 58 | 0 |
| all 31 `selectedReco` variables | | | 0 (max diff 3×10⁻⁶) |
| all 8 `selectedNu` variables (incl. `Pass_cut`, efficiency rasters) | | | 0 (max diff 3×10⁻⁷) |

The same comparison was repeated on 196 ICARUS **Run2** ReCAF2026 overlay
files (`Run2_ReCAF2026/flatcaf/000000/000000/`, ~10k events), covering the
Run2 fiducial-volume/run-period code path: 78/78 selected slices, 101/101
truth 1muNp interactions, all variables identical (max diff 6×10⁻⁶), zero
mismatches.

The residual ≲3×10⁻⁶ differences are pure floating-point rounding: CAFANA
computes in C++ double/float mixtures, cafpyana in numpy float64, both from
the same float32 CAF branches.

## Framework-inherent difference: how PID chi2 is computed

This is the one intentional difference between the two frameworks, and the
reason the evt dataframe carries **two chi2 flavors**:

1. **`*_cafana` (CAFANA-compat mode)** — replicates `chi2_ALG` from the
   MAPLE helper header exactly: it uses the **CAF-stored** `trk.calo[2].points.dedx`
   (i.e. whatever calibration was applied at CAF production time), the
   dE/dx-vs-residual-range templates from `dEdxrestemplates.root`, hit
   selection `0 < rr < 25` cm and `0.5 ≤ dE/dx ≤ 100` MeV/cm (first/last
   point dropped), resolution `σ = (0.04231 + 0.0001783·dEdx²)·dEdx`, and
   χ²/npt normalization. Implemented in `chi2pid_cafana.py`.
   The template ROOT files used by the two frameworks
   (`/exp/icarus/app/users/marterop/dev_areas/dEdxrestemplates.root` and the
   cvmfs `larsoft_data v1_02_02` copy) are **byte-identical** (same md5:
   `f06e5a0bc60dae2eae0dbb615a70702f`).

2. **`*_cafpyana` ("cafpyana" mode, physics default)** — the PID used by the
   GUMP analysis: dE/dx is **recomputed from the hit-level dQ/dx**
   (`integral/pitch`) with ICARUS gains, YZ / electron-lifetime / TPC-scale
   calibration databases, and ellipsoidal-modified-box recombination
   (`makedf/chi2pid.py: dedx(gain="ICARUS", calibrate="ICARUS")`), then χ²
   against the same templates. The hit selection also differs slightly from
   chi2_ALG: `rr < 25`, first/last point dropped, `dE/dx < 1000` (no 0.5
   MeV/cm lower bound, no upper 100 MeV/cm bound).

### Size of the difference

On the validation sample (cutting at the same GUMP-tuned working points
χ²µ ≤ 111 / χ²p ≥ 74 / χ²p-π boundary 92):

- muon-candidate count per slice: 8612 (cafana-PID) vs 8684 (cafpyana-PID),
  8606 in common;
- full 1muNp selection: 40 (cafana-PID) vs 40 (cafpyana-PID), **38 in common
  (2 lost + 2 gained, ~5% migration)**;
- on muon candidates, `χ²µ(cafpyana) − χ²µ(cafana)` has mean ≈ +1.2, rms ≈ 3.8.

### Workarounds / recommendations

- For an exact CAFANA replica, run the `maple_evt_cafanapid.py` config
  (or `pid_mode="cafana"`): validated identical above.
- For physics with the gump calibration treatment (and the calorimetric
  systematic variations, which are only defined for the recomputed dE/dx),
  use the default `pid_mode="cafpyana"` configs. If the ~5% selection migration
  matters, the χ² working points (111/74/92) could be re-tuned on the
  cafpyana-recomputed χ² distributions; both flavors are stored in every evt
  dataframe (`*_cafpyana` / `*_cafana` columns, plus `*_alt` cut booleans) so
  the migration can be measured in any production.

## Other implementation notes (no observable difference, but worth knowing)

- **NaN semantics.** All C++ cuts are written as "skip if outside", so NaN
  comparisons pass or fail in a specific direction. The pandas masks mirror
  this exactly (e.g. a NaN trackScore *passes* the `trackScore < 0.5` skip
  in `find_muon`; a NaN χ² falls through `id_pfp` to the shower block, not
  to Unknown).
- **`id_pfp` trackScore.** The trackScore branches in the helper's `id_pfp`
  are commented out in the CAFANA source: *every* PFP (track- or
  shower-like) goes through the χ²-proton decision first and falls through
  to the shower-energy block otherwise. The port reproduces this.
- **`std::min` with NaN.** The vertex-association distance
  `min(|start−vtx|, |end−vtx|)` keeps the C++ `std::min` NaN behavior
  (second argument NaN → first argument; first NaN → NaN).
- **MC vs data time windows.** `bar_flash` uses (0, 1.6) µs for MC and
  (−0.4, 1.5) µs for data, keyed off `rec.hdr.ismc`; the cryo-light window
  is (−0.6, 1.8) µs with PE > 3000/0.341 in both. Validation used MC.
- **Cutflow ordering.** `MaxCutPassed` (the `maxcut` column) uses the
  *cutflow* function's order (…FV → CRT veto → cryo-light → containment…),
  which differs from the boolean selection's evaluation order (FV →
  cryo-light → CRT veto → containment). Same selected set, different
  intermediate numbering; the port follows the cutflow ordering, validated
  against `Pass_cut`.
- **Efficiency rasters.** `Eff_raster_angle/nuscore/cryosel` record values
  from the last slice that *improves* the running max-cut (a loop-order
  dependent definition). Reproduced in `compare_cafana_cafpyana.py::rasters_for_nu`
  and validated; the same recipe should be used in notebooks.
- **`find_longest_proton`** requires proton length > 0 strictly (C++ starts
  the max at 0.0).
- **Smearing variations** (`smear5/13`, `sqsmear15`) use `np.random` and are
  not reproducible run-to-run (same as gump).
- **wgt psets.** The ReCAF2026 files do not carry every GENIE pset in
  cafpyana's default list; `make_maple_wgtdf` restricts to the psets present
  in each file (103 weight groups available in ReCAF2026).

## Datasets not validated

- Data (`ismc == 0`) samples: the data bar_flash/PE windows are implemented
  but the unblinded data lives on /pnfs/icarus or SAM, unreachable from the
  SBND gpvms where this validation ran.
- The `/pnfs/icarus/persistent/users/cfarnese/for_Nicola` Run4 concat files
  (the macro's currently-active dataset): same reachability caveat. The
  Run4_ReCAF2026 sample used here is the same production those files were
  concatenated from.
