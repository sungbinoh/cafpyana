# HNL nu_ee / nu_pi0 dataframe production

Two-step workflow for every sample: **submit grid jobs** to produce split
`.df` shards, then **concatenate** those shards into a single dataframe file
once the jobs finish.

## Cuts already baked in

Every config below calls `make_hnldf()` (in `makedf/make_hnldf.py`) with its
default `applyPreselection=True, savePfp=False` — none of the wrapper
functions (`make_hnldf_mevprtl`/`_wgt`, `make_hnldf_mcnu`/`_wgt`,
`make_hnldf_data`) override either default. That means **every sample
produced via this README's commands already has the following cuts applied
at production time**, before any downstream analysis-level selection:

- **Slice preselection** (`applyPreselection=True`):
  - Not clear cosmic (`slc.is_clear_cosmic == 0`)
  - Neutrino score > 0.5 (`slc.nu_score > 0.5`)
  - Standard-box fiducial volume (`InFV(..., det="SBND", inzback=0)`)
- **PFP sanity filter** (`savePfp=False`): drops PFPs with
  `trackScore <= 0`, `shw.bestplane_energy <= 0`, an invalid shower start
  position (`shw.start.{x,y,z} == -999`), `trk.len <= 0`, or
  `shw.bestplane_dEdx <= 0`.
- **Topology cuts** (`savePfp=False`, PFPs split into tracks/showers by
  `trackScore` vs. the `trackScore=0.51` threshold):
  - Exactly 0 tracks
  - Exactly 1 or 2 showers
  - The highest- and second-highest-energy showers become `primshw`/`secshw`
    (not a cut, just which PFP gets picked for each role)

Pass `applyPreselection=False`/`savePfp=True` explicitly (see
`make_hnldf_mcnu_nopreselect_savepfp`/`make_hnldf_mevprtl_nopreselect_savepfp`/
`make_hnldf_data_nopreselect_savepfp` in `makedf/make_hnldf.py`, not wired
into `submit_jobs.sh` by default) if you need the full, unselected
event/PFP-level dataframe instead.

## Mechanics

**Step 1 — submit.** Run from the `cafpyana` repo root, via `./submit_jobs.sh`
(uncomment the lines you need) or `run_df_maker.py` directly:

```
python run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/<config>.py -l <input.list> -ngrid <N> -o <output_name>
```

`-ngrid <N>` submits to the grid (`N` jobs); omit it to run locally instead.
Output shards land in `$CAFPYANA_GRID_OUT_DIR/dfs/<timestamp>__<output_name>/`,
one `<output_name>_<i>.df` file per job — the `<timestamp>` is stamped at
submission time (check the job's own printed output, or
`ls $CAFPYANA_GRID_OUT_DIR/dfs/ | grep <output_name>` to find it). Check
`condor_q`/that directory to confirm all jobs finished before moving to
Step 2.

**Step 2 — concatenate.** Combine every shard for a sample into one dataframe
file with `scripts/concat_hdf.py`:

```
python analysis_village/hnl_nuee_nupi0/scripts/concat_hdf.py \
    $CAFPYANA_GRID_OUT_DIR/dfs/<timestamp>__<output_name> \
    -o <final_output_dir>
```

- Positional `directory` is the per-job shard directory from Step 1.
- Omit the `keys` argument to concatenate every key found in the shards;
  pass specific keys (e.g. `nu slc`) to concatenate only those.
- `-o/--output-dir` defaults to the input directory if omitted.
- Each output row gets a `file_idx` index level identifying which source
  shard it came from (needed downstream to disambiguate `__ntuple`/`entry`
  numbering, which each shard resets independently).
- The result is itself split-formatted (`<key>_<i>` + a `split` key), sized
  to `--split`/`-s` GB per split (default 1.0 GB) — readable the same way as
  any other maker output.

## Config reference

| Sample | Config |
|---|---|
| MC MeVPrtl (HNL), no weights | `configs/hnl_mcmevprtl.py` |
| MC MeVPrtl (HNL), with BNB flux + G4 weights | `configs/hnl_mcmevprtl_wgt.py` |
| MC SM neutrinos, no weights | `configs/hnl_mcnu.py` |
| MC SM neutrinos, with GENIE/BNB flux/G4 weights | `configs/hnl_mcnu_wgt.py` |
| MC BNB detector variations (CV + light-yield DV) | `configs/hnl_mcnu_detvar.py` |
| Data (BNB on-beam / off-beam) | `configs/hnl_data.py` |

## Per-sample commands

Every submit command below is in `submit_jobs.sh` (uncomment to run); the
concat command underneath it is the matching Step 2 — just fill in the real
`<timestamp>` from Step 1's output directory. `$OUT` is shorthand for
`$CAFPYANA_GRID_OUT_DIR/dfs`.

### MC BNB nominal (full weighted dataframe, main analysis)

```bash
# submit
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcnu_wgt.py -l ./data/sample_lists/sbnd/mcbnb_cv.list -ngrid 2334 -o mcbnb_cv
# concat
python analysis_village/hnl_nuee_nupi0/scripts/concat_hdf.py $OUT/<timestamp>__mcbnb_cv -o <final_output_dir>
```

### MC BNB DetVar (CV + light-yield DV variations)

```bash
# submit
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcnu_detvar.py -l ./data/sample_lists/sbnd/mcbnb_cv.list -ngrid 937 -o detvar_cv_0
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcnu_detvar.py -l ./data/sample_lists/sbnd/mcbnb_0p94xLY.list -ngrid 312 -o detvar_0p94xly_0
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcnu_detvar.py -l ./data/sample_lists/sbnd/mcbnb_1p19xLY.list -ngrid 312 -o detvar_1p19xly_0
# concat (one per variation)
python analysis_village/hnl_nuee_nupi0/scripts/concat_hdf.py $OUT/<timestamp>__detvar_cv_0 -o <final_output_dir>
python analysis_village/hnl_nuee_nupi0/scripts/concat_hdf.py $OUT/<timestamp>__detvar_0p94xly_0 -o <final_output_dir>
python analysis_village/hnl_nuee_nupi0/scripts/concat_hdf.py $OUT/<timestamp>__detvar_1p19xly_0 -o <final_output_dir>
```

### MC HNL: nu_ee channel

`ngrid` retuned to ~709 files/job (see `submit_jobs.sh` header) — none of
these five submitted yet as of this writing.

```bash
# submit
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nuee_m35.list -ngrid 7 -o mchnl_nuee_m35
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nuee_m50.list -ngrid 7 -o mchnl_nuee_m50
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nuee_m75.list -ngrid 7 -o mchnl_nuee_m75
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nuee_m100.list -ngrid 7 -o mchnl_nuee_m100
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nuee_m125.list -ngrid 7 -o mchnl_nuee_m125
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nuee_m140.list -ngrid 7 -o mchnl_nuee_m140
# concat (one per mass point)
python analysis_village/hnl_nuee_nupi0/scripts/concat_hdf.py $OUT/<timestamp>__mchnl_nuee_m35 -o <final_output_dir>
python analysis_village/hnl_nuee_nupi0/scripts/concat_hdf.py $OUT/<timestamp>__mchnl_nuee_m50 -o <final_output_dir>
python analysis_village/hnl_nuee_nupi0/scripts/concat_hdf.py $OUT/<timestamp>__mchnl_nuee_m75 -o <final_output_dir>
python analysis_village/hnl_nuee_nupi0/scripts/concat_hdf.py $OUT/<timestamp>__mchnl_nuee_m100 -o <final_output_dir>
python analysis_village/hnl_nuee_nupi0/scripts/concat_hdf.py $OUT/<timestamp>__mchnl_nuee_m125 -o <final_output_dir>
python analysis_village/hnl_nuee_nupi0/scripts/concat_hdf.py $OUT/<timestamp>__mchnl_nuee_m140 -o <final_output_dir>
```

### MC HNL: nu_pi0 channel

```bash
# submit
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nupi0_m140.list -ngrid 7 -o mchnl_nupi0_m140
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nupi0_m165.list -ngrid 7 -o mchnl_nupi0_m165
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nupi0_m190.list -ngrid 7 -o mchnl_nupi0_m190
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nupi0_m215.list -ngrid 7 -o mchnl_nupi0_m215
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nupi0_m240.list -ngrid 7 -o mchnl_nupi0_m240
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_mcmevprtl_wgt.py -l ./data/sample_lists/sbnd/mchnl_nupi0_m260.list -ngrid 7 -o mchnl_nupi0_m260
# concat (one per mass point)
python analysis_village/hnl_nuee_nupi0/scripts/concat_hdf.py $OUT/<timestamp>__mchnl_nupi0_m140 -o <final_output_dir>
python analysis_village/hnl_nuee_nupi0/scripts/concat_hdf.py $OUT/<timestamp>__mchnl_nupi0_m165 -o <final_output_dir>
python analysis_village/hnl_nuee_nupi0/scripts/concat_hdf.py $OUT/<timestamp>__mchnl_nupi0_m190 -o <final_output_dir>
python analysis_village/hnl_nuee_nupi0/scripts/concat_hdf.py $OUT/<timestamp>__mchnl_nupi0_m215 -o <final_output_dir>
python analysis_village/hnl_nuee_nupi0/scripts/concat_hdf.py $OUT/<timestamp>__mchnl_nupi0_m240 -o <final_output_dir>
python analysis_village/hnl_nuee_nupi0/scripts/concat_hdf.py $OUT/<timestamp>__mchnl_nupi0_m260 -o <final_output_dir>
```

### Data (BNB + offbeam)

```bash
# submit
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_data.py -l ./data/sample_lists/sbnd/dtbnb_fix.list -ngrid 9 -o dtbnb_fix
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_data.py -l ./data/sample_lists/sbnd/dtbnb_roll.list -ngrid 7 -o dtbnb_roll
python3 run_df_maker.py -c ./analysis_village/hnl_nuee_nupi0/configs/hnl_data.py -l ./data/sample_lists/sbnd/dtoff.list -ngrid 216 -o dtoff
# concat
python analysis_village/hnl_nuee_nupi0/scripts/concat_hdf.py $OUT/<timestamp>__dtbnb_fix -o <final_output_dir>
python analysis_village/hnl_nuee_nupi0/scripts/concat_hdf.py $OUT/<timestamp>__dtbnb_roll -o <final_output_dir>
python analysis_village/hnl_nuee_nupi0/scripts/concat_hdf.py $OUT/<timestamp>__dtoff -o <final_output_dir>
```
