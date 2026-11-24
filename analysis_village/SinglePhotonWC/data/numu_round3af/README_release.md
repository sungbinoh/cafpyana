# SBND Wire-Cell νμ CC BDT: round 3 af (2026-11-23), the current νμ release

A retrained νμ tagger for SBND: the three small models (numu_1, numu_2, cosmict_10) and the final model, all trained on
SBND simulation whose charge (dQ/dx) inputs were corrected towards data. It replaces round 1 (`../numu_round1/`, final
model only), which stays in the repository for comparison.

Trained on SBND simulation reconstructed with Wire-Cell toolkit 0.38.0 (production WCLocal2026Nov16).

## Performance at 95 % purity

Measured on the SBND evaluation samples at Fermilab, which were not used in training. Their simulation was corrected in
the same way as the training simulation (below); data are never corrected. Data are beam-on Run 3 plus open Run 1
(9.57e18 POT).

| score | cut on `numu_score` | νμ CC efficiency | range | beam-on data / prediction |
|---|---:|---:|---|---:|
| shipped weights | 1.017 | 75.0 % | 74.7–75.3 | 1.058 ± 0.005 |
| round 1 | 1.011 | 80.1 % | 79.8–80.3 | 1.048 ± 0.005 |
| **round 3 af (this)** | **0.985** | **80.3 %** | 80.0–80.5 | **1.063 ± 0.005** |

- Signal is νμ / anti-νμ CC with the true vertex in the fiducial volume and deposited energy above 100 MeV.
- The efficiency range is where the efficiency moves when the purity curve is shifted by its own statistical error.
- Data / prediction is 40,222 selected against 37,854 predicted, at the level of the shipped weights.
- A version that retrained the small models on uncorrected simulation (round 2) lost about 7 % of the selected data
  relative to prediction. The data and simulation disagree on the muon dQ/dx, so that version's numu_1 cut where data and
  simulation differ. Correcting the simulation removes that loss: in the most muon-like tenth of numu_1 values, data /
  prediction is 0.91–0.97, against 0.11 for round 2.

## The simulation correction (part of the model)

Every prediction made with this model must use simulation corrected as it was trained (`correction.json`, also in
`manifest.json` under `simulation_correction`):

- **Charge:** q → 1.071 q on the seven dQ/dx inputs (MIP units): `numu_cc_1_medium_dQ_dx` (numu_1 segments) and
  `cosmict_{2,3,7}_dQ_dx_{end,front}`. Simulation only, filled values above 0 only. The factor was fitted by gpvm to
  data / prediction on these seven inputs.
- **Muon flag:** per numu_1 segment, `numu_cc_flag_1` stays on only if the corrected dQ/dx is still below the segment's
  `numu_cc_1_dQ_dx_cut`. The candidate's `numu_cc_flag` is recomputed from the segment flags, the numu_cc_2 flags and
  `numu_cc_flag_3`. About 1 % of flagged segments turn off.

`training/numu/make_corrected_tables.py MERGED OUT correction.json --samples <simulated samples> --update-flag1` writes a
corrected copy of the table builder's tables.

**Inside the chain, simulation must be corrected the same way.** The toolkit computes the inputs from the reconstructed
charge. Without the correction the number of selected BNB simulation candidates hardly changes (+0.09 %: 171,651 against
171,494), but 9,049 of 302,172 candidates move across the cut. The patch below applies the correction inside the chain, on
simulation only (see "Use it in Wire-Cell"). Without the patch, score simulation from corrected tables, as the evaluation
did. Data is never corrected and can be scored directly by the chain.

## Use it in Wire-Cell

`wirecell/uboone/weights/` holds four files under the shipped names, so no code or configuration change is needed. The SBND
configuration names them in `cfg/pgrapher/experiment/sbnd/clus.jsonnet`:

| file | model | trees | sha256 |
|---|---|---:|---|
| `numu_tagger1.weights.xml` | numu_1 | 692 | in `wirecell/chain_weights.json` |
| `numu_tagger2.weights.xml` | numu_2 | 115 | 〃 |
| `cos_tagger_10.weights.xml` | cosmict_10 | 28 | 〃 |
| `numu_scalars_scores_0923.xml` | final | 472 | 〃 |

`numu_tagger3.weights.xml` is not retrained; it still comes from `wire-cell-data`. Put `releases/numu_round3af/wirecell`
in front of `wire-cell-data` on `WIRECELL_PATH`. Every other weight file still comes from `wire-cell-data`.

**Select νμ CC with `numu_score > 0.985`** for 95 % purity. `numu_score` is on the shipped scale.

**Simulation: apply the charge correction in the chain.** Toolkit 0.38.0 needs `wirecell/patches/wct-0.38.0-numu_sim_dqdx_scale.patch`
(`git am` in wire-cell-toolkit, then rebuild):
- `UbooneNumuBDTScorer` gets the setting `numu_sim_dqdx_scale` (default 1 = off, output unchanged).
- It scales the seven dQ/dx inputs by the factor in the copies handed to the νμ BDTs, and switches `numu_cc_flag` off as
  `make_corrected_tables.py --update-flag1` does. The stored tagger values stay raw.
- `wct-pr.jsonnet` takes it as a TLA, and `sbnd/clus.jsonnet` passes it only when `reality == 'sim'`. Run simulation with
  `--tla-str reality=sim --tla-code numu_sim_dqdx_scale=1.071`, and data with `reality=data` (always 1).
- For the Aurora production workflow (`sbnd-wirecell-production`), `wirecell/patches/step2-numu_sim_dqdx_scale.patch` turns
  the environment variable `NUMU_SIM_DQDX_SCALE` into that TLA and passes it into the container.
- The job log shows `UbooneNumuBDTScorer: simulation dQ/dx scale 1.071 ...` on simulation.

How the files are written (`training/numu/make_numu_chain_weights.py`):
- **Small models:** each file holds the model's trees with the scale and offset of `manifest.json`. The chain then stores
  `numu_1_score` = float32(tanh(scale × max margin + offset)), the value the final model was trained on (min for
  cosmict_10). The toolkit's constants for a candidate with no object (−0.4, −0.1, 0.7) lie below every value a real object
  can give.
- **Final model:** written unscaled, so `numu_score` = 0.868589 × margin, as round 1.

Checks (`training/numu/check_numu_chain_xml.py`, results in `checks/`): every file was read back with an independent
TMVA evaluator that follows the chain's code (`UbooneNumuBDTScorer.cxx`), against the xgboost models, on every candidate
of the corrected BNB simulation, the NC sample and the beam-off training data (328,698 candidates).

| check | result |
|---|---|
| small-model values | max \|diff\| 5.5e-07 |
| final model, same small-model values | max \|diff\| 7.7e-06 on `numu_score` |
| full path, pass/fail at the cut | 4 of 328,698 candidates differ (all in BNB simulation; 171,494 pass) |

A last-bit float difference in a small-model value can cross a split of the final model. That moves about 1 % of
candidates by up to 0.3 in `numu_score`, but almost never across the cut.

**Full-chain test: passed (2026-11-23).** The chain was run with these files and the patch. Production step 1 was reused
and step 2 rerun with the production settings. Its stored scores were then compared with the release scored from the
builder's tables of the same output (`training/numu/compare_chain_vs_table.py`; simulation tables corrected first). Results
in `checks/chain/`:

| test | candidates | max \|diff\| `numu_score` | pass at 0.985: chain / table / differ |
|---|---:|---:|---|
| BNB simulation, ×1.071 in the chain (Aurora, 2 files) | 82 | 4.5e-06 | 44 / 44 / 0 |
| BNB simulation, ×1.071 in the chain (Fermilab, 2 files), against the evaluation's own scores | 31 | 6.0e-06 | 25 / 25 / 0 |
| beam-on data, no correction (Fermilab, 1 file) | 34 | 2.4e-06 (one candidate 0.083, far below the cut) | 22 / 22 / 0 |

- The small-model scores agree to within 1.8e-07.
- **The correction really acts:** against uncorrected tables, the same simulation output differs on 69 of 82 and 28 of 31
  candidates (up to 2.3 in `numu_score`).
- A candidate whose `numu_cc_flag` the correction switches off is reproduced as well.
- **Every other chain output is unchanged:** against the production files (shipped weights), all other candidate columns
  and all small-model input tables are identical, including the raw dQ/dx values.
- The one larger difference on data is the last-bit effect described above.

Not done: the toolkit's own standalone evaluators on table rows (wcgpu1). The full chain covers what they would test.

## Use it outside Wire-Cell

`l1/*.json` and `l2/numu_l2.json` are plain xgboost models (xgboost 3.4.1).

1. For each small model in `manifest.json` (`level1.models`):
   - score each object (inputs in the listed order, float32);
   - apply the protections listed under `protect` (numu_1: an infinite `numu_cc_1_dQ_dx_cut` becomes 10; cosmict_10: a
     missing `cosmict_10_angle_beam` becomes 0);
   - take the max margin (numu_1, numu_2) or the min (cosmict_10);
   - the value is float32(tanh(scale × margin + offset)), or the empty constant when there is no object;
   - it replaces `<model>_score`.
2. Final model: the 74 inputs in the order of `manifest.json`. Set the four protected inputs to 0 where they are missing.
   The margin is `predict(..., output_margin=True)`, `numu_score` = 0.868589 × margin, and the cut is margin 1.134.

From the table builder's tables (correct simulation first):

    python training/numu/make_corrected_tables.py MERGED CORRECTED releases/numu_round3af/correction.json \
        --samples wcl1116-bnb wcl1116-ncdr --update-flag1
    python training/numu/score_numu_release.py releases/numu_round3af CORRECTED SAMPLE OUT.parquet   # column d6_f0 = margin

## How it was trained

`training/numu/train_numu_l1.py` (small models), `training/numu/train_numu_l2.py` (final model) and
`training/numu/freeze_numu_l2.py`, all on the corrected tables.

- **Samples:** BNB simulation; NC events from the NC-enhanced (ncdr) sample as background only; the beam-off data for the
  cosmic background (not corrected).
- **Weights:** each sample is weighted to the same exposure (BNB 1, ncdr 0.010571, beam-off 4.91257).
- **Small models:** rule A labels (MicroBooNE's), as in round 2. Their outputs are mapped into the chain's range with
  wcgpu's `affine_for_range` recipe, and the final model is trained on those values.
- **Final model:** binary logistic xgboost, depth 6, learning rate 0.1, subsample 0.8, column sampling 0.5 per level and
  per node, gamma 0.2. Early stopping on a held-out development slice kept 472 trees.
- **Development slice:** Aurora's held-out simulation, corrected. Efficiency at 95 % purity is 86.2 ± 0.3 % here, against
  77.8 % for the shipped weights.

Variants trained beside this one (Fermilab evaluation, efficiency and data / prediction at 95 % purity) were all within
0.3 points of each other:

| variant | correction | numu_1 muon flag | efficiency | data / prediction |
|---|---|---|---:|---:|
| 3a | q → 1.071 q | not updated | 80.35 % | 1.061 |
| **3af (this)** | q → 1.071 q | updated | 80.33 % | 1.063 |
| 3b | a straight line in q | not updated | 80.38 % | 1.060 |
| 3bf | a straight line in q | updated | 80.63 % | 1.060 |

The line in 3b and 3bf is q → q (1.100 − 0.0245 clip(q, 0.5, 4)).

## Known limitations

- The data are 10–15 % above prediction just above the cut and 10–20 % below it at the highest scores. Every model,
  including the shipped weights, shows this.
- The single factor of 1.071 is an average. Data / simulation differs by angle and by dQ/dx (7–9 %). Variants with a line
  in dQ/dx (3b, 3bf) gave the same result within 0.3 points.
