# GUMP technote figures from the GUMPLE dataframes

This directory holds everything needed to regenerate the locally reproducible
figures of the GUMP technical note (Overleaf project `6973acd1e6f673f7a8b495e7`,
checked out beside this repo as `cafpyana/6973acd1e6f673f7a8b495e7/`) from a
GUMPLE production (`sbn-rewgted-25` as of 2026-10-24) and install them into the
Overleaf checkout.

| file | role |
|---|---|
| `make_technote_plots.sh` | master driver: runs every analysis step, checks the manifest, installs the figures, compiles the note |
| `figure_manifest.tsv` | one row per regenerated figure: `<plotbase-relative source>` TAB `<technote-relative destination>` |
| `technote_figures.py` | stdlib inventory tool: `list` / `unused` / `check` / `install` (parses `\includegraphics`, expands `\foreach` brace-aware, zero-argument macros such as `\FDFigDir`, and parameterized `\newcommand` figure macros at their call sites) |
| `sync_calibrations.py` | stdlib: compares the constants measured by the calibration steps (`flashfit`, `tracksplit`) with the literals in `../loaddf.py`; `--apply` writes them |
| `figure_placement.tsv` | snapshot of `technote_figures.py list --table` after the last install: every referenced image with its tex file, section, label and caption |
| `README.md` | this file |

The Claude skill `.claude/skills/technote-plots/SKILL.md` (repo root) describes
how to redo the whole exercise for a future iteration of the note.

## Prerequisites

* The GUMPLE dataframes: `--df-dir` must hold the 47 `.df` files of the
  production (`SBNDMCCV_*.df`, `ICARUSRun{2,4}_SpringMCOverlay_rewgt_*.df`, the
  detvar / dirt / beam-off / on-beam files, `SBNDIntimeMC.df`).
* The interpreter `/Users/gputnam-local/Work/fitter/env/bin/python` (+ its
  `jupyter`) -- the only local environment that imports cafpyana. Pass it to
  `run_gumple_serial.py` through `GUMPLE_PY` / `GUMPLE_JUPYTER`
  (`export GUMPLE_PY=/Users/gputnam-local/Work/fitter/env/bin/python
  GUMPLE_JUPYTER=/Users/gputnam-local/Work/fitter/env/bin/jupyter
  GUMPLE_KERNEL=python3`). `GUMPLE_KERNEL` overrides the notebooks' kernelspec
  (`venv_py310_cafpyana`, GPVM only); without it the `pid` steps fail at once
  with `NoSuchKernel`.
* Every analysis runs with cwd = `analysis_village/gump` (the serial runner does
  this for you).
* `latexmk` / `pdflatex` / `bibtex` on `PATH` (MacTeX) for the compile step.
* Disk: loaddf caches its preselected frames as `<hash>.h5` *inside the
  production directory*. A full first pass over a new production writes tens of
  GB there. Check `df -h` first.

## Quick start

```bash
cd analysis_village/gump/technote
./make_technote_plots.sh                       # everything, defaults below
./make_technote_plots.sh --only mcdata --only install --only compile
./make_technote_plots.sh --update-code           # write re-measured calibration constants into loaddf.py
./make_technote_plots.sh --only sync --plotbase <existing plotbase>   # just compare loaddf with the fits
DF_DIR=/path/to/sbn-rewgted-22/ PLOTBASE=$PWD/../../../plots-gumple-XX ./make_technote_plots.sh
```

Defaults (env or flag): `DF_DIR=/Users/gputnam/Work/osc/sbn-rewgted-25/`,
`PLOTBASE=<repo>/plots-gumple-2026-10-24-full`, `FLASHONLY_DF_DIR=<df-dir>-flashonly/`,
`TECHNOTE=<repo>/6973acd1e6f673f7a8b495e7`. Logs go to `<PLOTBASE>/logs/`
(`technote-serial.log`, one `<step>.log` per analysis, `technote-check.log`,
`technote-sync.log`, `technote-install.log`, `technote-unused.txt`, `technote-compile.log`,
`tex/main.log`). Steps are independent: a failure is reported in the summary
table and the chain continues -- except a calibration mismatch (phase 2 below),
which stops the chain before the downstream steps. Nothing is committed to either git repo.

## Phases: calibration steps come first

Three steps *measure* numbers that are hard-coded in `../loaddf.py` and used by
every later step, so the driver runs in phases (`run_gumple_serial.py` runs its
`--only` steps in its own order, so each phase is a separate call):

| phase | steps | measures | loaddf constant | used by |
|---|---|---|---|---|
| 1 calibrate | `flashfit` | all-trigger MC flash-PE scale `s_alltrig`, trigger syst. `trig_unc` = \|all-trigger - Contained\| (`flashpe_fit/trigger_scale_summary.tsv`) | `MC_PE_SCALE`, `MC_PE_SCALE_UNC` | every MC figure (MC flash PE is scaled at load; the flash-PE cut sees it), `syst`/`mcdata`/`sbsyst` `TRIG_PE_SCALE*` (Trig systematic) |
| | `tracksplit` | ICARUS split fraction per plane, Simple stage, Run2+Run4 rows of `tracksplit/tracksplit_summary.csv` | `SPLIT_FRAC` | `mcdata` + `sbsyst` track-split systematic; the `tracksplit` post-split figures themselves |
| | `sbndsplit` | SBND APA-gap split fraction | none (MC already over-splits; no SBND split applied) | report only |
| 2 sync | `sync_calibrations.py` | compares phase 1 with `loaddf.py` (3 decimals for the PE scale, 4 for split fractions) | -- | -- |
| 3 downstream | `mcdata constraint pid sbsyst corsika selection pot stability cvwgt flashpe nminus1` | -- | -- | -- |
| 4 | `check install compile` | -- | -- | -- |

On a phase-2 mismatch the driver **stops** unless `--update-code` (rewrite the
literals in `loaddf.py`, stamped `# sync_calibrations.py <date> from <plotbase>`,
then re-run `tracksplit` so its post-split MC uses the new fraction) or
`--allow-stale` (continue with the old constants) is given. Default is neither:
changing loaddf moves every MC figure, so it is a deliberate choice. The
tables quoting these numbers -- `tab:trigger-systematic` and the
track-splitting fractions in `Sections/Systematic Uncertainties.tex` -- are not
edited by the driver; take them from `flashpe_fit/trigger_scale_table.tex` and
`tracksplit/tracksplit_summary.txt`. `flashpe_fit_gumple.py` fits above
`gumple_cuts.FLASH_PE_THRESHOLD`, so a cut change re-measures `MC_PE_SCALE`.

## Steps

The technote-only steps `sbsyst`, `pot`, `stability`, `cvwgt`, `flashpe`,
`nminus1`, `sbndsplit` and `flashfit` are skipped by `run_gumple_serial.py`
unless named with `--only` (`TECHNOTE_ONLY`); the driver always names its
steps. Rows are in run order (phases above).

All analysis steps are invoked through `run_gumple_serial.py` (which passes
`--df-dir` / `--plotdir` to the scripts and `GUMPLE_DF_DIR` / `GUMPLE_PLOTBASE`
to the notebooks). Runtimes are from the 2026-09-02 run against
`sbn-rewgted-21` with warm loaddf caches.

| step | producer | outputs under `<PLOTBASE>` | technote destination | figures | runtime |
|---|---|---|---|---|---|
| `flashfit` | `flashpe_fit_gumple.py --schema gumple` (on `--df-dir`) and `--schema flashonly` (on `--flashonly-df-dir`, default `<df-dir>-flashonly/`) | `flashpe_fit/{gumple,flashonly}/{png,pdf}/`, `flashpe_fit/trigger_scale_{summary.tsv,table.tex}` | `Figures/Trigger-Systematic/{Contained,AllTriggers}/` | 6 | ~5 min; logs `flashfit-{gumple,flashonly}.log` |
| `tracksplit` | `TrackSplittingCorrection_GUMPLE.py` (all detectors, stages simple+gump, angular Run 2+4, 4 quantiles) | `tracksplit/{png,pdf}/` | `Figures/Track-Splitting/` | 36 | 0.6 min |
| `sbndsplit` | `TrackSplittingCorrection_SBND_GUMPLE.py` (SBND APA gap at z=250, y<90; both TPCs, East, West; simple+gump; angular) | `sbndtracksplit/{png,pdf}/` | `Figures/Track-Splitting/` | 10 (simple stage only) | ~20 min uncached (13 MC shards) |
| `sync` | `technote/sync_calibrations.py` (`--apply` with `--update-code`) | `logs/technote-sync.log` | `../loaddf.py` constants | -- | s |
| `mcdata` | `mcdata_comparison_gumple.py -d all --nproc 8` | `evtsel/{png,pdf}/` | `Figures/MC-Data-Comparisons/<Det>-Area-Norm/` | 43 x 3 | 7.2 min |
| `constraint` | `mcdata_constraint_gumple.py` (SBND only, `--proton-sel 1p`; inclusive + delta-p-sliced rounds) | `constraint/{png,pdf}/` | `Figures/MC-Data-Comparisons/SBND-Abs/` | 9 (of 26 written) | see `logs/constraint.log` |
| `pid` | `nb/MCDataComparisonPID-GUMPLE.ipynb` x {SBND, ICARUS Run2, ICARUS Run4}; `--chi2-flavor nominal --calo-model smearscan --cut-stages simplecosrej_trkqual_nop` | `pid/nominal/{png,pdf}/` | `Figures/PID-MC-Data/` | 62 | ~3 min |
| `sbsyst` | `signalbox_systematics_gumple.py` (every `SBNDMCCV_<i>.df`; replaces the notebook `signalbox` step, 2026-10-20) | `signalbox/{1p,mp,joint}/{png,pdf}/` | `Figures/Signal-Box-Systematics/{1p,mp,joint}/` | 41 | see `logs/sbsyst.log` |
| `corsika` | `offbeam_cosmicmc_gumple.py` | `corsika/{png,pdf}/` | `Figures/CORSIKA-Sytematic/` | 4 (of 40 written) | ~1 min |
| `selection` | `selection_plots_gumple.py` | `selection/{png,pdf}/` | `Figures/cut_stuff/` | 10 | ~5 min (mcnu + evtrec truth read for the efficiency denominator) |
| `pot` | `pot_quality_summary.py` | `pot/` (`<Det>_pot.pdf`, `bq_tables.tex`, `dq_tables.tex`, `pot_summary.txt`) | `Figures/DataQuality/`, `Sections/{bq,dq}_tables.tex` | 3 + 2 tables | see `logs/pot.log` |
| `stability` | `selection_stability_gumple.py --max-groups 15 --min-groups 5` | `stability/{png,pdf}/` | `Figures/DataQuality/` | 12 | see `logs/stability.log` |
| `cvwgt` | `cvweight_impact_gumple.py` | `cvwgt/{png,pdf}/` | `Figures/CV-Weights/` | 4 | see `logs/cvwgt.log` |
| `flashpe` | `flashpe_vs_trueE_gumple.py` | `flashpe_vs_trueE/{png,pdf}/` | `Figures/Trigger-Systematic/` | 3 | see `logs/flashpe.log` |
| `nminus1` | `nminus1_gumple.py -d all` (N-1 cut comparisons, area-normalized; mode breakdowns only, no muon track score, PFP count or far-shower length, in the note) | `evtsel_n-1/{png,pdf}/` | `Figures/MC-Data-Comparisons/N-1/` | 25 | see `logs/nminus1.log` |
| `check` | `technote_figures.py check` | `logs/technote-check.log` | -- | -- | s |
| `install` | `technote_figures.py install` + copy of `pot/{bq,dq}_tables.tex` | `logs/technote-install.log`, `logs/technote-unused.txt` | copies 355 files + 2 tables | -- | s |
| `compile` | `latexmk -pdf -bibtex` | `logs/tex/main.pdf`, `logs/tex/main.log` | -- | -- | ~1 min |

The `pid` flags encode the technote configuration: `sbn-rewgted-21` stores only
the `nominal` chi2 flavor; the note shows the dE/dx-smearing scan (0x / 1x / 2x)
inclusively for all three datasets and the drift-angle quintiles for SBND and
ICARUS Run 2 (the Run 4 quintiles are made too, from the full ReCAF2026 stream
since 2026-10-20, but are not yet in the note), all at the
`simplecosrej_trkqual_nop` stage.

## Where the figures go in the note

Generated from `technote_figures.py list --table` after the 2026-09-02 install
(275 unique images referenced by the compiled document; 244 regenerated here).
`<Det>` runs over `SBND`, `ICARUS-Run2`, `ICARUS-Run4` via `\foreach`.

| technote directory | section (`Sections/*.tex`) | figure labels | count | producer |
|---|---|---|---|---|
| `Figures/cut_stuff/{cosmic_cut,track_score,mu_len,mu_cand,p_cand}.png` | 3.7 Cuts (`Event Selection.tex`) | `fig:cosmic_cuts`, `fig:track_score` (x2, duplicate label), `fig:mu_cand`, `fig:p_cand` | 5 | `selection_plots_gumple.py` (`cut_breakdown_*.png`); `fig:efficnecy` is now the external `updated_eff.png` |
| `Figures/cut_stuff/{e_reco_ratio,del_p_ratio,e_reco_{low,med,hi}_delp}.png` | 3.8 Signal Box Distributions (`SignalBoxDistributions.tex`) | `fig:signalbox-dist-<var>` | 5 | `selection_plots_gumple.py` (`signalbox_*_near_far.png`) |
| `Figures/cut_stuff/1D_scan.png` | 3.7 Cuts | `fig:1Dscan` | 1 | **not regenerated** (rewgt10 cut optimization; FOM machinery not ported) |
| `Figures/Signal-Box-Systematics/{SBND,ICARUS}_{mu,prot}_chi2_of_{mu,prot}_cand_variations.pdf` | 4.2.4 Particle ID Variations (`Systematic Uncertainties.tex`) | `fig:SBNDPIDDetvar`, `fig:ICARUSPIDDetvar` | 8 | `SignalBoxSystematics-GUMPLE.ipynb` |
| `Figures/PID-MC-Data/<Det>_nominal_simplecosrej_trkqual_nop_*_bkdwnpdg_dedx_{nosmear,smear1x,smear2x}.pdf` and `..._nominal_{mu,p}thetadrift{0..4}.pdf` | 4.2.5 MC to Data Comparison PID Validation | `fig:<Det>-MCData-PID`, `fig:<Det>-MCData-PID-perangle-*` | 27 + 25 | `MCDataComparisonPID-GUMPLE.ipynb` |
| `Figures/Track-Splitting/ICARUS-Run{2,4}_simple_*_muend_zoom*.pdf`, `tracksplit_simple_{z0,west}_f_vs_theta_run2run4.pdf`, `ICARUS-Run4_simple_east_f_vs_theta.pdf`, `tracksplit_summary_simple.pdf` | 4.2.6 Track Splitting | `fig:track-split-*` | 29 | `TrackSplittingCorrection_GUMPLE.py` |
| `Figures/CORSIKA-Sytematic/SBND_pid1p_{mu_len,mu_costh,del_p,nu_E_calo}_wnorm.pdf` | 4.2.9 CORSIKA Modeling | `fig:CORSIKA-MC-data` | 4 | `offbeam_cosmicmc_gumple.py` |
| `Figures/Signal-Box-Systematics/{SBND,ICARUS}_altsamples.pdf` | 4.2.10 OOAV Background | `fig:AltSamples` | 2 | `SignalBoxSystematics-GUMPLE.ipynb` |
| `Figures/Signal-Box-Systematics/{flux,xsec,detector,g4}_correlation.pdf`, `{SBND,ICARUS}_signalbox_systematics{,_uncorr,_trueE_uncorr}.pdf`, `ratio_signalbox_systematics{,_trueE}.pdf`, `condconstraint_signalbox_systematics{,_trueE}.pdf` | 4.4 Impact on Signal Box | `fig:signal-box-*`, `fig:signalbox-ICARUS-condconstraint` | 14 | `SignalBoxSystematics-GUMPLE.ipynb` |
| `Figures/MC-Data-Comparisons/<Det>-Area-Norm/<Det>_{contained,twoprongcut,pid1p}_<var>.pdf` | 5 MC to Data Comparisons (`MCDataComparisons.tex`) | `fig:MCData-<Det>-{contained,twoprongcut,pidcuts}-*` | 43 x 3 | `mcdata_comparison_gumple.py` |
| `Figures/MC-Data-Comparisons/SBND-Abs/SBND_pid1p_{mu_p,nu_E_calo,nu_E_ccqe,nu_E_frac_diff}_absnorm.pdf`, `SBND_pid1p_nu_E_calo_constrained_by_nu_E_ccqe{,_lodp,_middp,_hidp}.pdf`, `SBND_signalbox_nu_E_calo_joint1pMp_absnorm.pdf` | 5.1 MC to Data Comparisons, SBND absolute normalization + CC-QE constraint (added by the author 2026-09-02) + joint 1p/Mp signal box (2026-10-24) | `fig:MCData-SBND-abs-nulep-kinematics`, `fig:MCData-SBND-abs-nu-ereco`, `fig:MCData-SBND-abs-nu-constrained`, `fig:MCData-SBND-abs-joint1pMp` | 9 | `mcdata_constraint_gumple.py` |
| `Figures/MC-Data-Comparisons/N-1/<Det>_n-1-{1p,mp}_<var>.pdf` | 5.3 N-1 Event Selection Comparisons (2026-10-24; installed from `plots-gumple-2026-10-23-rewgted25-trigunc`, which has the trigger systematic -- the 10-20 base's own `evtsel_n-1/` predates it, so rerun `--only nminus1` before a full install from 10-20) | `fig:MCData-nminus1-<Det>-{1p,mp}` | 25 (track-score, PFP-count and far-shower N-1 not shown) | `nminus1_gumple.py` |
| `Figures/Track-Splitting/SBND_simple_{east,west}_apa_muend_zoom_th{0..3}.pdf`, `SBND_tracksplit_simple_f_vs_theta_tpc.pdf`, `SBND_tracksplit_summary_tpc.pdf` | SBND Track Splitting at APA Gap (2026-10-24; installed from `plots-gumple-2026-10-24-rewgted25-sbnd-tracksplit`, i.e. the step's output without the `sbndtracksplit/` level) | `fig:track-split-MCdata-APA-angledep-SBND{East,West}`, `fig:track-split-f-vs-theta-SBND`, `fig:track-split-summary-SBND` | 10 | `TrackSplittingCorrection_SBND_GUMPLE.py` |

Figures referenced by the note that are **not** produced here (external inputs,
collaborators' uploads, or the PROfit sensitivity / fake-data studies) and are
left untouched (`technote_figures.py check` prints exactly this set as
"referenced images NOT covered by the manifest"; 419 images on 2026-10-24):
`Figures/chi2_particleID_profiles.pdf`, `icarus_sce.png`, `icarus_no_sce.png`,
`unblinding_full.png`, `Figures/SCE-Detvar/SBNDSpaceChargeMeasurement.png`,
`Figures/Proton-Energy-Resolution/*` (6), `Figures/cut_stuff/1D_scan.png`,
`Figures/cut_stuff/updated_eff.png` (replaced our `eff.png` on 2026-10-20; a
different efficiency definition, see its caption) and
`Figures/cut_stuff/cumulative_efficiencies_final.png` (Mp section, 2026-10-20)
-- no producer in this repo, so `selection`'s `efficiency.png` is no longer
installed --, `Figures/DENT_eff/*` (5), `Figures/AutoPlotter/*` (3),
`Figures/MAPLE_*` (MAPLE data/MC, fractional uncertainty panels, sensitivity),
`Figures/conservative_flux/*`, `Figures/SensitivityStudies/*`,
`Figures/FakeData/**` and `Figures/FakeDataGumple/*` (sections 6-7).

**Trigger scale-factor table must stay in sync.** `tab:trigger-systematic`
(`Sections/Systematic Uncertainties.tex`) quotes the fits drawn in the
`flashfit` figures: the all-trigger scale (the MC CV, `loaddf.MC_PE_SCALE`), the
Contained-stage scale, and the trigger systematic = |all-trigger - Contained|
(`loaddf.MC_PE_SCALE_UNC`, used by `mcdata_comparison_gumple`,
`signalbox_systematics_gumple` and `syst.v_flashscale`). Whenever the `flashfit`
figures are re-installed, update the table from
`<PLOTBASE>/flashpe_fit/trigger_scale_table.tex` and loaddf from
`trigger_scale_summary.tsv` (the fit script logs a WARNING when they disagree);
after a loaddf change re-run `flashpe` and `sbsyst` too (the MC flash PE moves).
The `flashonly` fit needs the ICARUS Run 4 10% unblind file in the flash-only
production: `make_unblind_subsample.py --input <flashonly>/ICARUS_SpringRun4BNB_ReCAF2026.df
--beamfile data/icarus_run4_full_fixed_fx.root --output <flashonly>/ICARUS_SpringRun4BNB_ReCAF2026_unblind.df`.

## What changed in the note on 2026-10-24 (trigger systematic)

* MC CV flash PE scale = all-trigger fit (SBND 0.647, ICARUS Run 2 0.574, Run 4
  0.369; was 0.642 / 0.632 / 0.358); trigger systematic = |all-trigger -
  Contained| (0.029 / 0.097 / 0.076; was |old CV - Contained| = 0.024 / 0.155 /
  0.065). Only the Trigger-section figures (`flashfit`, `flashpe`) and the 1p
  `*_Trig_variations` were regenerated at the new CV
  (`plots-gumple-2026-10-24-rewgted25`); every other figure in the note is still
  at the old CV scale until the next full refresh.
* `Figures/Trigger-Systematic/{SBND,ICARUS,ICARUSRun4}_flash_totalpe.pdf`
  replaced by `AllTriggers/` and `Contained/` (manifest rows added);
  `tab:trigger-systematic` gained a Contained column.

## What changed in the note on 2026-10-20 (sbn-rewgted-25)

* Every manifest figure regenerated from `sbn-rewgted-25` into
  `plots-gumple-2026-10-20-rewgted25` with the new beam- and data-quality code
  (`dataquality.py`): no LM875A/B/C requirement (good-beam definition of
  SBN-doc-49735-v1, `data/good_beam_lists.md`); corrected FOM from the 2026-10
  `*_fx` beam-quality files for every run set; run lists applied as provided
  (SBND 1e20 + post-1e20 "Good, SBND Run 1 total"; no "No End Time" rescue, no
  hand-added SBND run 18255); ICARUS Run 4 ReCAF2026 spills taken from the
  beam-quality ROOT file (`dq.EXTERNAL_BEAM`).
* ICARUS Run 4 MC/data comparisons (`mcdata`) use
  `ICARUS_SpringRun4BNB_ReCAF2026_unblind.df`, a random 10% of the CAF files of
  the full ReCAF2026 dataset with their own spills, made by
  `../make_unblind_subsample.py --input ICARUS_SpringRun4BNB_ReCAF2026.df
  --beamfile icarus_run4_full_fixed_fx.root --output
  ICARUS_SpringRun4BNB_ReCAF2026_unblind.df` (seed 20261020; 3.60e19 of
  3.61e20 POT). The full file is still used by `stability`, `pot` (POT plot,
  BQ/DQ tables) and `tracksplit` (`FULL_ONBEAM_DETECTORS`); the PID notebook
  reads the older `ICARUS_SpringRun4BNB_unblind.df`.
* Manifest: added CV-Weights (4), Signal-Box-Systematics `mp/` (10) and
  `joint/condconstraint_joint_*` (4), `Trigger-Systematic/*_trueE_vs_flashpe`
  (3); dropped rows whose destination the note no longer references
  (`Track-Splitting/ICARUS-Run2_gump_{z0,west}_f_vs_theta.pdf`, `cut_stuff/eff.png`).
* `Sections/{bq,dq}_tables.tex` regenerated by `pot_quality_summary.py`;
  `Sections/gump_recaf_sample_tables.tex` exposure cells recomputed from the -25
  DFs (raw hdr POT / ngenevt / off-beam gates; text untouched, generator not in
  this repo); `Sections/Samples.tex` BQ/DQ paragraphs updated.

## What changed in the note on 2026-09-02

* Pruned 483 image files (38.5 MB) that no compiled `\includegraphics` referenced
  (`technote_figures.py unused --keep-commented`). Kept: the 3 SBND SCE PNGs in a
  commented-out figure, the 9 `Figures/pid+angle/` files used by the orphan
  `Sections/SignalBoxDistributions_Han.tex` (never `\input`), and that file.
* `Sections/MCDataComparisons.tex`: GUMPLE names -- stage `pidcut` -> `pid1p`;
  `other_trk_length` -> `othr_pfp_length`, `other_shw_length` -> `n_shower`
  (the GUMPLE dataframes merge track and shower lengths into the longest
  non-candidate pfp, and carry shower / other counts). Removed a duplicated
  vertex-position figure that carried the `-pidcuts-muend` label, so the muon-
  and proton-endpoint figures now carry `-pidcuts-muend` / `-pidcuts-pend`
  (a stray `i` after the last label is gone too).
* `Sections/Systematic Uncertainties.tex`: CORSIKA figure paths `SBND_pidcut_*`
  -> `SBND_pid1p_*`.
* Old-name figures made unreferenced by those renames were pruned as well.
* 2026-09-02 (later): the author added three SBND absolute-normalization /
  CC-QE-constraint figures to section 5 (8 files under
  `Figures/MC-Data-Comparisons/SBND-Abs/`). They are the outputs of
  `mcdata_constraint_gumple.py`, so a `constraint` step was added to the driver
  and manifest and the figures regenerated on sbn-rewgted-21. The 12 sibling
  files uploaded alongside them but not referenced (`*_unconstrained*`,
  `*_constraintvar*`, the non-absnorm copies) were pruned.

## Caveats to carry into the text

* **No neutrino-score cut in production.** `gumple_cuts.cosmic_cut` is the
  opening-angle cut plus the CRT veto; `cut_breakdown_cosmic.png` still shows
  `nu_score` (marked "no cut in production"). Table `tab:cuts` and the cosmic-cut
  caption in section 3 describe a score cut.
* **The muon-length lower bound (40 cm) is applied in the production
  preselection**, so no candidates below it reach the dataframes; the
  `mu_len.png` figure shows only the upper (400 cm) cut acting.
* **CORSIKA normalization**: the `_wnorm` figures use the analysis' nominal
  10.7% normalization uncertainty (`mcdata_comparison_gumple.SBND_COSMIC_NORM`,
  `offbeam_cosmicmc_gumple.py --norm-unc nominal`, the default). The
  gate-normalized MC/data offset measured on sbn-rewgted-21 after the 1p
  selection is 6.4% (`--norm-unc measured`); it is printed in
  `<PLOTBASE>/corsika/SBND_offbeam_cosmicmc_summary.txt`.
* **Efficiency curves are now absolute.** The denominator is every generated
  neutrino interaction with true vertex in the analysis FV, taken from the
  `mcnu` table (reconstructed or not) joined to the GENIE event record
  (`evtrec`, status-1 final-state particles) for the leading muon / proton /
  charged-pion momenta; those reproduce the production's stored `true_*_p`
  exactly. Numerator and denominator both carry the POT scale and the aFF CV
  weight. Six panels: muon / proton / charged-pion vs momentum, and pi0
  events, numu CC QE events and numu CC 1u1p events (exactly one proton with
  T_p > 50 MeV, no charged pion, no pi0) vs true neutrino energy. The previous
  figure was post-selection over pre-selection with four panels, so the
  scale is not comparable (the old caption says "post-selection/pre-selection").
* **Signal-box systematics binning** (`nb/SignalBoxSystematics-GUMPLE.ipynb`,
  cell 54): reco `nu_E_calo` edges 0.3, 0.4, 0.45, ..., 0.95, 1.0, 1.25, 1.5 GeV
  (15 bins; the fine scheme, restored 2026-09-02) and true energy = the same
  edges wrapped with under/overflow (-2, ..., inf). The trueE plots drop their
  under/overflow bins; the reco plots drop edge bins only under the coarse
  0.0/0.4..1.0/1.5 scheme (`DROP_EDGE_BINS`).
* The `SBND-Abs` figures are absolutely normalized and unblind the SBND dev
  sample in muon momentum and neutrino energy. `mcdata_constraint_gumple.py`
  refuses any detector but SBND and any proton selection but 1p
  (`BlindedError`) -- do not extend the driver to ICARUS for this step.
* The sample tables (`gump_recaf_sample_tables.tex`) and the cut table were not
  regenerated.

## Maintaining the manifest

When the note changes (new figure, renamed variable, new section):

1. `python3 technote_figures.py list --technote <dir> --table` -- every image the
   compiled document references, with tex file / section / label / caption.
   Anything printed as `MISSING:` has no file on disk yet.
2. Find the producer of each new figure by its filename pattern (`grep -n
   savefig analysis_village/gump/*.py`, the `_save` helpers, the notebooks'
   `PLOTDIR + ...` calls) and add a `source<TAB>dest` row to `figure_manifest.tsv`.
   The source is relative to `PLOTBASE` (`evtsel/pdf/...`, `pid/nominal/pdf/...`,
   `tracksplit/pdf/...`, `signalbox/pdf/...`, `corsika/pdf/...`, `selection/png/...`).
3. `python3 technote_figures.py check ...` must report `problems: 0`; its
   "NOT covered" list should contain only the external figures.
4. `python3 technote_figures.py unused --technote <dir> --keep-commented` lists
   what to `git rm` in the technote checkout after an install.
