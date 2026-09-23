# GUMP technote figures from the GUMPLE dataframes

This directory holds everything needed to regenerate the locally reproducible
figures of the GUMP technical note (Overleaf project `6973acd1e6f673f7a8b495e7`,
checked out beside this repo as `cafpyana/6973acd1e6f673f7a8b495e7/`) from a
GUMPLE production (`sbn-rewgted-21` as of 2026-09-02) and install them into the
Overleaf checkout.

| file | role |
|---|---|
| `make_technote_plots.sh` | master driver: runs every analysis step, checks the manifest, installs the figures, compiles the note |
| `figure_manifest.tsv` | one row per regenerated figure: `<plotbase-relative source>` TAB `<technote-relative destination>` |
| `technote_figures.py` | stdlib inventory tool: `list` / `unused` / `check` / `install` (parses `\includegraphics`, expands `\foreach`) |
| `figure_placement.tsv` | snapshot of `technote_figures.py list --table` after the last install: every referenced image with its tex file, section, label and caption |
| `README.md` | this file |

The Claude skill `.claude/skills/technote-plots/SKILL.md` (repo root) describes
how to redo the whole exercise for a future iteration of the note.

## Prerequisites

* The GUMPLE dataframes: `--df-dir` must hold the 47 `.df` files of the
  production (`SBNDMCCV_*.df`, `ICARUSRun{2,4}_SpringMCOverlay_rewgt_*.df`, the
  detvar / dirt / beam-off / on-beam files, `SBNDIntimeMC.df`).
* The interpreter `/Users/gputnam-local/Work/fitter/env/bin/python` (+ its
  `jupyter`) -- the only local environment that imports cafpyana. It is
  hard-coded in `run_gumple_serial.py` and in the driver.
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
DF_DIR=/path/to/sbn-rewgted-22/ PLOTBASE=$PWD/../../../plots-gumple-XX ./make_technote_plots.sh
```

Defaults (env or flag): `DF_DIR=/Users/gputnam/Work/osc/sbn-rewgted-21/`,
`PLOTBASE=<repo>/plots-gumple-2026-09-02-rewgted21`,
`TECHNOTE=<repo>/6973acd1e6f673f7a8b495e7`. Logs go to `<PLOTBASE>/logs/`
(`technote-serial.log`, one `<step>.log` per analysis, `technote-check.log`,
`technote-install.log`, `technote-unused.txt`, `technote-compile.log`,
`tex/main.log`). Steps are independent: a failure is reported in the summary
table and the chain continues. Nothing is committed to either git repo.

## Steps

All analysis steps are invoked through `run_gumple_serial.py` (which passes
`--df-dir` / `--plotdir` to the scripts and `GUMPLE_DF_DIR` / `GUMPLE_PLOTBASE`
to the notebooks). Runtimes are from the 2026-09-02 run against
`sbn-rewgted-21` with warm loaddf caches.

| step | producer | outputs under `<PLOTBASE>` | technote destination | figures | runtime |
|---|---|---|---|---|---|
| `tracksplit` | `TrackSplittingCorrection_GUMPLE.py` (all detectors, stages simple+gump, angular Run 2, 4 quantiles) | `tracksplit/{png,pdf}/` | `Figures/Track-Splitting/` | 16 | 0.6 min |
| `mcdata` | `mcdata_comparison_gumple.py -d all --nproc 8` | `evtsel/{png,pdf}/` | `Figures/MC-Data-Comparisons/<Det>-Area-Norm/` | 43 x 3 | 7.2 min |
| `constraint` | `mcdata_constraint_gumple.py` (SBND only, `--proton-sel 1p`; inclusive + delta-p-sliced rounds) | `constraint/{png,pdf}/` | `Figures/MC-Data-Comparisons/SBND-Abs/` | 8 (of 26 written) | see `logs/constraint.log` |
| `pid` | `nb/MCDataComparisonPID-GUMPLE.ipynb` x {SBND, ICARUS Run2, ICARUS Run4}; `--chi2-flavor nominal --calo-model smearscan --no-angles "ICARUS Run4" --cut-stages simplecosrej_trkqual_nop` | `pid/nominal/{png,pdf}/` | `Figures/PID-MC-Data/` | 52 | ~3 min |
| `signalbox` | `nb/SignalBoxSystematics-GUMPLE.ipynb` | `signalbox/{png,pdf}/` | `Figures/Signal-Box-Systematics/` | 24 (of 68 written) | 4.7 min |
| `corsika` | `offbeam_cosmicmc_gumple.py` | `corsika/{png,pdf}/` | `Figures/CORSIKA-Sytematic/` | 4 (of 40 written) | ~1 min |
| `selection` | `selection_plots_gumple.py` | `selection/{png,pdf}/` | `Figures/cut_stuff/` | 11 | ~5 min (mcnu + evtrec truth read for the efficiency denominator) |
| `check` | `technote_figures.py check` | `logs/technote-check.log` | -- | -- | s |
| `install` | `technote_figures.py install` | `logs/technote-install.log`, `logs/technote-unused.txt` | copies 244 files | -- | s |
| `compile` | `latexmk -pdf -bibtex` | `logs/tex/main.pdf`, `logs/tex/main.log` | -- | -- | ~1 min |

The `pid` flags encode the technote configuration: `sbn-rewgted-21` stores only
the `nominal` chi2 flavor; the note shows the dE/dx-smearing scan (0x / 1x / 2x)
inclusively for all three datasets and the drift-angle quintiles for SBND and
ICARUS Run 2 only (Run 4 rides the prescaled stream), all at the
`simplecosrej_trkqual_nop` stage.

## Where the figures go in the note

Generated from `technote_figures.py list --table` after the 2026-09-02 install
(275 unique images referenced by the compiled document; 244 regenerated here).
`<Det>` runs over `SBND`, `ICARUS-Run2`, `ICARUS-Run4` via `\foreach`.

| technote directory | section (`Sections/*.tex`) | figure labels | count | producer |
|---|---|---|---|---|
| `Figures/cut_stuff/{cosmic_cut,track_score,mu_len,mu_cand,p_cand,eff}.png` | 3.7 Cuts (`Event Selection.tex`) | `fig:cosmic_cuts`, `fig:track_score` (x2, duplicate label), `fig:mu_cand`, `fig:p_cand`, `fig:efficnecy` | 6 | `selection_plots_gumple.py` (`cut_breakdown_*.png`, `efficiency.png`) |
| `Figures/cut_stuff/{e_reco_ratio,del_p_ratio,e_reco_{low,med,hi}_delp}.png` | 3.8 Signal Box Distributions (`SignalBoxDistributions.tex`) | `fig:signalbox-dist-<var>` | 5 | `selection_plots_gumple.py` (`signalbox_*_near_far.png`) |
| `Figures/cut_stuff/1D_scan.png` | 3.7 Cuts | `fig:1Dscan` | 1 | **not regenerated** (rewgt10 cut optimization; FOM machinery not ported) |
| `Figures/Signal-Box-Systematics/{SBND,ICARUS}_{mu,prot}_chi2_of_{mu,prot}_cand_variations.pdf` | 4.2.4 Particle ID Variations (`Systematic Uncertainties.tex`) | `fig:SBNDPIDDetvar`, `fig:ICARUSPIDDetvar` | 8 | `SignalBoxSystematics-GUMPLE.ipynb` |
| `Figures/PID-MC-Data/<Det>_nominal_simplecosrej_trkqual_nop_*_bkdwnpdg_dedx_{nosmear,smear1x,smear2x}.pdf` and `..._nominal_{mu,p}thetadrift{0..4}.pdf` | 4.2.5 MC to Data Comparison PID Validation | `fig:<Det>-MCData-PID`, `fig:<Det>-MCData-PID-perangle-*` | 27 + 25 | `MCDataComparisonPID-GUMPLE.ipynb` |
| `Figures/Track-Splitting/ICARUS-Run{2,4}_simple_*_muend_zoom*.pdf`, `ICARUS-Run2_gump_{z0,west}_f_vs_theta.pdf`, `tracksplit_summary_simple.pdf` | 4.2.6 Track Splitting | `fig:track-split-*`, `fig:placeholder` | 16 | `TrackSplittingCorrection_GUMPLE.py` |
| `Figures/CORSIKA-Sytematic/SBND_pid1p_{mu_len,mu_costh,del_p,nu_E_calo}_wnorm.pdf` | 4.2.9 CORSIKA Modeling | `fig:CORSIKA-MC-data` | 4 | `offbeam_cosmicmc_gumple.py` |
| `Figures/Signal-Box-Systematics/{SBND,ICARUS}_altsamples.pdf` | 4.2.10 OOAV Background | `fig:AltSamples` | 2 | `SignalBoxSystematics-GUMPLE.ipynb` |
| `Figures/Signal-Box-Systematics/{flux,xsec,detector,g4}_correlation.pdf`, `{SBND,ICARUS}_signalbox_systematics{,_uncorr,_trueE_uncorr}.pdf`, `ratio_signalbox_systematics{,_trueE}.pdf`, `condconstraint_signalbox_systematics{,_trueE}.pdf` | 4.4 Impact on Signal Box | `fig:signal-box-*`, `fig:signalbox-ICARUS-condconstraint` | 14 | `SignalBoxSystematics-GUMPLE.ipynb` |
| `Figures/MC-Data-Comparisons/<Det>-Area-Norm/<Det>_{contained,twoprongcut,pid1p}_<var>.pdf` | 5 MC to Data Comparisons (`MCDataComparisons.tex`) | `fig:MCData-<Det>-{contained,twoprongcut,pidcuts}-*` | 43 x 3 | `mcdata_comparison_gumple.py` |
| `Figures/MC-Data-Comparisons/SBND-Abs/SBND_pid1p_{mu_p,nu_E_calo,nu_E_ccqe,nu_E_frac_diff}_absnorm.pdf`, `SBND_pid1p_nu_E_calo_constrained_by_nu_E_ccqe{,_lodp,_middp,_hidp}.pdf` | 5 MC to Data Comparisons, SBND absolute normalization + CC-QE constraint (added by the author 2026-09-02) | `fig:MCData-SBND-abs-nulep-kinematics`, `fig:MCData-SBND-abs-nu-ereco`, `fig:MCData-SBND-abs-nu-constrained` | 8 | `mcdata_constraint_gumple.py` |

Figures referenced by the note that are **not** produced here (external inputs
or the PROfit sensitivity / fake-data studies) and were left untouched:
`Figures/chi2_particleID_profiles.pdf`, `icarus_sce.png`, `icarus_no_sce.png`,
`Figures/SCE-Detvar/SBNDSpaceChargeMeasurement.png`,
`Figures/Proton-Energy-Resolution/*` (6), `Figures/Trigger-Systematic/*_flash_totalpe.pdf` (3),
`Figures/cut_stuff/1D_scan.png`, `Figures/SensitivityStudies/*` (7 files, sections 6),
`Figures/FakeData/*` (10 files, section 7). `technote_figures.py check` prints
this list as "referenced images NOT covered by the manifest".

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
