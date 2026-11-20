---
name: technote-plots
argument-hint: "[update_code=true|false]  (default false)"
description: Regenerate and install the figures of a LaTeX technical note kept in an Overleaf git repo from this repository's plotting scripts and notebooks, after a new dataframe production. Inventories every \includegraphics (expanding \foreach), prunes unreferenced images, maps each figure to the script that makes it, runs the analyses through the master driver, installs the outputs by manifest, recompiles, and reports what to commit in both repos. Use when the user asks to update, refresh, or redo the plots in a technote / Overleaf document, to point the note at a new sbn-rewgted-NN production, or to prune unused figures from an Overleaf project.
---

# /technote-plots — refresh a technote's figures from the analysis scripts

Tooling lives in `analysis_village/gump/technote/` (`README.md` there is the
reference for the GUMP note specifically): `technote_figures.py` (inventory /
prune / check / install), `figure_manifest.tsv` (source → destination rows) and
`make_technote_plots.sh` (runs the analyses via `run_gumple_serial.py`, then
check → install → compile). Do **not** assume the note still has the sections,
figure names or `\foreach` lists it had last time — re-derive everything from
the current `.tex` files with the tool. Never commit or push; the user does that.

## Arguments

* `update_code` (default **false**; `update_code=true`, `update_code`, or
  `--update-code` in the args turn it on). Some analysis steps *measure*
  constants that are hard-coded in `analysis_village/gump/loaddf.py` and used by
  every later step:

  | calibration step | measures | loaddf constant |
  |---|---|---|
  | `flashfit` (trigger study) | all-trigger MC flash-PE scale; trigger syst. = \|all-trigger - Contained\| | `MC_PE_SCALE`, `MC_PE_SCALE_UNC` |
  | `tracksplit` (ICARUS track splitting) | split fraction per plane, Simple stage, Run2+Run4 | `SPLIT_FRAC` |
  | `sbndsplit` (SBND APA gap) | SBND split fraction | none (report only) |

  The driver always runs these first and then `technote/sync_calibrations.py`.
  * **false**: if the fits disagree with loaddf, the driver stops before the
    downstream steps. Report the table from `logs/technote-sync.log` and ask the
    user whether to update loaddf (re-run with `--update-code`), continue with the
    stale constants (`--allow-stale`), or stop. Do not edit loaddf yourself.
  * **true**: pass `--update-code`. The driver rewrites the literals in
    `loaddf.py`, re-draws `tracksplit` (its post-split MC uses `SPLIT_FRAC`) and
    continues. Report the `git diff analysis_village/gump/loaddf.py` and include it
    in what is staged.

  Either way, the tex that quotes these numbers (`tab:trigger-systematic`, the
  track-splitting fractions) is **not** edited unless the user asks; list the new
  values from `flashpe_fit/trigger_scale_table.tex` and
  `tracksplit/tracksplit_summary.txt` in the report.

## Step 1 — Baseline the Overleaf checkout

1. Locate the checkout (a directory with `main.tex`, `Sections/`, `Figures/`,
   `.git` whose remote is `git.overleaf.com`). If asked, `git pull` first; note
   the HEAD commit and confirm `git status` is clean before touching anything.
2. Compile from the project root into a scratch out-dir and record the
   baseline: exit code, `grep -c '^!' main.log` (must be 0), the sorted set of
   `Reference ... undefined` warnings, and the page count (`(N pages` in
   `main.log`). Every later compile must match this except where a figure was
   deliberately added or removed.
   ```bash
   latexmk -pdf -bibtex -interaction=nonstopmode -outdir=<scratch>/tex-before main.tex
   ```
3. Inventory: `technote_figures.py list --technote <dir> --table` gives every
   referenced image with tex file / section / label / caption; the stderr line
   gives references / unique images / missing-on-disk. Save it.

## Step 2 — Prune unreferenced figures (ask before being aggressive)

`technote_figures.py unused --technote <dir> --keep-commented` lists images no
compiled `\includegraphics` uses. `--keep-commented` protects files named only in
commented-out lines or in `.tex` files that are never `\input` (orphan drafts).
Ask the user whether to keep or drop those protected files and the orphan `.tex`
itself; default to keeping them. Remove with `git rm -q --` in batches (paths
may contain spaces or `$`), delete emptied directories, recompile, and confirm
the baseline numbers are unchanged.

## Step 3 — Map each referenced figure to its producer

Group the inventory by directory / filename pattern and find the script or
notebook that writes each pattern:

* `grep -n "savefig\|_save(" analysis_village/gump/*.py`, and for notebooks dump
  the `PLOTDIR + ...` / `savefig` lines from the `.ipynb` JSON.
* Output-name templates to expect: MC/data `<Det>_<cutstem>_<var>[suffix].pdf`
  (`mcdata_comparison_gumple.py`), PID
  `<Det>_<flavor>_<cuttag>_<var>_bkdwnpdg_<model>[angle].pdf` (PID notebook),
  track splitting `<Det>_<stage>_<plane>_muend_zoom[_th<i>].pdf` etc.,
  signal-box systematics `{SBND,ICARUS}_signalbox_systematics*.pdf`,
  `*_correlation.pdf`, `*_variations.pdf` (SignalBox notebook), CORSIKA
  `SBND_<cutstem>_<var>_wnorm.pdf` (`offbeam_cosmicmc_gumple.py`), SBND
  absolute-normalization and constraint `SBND_pid1p_<var>_absnorm.pdf`,
  `SBND_pid1p_<y>_constrained_by_<x>[_<dp>].pdf` (`mcdata_constraint_gumple.py`,
  SBND/1p only -- blinded elsewhere), cut breakdowns / efficiency / near-far
  (`selection_plots_gumple.py`).
* Authors often upload a producer's whole output directory alongside the few
  files they reference. Diff the note against the last refresh commit
  (`git diff <sha> HEAD -- Sections/`) to see what is actually new, add only
  the referenced files to the manifest, and let `unused` catch the rest.
* When a production renames a stage or variable (e.g. `pidcut` → `pid1p`,
  `other_trk_length` → `othr_pfp_length`), prefer editing the `.tex` to the new
  names over installing new plots under stale names — and tell the user which
  captions/tables may now disagree with the plots.
* Figures whose inputs are not on this machine (external PDFs, PROfit
  sensitivity / fake-data outputs, anything the user names) are left alone and
  listed explicitly in the final report. Collaborators also upload figures made
  by their own code straight into the note (e.g. `cut_stuff/updated_eff.png`,
  `DENT_eff/`, `AutoPlotter/` in Oct 2026): search the repo, all branches
  (`git log --all -S<name>`), before deciding; with no producer they are external,
  not manifest rows, even when they replace one of ours (drop our row).
* The inventory expands `\foreach` lists brace-aware (`\var/\cap in {a/{x, y}}`),
  zero-argument macros (`\providecommand{\FDFigDir}{...}`, also in `\input`) and
  parameterized `\newcommand` figure macros at their call sites. If `list` shows
  a path with `\name`, `#1` or caption fragments, extend the parser in
  `technote_figures.py` rather than working around it in the manifest.

Write / extend `figure_manifest.tsv` (`<source rel. to PLOTBASE>\t<dest rel. to
technote>`) so that `technote_figures.py check` reports `problems: 0` once the
analyses have run, and its "NOT covered" list is exactly the external set.

## Step 4 — Run the analyses

* Order matters: the driver runs phase 1 calibrate (`flashfit tracksplit
  sbndsplit`), phase 2 `sync`, phase 3 downstream (`mcdata constraint pid sbsyst
  corsika selection pot stability cvwgt flashpe nminus1`), phase 4 `check install
  compile`. Never run a downstream step on a plotbase whose calibration has not
  been synced (see Arguments). When running steps by hand with `--only`, keep the
  same order and run `--only sync` after any calibration step. A new step that
  measures a loaddf constant goes in `CALIB_STEPS` and `sync_calibrations.py`.
* Check disk first: loaddf writes `<hash>.h5` caches *inside* the production
  directory; a new production means tens of GB. Watch `df -h` during long steps.
* Run `make_technote_plots.sh --df-dir <production> --plotbase <new plots dir>
  --technote <checkout>` (or `--only <step>` for a subset). It runs the analyses
  through `run_gumple_serial.py` (cwd `analysis_village/gump`, the fitter-env
  interpreter, per-step logs in `<plotbase>/logs/`), then `check`, `install`,
  `compile`. Long steps belong in the background; tail the step logs.
* New or ported producers must keep the module-scope `FV` preselection name
  (loaddf's cache key hashes the preselection qualname) and take `--df-dir` /
  `--plotdir`; notebooks read `GUMPLE_DF_DIR` / `GUMPLE_PLOTBASE`. Add a new
  step to `run_gumple_serial.py` and to the driver rather than running ad hoc.
* Spot-check one figure per producer visually (Read the PNG) before installing.

## Step 5 — Install, recompile, report

**Table sync for the trigger figures:** the `flashpe_fit/` (`flashfit` step)
figures are quoted by `tab:trigger-systematic` in
`Sections/Systematic Uncertainties.tex` (scale factors and the
`gumple_cuts.FLASH_PE_THRESHOLD` values). Update it from
`<PLOTBASE>/flashpe_fit/trigger_scale_table.tex` when the user wants the tex
touched; otherwise flag every changed number in the report. loaddf itself is
handled by the `sync` phase (Arguments).


1. `technote_figures.py install` copies by manifest; then `unused
   --keep-commented` again lists figures orphaned by tex renames — `git rm` them.
2. Recompile; compare errors, undefined references and page count with the
   baseline and explain every difference.
3. Update `README.md` (steps, runtimes, placement table, caveats, tex edits) and
   the memory notes if conventions changed.
4. Final report to the user, standing on its own:
   * what was regenerated (counts per producer), what was deliberately not,
   * every `.tex` edit and any text/table now inconsistent with the new plots
     (e.g. cut values, quoted normalizations),
   * the exact file lists to `git add` / commit / push in the Overleaf repo and
     in this repo (`git status --short` in both), without committing yourself.
