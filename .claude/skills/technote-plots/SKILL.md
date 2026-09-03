---
name: technote-plots
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
  sensitivity / fake-data outputs, flash-PE fits, anything the user names) are
  left alone and listed explicitly in the final report.

Write / extend `figure_manifest.tsv` (`<source rel. to PLOTBASE>\t<dest rel. to
technote>`) so that `technote_figures.py check` reports `problems: 0` once the
analyses have run, and its "NOT covered" list is exactly the external set.

## Step 4 — Run the analyses

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
