#!/usr/bin/env bash
# make_technote_plots.sh -- regenerate every locally reproducible figure of the
# GUMP technote from a GUMPLE production and install them into the Overleaf
# checkout. See README.md (sibling) for what each step produces and where the
# figures land in the note.
#
#   usage: make_technote_plots.sh [--df-dir DIR] [--flashonly-df-dir DIR] [--plotbase DIR] [--technote DIR]
#                                 [--only STEP ...] [--update-code] [--allow-stale]
#                                 [--no-install] [--no-compile] [--no-check]
#
#   Steps run in four phases, because the calibration steps measure constants
#   (hard-coded in loaddf.py) that every later step uses:
#
#   1 calibrate
#     flashfit    flashpe_fit_gumple.py --schema {gumple,flashonly} -> <plotbase>/flashpe_fit/
#                 (all-trigger MC flash-PE scale + trigger systematic -> loaddf.MC_PE_SCALE{,_UNC})
#     tracksplit  TrackSplittingCorrection_GUMPLE.py           -> <plotbase>/tracksplit/
#                 (ICARUS Run2+Run4 split fractions, Simple stage -> loaddf.SPLIT_FRAC)
#     sbndsplit   TrackSplittingCorrection_SBND_GUMPLE.py (APA gap) -> <plotbase>/sbndtracksplit/
#                 (report only: no SBND split is applied)
#   2 sync        sync_calibrations.py: compare the measured constants with loaddf.py.
#                 On a mismatch: with --update-code, write them into loaddf.py and
#                 re-run tracksplit (its post-split figures draw the split MC at
#                 loaddf.SPLIT_FRAC); without it, stop before phase 3 unless
#                 --allow-stale is given.
#   3 downstream
#     mcdata      mcdata_comparison_gumple.py -d all           -> <plotbase>/evtsel/
#     constraint  mcdata_constraint_gumple.py (SBND abs-norm,  -> <plotbase>/constraint/
#                 CC-QE-energy constraint, incl. delta-p slices)
#     pid         nb/MCDataComparisonPID-GUMPLE.ipynb (x3 det) -> <plotbase>/pid/nominal/
#     sbsyst      signalbox_systematics_gumple.py              -> <plotbase>/signalbox/{1p,mp,joint}/
#     corsika     offbeam_cosmicmc_gumple.py                   -> <plotbase>/corsika/
#     selection   selection_plots_gumple.py                    -> <plotbase>/selection/
#     pot         pot_quality_summary.py (BQ/DQ POT, tables)   -> <plotbase>/pot/
#     stability   selection_stability_gumple.py                -> <plotbase>/stability/
#     cvwgt       cvweight_impact_gumple.py                    -> <plotbase>/cvwgt/
#     flashpe     flashpe_vs_trueE_gumple.py                   -> <plotbase>/flashpe_vs_trueE/
#     nminus1     nminus1_gumple.py -d all (N-1 cut comparisons) -> <plotbase>/evtsel_n-1/
#   4 check       technote_figures.py check   (manifest sources exist, dests referenced)
#     install     technote_figures.py install (copy into the technote checkout), plus
#                 <plotbase>/pot/{bq,dq}_tables.tex -> Sections/
#     compile     latexmk the technote into <plotbase>/logs/tex/ and report errors
#
#   The tex tables quoting the calibration results (tab:trigger-systematic, the
#   track-splitting fractions) are NOT edited: update them by hand from
#   <plotbase>/flashpe_fit/trigger_scale_table.tex and tracksplit/tracksplit_summary.txt.
#
# Analysis steps go through run_gumple_serial.py (same interpreter, cwd, env
# passthrough and per-step logs as the regular GUMPLE chain); a failing step is
# reported at the end but does not stop the later ones (only a calibration
# mismatch stops the chain, see phase 2). Nothing is committed.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GUMP_DIR="$(dirname "$HERE")"
REPO="$(cd "$GUMP_DIR/../.." && pwd)"
PY=python3

DF_DIR="${DF_DIR:-/Users/gputnam/Work/osc/sbn-rewgted-25/}"
FLASHONLY_DF_DIR="${FLASHONLY_DF_DIR:-}"   # empty -> run_gumple_serial default <df-dir>-flashonly/
PLOTBASE="${PLOTBASE:-$REPO/plots-gumple-2026-10-24-full}"
TECHNOTE="${TECHNOTE:-$REPO/6973acd1e6f673f7a8b495e7}"
MANIFEST="$HERE/figure_manifest.tsv"
ONLY=()
DO_INSTALL=1
DO_COMPILE=1
DO_CHECK=1
UPDATE_CODE=0
ALLOW_STALE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --df-dir)     DF_DIR="$2"; shift 2 ;;
    --flashonly-df-dir) FLASHONLY_DF_DIR="$2"; shift 2 ;;
    --plotbase)   PLOTBASE="$2"; shift 2 ;;
    --technote)   TECHNOTE="$2"; shift 2 ;;
    --manifest)   MANIFEST="$2"; shift 2 ;;
    --only)       ONLY+=("$2"); shift 2 ;;
    --no-install) DO_INSTALL=0; shift ;;
    --no-compile) DO_COMPILE=0; shift ;;
    --no-check)   DO_CHECK=0; shift ;;
    --update-code) UPDATE_CODE=1; shift ;;
    --allow-stale) ALLOW_STALE=1; shift ;;
    -h|--help)    sed -n '2,53p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

CALIB_STEPS=(flashfit tracksplit sbndsplit)
DOWNSTREAM_STEPS=(mcdata constraint pid sbsyst corsika selection pot stability cvwgt flashpe nminus1)
ALL_STEPS=("${CALIB_STEPS[@]}" sync "${DOWNSTREAM_STEPS[@]}" check install compile)
want() { # step in ONLY (or ONLY empty)
  [[ ${#ONLY[@]} -eq 0 ]] && return 0
  local s; for s in "${ONLY[@]}"; do [[ "$s" == "$1" ]] && return 0; done; return 1
}

mkdir -p "$PLOTBASE/logs"
echo "[cfg] df-dir    $DF_DIR"
echo "[cfg] plotbase  $PLOTBASE"
echo "[cfg] technote  $TECHNOTE"
echo "[cfg] manifest  $MANIFEST"
echo "[cfg] steps     ${ONLY[*]:-${ALL_STEPS[*]}}"
echo "[cfg] update-code $UPDATE_CODE  allow-stale $ALLOW_STALE"

declare -a RESULT_NAMES=() RESULT_CODES=()
record() { RESULT_NAMES+=("$1"); RESULT_CODES+=("$2"); }

# --- analysis steps, via the serial runner ---------------------------------
# The PID step is the technote configuration: the single 'nominal' chi2 flavor
# stored by sbn-rewgted-21, the dE/dx-smearing scan (0x/1x/2x) on the
# angle-inclusive plots, angular quintiles for every detector (ICARUS Run 4 reads
# the full ReCAF2026 stream since 2026-10-20), at the 'simplecosrej_trkqual_nop'
# cut stage. run_gumple_serial.py runs its --only steps in its own order, so
# each phase is a separate call.
serial() { # serial <label> <step>...   (the steps selected by --only, or all)
  local label="$1"; shift
  local steps=() s
  for s in "$@"; do want "$s" && steps+=("$s"); done
  [[ ${#steps[@]} -eq 0 ]] && return 0
  serial_run "$label" "${steps[@]}"
}
serial_run() { # serial_run <label> <step>...  (unconditionally)
  local label="$1"; shift
  local args=() s
  for s in "$@"; do args+=(--only "$s"); done
  local fo=(); [[ -n "$FLASHONLY_DF_DIR" ]] && fo=(--flashonly-df-dir "$FLASHONLY_DF_DIR")
  ( cd "$GUMP_DIR" && "$PY" run_gumple_serial.py \
      --df-dir "$DF_DIR" "${fo[@]+"${fo[@]}"}" --skip-wait --plotbase "$PLOTBASE" \
      "${args[@]}" \
      --chi2-flavor nominal --calo-model smearscan \
      --cut-stages simplecosrej_trkqual_nop \
      2>&1 | tee -a "$PLOTBASE/logs/technote-serial.log" )
  record "$label (run_gumple_serial.py)" "${PIPESTATUS[0]}"
}

# phase 1: calibration steps
serial "calibrate" "${CALIB_STEPS[@]}"

# phase 2: are loaddf's constants what phase 1 measured?
STALE=0
if want sync; then
  python3 "$HERE/sync_calibrations.py" --plotbase "$PLOTBASE" 2>&1 | tee "$PLOTBASE/logs/technote-sync.log"
  rc=${PIPESTATUS[0]}
  if [[ $rc -eq 3 && $UPDATE_CODE -eq 1 ]]; then
    python3 "$HERE/sync_calibrations.py" --plotbase "$PLOTBASE" --apply 2>&1 | tee -a "$PLOTBASE/logs/technote-sync.log"
    rc=${PIPESTATUS[0]}
    record "sync (loaddf.py updated)" "$rc"
    # the post-split track-splitting figures draw the split MC at loaddf.SPLIT_FRAC
    [[ $rc -eq 0 ]] && serial_run "re-draw tracksplit" tracksplit
  elif [[ $rc -eq 3 ]]; then
    STALE=1
    record "sync (loaddf.py differs from the fits)" 3
    echo "[sync] loaddf.py constants differ from the calibration fits (see above)."
    echo "[sync] Re-run with --update-code to write them, or --allow-stale to continue anyway."
  else
    record "sync" "$rc"
  fi
fi

# phase 3: everything that uses the calibrated constants
STOP=0
if [[ $STALE -eq 1 && $ALLOW_STALE -eq 0 ]]; then
  echo "[sync] STOPPING before the downstream steps."
  STOP=1
else
  serial "downstream" "${DOWNSTREAM_STEPS[@]}"
fi

# --- manifest check / install ------------------------------------------------
if want check && [[ $DO_CHECK -eq 1 && $STOP -eq 0 ]]; then
  python3 "$HERE/technote_figures.py" check --technote "$TECHNOTE" \
      --plotbase "$PLOTBASE" --manifest "$MANIFEST" 2>&1 | tee "$PLOTBASE/logs/technote-check.log"
  record "check" "${PIPESTATUS[0]}"
fi

if want install && [[ $DO_INSTALL -eq 1 && $STOP -eq 0 ]]; then
  python3 "$HERE/technote_figures.py" install --technote "$TECHNOTE" \
      --plotbase "$PLOTBASE" --manifest "$MANIFEST" 2>&1 | tee "$PLOTBASE/logs/technote-install.log"
  record "install" "${PIPESTATUS[0]}"
  # the beam-/data-quality tables are \input by Sections/Samples.tex
  for t in bq_tables.tex dq_tables.tex; do
    if [[ -f "$PLOTBASE/pot/$t" ]]; then
      cp "$PLOTBASE/pot/$t" "$TECHNOTE/Sections/$t" && echo "[install] pot/$t -> Sections/$t"
    else
      echo "[install] MISSING $PLOTBASE/pot/$t"
    fi
  done
  # figures the note no longer references (renamed stages/variables) are listed,
  # not deleted -- prune them by hand with `git rm` in the technote checkout:
  python3 "$HERE/technote_figures.py" unused --technote "$TECHNOTE" --keep-commented \
      > "$PLOTBASE/logs/technote-unused.txt" 2>&1
  echo "[info] $(grep -vc '^#' "$PLOTBASE/logs/technote-unused.txt") unreferenced image(s) listed in $PLOTBASE/logs/technote-unused.txt"
fi

# --- compile -----------------------------------------------------------------
if want compile && [[ $DO_COMPILE -eq 1 && $STOP -eq 0 ]]; then
  OUT="$PLOTBASE/logs/tex"; mkdir -p "$OUT"
  ( cd "$TECHNOTE" && latexmk -pdf -bibtex -interaction=nonstopmode -outdir="$OUT" main.tex \
      > "$PLOTBASE/logs/technote-compile.log" 2>&1 )
  rc=$?
  record "compile" "$rc"
  echo "[compile] exit $rc; errors: $(grep -c '^!' "$OUT/main.log" 2>/dev/null || true)"
  grep -h "Reference .* undefined\|Output written" "$OUT/main.log" 2>/dev/null | sort -u | sed 's/^/[compile] /'
fi

echo
echo "============================================================"
echo "TECHNOTE PLOT SUMMARY  (plotbase $PLOTBASE)"
fail=0
for i in "${!RESULT_NAMES[@]}"; do
  rc=${RESULT_CODES[$i]}
  printf "  %-34s %s\n" "${RESULT_NAMES[$i]}" "$([[ $rc -eq 0 ]] && echo OK || echo "FAILED (exit $rc)")"
  [[ $rc -ne 0 ]] && fail=1
done
echo "============================================================"
exit $fail
