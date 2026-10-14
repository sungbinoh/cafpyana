#!/usr/bin/env bash
# make_technote_plots.sh -- regenerate every locally reproducible figure of the
# GUMP technote from a GUMPLE production and install them into the Overleaf
# checkout. See README.md (sibling) for what each step produces and where the
# figures land in the note.
#
#   usage: make_technote_plots.sh [--df-dir DIR] [--plotbase DIR] [--technote DIR]
#                                 [--only STEP ...] [--no-install] [--no-compile] [--no-check]
#
#   steps (in order): tracksplit mcdata constraint pid signalbox corsika selection check install compile
#     tracksplit  TrackSplittingCorrection_GUMPLE.py           -> <plotbase>/tracksplit/
#     mcdata      mcdata_comparison_gumple.py -d all           -> <plotbase>/evtsel/
#     constraint  mcdata_constraint_gumple.py (SBND abs-norm,  -> <plotbase>/constraint/
#                 CC-QE-energy constraint, incl. delta-p slices)
#     pid         nb/MCDataComparisonPID-GUMPLE.ipynb (x3 det) -> <plotbase>/pid/nominal/
#     signalbox   nb/SignalBoxSystematics-GUMPLE.ipynb         -> <plotbase>/signalbox/
#     corsika     offbeam_cosmicmc_gumple.py                   -> <plotbase>/corsika/
#     selection   selection_plots_gumple.py                    -> <plotbase>/selection/
#     check       technote_figures.py check   (manifest sources exist, dests referenced)
#     install     technote_figures.py install (copy into the technote checkout)
#     compile     latexmk the technote into <plotbase>/logs/tex/ and report errors
#
# Analysis steps go through run_gumple_serial.py (same interpreter, cwd, env
# passthrough and per-step logs as the regular GUMPLE chain); a failing step is
# reported at the end but does not stop the later ones. Nothing is committed.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GUMP_DIR="$(dirname "$HERE")"
REPO="$(cd "$GUMP_DIR/../.." && pwd)"
PY=python3

DF_DIR="${DF_DIR:-/exp/sbnd/data/users/gputnam/GUMPLE/sbn-rewgted-21/}"
PLOTBASE="${PLOTBASE:-$REPO/plots-gumple-2026-09-02-rewgted21}"
TECHNOTE="${TECHNOTE:-$REPO/6973acd1e6f673f7a8b495e7}"
MANIFEST="$HERE/figure_manifest.tsv"
ONLY=()
DO_INSTALL=1
DO_COMPILE=1
DO_CHECK=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --df-dir)     DF_DIR="$2"; shift 2 ;;
    --plotbase)   PLOTBASE="$2"; shift 2 ;;
    --technote)   TECHNOTE="$2"; shift 2 ;;
    --manifest)   MANIFEST="$2"; shift 2 ;;
    --only)       ONLY+=("$2"); shift 2 ;;
    --no-install) DO_INSTALL=0; shift ;;
    --no-compile) DO_COMPILE=0; shift ;;
    --no-check)   DO_CHECK=0; shift ;;
    -h|--help)    sed -n '2,25p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

ALL_STEPS=(tracksplit mcdata constraint pid signalbox corsika selection check install compile)
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

declare -a RESULT_NAMES=() RESULT_CODES=()
record() { RESULT_NAMES+=("$1"); RESULT_CODES+=("$2"); }

# --- analysis steps, via the serial runner ---------------------------------
# The PID step is the technote configuration: the single 'nominal' chi2 flavor
# stored by sbn-rewgted-21, the dE/dx-smearing scan (0x/1x/2x) on the
# angle-inclusive plots, angular quintiles for SBND and ICARUS Run 2 only (Run 4
# rides the prescaled stream), at the 'simplecosrej_trkqual_nop' cut stage.
SERIAL_STEPS=()
for s in tracksplit mcdata constraint pid signalbox corsika selection; do
  want "$s" && SERIAL_STEPS+=(--only "$s")
done
if [[ ${#SERIAL_STEPS[@]} -gt 0 ]]; then
  ( cd "$GUMP_DIR" && "$PY" run_gumple_serial.py \
      --df-dir "$DF_DIR" --skip-wait --plotbase "$PLOTBASE" \
      "${SERIAL_STEPS[@]}" \
      --chi2-flavor nominal --calo-model smearscan \
      --no-angles "ICARUS Run4" --cut-stages simplecosrej_trkqual_nop \
      2>&1 | tee "$PLOTBASE/logs/technote-serial.log" )
  record "analysis (run_gumple_serial.py)" "${PIPESTATUS[0]}"
fi

# --- manifest check / install ------------------------------------------------
if want check && [[ $DO_CHECK -eq 1 ]]; then
  python3 "$HERE/technote_figures.py" check --technote "$TECHNOTE" \
      --plotbase "$PLOTBASE" --manifest "$MANIFEST" 2>&1 | tee "$PLOTBASE/logs/technote-check.log"
  record "check" "${PIPESTATUS[0]}"
fi

if want install && [[ $DO_INSTALL -eq 1 ]]; then
  python3 "$HERE/technote_figures.py" install --technote "$TECHNOTE" \
      --plotbase "$PLOTBASE" --manifest "$MANIFEST" 2>&1 | tee "$PLOTBASE/logs/technote-install.log"
  record "install" "${PIPESTATUS[0]}"
  # figures the note no longer references (renamed stages/variables) are listed,
  # not deleted -- prune them by hand with `git rm` in the technote checkout:
  python3 "$HERE/technote_figures.py" unused --technote "$TECHNOTE" --keep-commented \
      > "$PLOTBASE/logs/technote-unused.txt" 2>&1
  echo "[info] $(grep -vc '^#' "$PLOTBASE/logs/technote-unused.txt") unreferenced image(s) listed in $PLOTBASE/logs/technote-unused.txt"
fi

# --- compile -----------------------------------------------------------------
if want compile && [[ $DO_COMPILE -eq 1 ]]; then
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
