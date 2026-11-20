#!/bin/bash
# Local MAPLE dataframe production over ALL available samples.
#
# Usage:
#   ./run_local_all.sh [config] [ncpu] [outdir]
#
#   [config]: config name in analysis_village/maple/configs (default maple_evt)
#   [ncpu]:   number of processes per sample (default 8)
#   [outdir]: output directory (default $CAFPYANA_WD/maple-out)
#
# Samples are taken from every list in production/file-lists/ (regenerate
# with make_filelists.sh). A sample is skipped if its output .df already
# exists, so the script can be re-run to resume.
#
# Run from the cafpyana root AFTER `source setup.sh`.

set -e

CONFIG=${1:-maple_evt}
NCPU=${2:-8}

HERE="$(cd "$(dirname "$0")" && pwd)"
CAFPYANA=${CAFPYANA_WD:?"source setup.sh in the cafpyana root first"}
OUTDIR=${3:-$CAFPYANA/maple-out}
mkdir -p "$OUTDIR"

LISTS=("$HERE"/file-lists/*.list)
[ -e "${LISTS[0]}" ] || { echo "no file lists in $HERE/file-lists (run make_filelists.sh first)"; exit 1; }

echo "Producing ${#LISTS[@]} samples with config $CONFIG, $NCPU processes each"
for LIST in "${LISTS[@]}"; do
  NAME=$(basename "$LIST" .list)
  OUT=$OUTDIR/${NAME}_${CONFIG}.df
  if [ -f "$OUT" ]; then
    echo "== $NAME: $OUT exists, skipping"
    continue
  fi
  echo "== $NAME: $(wc -l < "$LIST") files"
  "$HERE"/run_local.sh "$NAME" "$CONFIG" "$NCPU" "$OUTDIR"
done

echo "All done. Outputs in $OUTDIR"
