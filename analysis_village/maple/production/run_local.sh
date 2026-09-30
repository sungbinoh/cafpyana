#!/bin/bash
# Local (multiprocessing) MAPLE dataframe production for one sample.
#
# Usage:
#   ./run_local.sh <sample> [config] [ncpu] [outdir]
#
#   <sample>: name of a list in production/file-lists/ (e.g. Run4_MC), or a
#             path to a file list
#   [config]: config name in analysis_village/maple/configs (default maple_evt)
#   [ncpu]:   number of processes (default 8)
#   [outdir]: output directory (default $CAFPYANA_WD/maple-out)
#
# Run from the cafpyana root AFTER `source setup.sh`.

set -e

SAMPLE=${1:?"usage: $0 <sample> [config] [ncpu] [outdir]"}
CONFIG=${2:-maple_evt}
NCPU=${3:-8}

HERE="$(cd "$(dirname "$0")" && pwd)"
CAFPYANA=${CAFPYANA_WD:?"source setup.sh in the cafpyana root first"}
OUTDIR=${4:-$CAFPYANA/maple-out}
mkdir -p "$OUTDIR"

if [ -f "$SAMPLE" ]; then
  LIST=$SAMPLE
  NAME=$(basename "$SAMPLE" .list)
else
  LIST=$HERE/file-lists/$SAMPLE.list
  NAME=$SAMPLE
fi
[ -f "$LIST" ] || { echo "file list $LIST not found (run make_filelists.sh?)"; exit 1; }

cd "$CAFPYANA"
echo "Producing $NAME with $CONFIG on $(wc -l < "$LIST") files, $NCPU processes"
python run_df_maker.py \
  -c analysis_village/maple/configs/$CONFIG.py \
  -l "$LIST" \
  -o "$OUTDIR/${NAME}_${CONFIG}" \
  -ncpu $NCPU
echo "Output: $OUTDIR/${NAME}_${CONFIG}.df"
