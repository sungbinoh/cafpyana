#!/bin/bash
# Grid (jobsub) MAPLE dataframe production for one sample.
#
# !! PREREQUISITE: this MAPLE code must be pushed to the GUMP branch of
# !! github.com/gputnam/cafpyana -- grid jobs clone the repo and check out
# !! origin/GUMP at run time
# !! (see bin/grid_executable.sh). DO NOT submit before pushing.
#
# Usage:
#   ./run_grid.sh <sample> [config] [nfileperjob]
#
#   <sample>:      name of a list in production/file-lists/ (e.g. Run4_MC)
#   [config]:      config name in analysis_village/maple/configs (default maple_evt)
#   [nfileperjob]: input files per grid job (default 200)
#
# Run from the cafpyana root AFTER `source setup.sh`, with a valid token
# (htgettoken -a htvaultprod.fnal.gov -i sbnd).
# Output lands in $CAFPYANA_GRID_OUT_DIR/dfs4/<timestamp>_<name>/ ;
# merge with dfadd.py and/or filter with preselct.py afterwards.

set -e

SAMPLE=${1:?"usage: $0 <sample> [config] [nfileperjob]"}
CONFIG=${2:-maple_evt}
NFILEPERJOB=${3:-200}

HERE="$(cd "$(dirname "$0")" && pwd)"
CAFPYANA=${CAFPYANA_WD:?"source setup.sh in the cafpyana root first"}

if [ -f "$SAMPLE" ]; then
  LIST=$SAMPLE
  NAME=$(basename "$SAMPLE" .list)
else
  LIST=$HERE/file-lists/$SAMPLE.list
  NAME=$SAMPLE
fi
[ -f "$LIST" ] || { echo "file list $LIST not found (run make_filelists.sh?)"; exit 1; }

NFILES=$(wc -l < "$LIST")
NGRID=$(( (NFILES + NFILEPERJOB - 1) / NFILEPERJOB ))

cd "$CAFPYANA"
htgettoken -a htvaultprod.fnal.gov -i sbnd

echo "Submitting $NAME with $CONFIG: $NFILES files in $NGRID jobs"
python run_df_maker.py \
  -c analysis_village/maple/configs/$CONFIG.py \
  -l "$LIST" \
  -o "${NAME}_${CONFIG}" \
  -ngrid $NGRID
