#!/bin/bash
# Grid (jobsub) MAPLE dataframe production over ALL available samples.
#
# !! PREREQUISITE: this MAPLE code must be pushed to the GUMP branch of
# !! github.com/gputnam/cafpyana -- grid jobs clone the repo and check out
# !! origin/GUMP at run time
# !! (see bin/grid_executable.sh). DO NOT submit before pushing.
#
# Usage:
#   ./run_grid_all.sh [config] [nfileperjob]
#
#   [config]:      config name in analysis_village/maple/configs (default maple_evt)
#   [nfileperjob]: input files per grid job (default 200)
#
# Samples are taken from every list in production/file-lists/ (regenerate
# with make_filelists.sh).
#
# Run from the cafpyana root AFTER `source setup.sh`, with a valid token.
# Output lands in $CAFPYANA_GRID_OUT_DIR/dfs4/<timestamp>_<name>/ ;
# merge with dfadd.py and/or filter with preselct.py afterwards.

set -e

CONFIG=${1:-maple_evt}
NFILEPERJOB=${2:-200}

HERE="$(cd "$(dirname "$0")" && pwd)"
CAFPYANA=${CAFPYANA_WD:?"source setup.sh in the cafpyana root first"}

LISTS=("$HERE"/file-lists/*.list)
[ -e "${LISTS[0]}" ] || { echo "no file lists in $HERE/file-lists (run make_filelists.sh first)"; exit 1; }

htgettoken -a htvaultprod.fnal.gov -i sbnd

echo "Submitting ${#LISTS[@]} samples with config $CONFIG, $NFILEPERJOB files/job"
for LIST in "${LISTS[@]}"; do
  NAME=$(basename "$LIST" .list)
  echo "== $NAME: $(wc -l < "$LIST") files"
  "$HERE"/run_grid.sh "$NAME" "$CONFIG" "$NFILEPERJOB"
done

echo "All submissions done. Monitor with: jobsub_q --user \$USER"
