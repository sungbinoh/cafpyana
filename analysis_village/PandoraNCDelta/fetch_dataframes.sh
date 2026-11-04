#!/usr/bin/env bash
# Copy NCDelta .df files from a Fermilab GPVM (where /pnfs is mounted) to a
# local directory. Needs a valid Kerberos ticket: kinit <user>@FNAL.GOV
#
# Usage: ./fetch_dataframes.sh [remote_dir] [local_dir]
#   REMOTE_HOST / REMOTE_USER environment variables override the defaults.
set -euo pipefail

REMOTE_HOST="${REMOTE_HOST:-sbndgpvm01.fnal.gov}"
REMOTE_USER="${REMOTE_USER:-nrowe}"
REMOTE_DIR="${1:-/pnfs/sbnd/scratch/users/nrowe/cafpyana_out/dfs/2026_10_17_164529__ncdelta.df/}"
REMOTE_DIR="${REMOTE_DIR%/}"
LOCAL_DIR="${2:-data/$(basename "$REMOTE_DIR")}"

mkdir -p "$LOCAL_DIR"
# Only the .df files; the grid logs next to them are not needed.
rsync -av --progress --include='*.df' --exclude='*' \
    "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_DIR}/" "$LOCAL_DIR/"
echo "Copied $(ls "$LOCAL_DIR"/*.df | wc -l) file(s) to $LOCAL_DIR"
