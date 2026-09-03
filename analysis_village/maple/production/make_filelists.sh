#!/bin/bash
# Build file lists for the MAPLE samples reachable from the SBND gpvms.
#
# These are the twester ReCAF2026 flat CAF productions referenced by the
# CAFANA-MAPLE macro (make_tree_eff_sel.C). The macro's other datasets live
# on /pnfs/icarus or CNAF storage, which are NOT mounted on SBND machines --
# see the note at the bottom.

HERE="$(cd "$(dirname "$0")" && pwd)"
OUT=$HERE/file-lists
mkdir -p $OUT

BASE=/pnfs/sbn/scratch/users/twester

declare -A SAMPLES=(
  [Run2_MC]="$BASE/Run2_ReCAF2026/flatcaf"
  [Run4_MC]="$BASE/Run4_ReCAF2026/flatcaf"
  [Run4_Dirt]="$BASE/Run4_Dirt_ReCAF2026/flatcaf"
  [Run2_WM]="$BASE/Run2_WM_ReCAF2026/flatcaf"
  [Run2_WM_YZ]="$BASE/Run2_WM_YZ_ReCAF2026/flatcaf"
)

for name in "${!SAMPLES[@]}"; do
  dir=${SAMPLES[$name]}
  if [ ! -d "$dir" ]; then
    echo "WARNING: $dir not found, skipping $name"
    continue
  fi
  echo "Listing $name from $dir ..."
  find "$dir" -name '*.flat.caf.root' | sort > "$OUT/$name.list"
  echo "  $(wc -l < "$OUT/$name.list") files -> $OUT/$name.list"
done

cat <<'EOF'

NOTE: the following CAFANA-MAPLE datasets are NOT reachable from SBND
machines and are not listed here:
  - /pnfs/icarus/persistent/users/cfarnese/for_Nicola/*         (Run4 concat + offbeam)
  - /storage/gpfs_data/icarus/... (CNAF)
  - SAM definitions (Icaruspro_2025_..._unblind data / offbeam)
Run those from an ICARUS gpvm, or via SAM with the icarus VO.
EOF
