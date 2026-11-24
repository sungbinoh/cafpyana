#!/bin/bash
# xrootd URLs of every .root file under a /pnfs directory, found with `find` (no SAM), one per line, for
# dump_waveforms.py -i <list>. Much faster than looking files up in SAM one by one.
#
#   find_reco1_urls.sh <top directory | a URL or path of one file in it> <output list> [name pattern]
#
#   <top>     a /pnfs directory (e.g. .../Gen2_2026/CV/reco1), or a root:// URL / /pnfs path of one
#             file: the directory is then taken up to and including its "reco1" component
#   pattern   find -name pattern (default '*.root')
# The top directory's subdirectories are searched in parallel (JOBS, default 8). /pnfs/<area>/... becomes
# root://fndcadoor.fnal.gov:1094//pnfs/fnal.gov/usr/<area>/... (both /pnfs/<area> and
# /pnfs/fnal.gov/usr/<area> inputs work). Files outside the SAM definition (e.g. from failed jobs) are
# included: the directory, not the definition, decides.
#
#   MC:   find_reco1_urls.sh /pnfs/sbn/data_add/sbn_nd/aurora/mc/v10_14_02_03/prodgenie_corsika_proton_rockbox0p1_sbnd/Gen2_2026/CV/reco1 reco1_urls_mc_all.txt
#   data: find_reco1_urls.sh "$(head -1 reco1_urls_data.txt)" reco1_urls_data_all.txt
set -u
if [ $# -lt 2 ]; then sed -n '2,19p' "$0"; exit 1; fi
top=$1; out=$2; pat=${3:-*.root}
JOBS=${JOBS:-8}

# a URL or file path -> its local /pnfs directory, cut after the "reco1" component
top=$(echo "$top" | sed -e 's|^root://[^/]*/*|/|' -e 's|^/*pnfs/fnal.gov/usr/|/pnfs/|')
if [ -f "$top" ] || [[ "$top" == *.root ]]; then
    case "$top" in
        */reco1/*) top="${top%%/reco1/*}/reco1" ;;
        *) top=$(dirname "$top") ;;
    esac
fi
if [ ! -d "$top" ]; then echo "not a directory: $top (is /pnfs mounted here?)" >&2; exit 1; fi
echo "searching $top for '$pat' with $JOBS parallel finds ..." >&2

tmp=$(mktemp -d)
find "$top" -maxdepth 1 -type f -name "$pat" > "$tmp/top"                # files directly in the top
# every subdirectory in parallel, each find into its own file (no interleaved lines)
find "$top" -mindepth 1 -maxdepth 1 -type d -print0 |
    xargs -0 -P "$JOBS" -I{} sh -c 'find "$1" -type f -name "$2" > "$3/$(basename "$1").$$"' _ {} "$pat" "$tmp"
cat "$tmp"/* | sed -e 's|^/pnfs/fnal.gov/usr/|/pnfs/|' -e 's|^/pnfs/|root://fndcadoor.fnal.gov:1094//pnfs/fnal.gov/usr/|' |
    sort -u > "$out"
rm -rf "$tmp"
echo "$(wc -l < "$out") files -> $out" >&2
