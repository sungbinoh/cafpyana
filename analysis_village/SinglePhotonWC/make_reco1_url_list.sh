#!/bin/bash
# xrootd URLs of reco1 files from a SAM definition, one per line, for dump_waveforms.py -i <list>.
#
#   make_reco1_url_list.sh <defname> <n files | all> <output list> [file names]
#
#   <n files>  the first n files of the definition (`with limit n`, quick)
#   all        every file of the definition (list-definition-files: slow for a big definition,
#              ~749k files for the Gen2 MC; the names are kept in <output list>.names)
#   [file names]  a list of names to look up instead (e.g. reco1_files.txt from
#              select_waveform_channels.py); <defname> is then not used
#
# Each name goes through `samweb get-file-access-url --schema=root`; if that gives nothing, the
# `locate-file` location becomes an xrootd URL (/pnfs/<area>/... -> root://fndcadoor.fnal.gov:1094/
# /pnfs/fnal.gov/usr/<area>/...). JOBS lookups run in parallel (default 8). Names without a location go
# to <output list>.missing. Rerunning with the same output resumes: names already in the list are skipped.
#
#   JOBS=8 ./make_reco1_url_list.sh aurora_SBND2026A_gen2_BNBLight_prodgenie_corsika_proton_rockbox0p1_sbnd_CV_v10_14_02_03_reco1_sbnd all reco1_urls_mc_all.txt
#   JOBS=8 ./make_reco1_url_list.sh data_MCP2025C_Fall25-Run1_BNB_RollingDev_bnblight_v10_14_02_reco1_sbnd all reco1_urls_data_all.txt
#   (for days-long runs: nohup ... > log 2>&1 &)
set -u
if [ $# -lt 3 ]; then sed -n '2,22p' "$0"; exit 1; fi
def=$1; n=$2; out=$3; names_in=${4:-}
JOBS=${JOBS:-8}
names=$out.names

# 1. the file names
if [ -n "$names_in" ]; then
    grep -v '^#' "$names_in" | grep '\.root$' > "$names"
elif [ ! -s "$names" ]; then
    if [ "$n" = "all" ]; then
        echo "listing every file of $def (slow) ..." >&2
        samweb -e sbnd list-definition-files "$def" > "$names.tmp" || { echo "list-definition-files failed" >&2; exit 1; }
    else
        echo "listing $n files of $def ..." >&2
        samweb -e sbnd list-files "defname: $def with limit $n" > "$names.tmp" || { echo "list-files failed" >&2; exit 1; }
    fi
    grep '\.root$' "$names.tmp" > "$names"; rm -f "$names.tmp"
else
    echo "reusing the names in $names" >&2
fi
total=$(wc -l < "$names")
if [ "$total" -eq 0 ]; then echo "no .root file names: check the definition name" >&2; exit 1; fi

# 2. names -> URLs, skipping those already resolved (resume)
touch "$out" "$out.missing"
todo=$(mktemp)
awk -F/ 'FILENAME == ARGV[1] {done[$NF] = 1; next} !($0 in done)' "$out" "$names" | grep -vxF -f "$out.missing" > "$todo" || true
echo "$total names, $(wc -l < "$todo") to look up with $JOBS parallel jobs" >&2

lookup() {
    f=$1
    # only a root:// line is a URL (samweb can print an error message on stdout)
    url=$(timeout 120 samweb -e sbnd get-file-access-url --schema=root "$f" 2>/dev/null | grep -m1 '^root://')
    if [ -z "$url" ]; then
        loc=$(timeout 120 samweb -e sbnd locate-file "$f" 2>/dev/null | grep -m1 -E '^[a-z]+:/pnfs/')
        if [ -n "$loc" ]; then
            dir=$(echo "$loc" | sed -e 's/^[a-z]*://' -e 's/(.*)$//')
            url="root://fndcadoor.fnal.gov:1094/$(echo "$dir" | sed 's|^/pnfs/\([^/]*\)/|/pnfs/fnal.gov/usr/\1/|')/$f"
        fi
    fi
    if [ -n "$url" ]; then echo "URL $url"; else echo "MISSING $f"; fi
}
export -f lookup
xargs -d '\n' -P "$JOBS" -I{} bash -c 'lookup "$1"' _ {} < "$todo" | while read -r tag val; do
    if [ "$tag" = URL ]; then echo "$val" >> "$out"; else echo "$val" >> "$out.missing"; fi
done
rm -f "$todo"
echo "$(wc -l < "$out") URLs in $out, $(wc -l < "$out.missing") without a location in $out.missing (of $total)" >&2
