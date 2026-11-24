#!/bin/bash
# xrootd URLs of reco1 files from SAM, one per line, for dump_waveforms.py -i <list>.
#
#   make_reco1_url_list.sh <defname> <n files> <output list> [file names]
#
# With a file-names list (e.g. reco1_files.txt from select_waveform_channels.py) those files are looked
# up; otherwise the first <n files> of the definition (`with limit`: never list a big definition whole).
# Each name goes through `samweb get-file-access-url --schema=root`; if that gives nothing, the
# `locate-file` location is turned into an xrootd URL (/pnfs/sbnd/... -> root://fndcadoor.fnal.gov:1094/
# /pnfs/fnal.gov/usr/sbnd/...). Names without any location are listed in <output list>.missing.
#
# MC:   aurora_SBND2026A_gen2_BNBLight_prodgenie_corsika_proton_rockbox0p1_sbnd_CV_v10_14_02_03_reco1_sbnd
# data: data_MCP2025C_Fall25-Run1_BNB_RollingDev_bnblight_v10_14_02_reco1_sbnd
set -u
def=$1; n=$2; out=$3; names=${4:-}
: > "$out"; : > "$out.missing"
if [ -n "$names" ]; then
    list=$(grep -v '^#' "$names")
else
    echo "listing $n files of $def ..." >&2
    list=$(timeout 1800 samweb -e sbnd list-files "defname: $def with limit $n")
fi
for f in $list; do
    url=$(timeout 120 samweb -e sbnd get-file-access-url --schema=root "$f" 2>/dev/null | head -1)
    if [ -z "$url" ]; then
        loc=$(timeout 120 samweb -e sbnd locate-file "$f" 2>/dev/null | head -1)
        if [ -n "$loc" ]; then
            dir=$(echo "$loc" | sed -e 's/^[a-z]*://' -e 's/(.*)$//')
            url="root://fndcadoor.fnal.gov:1094/${dir/\/pnfs\/sbnd//pnfs/fnal.gov/usr/sbnd}/$f"
        fi
    fi
    if [ -n "$url" ]; then echo "$url" >> "$out"; else echo "$f" >> "$out.missing"; fi
done
echo "$(wc -l < "$out") URLs in $out, $(wc -l < "$out.missing") without a location in $out.missing" >&2
