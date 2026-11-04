#!/bin/bash
# Fetch ROOT files via direct /pnfs/ POSIX path using samweb locate-file

DEST_DIR='/exp/sbnd/data/users/nrowe/shw_failures'
MAX_FILES=50
FILE_COUNT=0
mkdir -p "$DEST_DIR"

DEFNAME='aurora_SBND2026A_gen2_BNBLight_prodgenie_corsika_proton_rockbox0p1_sbnd_CV_v10_14_02_03_reco1_sbnd'

# NC pi0 background slices where the second decay photon is GENUINELY MISSING:
#   - it deposited >60 MeV in the TPC and was contained
#   - conversion points >30 cm apart, opening angle >40 deg (resolvable by eye)
#   - no second shower PFP in the slice, and no matching shower in any other slice
#   - reconstructed shower energy agrees with the LEADING photon alone to 15%,
#     and disagrees with the sum of both by >25% -- so the charge is not merged in
# Top 50 by missed visible energy (65-174 MeV). Format run:subrun:evt
TARGET_ENTRIES=(
    "498:99:26"         # 1g0p balanced slc1, missed gamma 87 MeV, sep 80 cm, 64 deg
    "502:75:16"         # 1g1p balanced slc0, missed gamma 73 MeV, sep 52 cm, 118 deg
    "525:48:6"          # 1g1p high_purity slc1, missed gamma 113 MeV, sep 44 cm, 64 deg
    "534:85:41"         # 1g1p balanced slc0, missed gamma 65 MeV, sep 59 cm, 96 deg
    "730:82:40"         # 1g1p balanced slc1, missed gamma 77 MeV, sep 67 cm, 85 deg
    "738:86:12"         # 1g0p high_purity slc2, missed gamma 67 MeV, sep 40 cm, 84 deg
    "784:25:40"         # 1g0p balanced slc4, missed gamma 77 MeV, sep 101 cm, 62 deg
    "915:14:29"         # 1g0p balanced slc1, missed gamma 71 MeV, sep 43 cm, 69 deg
    "1025:37:34"        # 1g0p high_purity slc1, missed gamma 105 MeV, sep 31 cm, 44 deg
    "1048:75:29"        # 1g1p balanced slc1, missed gamma 73 MeV, sep 51 cm, 89 deg
    "1173:91:18"        # 1g0p balanced slc2, missed gamma 82 MeV, sep 66 cm, 65 deg
    "1198:94:23"        # 1g0p high_purity slc0, missed gamma 174 MeV, sep 48 cm, 45 deg
    "1199:61:24"        # 1g0p high_purity slc0, missed gamma 108 MeV, sep 49 cm, 48 deg
    "1223:24:49"        # 1g1p balanced slc0, missed gamma 78 MeV, sep 61 cm, 67 deg
    "1287:38:11"        # 1g1p balanced slc1, missed gamma 69 MeV, sep 56 cm, 92 deg
    "1301:38:17"        # 1g0p balanced slc3, missed gamma 81 MeV, sep 31 cm, 60 deg
    "1344:42:16"        # 1g0p high_purity slc2, missed gamma 131 MeV, sep 57 cm, 49 deg
    "1385:8:50"         # 1g1p balanced slc1, missed gamma 66 MeV, sep 36 cm, 108 deg
    "1410:70:43"        # 1g0p high_purity slc1, missed gamma 102 MeV, sep 42 cm, 59 deg
    "1554:0:26"         # 1g1p high_purity slc0, missed gamma 136 MeV, sep 31 cm, 53 deg
    "1600:27:15"        # 1g0p balanced slc0, missed gamma 140 MeV, sep 41 cm, 43 deg
    "1803:1:25"         # 1g0p balanced slc2, missed gamma 115 MeV, sep 43 cm, 46 deg
    "1972:1:13"         # 1g1p balanced slc2, missed gamma 72 MeV, sep 106 cm, 106 deg
    "1983:1:31"         # 1g0p high_purity slc2, missed gamma 91 MeV, sep 70 cm, 63 deg
    "1999:19:4"         # 1g0p high_purity slc2, missed gamma 76 MeV, sep 37 cm, 66 deg
    "2004:63:26"        # 1g0p balanced slc2, missed gamma 67 MeV, sep 134 cm, 91 deg
    "2025:4:49"         # 1g1p high_purity slc1, missed gamma 109 MeV, sep 54 cm, 69 deg
    "2096:32:41"        # 1g0p high_purity slc2, missed gamma 98 MeV, sep 42 cm, 56 deg
    "2134:91:11"        # 1g0p high_purity slc2, missed gamma 96 MeV, sep 37 cm, 62 deg
    "2141:83:35"        # 1g0p high_purity slc0, missed gamma 168 MeV, sep 50 cm, 47 deg
    "2194:75:39"        # 1g1p balanced slc1, missed gamma 77 MeV, sep 35 cm, 114 deg
    "2224:38:31"        # 1g1p balanced slc0, missed gamma 73 MeV, sep 33 cm, 125 deg
    "2238:97:12"        # 1g0p balanced slc1, missed gamma 67 MeV, sep 57 cm, 76 deg
    "2254:0:11"         # 1g0p balanced slc2, missed gamma 80 MeV, sep 40 cm, 58 deg
    "2286:12:38"        # 1g0p balanced slc3, missed gamma 123 MeV, sep 58 cm, 45 deg
    "2299:8:26"         # 1g1p high_purity slc1, missed gamma 111 MeV, sep 63 cm, 45 deg
    "2304:97:31"        # 1g0p balanced slc0, missed gamma 73 MeV, sep 66 cm, 69 deg
    "2314:95:28"        # 1g0p balanced slc0, missed gamma 67 MeV, sep 117 cm, 81 deg
    "2405:72:44"        # 1g0p high_purity slc2, missed gamma 73 MeV, sep 122 cm, 64 deg
    "2410:63:14"        # 1g0p high_purity slc0, missed gamma 102 MeV, sep 134 cm, 74 deg
    "2521:12:7"         # 1g1p balanced slc2, missed gamma 92 MeV, sep 69 cm, 48 deg
    "2527:61:23"        # 1g0p balanced slc1, missed gamma 82 MeV, sep 51 cm, 77 deg
    "2576:4:29"         # 1g0p high_purity slc5, missed gamma 65 MeV, sep 48 cm, 94 deg
    "2667:30:32"        # 1g0p balanced slc3, missed gamma 69 MeV, sep 64 cm, 83 deg
    "2670:5:44"         # 1g1p balanced slc1, missed gamma 65 MeV, sep 58 cm, 88 deg
    "2719:93:29"        # 1g1p balanced slc3, missed gamma 83 MeV, sep 35 cm, 83 deg
    "2735:11:20"        # 1g1p balanced slc1, missed gamma 78 MeV, sep 40 cm, 69 deg
    "2760:88:16"        # 1g1p balanced slc2, missed gamma 65 MeV, sep 42 cm, 106 deg
    "2766:56:1"         # 1g0p balanced slc1, missed gamma 140 MeV, sep 46 cm, 57 deg
    "2791:1:49"         # 1g0p balanced slc2, missed gamma 102 MeV, sep 61 cm, 56 deg
)

for entry in "${TARGET_ENTRIES[@]}"; do
    if [ "$FILE_COUNT" -ge "$MAX_FILES" ]; then
        echo "Reached max limit of ${MAX_FILES} files. Stopping download."
        break
    fi
    
    IFS=":" read -r run subrun evts <<< "$entry"
    target="${run}.${subrun}"
    prefix="run${run}_subrun${subrun}_evt${evts}"
    
    echo "=== Processing run.subrun: ${target} (Events: ${evts}) ==="
    FILES=$(samweb -e sbnd list-files "defname:${DEFNAME} and run_number ${target}")
    
    if [ -z "$FILES" ]; then
        echo "WARNING: No files returned for run.subrun ${target}"
        continue
    fi
    
    for fname in $FILES; do
        if [ "$FILE_COUNT" -ge "$MAX_FILES" ]; then
            break
        fi
        
        out_path="${DEST_DIR}/${prefix}_${fname}"
        
        if [ -s "$out_path" ]; then
            echo "  File ${out_path} already exists. Skipping."
            FILE_COUNT=$((FILE_COUNT + 1))
            continue
        fi
        
        # Extract physical /pnfs/ directory using samweb locate-file
        LOC_DIR=$(samweb -e sbnd locate-file "$fname" 2>/dev/null | grep -o "/pnfs/[^:]*" | head -n 1)
        
        if [ -n "$LOC_DIR" ]; then
            pnfs_path="${LOC_DIR}/${fname}"
            echo "  Copying POSIX: ${pnfs_path} -> ${out_path}"
            cp "$pnfs_path" "$out_path"
            
            if [ -s "$out_path" ]; then
                FILE_COUNT=$((FILE_COUNT + 1))
                echo "  SUCCESS: Staged to ${out_path} (Total: ${FILE_COUNT}/${MAX_FILES})"
            else
                echo "  ERROR: cp failed or produced empty file for ${fname}."
                rm -f "$out_path"
            fi
        else
            echo "  ERROR: Could not locate /pnfs directory for ${fname}"
        fi
    done
done

echo "Finished! Total files staged: ${FILE_COUNT}"
