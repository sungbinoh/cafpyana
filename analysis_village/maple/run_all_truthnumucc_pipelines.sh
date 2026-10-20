#!/bin/bash
#
# TRUTH numu-CC pipeline: same structure as run_all_gump_pipelines.sh /
# run_all_maple{np,mp}_pipelines.sh, but the selection is the truth-based
# gmpl.all_truth_numucc_cuts (slices matched to a true numu CC interaction with
# the true vertex in the FV) instead of a reconstructed gump/maple selection.
#
# MC samples only: data (on-beam / off-beam) carries no truth, so it is skipped.
# Dirt is not run either.
#
# Switch the selection by editing `selection` below (or SELECTION=... in the
# environment):
#   gmpl.all_truth_numucc_cuts       numu / numubar CC, true vertex in the FV
#   gmpl.all_truth_numucc_nofv_cuts  numu / numubar CC, no FV requirement
#
# Other overrides (environment): GRAY_PREFIX (input .df directory), OUTPUT
# (output directory), MAX_JOBS.

# The pipeline scripts use paths relative to analysis_village/<dir>/ (../gumple/,
# ../gump/, ...), so always run from this script's own directory.
cd "$(dirname "$(readlink -f "$0")")" || exit 1

gray_prefix="${GRAY_PREFIX:-/exp/sbnd/data/users/gputnam/GUMPLE/sbn-rewgted-24/}"
gumple_prefix='../gumple/'
output="${OUTPUT:-/exp/sbnd/data/users/nrowe/sbn-rewgted-24-refom-light/TRUTHNUMUCC/}"
MAX_JOBS="${MAX_JOBS:-8}"

selection="${SELECTION:-gmpl.all_truth_numucc_cuts}"
splinedir="${selection#*.}"

mkdir -p "${output}"

echo "========================================================"
echo " Starting TRUTH numuCC TTree Processing Batch Run...    "
echo " selection: ${selection}"
echo " input:     ${gray_prefix}"
echo " output:    ${output}"
echo "========================================================"

echo "Remaking det var maps..."
python3 ${gumple_prefix}rwt_map.py -t -s ${selection} -o ${splinedir} --treeoutdir ${output} -d ${gray_prefix} -j 2 -b "2D" \
    || { echo "rwt_map.py failed, aborting"; exit 1; }

pids=()

# Block until fewer than MAX_JOBS background jobs are running.
throttle() {
    while [ "$(jobs -rp | wc -l)" -ge "${MAX_JOBS}" ]; do
        sleep 5
    done
}

# run_mc <with_weights 0|1> <input .df>
#   with_weights=1 adds -w (systematic weights) and -f (detvar splines made above),
#   exactly as the gump/maple scripts do for the nominal MC samples.
#   Output is <output>/<input basename>_sbruce.root.
run_mc() {
    local weights="$1" file="$2"
    local base
    base="$(basename "${file}" .df)"
    [ -e "${file}" ] || { echo "WARNING: ${file} not found, skipping"; return; }
    throttle
    local extra=()
    [ "${weights}" = "1" ] && extra=(-w -f "${splinedir}")
    echo "Launching ${base}"
    python3 ${gumple_prefix}/run_gumple_pipeline.py \
        -c mc \
        "${extra[@]}" \
        -s ${selection} \
        -i "${file}" \
        -o "${output}${base}_sbruce.root" &
    pids+=($!)
}

### Nominal MC (with weights)
echo "--> Staging SBND MC files..."
for file in ${gray_prefix}SBNDMCCV_*.df; do
    run_mc 1 "${file}"
done

echo "--> Staging ICARUS Run 4 MC files..."
for file in ${gray_prefix}ICARUSRun4_SpringMCOverlay_rewgt_*.df; do
    run_mc 1 "${file}"
done

echo "--> Staging ICARUS Run 2 MC files..."
for file in ${gray_prefix}ICARUSRun2_SpringMCOverlay_rewgt_*.df; do
    run_mc 1 "${file}"
done

# NB: the gump/maple scripts end by calling $gumple_prefix/run_all.sh, which makes
# fake-data-weighted copies of the sbruce files. That is not part of this pipeline.

### Wait for everything and report failures (the other scripts would exit 0 regardless)
failed=0
for pid in "${pids[@]}"; do
    wait "${pid}" || failed=$((failed + 1))
done
echo "========================================================"
if [ "${failed}" -eq 0 ]; then
    echo " Done: ${#pids[@]} jobs launched, all succeeded."
else
    echo " Done: ${#pids[@]} jobs launched, ${failed} FAILED."
fi
echo "========================================================"
[ "${failed}" -eq 0 ]
