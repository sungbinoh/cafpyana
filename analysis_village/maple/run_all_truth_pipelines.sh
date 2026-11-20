#!/bin/bash
#
# TRUTH pipeline: same structure as run_all_gump_pipelines.sh / run_all_maple{np,mp}_pipelines.sh,
# but the selection is a truth-based one from gumple_cuts.py instead of a reconstructed gump/maple
# selection.  Default: gmpl.all_truth_maplenp_cuts -- the true "np" signal: a muon and any number
# (>= 1) of protons, nothing else (Reco_class 1 or 2).
#
# MC only (data carries no truth), no dirt, and NO systematic weights by default.
#
# Selection -- edit `selection` below or set SELECTION=... in the environment.  The output directory
# is derived from it, so different selections never share a directory:
#   gmpl.all_truth_maplenp_cuts       muon + any number of protons       -> .../TRUTHMAPLENP/
#   gmpl.all_truth_numucc_cuts        any numu / numubar CC, true FV     -> .../TRUTHNUMUCC/
#   gmpl.all_truth_numucc_nofv_cuts   any numu / numubar CC, no FV cut   -> .../TRUTHNUMUCCNOFV/
#
# Environment overrides:
#   WEIGHTS=1    also make the detector-variation maps (rwt_map.py) and run the nominal MC with -w -f
#                (systematic weights).  Default 0: no weights, no rwt_map step.
#   GRAY_PREFIX  input .df directory          OUTPUT  output directory          MAX_JOBS  parallel jobs

# The pipeline scripts use paths relative to analysis_village/<dir>/ (../gumple/, ../gump/, ...), so
# always run from this script's own directory.
cd "$(dirname "$(readlink -f "$0")")" || exit 1

gray_prefix="${GRAY_PREFIX:-/exp/sbnd/data/users/gputnam/GUMPLE/sbn-rewgted-24/}"
gumple_prefix='../gumple/'
MAX_JOBS="${MAX_JOBS:-8}"
WEIGHTS="${WEIGHTS:-0}"

selection="${SELECTION:-gmpl.all_truth_maplenp_cuts}"
splinedir="${selection#*.}"
case "${splinedir}" in
    all_truth_maplenp_cuts)      outname=TRUTHMAPLENP ;;
    all_truth_numucc_cuts)       outname=TRUTHNUMUCC ;;
    all_truth_numucc_nofv_cuts)  outname=TRUTHNUMUCCNOFV ;;
    *)                           outname="TRUTH_${splinedir}" ;;
esac
output="${OUTPUT:-/exp/sbnd/data/users/nrowe/sbn-rewgted-24-refom-light/${outname}/}"

mkdir -p "${output}" || { echo "cannot create ${output}"; exit 1; }

echo "========================================================"
echo " Starting TRUTH TTree Processing Batch Run...           "
echo " selection: ${selection}"
echo " weights:   ${WEIGHTS}"
echo " input:     ${gray_prefix}"
echo " output:    ${output}"
echo "========================================================"

if [ "${WEIGHTS}" = "1" ]; then
    echo "Remaking det var maps..."
    python3 ${gumple_prefix}rwt_map.py -t -s ${selection} -o ${splinedir} --treeoutdir ${output} -d ${gray_prefix} -j 2 -b "2D" \
        || { echo "rwt_map.py failed, aborting"; exit 1; }
fi

pids=()

# Block until fewer than MAX_JOBS background jobs are running.
throttle() {
    while [ "$(jobs -rp | wc -l)" -ge "${MAX_JOBS}" ]; do
        sleep 5
    done
}

# run_mc <with_weights 0|1> <input .df>
#   with_weights=1 adds -w (systematic weights) and -f (detvar splines made above), as the gump/maple
#   scripts do for the nominal MC samples.  Output is <output>/<input basename>_sbruce.root.
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

### Nominal MC (weights only if WEIGHTS=1)
echo "--> Staging SBND MC files..."
for file in ${gray_prefix}SBNDMCCV_*.df; do
    run_mc "${WEIGHTS}" "${file}"
done

echo "--> Staging ICARUS Run 4 MC files..."
for file in ${gray_prefix}ICARUSRun4_SpringMCOverlay_rewgt_*.df; do
    run_mc "${WEIGHTS}" "${file}"
done

echo "--> Staging ICARUS Run 2 MC files..."
for file in ${gray_prefix}ICARUSRun2_SpringMCOverlay_rewgt_*.df; do
    run_mc "${WEIGHTS}" "${file}"
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
