#!/bin/bash

# Define the absolute input storage directories
gray_prefix='/exp/sbnd/data/users/gputnam/GUMPLE/sbn-rewgted-21/'
gumple_prefix='../gumple/'
output='/exp/sbnd/data/users/nrowe/MAPLE/sbn-rewgted-21/'
MAX_JOBS=8

# Navigate to the working directory context
echo "========================================================"
echo " Starting GUMP TTree Processing Batch Run...            "
echo "========================================================"

echo "Remaking det var maps..."
selection="gmpl.all_maplemp_cuts"
splinedir="${selection#*.}"

python3 ${gumple_prefix}rwt_map.py -s ${selection} -o ${splinedir} -d ${gray_prefix} -b "1D"

### 1. SBND MC (20 files, 0 to 9)
echo "--> Staging SBND Spring MC Files..."
for file in ${gray_prefix}SBNDMCCV_*.df
do
    while [ $(jobs -rp | wc -l) -ge $MAX_JOBS ]; do
        sleep 10 # Check every 2 seconds
    done

    filename=$(basename "$file")
    i="${filename#SBNDMCCV_}"
    i="${i%.df}"

    echo "Launching SBND MC Step $i"
    python3 ${gumple_prefix}/run_gumple_pipeline.py \
        -c mc \
        -w \
        -f ${splinedir} \
	-s ${selection} \
        -i ${gray_prefix}SBNDMCCV_${i}.df \
        -o ${output}SBNDMCCV_${i}_sbruce.root &
done

### 4. ICARUS Run 2 OffBeam
echo "--> Launching ICARUS Run 2 OffBeam Data..."
python3 ${gumple_prefix}/run_gumple_pipeline.py \
    -c data \
    -f ${splinedir} \
    -s ${selection} \
    -i ${gray_prefix}ICARUS_SpringRun2BNBOff_unblind.df \
    -o ${output}ICARUS_SpringRun2BNBOff_unblind_sbruce.root &

while [ $(jobs -rp | wc -l) -ge $MAX_JOBS ]; do
    sleep 10 # Check every 2 seconds
done

### 5. ICARUS Run 4 OffBeam
echo "--> Launching ICARUS Run 4 OffBeam Data..."
python3 ${gumple_prefix}/run_gumple_pipeline.py \
    -c data \
    -s ${selection} \
    -i ${gray_prefix}ICARUS_SpringRun4BNBOff_unblind.df \
    -o ${output}ICARUS_SpringRun4BNBOff_unblind_sbruce.root &

while [ $(jobs -rp | wc -l) -ge $MAX_JOBS ]; do
    sleep 10 # Check every 2 seconds
done
### 2. ICARUS Run 4 MC (10 files, 0 to 9)
echo "--> Staging ICARUS Run 4 MC Files..."
for file in ${gray_prefix}ICARUSRun4_SpringMCOverlay_rewgt_*.df
do
    while [ $(jobs -rp | wc -l) -ge $MAX_JOBS ]; do
        sleep 10 # Check every 2 seconds
    done

    filename=$(basename "$file")
    i="${filename#ICARUSRun4_SpringMCOverlay_rewgt_}"
    i="${i%.df}"

    echo "Launching ICARUS Run 4 MC Step $i"
    python3 ${gumple_prefix}/run_gumple_pipeline.py \
        -c mc \
        -w \
        -f ${splinedir} \
	-s ${selection} \
        -i ${gray_prefix}ICARUSRun4_SpringMCOverlay_rewgt_${i}.df \
        -o ${output}ICARUSRun4_SpringMCOverlay_rewgt_${i}_sbruce.root &
done

while [ $(jobs -rp | wc -l) -ge $MAX_JOBS ]; do
    sleep 10 # Check every 2 seconds
done

### 8. ICARUS Run 4 Dirt
echo "--> Launching ICARUS Run 4 Dirt..."
python3 ${gumple_prefix}/run_gumple_pipeline.py \
    -c data \
    -s ${selection} \
    -i ${gray_prefix}ICARUSRun4_Spring_Overlay_Dirt.df \
    -o ${output}ICARUSRun4_Spring_Overlay_Dirt_sbruce.root &

while [ $(jobs -rp | wc -l) -ge $MAX_JOBS ]; do
    sleep 10 # Check every 2 seconds
done

### 9. SBND Dirt
echo "--> Launching SBND Dirt..."
python3 ${gumple_prefix}/run_gumple_pipeline.py \
    -c mc \
    -s ${selection} \
    -i ${gray_prefix}SBND_SpringLowEMC.df \
    -o ${output}SBND_SpringLowEMC_sbruce.root &

### 2. ICARUS Run 2 MC (files, 0 to 3)
echo "--> Staging ICARUS Run 2 MC Files..."
for file in ${gray_prefix}ICARUSRun2_SpringMCOverlay_rewgt_*.df
do
    while [ $(jobs -rp | wc -l) -ge $MAX_JOBS ]; do
        sleep 10 # Check every 2 seconds
    done

    filename=$(basename "$file")
    i="${filename#ICARUSRun2_SpringMCOverlay_rewgt_}"
    i="${i%.df}"

    echo "Launching ICARUS Run 4 MC Step $i"
    python3 ${gumple_prefix}/run_gumple_pipeline.py \
        -c mc \
        -w \
        -f ${splinedir} \
	-s ${selection} \
        -i ${gray_prefix}ICARUSRun2_SpringMCOverlay_rewgt_${i}.df \
        -o ${output}ICARUSRun2_SpringMCOverlay_rewgt_${i}_sbruce.root &
done

while [ $(jobs -rp | wc -l) -ge $MAX_JOBS ]; do
    sleep 10 # Check every 2 seconds
done

### 6. SBND OffBeam
echo "--> Launching SBND OffBeam Data..."
python3 ${gumple_prefix}/run_gumple_pipeline.py \
    -c data \
    -s ${selection} \
    -i ${gray_prefix}SBND_SpringBNBOffData.df \
    -o ${output}SBND_SpringBNBOffData_sbruce.root &

while [ $(jobs -rp | wc -l) -ge $MAX_JOBS ]; do
    sleep 10 # Check every 2 seconds
done

### 7. ICARUS Run 2 Dirt
echo "--> Launching ICARUS Run 2 Dirt..."
python3 ${gumple_prefix}/run_gumple_pipeline.py \
    -c mc \
    -s ${selection} \
    -i ${gray_prefix}ICARUSRun2_Spring_Overlay_Dirt.df \
    -o ${output}ICARUSRun2_Spring_Overlay_Dirt_sbruce.root &

while [ $(jobs -rp | wc -l) -ge $MAX_JOBS ]; do
    sleep 10 # Check every 2 seconds
done
