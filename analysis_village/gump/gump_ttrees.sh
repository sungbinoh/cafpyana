gray_prefix='/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-6/'
nate_prefix='/exp/sbnd/data/users/nrowe/GUMP/sbn-rewgted-6/'
output='/exp/sbnd/data/users/nrowe/GUMP/sbn-rewgted-6/'

cd $CAFPYANA_WD 

### SBND MC
for i in {0..19}
do
    echo Running step $i
    python3 run_ttree_maker.py -c analysis_village/gump/configs/gump_ttree_mc.py -i ${gray_prefix}SBND_SpringMC_rewgt_5_${i}.df -o ${output}SBND_SpringMC_rewgt_5_${i}_ttrees.root &
done

### ICARUS Run 4 MC
for i in {0..9}
do
    echo Running step $i
    python3 run_ttree_maker.py -c analysis_village/gump/configs/gump_ttree_mc.py -i ${gray_prefix}ICARUSRun4_SpringMCOverlay_rewgt_${i}.df -o ${output}ICARUSRun4_SpringMCOverlay_rewgt_${i}_ttrees.root &
done

### ICARUS Run 2 MC
python3 run_ttree_maker.py -c analysis_village/gump/configs/gump_ttree_mc.py -i  ${nate_prefix}ICARUSRun2_SpringMCOverlay_rewgt.df -o ${output}ICARUSRun2_SpringMCOverlay_rewgt_ttrees.root &

### ICARUS Run 2 OffBeam
python3 run_ttree_maker.py -c analysis_village/gump/configs/gump_ttree_data.py -i  ${nate_prefix}ICARUSRun2_SpringRun2BNBOff_unblind_prescaled.df -o ${output}ICARUSRun2_SpringRun2BNBOff_unblind_prescaled_ttrees.root &

### ICARUS Run 4 OffBeam
python3 run_ttree_maker.py -c analysis_village/gump/configs/gump_ttree_data.py -i  ${nate_prefix}ICARUSRun4_SpringRun2BNBOff_unblind_prescaled.df -o ${output}ICARUSRun4_SpringRun2BNBOff_unblind_prescaled_ttrees.root &

### SBND OffBeam
python3 run_ttree_maker.py -c analysis_village/gump/configs/gump_ttree_data.py -i  ${gray_prefix}SBND_SpringBNBOffData_5000.df -o ${output}SBND_SpringBNBOffData_5000_ttrees.root &

### ICARUS Run 2 Dirt
python3 run_ttree_maker.py -c analysis_village/gump/configs/gump_ttree_data.py -i ${nate_prefix}ICARUSRun2_Spring_Overlay_Dirt.df -o ${output}ICARUSRun2_Spring_Overlay_Dirt_ttrees.root &

### ICARUS Run 4 Dirt
python3 run_ttree_maker.py -c analysis_village/gump/configs/gump_ttree_data.py -i ${nate_prefix}ICARUSRun4_Spring_Overlay_Dirt.df -o ${output}ICARUSRun4_Spring_Overlay_Dirt_ttrees.root &

### SBND Dirt
python3 run_ttree_maker.py -c analysis_village/gump/configs/gump_ttree_data.py -i ${gray_prefix}SBND_SpringLowEMC.df -o ${output}SBND_SpringLowEMC_ttrees.root & 

cd analysis_village/gump
