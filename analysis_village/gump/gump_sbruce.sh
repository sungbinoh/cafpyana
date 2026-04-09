prefix='/exp/sbnd/data/users/nrowe/GUMP/sbn-rewgted-6/'

for i in {0..19}
do
    echo Running step $i
    root -q '../PROfit/MakesBruce.C("'$prefix'SBND_SpringMC_rewgt_5_'$i'_ttrees.root", "'$prefix'SBND_SpringMC_rewgt_5_'$i'_sbruce.root")' &
done

### ICARUS Run 4 MC
for i in {0..9}
do
    echo Running step $i
    root -q '../PROfit/MakesBruce.C("'$prefix'ICARUSRun4_SpringMCOverlay_rewgt_'$i'_ttrees.root", "'$prefix'ICARUSRun4_SpringMCOverlay_rewgt_'$i'_sbruce.root")' &
done

### ICARUS Run 2 MC
root -q '../PROfit/MakesBruce.C("'$prefix'ICARUSRun2_SpringMCOverlay_rewgt_ttrees.root", "'$prefix'ICARUSRun2_SpringMCOverlay_rewgt_sbruce.root")' &

### ICARUS Run 2 Offbeam
root -q '../PROfit/MakesBruce.C("'$prefix'ICARUSRun2_SpringRun2BNBOff_unblind_prescaled_ttrees.root", "'$prefix'ICARUSRun2_SpringRun2BNBOff_unblind_prescaled_sbruce.root")' &

### ICARUS Run 4 Offbeam
root -q '../PROfit/MakesBruce.C("'$prefix'ICARUSRun4_SpringRun2BNBOff_unblind_prescaled_ttrees.root", "'$prefix'ICARUSRun4_SpringRun2BNBOff_unblind_prescaled_sbruce.root")' &

### SBND OffBeam
root -q '../PROfit/MakesBruce.C("'$prefix'SBND_SpringBNBOffData_5000_ttrees.root", "'$prefix'SBND_SpringBNBOffData_5000_sbruce.root")' &

### ICARUS Run 2 Dirt
root -q '../PROfit/MakesBruce.C("'$prefix'ICARUSRun2_Spring_Overlay_Dirt_ttrees.root", "'$prefix'ICARUSRun2_Spring_Overlay_Dirt_sbruce.root")' &

### ICARUS Run 4 Dirt
root -q '../PROfit/MakesBruce.C("'$prefix'ICARUSRun4_Spring_Overlay_Dirt_ttrees.root", "'$prefix'ICARUSRun4_Spring_Overlay_Dirt_sbruce.root")' &

### SBND Dirt
root -q '../PROfit/MakesBruce.C("'$prefix'SBND_SpringLowEMC_ttrees.root", "'$prefix'SBND_SpringLowEMC_sbruce.root")' &
