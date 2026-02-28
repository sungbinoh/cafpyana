prefix='/exp/sbnd/data/users/nrowe/GUMP/sbn-rewgted-5.3/'

for i in {0..19}
do
    echo Running step $i
    root -q '../PROfit/MakesBruce.C("'$prefix'SBND_SpringMC_rewgt_'$i'_ttrees.root", "'$prefix'SBND_SpringMC_rewgt_'$i'_sbruce.root")'
done

root -q '../PROfit/MakesBruce.C("'$prefix'ICARUS_SpringMCOverlay_rewgt_ttrees.root", "'$prefix'ICARUS_SpringMCOverlay_rewgt_sbruce.root")'
root -q '../PROfit/MakesBruce.C("'$prefix'ICARUS_SpringRun2BNBOff_unblind_prescaled_ttrees.root", "'$prefix'ICARUS_SpringRun2BNBOff_unblind_prescaled_sbruce.root")'
root -q '../PROfit/MakesBruce.C("'$prefix'SBND_SpringBNBOffData_5000_ttrees.root", "'$prefix'SBND_SpringBNBOffData_5000_sbruce.root")'
root -q '../PROfit/MakesBruce.C("'$prefix'ICARUS_SpringMCDirt_slimwgt_ttrees.root", "'$prefix'ICARUS_SpringMCDirt_slimwgt_sbruce.root")'
root -q '../PROfit/MakesBruce.C("'$prefix'SBND_SpringLowEMC_ttrees.root", "'$prefix'SBND_SpringLowEMC_sbruce.root")'
