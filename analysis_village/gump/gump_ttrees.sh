prefix='/exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-5/'
output='/exp/sbnd/data/users/nrowe/GUMP/sbn-rewgted-5.3/'

for i in {0..19}
do
    echo Running step $i
    python3 run_ttree_maker.py -c analysis_village/gump/configs/gump_ttree_mc.py -i ${prefix}SBND_SpringMC_rewgt_${i}.df -o ${output}SBND_SpringMC_rewgt_${i}_ttrees.root
done

python3 run_ttree_maker.py -c analysis_village/gump/configs/gump_ttree_mc.py -i ${prefix}ICARUS_SpringMCOverlay_rewgt.df -o ${output}ICARUS_SpringMCOverlay_rewgt_ttrees.root 

python3 run_ttree_maker.py -c analysis_village/gump/configs/gump_ttree_data.py -i  ${prefix}ICARUS_SpringRun2BNBOff_unblind_prescaled.df -o ${output}ICARUS_SpringRun2BNBOff_unblind_prescaled_ttrees.root

python3 run_ttree_maker.py -c analysis_village/gump/configs/gump_ttree_data.py -i  ${prefix}SBND_SpringBNBOffData_5000.df -o ${output}SBND_SpringBNBOffData_5000_ttrees.root

python3 run_ttree_maker.py -c analysis_village/gump/configs/gump_ttree_data.py -i ${prefix}ICARUS_SpringMCDirt_slimwgt.df -o ${output}ICARUS_SpringMCDirt_slimwgt_ttrees.root 

python3 run_ttree_maker.py -c analysis_village/gump/configs/gump_ttree_data.py -i ${prefix}SBND_SpringLowEMC.df -o ${output}SBND_SpringLowEMC_ttrees.root 
