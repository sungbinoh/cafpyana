for i in {1..7}
do
    echo Running step $i
    python3 run_ttree_maker.py -c configs/gump_ttree_mc.py -i /exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-3-fix-zexp-again-again/SBND_SpringMC_rewgt_${i}.df -o /exp/sbnd/data/users/nrowe/GUMP/sbn-rewgted-3-fix-zexp-again-again/SBND_SpringMC_rewgt_${i}_ttrees.root
done


python3 run_ttree_maker.py -c configs/gump_ttree_mc.py -i /exp/sbnd/data/users/gputnam/GUMP/sbn-rewgted-3-fix-zexp-again-again/ICARUS_SpringMC_Dev_rewgt.df -o /exp/sbnd/data/users/nrowe/GUMP/sbn-rewgted-3-fix-zexp-again-again/ICARUS_SpringMC_Dev_rewgt_ttrees.root 

python3 run_ttree_maker.py -c configs/gump_ttree_data.py -i  /exp/sbnd/data/users/gputnam/GUMP/sbn-wgted/ICARUS_Run2_BNBoff_uncalo_prescaled.df -o /exp/sbnd/data/users/nrowe/GUMP/sbn-rewgted-3-fix-zexp-again-again/ICARUS_Run2_BNBoff_uncalo_prescaled_ttrees.root
python3 run_ttree_maker.py -c configs/gump_ttree_data.py -i  /exp/sbnd/data/users/gputnam/GUMP/sbn-wgted/SBND_SPINE_SpringBNBOffData.df -o /exp/sbnd/data/users/nrowe/GUMP/sbn-rewgted-3-fix-zexp-again-again/SBND_SPINE_SpringBNBOffData_ttrees.root

