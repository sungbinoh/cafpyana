python3 render_config.py \
    --hdf-dir /exp/sbnd/data/users/gputnam/GUMPLE/sbn-rewgted-24/ \
    -t TruthNumuCCTemplate.xml.j2 \
    -d /exp/sbnd/data/users/nrowe/sbn-rewgted-24-refom-light \
    -o /exp/sbnd/data/users/nrowe/sbn-rewgted-24-refom-light/TRUTHNUMUCC/TruthNumuCC.xml
    
python3 render_config.py \
    --hdf-dir /exp/sbnd/data/users/gputnam/GUMPLE/sbn-rewgted-24/ \
    -t TruthMapleNPTemplate.xml.j2 \
    -d /flare/neutrinoGPU/SBN_PROfit/sBruce/Oct7/sbn-rewgted-24-refom-light \
    -o /exp/sbnd/data/users/nrowe/sbn-rewgted-24-refom-light/TRUTHMAPLENP/TruthMapleNPANL.xml
