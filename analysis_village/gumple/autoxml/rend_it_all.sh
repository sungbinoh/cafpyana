#!/bin/bash
### One invocation: every sample's POT is read once and reused for all 18 configs.
### -t, -d and -o are positionally paired, so the three lists must stay in step.

HDF=/exp/sbnd/data/users/gputnam/GUMPLE/sbn-rewgted-24/
EXP=/exp/sbnd/data/users/nrowe/sbn-rewgted-24-refom
ANL=/flare/neutrinoGPU/SBN_PROfit/sBruce/Oct1/sbn-rewgted-24-refom
OUT=/exp/sbnd/data/users/nrowe/sbn-rewgted-24-refom

python3 render_config.py --hdf-dir $HDF \
    -t GumpleTemplate.xml.j2 \
       GumpTemplate.xml.j2 \
       MapleNPTemplate.xml.j2 \
       GumpleTemplate.xml.j2 \
       GumpTemplate.xml.j2 \
       MapleNPTemplate.xml.j2 \
       DataMC/sbnd_bash_rwt24_GUMP.xml.j2 \
       DataMC/sbnd_data_bash_GUMP.xml.j2 \
       DataMC/icarus2_bash_rwt24_GUMP.xml.j2 \
       DataMC/icarus4_bash_rwt24_GUMP.xml.j2 \
       DataMC/icarus_data2_bash_GUMP.xml.j2 \
       DataMC/icarus4_data_bash_GUMP.xml.j2 \
       DataMC/sbnd_bash_rwt24_MAPLEMP.xml.j2 \
       DataMC/sbnd_data_bash_MAPLEMP.xml.j2 \
       DataMC/icarus2_bash_rwt24_MAPLEMP.xml.j2 \
       DataMC/icarus4_bash_rwt24_MAPLEMP.xml.j2 \
       DataMC/icarus_data2_bash_MAPLEMP.xml.j2 \
       DataMC/icarus4_data_bash_MAPLEMP.xml.j2 \
       DataMC/sbnd_bash_rwt24_MAPLENP.xml.j2 \
       DataMC/sbnd_data_bash_MAPLENP.xml.j2 \
       DataMC/icarus2_bash_rwt24_MAPLENP.xml.j2 \
       DataMC/icarus4_bash_rwt24_MAPLENP.xml.j2 \
       DataMC/icarus_data2_bash_MAPLENP.xml.j2 \
       DataMC/icarus4_data_bash_MAPLENP.xml.j2 \
    -d $EXP \
       $EXP \
       $EXP \
       $ANL \
       $ANL \
       $ANL \
       $EXP \
       $EXP \
       $EXP \
       $EXP \
       $EXP \
       $EXP \
       $EXP \
       $EXP \
       $EXP \
       $EXP \
       $EXP \
       $EXP \
       $EXP \
       $EXP \
       $EXP \
       $EXP \
       $EXP \
       $EXP \
    -o $OUT/MAPLEMP/GumpleOct1.xml \
       $OUT/GUMP/GumpOct1.xml \
       $OUT/MAPLENP/MapleNPOct1.xml \
       $OUT/MAPLEMP/GumpleANLOct1.xml \
       $OUT/GUMP/GumpANLOct1.xml \
       $OUT/MAPLENP/MapleNPANLOct1.xml \
       $OUT/GUMP/sbnd_bash_rwt24_GUMP.xml \
       $OUT/GUMP/sbnd_data_bash_GUMP.xml \
       $OUT/GUMP/icarus2_bash_rwt24_GUMP.xml \
       $OUT/GUMP/icarus4_bash_rwt24_GUMP.xml \
       $OUT/GUMP/icarus_data2_bash_GUMP.xml \
       $OUT/GUMP/icarus4_data_bash_GUMP.xml \
       $OUT/MAPLEMP/sbnd_bash_rwt24_MAPLEMP.xml \
       $OUT/MAPLEMP/sbnd_data_bash_MAPLEMP.xml \
       $OUT/MAPLEMP/icarus2_bash_rwt24_MAPLEMP.xml \
       $OUT/MAPLEMP/icarus4_bash_rwt24_MAPLEMP.xml \
       $OUT/MAPLEMP/icarus_data2_bash_MAPLEMP.xml \
       $OUT/MAPLEMP/icarus4_data_bash_MAPLEMP.xml \
       $OUT/MAPLENP/sbnd_bash_rwt24_MAPLENP.xml \
       $OUT/MAPLENP/sbnd_data_bash_MAPLENP.xml \
       $OUT/MAPLENP/icarus2_bash_rwt24_MAPLENP.xml \
       $OUT/MAPLENP/icarus4_bash_rwt24_MAPLENP.xml \
       $OUT/MAPLENP/icarus_data2_bash_MAPLENP.xml \
       $OUT/MAPLENP/icarus4_data_bash_MAPLENP.xml
