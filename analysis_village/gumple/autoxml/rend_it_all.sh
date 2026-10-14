#!/bin/bash
### One invocation: every sample's POT is read once and reused for all 18 configs.
### -t, -d and -o are positionally paired, so the three lists must stay in step.

HDF=/exp/sbnd/data/users/gputnam/GUMPLE/sbn-rewgted-25/
EXP=/exp/sbnd/data/users/nrowe/sbn-rewgted-25
ANL=/flare/neutrinoGPU/SBN_PROfit/sBruce/Oct13/sbn-rewgted-25
OUT=/exp/sbnd/data/users/nrowe/sbn-rewgted-25

python3 render_config.py --hdf-dir $HDF \
    -t GumpleTemplate.xml.j2 \
       GumpTemplate.xml.j2 \
       MapleNPTemplate.xml.j2 \
       GumpleR4Only_mc.xml.j2 \
       GumpleAr25_data.xml.j2 \
       GumpleNoCVWgtTemplate.xml.j2 \
       GumpleTemplate.xml.j2 \
       GumpTemplate.xml.j2 \
       MapleNPTemplate.xml.j2 \
       GumpleR4Only_mc.xml.j2 \
       GumpleAr25_data.xml.j2 \
       GumpleNoCVWgtTemplate.xml.j2 \
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
       $EXP \
       $EXP \
       $EXP \
       $ANL \
       $ANL \
       $ANL \
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
    -o $OUT/MAPLEMP/GumpleOct13.xml \
       $OUT/GUMP/GumpOct13.xml \
       $OUT/MAPLENP/MapleNPOct13.xml \
       $OUT/MAPLEMP/GumpleR4OnlyOct13.xml \
       $OUT/MAPLEMP/GumpleAr25Oct13.xml \
       $OUT/MAPLEMP/GumpleNoCVWgtOct13.xml \
       $OUT/MAPLEMP/GumpleANLOct13.xml \
       $OUT/GUMP/GumpANLOct13.xml \
       $OUT/MAPLENP/MapleNPANLOct13.xml \
       $OUT/MAPLEMP/GumpleANLR4OnlyOct13.xml \
       $OUT/MAPLEMP/GumpleANLAr25Oct13.xml \
       $OUT/MAPLEMP/GumpleANLNoCVWgtOct13.xml \
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
