#!/bin/bash
### One invocation: every sample's POT is read once and reused for all 18 configs.
### -t, -d and -o are positionally paired, so the three lists must stay in step.

### Date label used in the ANL path, output names and plot tags
LABEL=${LABEL:-Oct19}

HDF=/exp/sbnd/data/users/gputnam/GUMPLE/sbn-rewgted-25/
EXP=/exp/sbnd/data/users/nrowe/sbn-rewgted-25-fresh
ANL=/flare/neutrinoGPU/SBN_PROfit/sBruce/$LABEL/sbn-rewgted-25-fresh
OUT=/exp/sbnd/data/users/nrowe/sbn-rewgted-25-fresh

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
       DataMC/sbnd_bash_rwt_GUMP.xml.j2 \
       DataMC/sbnd_data_bash_GUMP.xml.j2 \
       DataMC/icarus2_bash_rwt_GUMP.xml.j2 \
       DataMC/icarus4_bash_rwt_GUMP.xml.j2 \
       DataMC/icarus_data2_bash_GUMP.xml.j2 \
       DataMC/icarus4_data_bash_GUMP.xml.j2 \
       DataMC/sbnd_bash_rwt_MAPLEMP.xml.j2 \
       DataMC/sbnd_data_bash_MAPLEMP.xml.j2 \
       DataMC/icarus2_bash_rwt_MAPLEMP.xml.j2 \
       DataMC/icarus4_bash_rwt_MAPLEMP.xml.j2 \
       DataMC/icarus_data2_bash_MAPLEMP.xml.j2 \
       DataMC/icarus4_data_bash_MAPLEMP.xml.j2 \
       DataMC/sbnd_bash_rwt_MAPLENP.xml.j2 \
       DataMC/sbnd_data_bash_MAPLENP.xml.j2 \
       DataMC/icarus2_bash_rwt_MAPLENP.xml.j2 \
       DataMC/icarus4_bash_rwt_MAPLENP.xml.j2 \
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
    -o $OUT/MAPLEMP/Gumple$LABEL.xml \
       $OUT/GUMP/Gump$LABEL.xml \
       $OUT/MAPLENP/MapleNP$LABEL.xml \
       $OUT/MAPLEMP/GumpleR4Only$LABEL.xml \
       $OUT/MAPLEMP/GumpleAr25$LABEL.xml \
       $OUT/MAPLEMP/GumpleNoCVWgt$LABEL.xml \
       $OUT/MAPLEMP/GumpleANL$LABEL.xml \
       $OUT/GUMP/GumpANL$LABEL.xml \
       $OUT/MAPLENP/MapleNPANL$LABEL.xml \
       $OUT/MAPLEMP/GumpleANLR4Only$LABEL.xml \
       $OUT/MAPLEMP/GumpleANLAr25$LABEL.xml \
       $OUT/MAPLEMP/GumpleANLNoCVWgt$LABEL.xml \
       $OUT/GUMP/sbnd_bash_rwt_GUMP.xml \
       $OUT/GUMP/sbnd_data_bash_GUMP.xml \
       $OUT/GUMP/icarus2_bash_rwt_GUMP.xml \
       $OUT/GUMP/icarus4_bash_rwt_GUMP.xml \
       $OUT/GUMP/icarus_data2_bash_GUMP.xml \
       $OUT/GUMP/icarus4_data_bash_GUMP.xml \
       $OUT/MAPLEMP/sbnd_bash_rwt_MAPLEMP.xml \
       $OUT/MAPLEMP/sbnd_data_bash_MAPLEMP.xml \
       $OUT/MAPLEMP/icarus2_bash_rwt_MAPLEMP.xml \
       $OUT/MAPLEMP/icarus4_bash_rwt_MAPLEMP.xml \
       $OUT/MAPLEMP/icarus_data2_bash_MAPLEMP.xml \
       $OUT/MAPLEMP/icarus4_data_bash_MAPLEMP.xml \
       $OUT/MAPLENP/sbnd_bash_rwt_MAPLENP.xml \
       $OUT/MAPLENP/sbnd_data_bash_MAPLENP.xml \
       $OUT/MAPLENP/icarus2_bash_rwt_MAPLENP.xml \
       $OUT/MAPLENP/icarus4_bash_rwt_MAPLENP.xml \
       $OUT/MAPLENP/icarus_data2_bash_MAPLENP.xml \
       $OUT/MAPLENP/icarus4_data_bash_MAPLENP.xml
RENDER_RC=$?
if [ $RENDER_RC -ne 0 ]; then
    echo "render_config.py failed (exit $RENDER_RC), skipping plots" >&2
    exit $RENDER_RC
fi

### ---------------------------------------------------------------------------
### PROfit plots for the Gump, Gumple and MapleNP XMLs
### ---------------------------------------------------------------------------
PROFIT=/exp/sbnd/app/users/nrowe/xml_tests/Elephant_Vanishes/build/bin/PROfit

### save python for later...
PY=$(command -v python3); PY_PATH=$PYTHONPATH

### PROfit dependencies (AL9 setup from the PROfit README)
source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
spack load gcc@12.2.0
spack load cmake@3.27.7
spack load boost/iylabz2
spack load hdf5@1.14.3
spack load root@6.28.12
spack load ifdhc@2.7.2

### PROfit writes its cache and plots into the current directory
mkdir -p $OUT/plots && cd $OUT/plots

$PROFIT -x $OUT/GUMP/Gump$LABEL.xml       -t Gump$LABEL    -n 8 --log Gump$LABEL.log    plot
$PROFIT -x $OUT/MAPLEMP/Gumple$LABEL.xml  -t Gumple$LABEL  -n 8 --log Gumple$LABEL.log  plot
$PROFIT -x $OUT/MAPLENP/MapleNP$LABEL.xml -t MapleNP$LABEL -n 8 --log MapleNP$LABEL.log plot

### ---------------------------------------------------------------------------
### Sensitivity surfaces (sin^2 2theta_mumu vs dm^2) and a Gumple/Gump/MapleNP comparison
### ---------------------------------------------------------------------------
$PROFIT -x $OUT/MAPLEMP/Gumple$LABEL.xml  -t Gumple$LABEL  -o surf -n 20 --log Gumple${LABEL}_surf.log  surface
$PROFIT -x $OUT/GUMP/Gump$LABEL.xml       -t Gump$LABEL    -o surf -n 20 --log Gump${LABEL}_surf.log    surface
$PROFIT -x $OUT/MAPLENP/MapleNP$LABEL.xml -t MapleNP$LABEL -o surf -n 20 --log MapleNP${LABEL}_surf.log surface

cd - > /dev/null
PLOTTER=../plotter.py
PYTHONPATH=$PY_PATH $PY $PLOTTER -f $OUT/plots/Gumple${LABEL}_surf_surface.txt $OUT/plots/Gump${LABEL}_surf_surface.txt $OUT/plots/MapleNP${LABEL}_surf_surface.txt \
	                 -l Gumple Gump MapleNP -o $OUT/plots/surface_comparison_$LABEL.pdf

## Post the plots to Slack (slack_upload.sh sits next to this script)
./slack_upload.sh "PROfit plots for $LABEL" \
    $OUT/plots/Gumple${LABEL}_v1_PROplot_Variable_0_ErrorBand.pdf \
    $OUT/plots/Gump${LABEL}_v1_PROplot_Variable_0_ErrorBand.pdf \
    $OUT/plots/MapleNP${LABEL}_v1_PROplot_Variable_0_ErrorBand.pdf \
    $OUT/plots/surface_comparison_$LABEL.pdf
