#!/bin/bash
# Grid (jobsub) MAPLE dataframe production over the GUMP reCAF file-lists,
# one pass per list, with output names matching the GUMP sbn-rewgted-14
# production so the two DF sets compare 1:1.
#
# !! PREREQUISITE: this MAPLE code must be pushed to the GUMP branch of
# !! github.com/gputnam/cafpyana -- grid jobs clone the repo and check out
# !! origin/GUMP at run time (see bin/grid_executable.sh).
# !! DO NOT submit before pushing.
#
# Run from the cafpyana root AFTER `source setup.sh`.
# Grid outputs land in $CAFPYANA_GRID_OUT_DIR/dfs4/<timestamp>__<NAME>/ ;
# merge into $OUTDIR afterwards with production/dfadd_all.sh.
#
# Config mapping (one production per GUMP file-list): only the PID-free MAPLE
# PRESELECTION is applied at DF-build time (NOT the full selection); the full,
# calo-varied selection is applied post-hoc from the preselected slices.
#   MC CV lists   -> maple_evtrec_presel  (presel + GUMP CV wgts + evtrec)
#   other MC      -> maple_evt_presel     (presel)
#   data lists    -> maple_evt_data_presel (presel, no mcnu/wgt)
# The GUMP duplicate no-wgt passes over the CV lists (SBNDMCCV_Nom,
# ICARUSRun2/4_Spring_Overlay_Nom) are intentionally skipped.

set -e

OUTDIR=/exp/sbnd/data/users/gputnam/MAPLE/sbn-rewgted-14
FILELISTS=/exp/sbnd/app/users/gputnam/osc2/cafpyana/analysis_village/gump/file-lists-reCAF
CFGDIR=analysis_village/maple/configs

CAFPYANA=${CAFPYANA_WD:?"source setup.sh in the cafpyana root first"}
cd "$CAFPYANA"
mkdir -p "$OUTDIR"

# submit <name> <list> <config> <nfileperjob>
submit() {
    NAME=$1
    LIST=$FILELISTS/$2.list
    CONFIG=$3
    NFILEPERJOB=$4

    if [ -f "$OUTDIR/$NAME.df" ] || [ -f "$OUTDIR/${NAME}_0.df" ]; then
        echo "$NAME: already merged in $OUTDIR, skipping"
        return 0
    fi
    [ -f "$LIST" ] || { echo "$NAME: file list $LIST not found"; exit 1; }

    NFILES=$(wc -l < "$LIST")
    NGRID=$(( (NFILES + NFILEPERJOB - 1) / NFILEPERJOB ))

    htgettoken -a htvaultprod.fnal.gov -i sbnd

    echo "Submitting $NAME ($CONFIG): $NFILES files in $NGRID jobs"
    python run_df_maker.py \
        -c $CFGDIR/$CONFIG.py \
        -l "$LIST" \
        -o "$NAME" \
        -ngrid $NGRID
}

### MC CV (full sel + weights + GENIE event record) ###
submit SBNDMCCV                          SBNDMCCV        maple_evtrec_presel 100
submit ICARUSRun2_SpringMCOverlay_rewgt  ICARUSRun2MCCV  maple_evtrec_presel 100
submit ICARUSRun4_SpringMCOverlay_rewgt  ICARUSRun4MCCV  maple_evtrec_presel 100

### SBND MC variations (full sel) ###
submit SBND_SpringMC_Nom          SBNDMCNominal   maple_evt_presel 200
submit SBND_SpringMC_WMNom        SBNDMCWMNominal maple_evt_presel 200
submit SBND_SpringMC_0xSCE        SBNDMC0xSCE     maple_evt_presel 200
submit SBND_SpringMC_2xSCE        SBNDMC2xSCE     maple_evt_presel 200
submit SBND_SpringMC_WMXThetaXW   SBNDMCWMXThXW   maple_evt_presel 200
submit SBND_SpringMC_WMYZ         SBNDMCWMYZ      maple_evt_presel 200
submit SBND_SpringMC_DENT         SBNDMCDENT      maple_evt_presel 200
submit SBND_SpringLowEMC          SBNDMCDirt      maple_evt_presel 200
submit SBNDIntimeMC               SBNDIntimeMC    maple_evt_presel 200
submit SBNDAr25                   SBNDAr25        maple_evt_presel 200

### ICARUS MC variations (full sel) ###
submit ICARUSRun2_Spring_Overlay_WMXThXW  ICARUSRun2WMXThXW maple_evt_presel 200
submit ICARUSRun2_Spring_Overlay_WMYZ     ICARUSRun2WMYZ    maple_evt_presel 200
submit ICARUSRun2_Spring_Overlay_SCE      ICARUSRun2MCSCE   maple_evt_presel 200
submit ICARUSRun4_Spring_Overlay_WMXThXW  ICARUSRun4WMXThXW maple_evt_presel 200
submit ICARUSRun4_Spring_Overlay_WMYZ     ICARUSRun4WMYZ    maple_evt_presel 200
submit ICARUSRun4_Spring_Overlay_SCE      ICARUSRun4MCSCE   maple_evt_presel 200
submit ICARUSRun4_Spring_Overlay_Dirt     ICARUSRun4MCDirt  maple_evt_presel 200
submit ICARUS_Spring_Overlay_Dirt         ICARUSRun2MCDirt  maple_evt_presel 200

### Data (full sel, no truth) ###
submit SBND_SpringBNBData_FixedDev    SBNDDataFixedDev    maple_evt_data_presel 200
submit SBND_SpringBNBData_RollingDev  SBNDDataRollingDev  maple_evt_data_presel 200
submit SBND_SpringBNBOffData          SBNDDataOffbeam     maple_evt_data_presel 200
submit ICARUS_SpringRun2BNB_unblind     ICARUSRun2BeamOn  maple_evt_data_presel 200
submit ICARUS_SpringRun2BNBOff_unblind  ICARUSRun2BeamOff maple_evt_data_presel 200
submit ICARUS_SpringRun4BNB_unblind     ICARUSRun4BeamOn  maple_evt_data_presel 200
submit ICARUS_SpringRun4BNBOff_unblind  ICARUSRun4BeamOff maple_evt_data_presel 200
