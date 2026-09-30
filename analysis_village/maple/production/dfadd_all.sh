#!/bin/bash
# Post-process the MAPLE grid production over the GUMP reCAF file-lists:
# for each sample, find the newest grid output directory in
# $CAFPYANA_GRID_OUT_DIR/dfs4/<timestamp>__<NAME>/ and dfadd the per-job
# .df files into $OUTDIR/<NAME>.df (rolling to <NAME>_0.df, <NAME>_1.df, ...
# above the dfadd 30 GB file-split).
#
# A sample is skipped when the merged output already exists (idempotent --
# rerun until everything is merged), or when the grid production is
# incomplete (fewer job outputs than submitted jobs) unless -f is given.
#
# Usage: ./dfadd_all.sh [-f] [sample ...]
#   -f          merge even if some grid jobs are missing
#   [sample]    restrict to the given output names (default: all)
#
# Run from the cafpyana root AFTER `source setup.sh`.

OUTDIR=/exp/sbnd/data/users/gputnam/MAPLE/sbn-rewgted-14
# Only consider grid-output dirs from this MAPLE campaign onward: the GUMP
# sbn-rewgted-13 campaign used the same scratch area and the SAME sample
# names, so older complete GUMP dirs must never be picked up.
MINTS=2026_08_12_1149

CAFPYANA=${CAFPYANA_WD:?"source setup.sh in the cafpyana root first"}
GRIDOUT=${CAFPYANA_GRID_OUT_DIR:?"CAFPYANA_GRID_OUT_DIR not set; source setup.sh"}
cd "$CAFPYANA"
mkdir -p "$OUTDIR"

FORCE=0
if [ "$1" == "-f" ]; then FORCE=1; shift; fi
ONLY=("$@")

ALL_SAMPLES="SBNDMCCV ICARUSRun2_SpringMCOverlay_rewgt ICARUSRun4_SpringMCOverlay_rewgt
SBND_SpringMC_Nom SBND_SpringMC_WMNom SBND_SpringMC_0xSCE SBND_SpringMC_2xSCE
SBND_SpringMC_WMXThetaXW SBND_SpringMC_WMYZ SBND_SpringMC_DENT SBND_SpringLowEMC SBNDIntimeMC SBNDAr25
ICARUSRun2_Spring_Overlay_WMXThXW ICARUSRun2_Spring_Overlay_WMYZ ICARUSRun2_Spring_Overlay_SCE
ICARUSRun4_Spring_Overlay_WMXThXW ICARUSRun4_Spring_Overlay_WMYZ ICARUSRun4_Spring_Overlay_SCE
ICARUSRun4_Spring_Overlay_Dirt ICARUS_Spring_Overlay_Dirt
SBND_SpringBNBData_FixedDev SBND_SpringBNBData_RollingDev SBND_SpringBNBOffData
ICARUS_SpringRun2BNB_unblind ICARUS_SpringRun2BNBOff_unblind
ICARUS_SpringRun4BNB_unblind ICARUS_SpringRun4BNBOff_unblind"

NDONE=0; NSKIP=0; NMERGED=0; NFAIL=0

for NAME in $ALL_SAMPLES; do
    if [ ${#ONLY[@]} -gt 0 ]; then
        case " ${ONLY[*]} " in *" $NAME "*) ;; *) continue ;; esac
    fi

    if [ -f "$OUTDIR/$NAME.df" ] || [ -f "$OUTDIR/${NAME}_0.df" ]; then
        echo "$NAME: already merged"
        NDONE=$((NDONE+1)); continue
    fi

    # newest production dir whose basename ends exactly in __<NAME>
    # (anchored match: SBNDMCCV must not pick up e.g. ..._rewgt dirs)
    # and whose timestamp is within this campaign (>= MINTS)
    DFDIR=""
    for d in "$GRIDOUT"/dfs4/*__"$NAME"; do
        [ -d "$d" ] || continue
        TS=$(basename "$d"); TS=${TS%__$NAME}
        [[ "$TS" < "$MINTS" ]] && continue
        case "$(basename "$d")" in *__"$NAME") DFDIR=$d ;; esac
    done
    if [ -z "$DFDIR" ]; then
        echo "$NAME: no grid output dir found, skipping"
        NSKIP=$((NSKIP+1)); continue
    fi

    TS=$(basename "$DFDIR"); TS=${TS%__$NAME}
    LOGDIR=$GRIDOUT/logs/${TS}__${NAME}_log
    NJOBS=$(ls "$LOGDIR"/run_*.sh 2>/dev/null | wc -l)
    NOUT=$(ls "$DFDIR/$NAME"_*.df 2>/dev/null | wc -l)

    if [ "$NOUT" -eq 0 ]; then
        echo "$NAME: no job outputs yet in $DFDIR, skipping"
        NSKIP=$((NSKIP+1)); continue
    fi
    if [ "$NOUT" -lt "$NJOBS" ] && [ $FORCE -eq 0 ]; then
        echo "$NAME: MISSING $((NJOBS-NOUT))/$NJOBS job outputs in $DFDIR, skipping (use -f to force)"
        NSKIP=$((NSKIP+1)); continue
    fi

    echo "$NAME: merging $NOUT/$NJOBS job outputs from $DFDIR"
    if python dfadd.py -o "$OUTDIR/$NAME.df" "$DFDIR/$NAME"_*.df; then
        ls -lh "$OUTDIR/$NAME"*.df
        NMERGED=$((NMERGED+1))
    else
        echo "$NAME: dfadd FAILED"
        rm -f "$OUTDIR/$NAME.df" "$OUTDIR/${NAME}"_[0-9]*.df
        NFAIL=$((NFAIL+1))
    fi
done

echo "dfadd_all: $NDONE already done, $NMERGED merged, $NSKIP skipped, $NFAIL failed"
[ $NFAIL -eq 0 ]
