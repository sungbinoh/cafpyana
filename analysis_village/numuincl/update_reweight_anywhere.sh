#!/bin/bash

[[ -d '/scratch' ]] && export APPTAINER_BINDPATH+="${APPTAINER_BINDPATH:+,}/scratch"
[[ -d "/run/user/$(id -u)" ]] && export APPTAINER_BINDPATH+="${APPTAINER_BINDPATH:+,}/run/user/$(id -u)"

TEMPFILE=$(mktemp --suffix=".list")
echo $1 > $TEMPFILE
/cvmfs/oasis.opensciencegrid.org/mis/apptainer/current/bin/apptainer exec \
--pid --ipc \
-B /etc/hosts,/tmp,/opt,/cvmfs,/exp/,/pnfs/,/nashome \
-B /etc/profile.d/jobsub_lite.sh,/etc/condor/,/etc/grid-security/ \
--home ~/:${HOME} --pwd ${PWD} \
/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-dev-sl7:jsl \
/bin/bash -c  "cd /exp/sbnd/app/users/gputnam/Ar23-knobs && source nusyst_env.sh && UpdateReweight -c nusystematics/fcl/Ar23_new_knobs_local.ParamHeader.fcl -i $TEMPFILE -o $2"
rm $TEMPFILE
