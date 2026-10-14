#!/bin/bash

echo "BEARER_TOKEN_FILE is set to: $BEARER_TOKEN_FILE"

######################################################
#### ROOT + GENIE for the raw GENIE event-record (evtrec) reader
######################################################
# The SBND flat CAFs store GenieEvtRecTree as a raw genie::NtpMCEventRecord object
# that uproot cannot deserialize; make_genie_evtrec_df (evtrec DF) then falls back to
# pyROOT + GENIE 3.06.02 (makedf/genie_evtrec.py). Mirror setup.sh: load root from the
# spack-fnal-v1.0.0 instance (this sets python to 3.10.16, so the venv built below has a
# matching pyROOT) and genie@3.06.02, exporting $GENIE so GENIE's Messenger finds its XML
# config. Best-effort + guarded: if unavailable the ICARUS (uproot StdHep) evtrec path and
# all non-evtrec DFs still work; only the SBND raw-object path needs this.
SPACK_FNAL="/cvmfs/larsoft.opensciencegrid.org/spack-fnal-v1.0.0/setup-env.sh"
if [ -f "$SPACK_FNAL" ]; then
  source "$SPACK_FNAL"
  # hdf5 + xrootd + root, matching setup.sh, so the venv is built on spack-root's
  # python 3.10.16 and the shipped py3.10 XRootD build has its libs on LD_LIBRARY_PATH.
  spack load hdf5@1.12.2%gcc@12.5.0 arch=linux-almalinux9-x86_64_v2 2>/dev/null
  spack load xrootd@5.6.9%gcc@12.5.0 arch=linux-almalinux9-x86_64_v2 2>/dev/null
  spack load root arch=linux-almalinux9-x86_64_v2
  echo "@@ root set up: $(command -v root-config) $(root-config --version 2>/dev/null)"
  echo "@@ python for venv: $(command -v python) $(python --version 2>&1)"
  if spack load genie@3.06.02 2>/dev/null; then
    export ROOT_INCLUDE_PATH=${GENIE}/include/GENIE:${ROOT_INCLUDE_PATH}
    echo "@@ GENIE set up ($GENIE)"
  else
    echo "@@ Warning: genie@3.06.02 not loaded; SBND raw GENIE event-record reading unavailable"
  fi
else
  echo "@@ Warning: $SPACK_FNAL not found; ROOT/GENIE (SBND evtrec) unavailable"
fi

######################################################
#### setup virtual python env if it is not already set
######################################################
# Define the Python version and virtual environment name
#!/bin/bash
cd envs
# Define the Python version and virtual environment name
PYTHON_VERSION=3.10.16
VENV_NAME=venv_py310_cafpyana

# Check if virtual environment already exists
if [ -d "$VENV_NAME" ]; then
    echo "Virtual environment '$VENV_NAME' already exists. Activating it."
else
    # Create the virtual environment
    python -m venv $VENV_NAME
    echo "Virtual environment '$VENV_NAME' created."
fi

# Activate the virtual environment
source $VENV_NAME/bin/activate

# Upgrade pip
pip install --upgrade pip
pip install wheel setuptools
cat pip_requirements.txt
pip install -r pip_requirements.txt

# Deactivate virtual environment
echo "Virtual environment '$VENV_NAME' set up successfully with Python $PYTHON_VERSION and required packages installed."
echo "If you want to exit from this virtual env, do $ deactivate"
echo "If you want to activate this virtual evn again, do $ source '$VENV_NAME'/bin/activate"
cd ..
######################################################

export PYTHONPATH=$PYTHONPATH:$PWD/..
export CAFPYANA_WD=`pwd`
echo $CAFPYANA_WD
echo "BEARER_TOKEN_FILE is set to: $BEARER_TOKEN_FILE"
