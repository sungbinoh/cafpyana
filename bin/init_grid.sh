#!/bin/bash 

echo "BEARER_TOKEN_FILE is set to: $BEARER_TOKEN_FILE"
#export machine=${HOSTNAME}
#if [[ $machine == *sbnd* ]]; then
#  echo "working on a sbnd machine"
#  source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
#  export CAFPYANA_GRID_OUT_DIR="/pnfs/sbnd/scratch/users/$USER/cafpyana_out"
#  mkdir -p $CAFPYANA_GRID_OUT_DIR
#fi
#if [[ $machine == *icarus* ]]; then
#  echo "working on a icarus machine"
#  source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
#  export CAFPYANA_GRID_OUT_DIR="/pnfs/icarus/scratch/users/$USER/cafpyana_out"
#  mkdir -p $CAFPYANA_GRID_OUT_DIR
#fi
#spack load hdf5@1.14.3
#spack load xrootd@5.6.1

######################################################
#### setup virtual python env if it is not already set
######################################################
# Define the Python version and virtual environment name
#!/bin/bash
cd envs
# Define the Python version and virtual environment name
PYTHON_VERSION=3.9.21
VENV_NAME=venv_py39_cafpyana

# Check if virtual environment already exists
if [ -d "$VENV_NAME" ]; then
    echo "Virtual environment '$VENV_NAME' already exists. Activating it."
else
    # Create the virtual environment
    python -m venv $VENV_NAME
    echo "Virtual environment '$VENV_NAME' created."
fi

echo "@@ before pip test:"
xrdcp root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/sbn/data_add/sbn_nd/poms_production/mc/MCP2025C_1e20_v10_06_00_09/v10_06_00_09/prodgenie_corsika_proton_rockbox_sbnd/CV/caf/d2/64/caf.flat.caf-3750ae6a-8d07-44db-b813-1a16c0407cf4.root before_pip.root
ls -alh

# Activate the virtual environment
source $VENV_NAME/bin/activate

# Upgrade pip
pip install --upgrade pip
pip install wheel setuptools
pip install -r pip_requirements.txt
echo "@@ after pip test:"
xrdcp root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/sbn/data_add/sbn_nd/poms_production/mc/MCP2025C_1e20_v10_06_00_09/v10_06_00_09/prodgenie_corsika_proton_rockbox_sbnd/CV/caf/d2/64/caf.flat.caf-3750ae6a-8d07-44db-b813-1a16c0407cf4.root after_pip.root
ls -alh

#htgettoken -a htvaultprod.fnal.gov -i sbnd

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
