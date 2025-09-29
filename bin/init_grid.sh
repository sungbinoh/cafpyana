#!/bin/bash 

export machine=${HOSTNAME}
echo "working on machine $HOSTNAME"

source /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh
export CAFPYANA_GRID_OUT_DIR="/pnfs/sbn/scratch/users/$USER/cafpyana_out"
mkdir -p $CAFPYANA_GRID_OUT_DIR
spack load hdf5@1.14.3%gcc@12.2.0 arch=linux-almalinux9-x86_64_v3
spack load xrootd@5.6.9%gcc@12.2.0
spack load ifdhc@2.7.2
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

# Activate the virtual environment
source $VENV_NAME/bin/activate

# Upgrade pip
pip install --upgrade pip
pip install wheel setuptools
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
