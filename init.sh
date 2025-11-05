#!/bin/bash 

export machine=${HOSTNAME}

# try: get SAM Web from FNAL tools env. Do this first, since next
# spack will override this one
SPACK_ROOT_FNAL_TOOLS="/cvmfs/fermilab.opensciencegrid.org/packages/common/setup-env.sh"
source "${SPACK_ROOT_FNAL_TOOLS}"
spack load sam-web-client@3.6%gcc@11.4.1 arch=linux-almalinux9-x86_64_v2

which samweb > /dev/null 2>&1 && echo "SAMWeb client set up" || echo "Warning: SAMWeb client not set up"

SPACK_ROOT="/cvmfs/larsoft.opensciencegrid.org/spack-fnal-v1.0.0/setup-env.sh"
if [[ $machine == *sbnd* || $machine == *jupyter* ]]; then
  echo "working on a sbnd machine"
  source "${SPACK_ROOT}"
  export CAFPYANA_GRID_OUT_DIR="/pnfs/sbnd/scratch/users/$USER/cafpyana_out"
  mkdir -p $CAFPYANA_GRID_OUT_DIR
fi
if [[ $machine == *icarus* ]]; then
  echo "working on a icarus machine"
  source "${SPACK_ROOT}"
  export CAFPYANA_GRID_OUT_DIR="/pnfs/icarus/scratch/users/$USER/cafpyana_out"
  mkdir -p $CAFPYANA_GRID_OUT_DIR
fi
spack load hdf5@1.12.2%gcc@12.5.0 arch=linux-almalinux9-x86_64_v2
spack load xrootd@5.6.9%gcc@12.5.0 arch=linux-almalinux9-x86_64_v2

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

######################################################
# Install DUNE plot style, https://github.com/DUNE/dune_plot_style/blob/main/README.md#standalone-python-setup
######################################################
#cd envs
#export DUNE_PLOT_STYLE_LATEST_TAG=`curl --silent "https://api.github.com/repos/DUNE/dune_plot_style/releases" | jq -r 'map(select(.prerelease == false)) | first | .tag_name'`
#wget --no-check-certificate https://github.com/DUNE/dune_plot_style/archive/refs/tags/${DUNE_PLOT_STYLE_LATEST_TAG}.tar.gz -O dune_plot_style.tar.gz
#tar -xvzf dune_plot_style.tar.gz
#cd dune_plot_style-01_01/
#python3 -m pip install .
#cd ../..

######################################################
# need to install uuid in the EAF
######################################################
if [[ $machine == *jupyter* ]]; then
    echo "Installing uuid for since you are in EAF"
    wget https://www.kernel.org/pub/linux/utils/util-linux/v2.39/util-linux-2.39.3.tar.xz
    tar xvf util-linux-2.39.3.tar.xz
    rm util-linux-2.39.3.tar.xz
    cd util-linux-2.39.3
    ./configure --prefix="$(pwd)/local" --disable-all-programs --enable-libuuid
    make -j$(nproc)
    make install
    export C_INCLUDE_PATH="$(pwd)/local/include:$C_INCLUDE_PATH"
    export CPLUS_INCLUDE_PATH="$(pwd)/local/include:$CPLUS_INCLUDE_PATH"
    export LD_LIBRARY_PATH="$(pwd)/local/lib:$LD_LIBRARY_PATH"
    export PKG_CONFIG_PATH="$(pwd)/local/lib/pkgconfig:$PKG_CONFIG_PATH"
    cd ../
fi

######################################################
# Needed to install xrootd -- which, by the way, is super annoying
###################################################### 
OLDPATH=$PATH
PATH=$PATH:$PWD
cd envs
ln -s /cvmfs/larsoft.opensciencegrid.org/products/cmake/v3_22_2/Linux64bit+3.10-2.17/bin/cmake cmake3
which cmake3
wget https://files.pythonhosted.org/packages/96/e9/32107ac154c33c6bafd53a5f8444290938c3557210276e5deabb82f74b8f/xrootd-5.6.1.tar.gz
tar -zxvf xrootd-5.6.1.tar.gz
rm xrootd-5.6.1.tar.gz
cd xrootd-5.6.1
python setup.py install
cd ../..
PATH=$OLDPATH

###################################################### 

export PYTHONPATH=$PYTHONPATH:$PWD/..
export LD_LIBRARY_PATH=$VIRTUAL_ENV/lib/python3.9/site-packages/xrootd-5.6.1-py3.9-linux-x86_64.egg/pyxrootd:$LD_LIBRARY_PATH
export CAFPYANA_WD=`pwd`
