#!/bin/bash 

export machine=${HOSTNAME}
export CAFPYANA_DIR=$(pwd)

# try: get SAM Web from FNAL tools env. Do this first, since next
# spack will override this one
SPACK_ROOT_FNAL_TOOLS="/cvmfs/fermilab.opensciencegrid.org/packages/common/setup-env.sh"
source "${SPACK_ROOT_FNAL_TOOLS}" > /dev/null 2>&1
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

# this sets Python to 3.10.16
spack load root arch=linux-almalinux9-x86_64_v2

# EAF has different library needs
EAF=0
if [[ $machine == *jupyter* ]]; then
    EAF=1
fi

######################################################
#### setup virtual python env if it is not already set
######################################################
# Define the Python version and virtual environment name
ENV_DIR=${CAFPYANA_DIR}/envs
LOGDIR=${ENV_DIR}/logs
mkdir -p ${LOGDIR}

# Define the Python version and virtual environment name
PYTHON_VERSION=3.10.16
VENV_NAME=venv_py310_cafpyana
VENV_DIR="${ENV_DIR}/${VENV_NAME}"

# Check if virtual environment already exists
if [ -d "${VENV_DIR}" ]; then
    echo "Virtual environment '$VENV_NAME' already exists. Activating it."
else
    # Create the virtual environment
    python -m venv ${VENV_DIR}
    echo "Virtual environment '$VENV_NAME' created."
fi

# Activate the virtual environment
source ${VENV_DIR}/bin/activate

# check requirements are installed
if pip freeze -r ${ENV_DIR}/pip_requirements.txt 2>&1 | grep -q "not installed"; then
    echo "Updating virtual env installation"

    PIPLOG=${LOGDIR}/init_pip.log
    echo $(date) >> ${PIPLOG}
    pip install --upgrade pip | tee -a ${PIPLOG}
    pip install wheel setuptools | tee -a ${PIPLOG}
    pip install -r ${ENV_DIR}/pip_requirements.txt | tee -a ${PIPLOG}
fi

# Deactivate virtual environment
echo "Virtual environment '$VENV_NAME' set up successfully with Python $PYTHON_VERSION and required packages installed."
echo "If you want to exit from this virtual env, do $ deactivate"
echo "If you want to activate this virtual evn again, do $ source '$VENV_NAME'/bin/activate"

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
if [ $EAF -eq 1 ]; then
    export C_INCLUDE_PATH="${ENV_DIR}/local/include:$C_INCLUDE_PATH"
    export CPLUS_INCLUDE_PATH="${ENV_DIR}/local/include:$CPLUS_INCLUDE_PATH"
    export LD_LIBRARY_PATH="${ENV_DIR}/local/lib:$LD_LIBRARY_PATH"
    export PKG_CONFIG_PATH="${ENV_DIR}/local/lib/pkgconfig:$PKG_CONFIG_PATH"

    # check for UUID
    if [ ! -f "${ENV_DIR}/local/lib/libuuid.so" ]; then
       cd ${ENV_DIR}
       echo "Installing uuid for since you are in EAF"
       UUIDLOG=${LOGDIR}/init_uuid.log
       echo $(date) >> ${UUIDLOG}
       wget https://www.kernel.org/pub/linux/utils/util-linux/v2.39/util-linux-2.39.3.tar.xz
       tar xvf util-linux-2.39.3.tar.xz
       rm util-linux-2.39.3.tar.xz
       cd util-linux-2.39.3
       ./configure --prefix="${ENV_DIR}/local" --disable-all-programs --enable-libuuid | tee -a ${UUIDLOG}
       make -j$(nproc) | tee -a ${UUIDLOG}
       make install | tee -a ${UUIDLOG}
       cd ${CAFPYANA_DIR}
    fi
fi

######################################################
# Needed to install xrootd -- which, by the way, is super annoying
###################################################### 
python -c "import XRootD" > /dev/null 2>&1 || {
    echo "Could not import XRootD! Attempting to install from source..."
    cd ${ENV_DIR}
    XROOTLOG=${LOGDIR}/init_xroot.log
    echo $(date) >> ${XROOTLOG}

    OLDPATH=$PATH
    PATH=$PATH:${ENV_DIR}
    ln -s /cvmfs/larsoft.opensciencegrid.org/products/cmake/v3_22_2/Linux64bit+3.10-2.17/bin/cmake cmake3
    echo "Using cmake at $(which cmake3)" | tee -a ${XROOTLOG}
    wget https://files.pythonhosted.org/packages/fd/4f/419b8caec575ab4133f41c37c39bd742251a4dc6a208a97a3fd772031fe7/xrootd-5.6.9.tar.gz
    tar -zxvf xrootd-5.6.9.tar.gz
    rm xrootd-5.6.9.tar.gz
    cd xrootd-5.6.9
    # this seems to be needed depending on SSL version -- GPVMs need it while EAF doesn't
    if [ $EAF -eq 0 ]; then
        sed -i 's/SSL_CTX_flush_sessions/SSL_CTX_flush_sessions_ex/g' src/XrdTls/XrdTlsContext.cc
    fi
    python setup.py install 2>&1 | tee -a ${XROOTLOG}
    cd ${CAFPYANA_DIR}
    PATH=$OLDPATH
}

###################################################### 

export PYTHONPATH=$PYTHONPATH:${CAFPYANA_DIR}
export LD_LIBRARY_PATH=${VENV_DIR}/lib/python3.9/site-packages/xrootd-5.6.9-py3.9-linux-x86_64.egg/pyxrootd:$LD_LIBRARY_PATH
export CAFPYANA_WD=${CAFPYANA_DIR}

htgettoken -a htvaultprod.fnal.gov -i sbnd
