Order of operations for using these scripts:
1. Use flat caf files and an input for rundf.py
2. Using the dataframes from this, run the GUMP selection using nb/gump_selection.py
3. You can then convert these outputs to sBruce trees to use with ProFit using PyRootConvert.py.

For PyRootConvert set up, use (from Sl7 container) 

source /cvmfs/larsoft.opensciencegrid.org/setup_larsoft.sh
setup root v6_28_12 -q e26:p3915:prof
