# Config for preselct.py: post-hoc filtering of a "no selection" gump .df
# to the PID-free MAPLE preselection.
#   python preselct.py -c analysis_village/gump/configs/preselect_gump.py -o out.df in.df
CANDIDATE = "evt"
CUT = lambda df: df.gump_presel
KEEP_ALL = ["hdr", "trig", "bnb", "mcnu", "wgt", "evtrec", "histpotdf", "histgenevtdf"]
DROP = []
