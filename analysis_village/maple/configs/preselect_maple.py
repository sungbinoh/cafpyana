# Config for preselct.py: post-hoc filtering of a "no selection" maple .df
# to the PID-free MAPLE preselection.
#   python preselct.py -c analysis_village/maple/configs/preselect_maple.py -o out.df in.df
CANDIDATE = "evt"
CUT = lambda df: df.maple_presel
KEEP_ALL = ["hdr", "trig", "bnb", "mcnu", "wgt", "evtrec", "histpotdf", "histgenevtdf"]
DROP = []
