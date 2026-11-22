from analysis_village.SinglePhotonWC.configs.sp_mc import DFS, NAMES
from analysis_village.SinglePhotonWC.make_weights import make_sp_wgtdf

# sp_mc.py plus the systematic weights per true neutrino ("wgt", make_weights.py;
# as gump's gump_no_cuts_wgt.py), index as mcnu.
DFS = DFS + [make_sp_wgtdf]
NAMES = NAMES + ["wgt"]
