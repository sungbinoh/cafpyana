from analysis_village.SinglePhotonWC.configs.sp_mc import DFS, NAMES
from analysis_village.SinglePhotonWC.make_weights import make_sp_slimwgtdf

# sp_mc.py plus the systematic weights per true neutrino ("slimwgt", make_weights.py;
# as gump's gump_no_cuts_slimwgt.py), index as mcnu.
DFS = DFS + [make_sp_slimwgtdf]
NAMES = NAMES + ["slimwgt"]
