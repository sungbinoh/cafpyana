# MAPLE evt production, data, with the muon-candidate end-hit table: same as
# maple_evt_data_presel.py, plus "muendhit" -- the x/y/z (and rr) of every
# collection-plane calo point of the muon candidate within its last 25 cm
# (rr <= MU_ENDHIT_RR_MAX), indexed (entry, slc, hit) so it joins onto evt:
#   evt.join(muendhit, how="inner")
# make_maple_mu_endhit_df reads a table cached by the evt builder, so it must
# come after make_maple_evt_presel_data_muendhit_df in DFS.
from analysis_village.maple.makedf import *

DFS = [make_maple_evt_presel_data_muendhit_df, make_maple_mu_endhit_df, make_hdrdf, make_triggerdf, make_maple_bnbdf]
NAMES = ["evt", "muendhit", "hdr", "trig", "bnb"]
