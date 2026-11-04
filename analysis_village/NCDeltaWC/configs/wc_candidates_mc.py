from analysis_village.NCDeltaWC.make_wc_candidates import make_wc_candidate_df
from analysis_village.NCDeltaWC.make_truth import make_ncdeltawc_truthkinedf

# wc_candidates.py plus mcnu: every true neutrino with its final-state
# kinematics (make_truth.py), for efficiency denominators and truth studies.
# The wc table already carries the same true_* columns for its matched
# neutrino; mcnu adds the neutrinos without a candidate.
DFS = [make_wc_candidate_df, make_ncdeltawc_truthkinedf]
NAMES = ["wc", "mcnu"]
