from analysis_village.NCDeltaWC.make_wc_candidates import make_wc_candidate_df
from analysis_village.NCDeltaWC.make_truth import make_ncdeltawc_truthkinedf

# wc_candidates.py plus mcnu: every true neutrino with its final-state
# kinematics (make_truth.py), for truth studies and efficiency denominators.
# Candidates join to it on (__ntuple, entry, tmatch_idx).
DFS = [make_wc_candidate_df, make_ncdeltawc_truthkinedf]
NAMES = ["wc", "mcnu"]
