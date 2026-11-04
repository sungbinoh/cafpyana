from analysis_village.NCDeltaWC.make_wc_candidates import make_wc_candidate_df
from analysis_village.SinglePhotonWC.make_truth import make_sp_truthdf

# wc: one row per Wire-Cell candidate, reco only (the candidate table is
# analysis-independent, so it is shared with NCDeltaWC).
# mcnu: every true neutrino with the single-photon truth (make_truth.py).
# Join on (__ntuple, entry, tmatch_idx). POT from histpotdf / dfio totals.
DFS = [make_wc_candidate_df, make_sp_truthdf]
NAMES = ["wc", "mcnu"]
