from analysis_village.NCDeltaWC.make_wc_candidates import make_wc_candidate_df
from analysis_village.SinglePhotonWC.make_truth import make_sp_truthdf
from analysis_village.SinglePhotonWC.make_wc_tagger import make_wc_tagger_df

# sp_mc.py plus wctag: every scalar T_tagger leaf (tag_*), same index as wc,
# for BDT inputs. About 1 MB per merged CAF more than sp_mc.py.
DFS = [make_wc_candidate_df, make_wc_tagger_df, make_sp_truthdf]
NAMES = ["wc", "wctag", "mcnu"]
