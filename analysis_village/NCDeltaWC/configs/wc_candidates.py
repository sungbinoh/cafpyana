from analysis_village.NCDeltaWC.make_wc_candidates import make_wc_candidate_df

# POT: run_df_maker adds histpotdf (TotalPOT per CAF) to every output; use
# that, or dfio.load()'s totals. Summing rec.hdr.pot undercounts on these
# merged CAFs (~5%): subruns with no event left in recTree drop out.
DFS = [make_wc_candidate_df]
NAMES = ["wc"]
