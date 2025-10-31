from makedf.makedf import *
from analysis_village.cc2p.makedf.make_cc2pdf import *

# DFS = [make_cc2pdf, make_hdrdf, make_potdf_bnb, make_mcprimdf,make_mcnuwgtdf]
# NAMES = ["cc2p", "hdr", "pot" , "mcprim","mcnuwgt"]

DFS = [make_pandora_df, make_hdrdf, make_potdf_bnb, make_mcnuwgtdf]
NAMES = ["evt", "hdr", "pot", "mcnuwgt" ]