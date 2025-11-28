from makedf.makedf import make_hdrdf, make_mcprimdf
from analysis_village.kaon.makedf import make_kaon_mcdf

DFS = [make_kaon_mcdf, make_hdrdf, make_mcprimdf]
NAMES = ["kmc", "hdr", "mcprim"]
