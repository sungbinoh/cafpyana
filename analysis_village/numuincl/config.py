from makedf.makedf import make_hdrdf, make_mcnudf
from analysis_village.numuincl.makedf1muX import make_spine_evtdf_wgt, make_pandora_evtdf_wgt

DFS = [
      make_hdrdf,
       make_mcnudf,
    #    make_pfpdf,
    #    make_mcdf,
    #    make_slicedf,
       make_spine_evtdf_wgt,
    #    make_spinepartdf,
    #    make_spineinterdf,
       make_pandora_evtdf_wgt,
       ]
NAMES = [
  'hdr',
  'mcnu',
#   'pfp',
#   'mcprim',
#   'slice',
  'evt',
#   'particle',
#   'interaction',
  'evt_pand'
  ]