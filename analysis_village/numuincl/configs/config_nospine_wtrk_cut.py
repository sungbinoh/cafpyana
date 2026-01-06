from makedf.makedf import make_hdrdf, make_mcnudf, make_trkdf
from analysis_village.numuincl.makedf1muX import make_pandora_evtdf_processed_signal_cut, make_pandora_evtdf_processed_selected_cut

DFS = [
      make_hdrdf,
       make_mcnudf,
        #make_custom_trkdf,
    #    make_mcdf,
    #    make_slicedf,
       #make_spine_evtdf_wgt,
    #    make_spinepartdf,
    #    make_spineinterdf,
       make_pandora_evtdf_processed_signal_cut,
       make_pandora_evtdf_processed_selected_cut

       ]
NAMES = [
  'hdr',
  'mcnu',
   #'trk',
#   'mcprim',
#   'slice',
  #'evt',
#   'particle',
#   'interaction',
  'evt_pand_signal',
  'evt_pand_selected'
  ]