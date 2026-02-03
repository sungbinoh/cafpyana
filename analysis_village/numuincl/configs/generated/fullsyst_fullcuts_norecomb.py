from importlib import import_module
from analysis_village.numuincl import makedf1muX as maker

# Auto-generated for job 'fullsyst_fullcuts_norecomb'
maker.INCLUDE_WEIGHTS = True
maker.SLIM = False
maker.set_update_recomb(False)
maker.VERBOSE = False
maker.ADD_STAT_UNC = False
maker.apply_setting_dependencies()

_base = import_module('analysis_village.numuincl.configs.base')
DFS = [
    getattr(import_module('makedf.makedf'), 'make_hdrdf'),
    getattr(import_module('analysis_village.numuincl.makedf1muX'), 'make_mcnu_processed'),
    getattr(import_module('analysis_village.numuincl.makedf1muX'), 'make_pandora_evtdf_processed_signal_cut'),
    getattr(import_module('analysis_village.numuincl.makedf1muX'), 'make_pandora_evtdf_processed_selected_cut'),
]
NAMES = ['hdr', 'mcnu', 'evt_pand_signal', 'evt_pand_selected']
