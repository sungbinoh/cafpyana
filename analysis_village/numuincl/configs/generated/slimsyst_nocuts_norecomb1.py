from importlib import import_module
from analysis_village.numuincl import makedf1muX as maker

# Auto-generated for job 'slimsyst_nocuts_norecomb1'
maker.INCLUDE_WEIGHTS = True
maker.SLIM = True
maker.set_update_recomb(False)
maker.VERBOSE = True
maker.ADD_STAT_UNC = False
maker.apply_setting_dependencies()

_base = import_module('analysis_village.numuincl.configs.base')
DFS = [
    getattr(import_module('makedf.makedf'), 'make_hdrdf'),
    getattr(import_module('makedf.makedf'), 'make_mcnudf'),
    getattr(import_module('analysis_village.numuincl.makedf1muX'), 'make_pandora_evtdf_processed'),
]
NAMES = ['hdr', 'mcnu', 'evt_pand']
