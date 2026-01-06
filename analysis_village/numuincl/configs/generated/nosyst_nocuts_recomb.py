from importlib import import_module
from analysis_village.numuincl import makedf1muX as maker

# Auto-generated for job 'nosyst_nocuts_recomb'
maker.INCLUDE_WEIGHTS = False
maker.SLIM = False
maker.set_update_recomb(True)
maker.VERBOSE = False
maker.apply_setting_dependencies()

_base = import_module('analysis_village.numuincl.configs.base')
DFS = [
    getattr(import_module('makedf.makedf'), 'make_hdrdf'),
    getattr(import_module('makedf.makedf'), 'make_mcnudf'),
    getattr(import_module('analysis_village.numuincl.makedf1muX'), 'make_pandora_evtdf_processed'),
    getattr(import_module('analysis_village.numuincl.makedf1muX'), 'make_custom_trkdf'),
]
NAMES = ['hdr', 'mcnu', 'evt_pand', 'trk']
if hasattr(_base, 'PREPROCESS'):
    PREPROCESS = _base.PREPROCESS
