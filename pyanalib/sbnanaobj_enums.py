"""
Load enums from sbnanaobj header file
"""

import pathlib

_SBNANAOBJ_VERSION = 'v10_00_14'
_SBNANAOBJ_PATH = pathlib.Path(f'/cvmfs/sbn.opensciencegrid.org/products/sbn/sbnanaobj/{_SBNANAOBJ_VERSION}/include/sbnanaobj/StandardRecord/SREnums.h')
_CAF_NAMESPACE = 'caf'
_ROOT_IGNORE_ATTR = ['as_integer_ratio', 'bit_count', 'bit_length', 'conjugate', 'denominator', 'from_bytes', 'imag', 'numerator', 'real', 'to_bytes']


def _load_enums():
    '''Use the ROOT interpreter to load enums.'''
    try:
        import ROOT
    except ImportError:
        print('Error: Could not import ROOT. sbnanaobj enum module requires ROOT.')
        return

    if not _SBNANAOBJ_PATH.is_file():
        print(f'Error: Could not load sbnanaobj header: {_SBNANAOBJ_PATH}. Does your machine have access to CVMFS?')
        return

    print(f'Enums will be loaded from sbnanaobj version {_SBNANAOBJ_VERSION}')

    with open(_SBNANAOBJ_PATH) as f:
        for l in f.readlines():
            ROOT.gInterpreter.ProcessLine(l)

    caf = getattr(ROOT, _CAF_NAMESPACE)
    enum_names = [name for name in dir(caf) if not name.startswith('__')]

    for name in enum_names:
        e = getattr(caf, name)
        field_names = [f for f in dir(e) if not (f in _ROOT_IGNORE_ATTR or f.startswith('__'))]

        # create forward & backward enums
        # sort by value
        py_enum_fwd = {f: getattr(e, f) for f in field_names}
        py_enum_fwd = dict(sorted(py_enum_fwd.items(), key=lambda item: item[1]))

        # sort by key
        py_enum_rev = {getattr(e, f): f for f in field_names}
        py_enum_rev = dict(sorted(py_enum_rev.items()))

        globals()[name] = py_enum_fwd
        globals()[f'inv_{name}'] = py_enum_rev

_load_enums()
