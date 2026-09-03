"""
Read the raw GENIE event record (the GHEP particle stack) out of the
``GenieEvtRecTree`` of a flat CAF file.

Some productions (e.g. the ICARUS-style flat CAFs) split ``GenieEvtRecTree`` into
flat ``StdHep*`` branches that uproot can read directly.  Other productions (e.g. the
SBND2026A flat CAFs) instead store the *raw* ``genie::NtpMCEventRecord`` C++ object
(-> ``genie::EventRecord`` -> ``GHepRecord`` -> ``TClonesArray<GHepParticle>``), which
uproot cannot deserialize because it relies on GENIE's custom ROOT streamers.

This module handles the second case: it loads the GENIE dictionaries through pyROOT,
JIT-compiles a tiny C++ reader (via ``gInterpreter``) that walks the GHEP stack with a
``TTreeReader``, and returns the per-particle / per-event quantities as numpy arrays.

The only hard requirement is that the GENIE 3.06.02 libraries are available and that the
``$GENIE`` environment variable is set *before* the first GENIE object is constructed
(otherwise GENIE's Messenger fails while looking up its XML config).  ``setup.sh`` sets
this up via ``spack load genie@3.06.02``; if it is missing we try to locate GENIE and set
it up ourselves, and otherwise raise a clear error.
"""

import os
import glob
import subprocess

import numpy as np

# GENIE framework libraries, in dependency-load order.
_GENIE_LIBS = [
    "libGFwMsg", "libGFwReg", "libGFwAlg", "libGFwNum", "libGFwUtl",
    "libGFwParDat", "libGFwInt", "libGFwGHEP", "libGFwEG", "libGFwNtp",
]

# Known cvmfs spack locations to fall back on if $GENIE is not already set.
_GENIE_FALLBACK_GLOBS = [
    "/cvmfs/larsoft.opensciencegrid.org/spack-fnal-v1.0.0/opt/spack/*/genie-3.06.02-*",
    "/cvmfs/larsoft.opensciencegrid.org/spack-fnal-v1.1.1/opt/spack/*/genie-3.06.02-*",
]

# The JIT-compiled C++ reader fills these module-level std::vectors on every call.
# Per-particle vectors are parallel (one entry per GHEP particle); per-event vectors
# carry one entry per tree entry (neutrino interaction).
_CXX_READER = r'''
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include <vector>
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"

// --- per-particle columns ---
std::vector<int>    g_entry, g_pindex, g_pdg, g_status, g_rescat;
std::vector<int>    g_fmother, g_lmother, g_fdaughter, g_ldaughter;
std::vector<double> g_px, g_py, g_pz, g_E, g_vx, g_vy, g_vz, g_vt;
// --- per-event columns ---
std::vector<int>    g_evt_entry, g_nstdhep, g_evtnum;
std::vector<double> g_xsec, g_diff_xsec, g_weight, g_prob;
std::vector<double> g_evtvtx_x, g_evtvtx_y, g_evtvtx_z, g_evtvtx_t;
long long g_nbad = 0;

void _cafpyana_clear_genie() {
    g_entry.clear(); g_pindex.clear(); g_pdg.clear(); g_status.clear(); g_rescat.clear();
    g_fmother.clear(); g_lmother.clear(); g_fdaughter.clear(); g_ldaughter.clear();
    g_px.clear(); g_py.clear(); g_pz.clear(); g_E.clear();
    g_vx.clear(); g_vy.clear(); g_vz.clear(); g_vt.clear();
    g_evt_entry.clear(); g_nstdhep.clear(); g_evtnum.clear();
    g_xsec.clear(); g_diff_xsec.clear(); g_weight.clear(); g_prob.clear();
    g_evtvtx_x.clear(); g_evtvtx_y.clear(); g_evtvtx_z.clear(); g_evtvtx_t.clear();
    g_nbad = 0;
}

long long cafpyana_read_genie_evtrec(const char* fn) {
    _cafpyana_clear_genie();
    TFile* f = TFile::Open(fn);
    if (!f || f->IsZombie()) return -1;
    TTreeReader reader("GenieEvtRecTree", f);
    TTreeReaderValue<genie::NtpMCEventRecord> ntp(reader, "GenieEvtRec");
    int ie = 0;
    while (reader.Next()) {
        try {
            genie::EventRecord* ev = ntp->event;
            if (!ev) { ie++; continue; }
            int np = ev->GetEntries();
            g_evt_entry.push_back(ie);
            g_nstdhep.push_back(np);
            g_evtnum.push_back((int) ntp->hdr.ievent);
            g_xsec.push_back(ev->XSec());
            g_diff_xsec.push_back(ev->DiffXSec());
            g_weight.push_back(ev->Weight());
            g_prob.push_back(ev->Probability());
            TLorentzVector* vtx = ev->Vertex();
            g_evtvtx_x.push_back(vtx ? vtx->X() : 0.);
            g_evtvtx_y.push_back(vtx ? vtx->Y() : 0.);
            g_evtvtx_z.push_back(vtx ? vtx->Z() : 0.);
            g_evtvtx_t.push_back(vtx ? vtx->T() : 0.);
            for (int j = 0; j < np; j++) {
                genie::GHepParticle* p = ev->Particle(j);
                if (!p) continue;
                g_entry.push_back(ie);         g_pindex.push_back(j);
                g_pdg.push_back(p->Pdg());     g_status.push_back(p->Status());
                g_rescat.push_back(p->RescatterCode());
                g_px.push_back(p->Px());       g_py.push_back(p->Py());
                g_pz.push_back(p->Pz());       g_E.push_back(p->E());
                g_vx.push_back(p->Vx());       g_vy.push_back(p->Vy());
                g_vz.push_back(p->Vz());       g_vt.push_back(p->Vt());
                g_fmother.push_back(p->FirstMother());   g_lmother.push_back(p->LastMother());
                g_fdaughter.push_back(p->FirstDaughter()); g_ldaughter.push_back(p->LastDaughter());
            }
        } catch (std::exception& e) {
            g_nbad++;
        }
        ie++;
    }
    f->Close();
    return (long long) g_entry.size();
}
'''

_loaded = False


def _find_genie_prefix():
    """Return the GENIE install prefix, or None if it can't be found."""
    p = os.environ.get("GENIE")
    if p and os.path.isdir(p):
        return p
    # try spack
    try:
        out = subprocess.check_output(
            ["spack", "location", "-i", "genie@3.06.02"],
            stderr=subprocess.DEVNULL, text=True).strip()
        if out and os.path.isdir(out):
            return out
    except Exception:
        pass
    # try known cvmfs locations
    for pattern in _GENIE_FALLBACK_GLOBS:
        hits = sorted(glob.glob(pattern))
        if hits:
            return hits[0]
    return None


def _ensure_genie_loaded():
    """Load the GENIE dictionaries and JIT-compile the C++ reader (idempotent)."""
    global _loaded
    if _loaded:
        return
    import ROOT

    prefix = _find_genie_prefix()
    if prefix is None:
        raise RuntimeError(
            "Could not locate the GENIE installation needed to read the raw "
            "GENIE event record from this flat CAF. Please run "
            "`spack load genie@3.06.02` (this is done for you by setup.sh).")

    # $GENIE must be set before any GENIE object is constructed, otherwise GENIE's
    # Messenger throws 'basic_string: construction from null' looking up its XML config.
    os.environ["GENIE"] = prefix
    inc = os.path.join(prefix, "include", "GENIE")
    lib = os.path.join(prefix, "lib")
    if inc not in os.environ.get("ROOT_INCLUDE_PATH", "").split(":"):
        os.environ["ROOT_INCLUDE_PATH"] = inc + ":" + os.environ.get("ROOT_INCLUDE_PATH", "")
    if lib not in os.environ.get("LD_LIBRARY_PATH", "").split(":"):
        os.environ["LD_LIBRARY_PATH"] = lib + ":" + os.environ.get("LD_LIBRARY_PATH", "")
    ROOT.gInterpreter.AddIncludePath(inc)

    for name in _GENIE_LIBS:
        # load by full path so this works even if LD_LIBRARY_PATH wasn't set by spack
        full = os.path.join(lib, name + ".so")
        ROOT.gSystem.Load(full if os.path.exists(full) else name)
    # force the dictionary to initialize before we compile against the headers
    ROOT.TClass.GetClass("genie::NtpMCEventRecord")

    if not ROOT.gInterpreter.Declare(_CXX_READER):
        raise RuntimeError("Failed to JIT-compile the GENIE event-record reader.")
    _loaded = True


# columns filled by the C++ reader: (python name -> global vector name, numpy dtype)
_PART_COLS = {
    "pindex": ("g_pindex", np.int32),
    "pdg": ("g_pdg", np.int32),
    "status": ("g_status", np.int32),
    "rescat": ("g_rescat", np.int32),
    "px": ("g_px", np.float64), "py": ("g_py", np.float64),
    "pz": ("g_pz", np.float64), "E": ("g_E", np.float64),
    "vx": ("g_vx", np.float64), "vy": ("g_vy", np.float64),
    "vz": ("g_vz", np.float64), "vt": ("g_vt", np.float64),
    "fmother": ("g_fmother", np.int32), "lmother": ("g_lmother", np.int32),
    "fdaughter": ("g_fdaughter", np.int32), "ldaughter": ("g_ldaughter", np.int32),
}
_EVT_COLS = {
    "nstdhep": ("g_nstdhep", np.int32),
    "evtnum": ("g_evtnum", np.int32),
    "xsec": ("g_xsec", np.float64),
    "diff_xsec": ("g_diff_xsec", np.float64),
    "weight": ("g_weight", np.float64),
    "prob": ("g_prob", np.float64),
    "vtx_x": ("g_evtvtx_x", np.float64), "vtx_y": ("g_evtvtx_y", np.float64),
    "vtx_z": ("g_evtvtx_z", np.float64), "vtx_t": ("g_evtvtx_t", np.float64),
}


def _v2np(vec, dtype):
    """Copy a cppyy std::vector into a numpy array."""
    a = np.array(vec, dtype=dtype)
    return a if a.size else np.empty(0, dtype=dtype)


def read_genie_evtrec(path):
    """Read the raw GENIE event record from ``path`` (a local or xrootd file path).

    Returns a dict with:
      - ``entry``: per-particle tree-entry index (int array)
      - the per-particle columns in ``_PART_COLS``
      - ``event``: a dict of per-event (one row per tree entry) columns, keyed by the
        entries in ``_EVT_COLS`` plus ``entry`` (the tree-entry index)
      - ``nbad``: number of entries that failed to read

    All momenta/energies are in GeV; the GHEP particle positions (``vx..vt``) are in the
    native GENIE GHEP units (fermi / yoctoseconds), while the event vertex (``vtx_*``) is
    the interaction vertex as stored by GENIE.
    """
    _ensure_genie_loaded()
    import ROOT

    n = ROOT.cafpyana_read_genie_evtrec(str(path))
    if n < 0:
        raise IOError("Could not open GenieEvtRecTree in file: %s" % path)

    out = {"entry": _v2np(ROOT.g_entry, np.int32)}
    for name, (gname, dt) in _PART_COLS.items():
        out[name] = _v2np(getattr(ROOT, gname), dt)

    event = {"entry": _v2np(ROOT.g_evt_entry, np.int32)}
    for name, (gname, dt) in _EVT_COLS.items():
        event[name] = _v2np(getattr(ROOT, gname), dt)
    out["event"] = event
    out["nbad"] = int(ROOT.g_nbad)
    return out
