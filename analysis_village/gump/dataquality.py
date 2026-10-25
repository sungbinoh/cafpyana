"""Beam-quality and data-quality (good-run) cuts for the GUMPLE data samples.

Beam quality (per spill) follows the good-beam definition of data/good_beam_lists.md
(SBN-doc-49735-v1):

  1. beam intensity: TOR860 > 1e11 and TOR875 > 1e11 protons
  2. horn current:   173 <= THCURR <= 175 kA
  3. figure of merit: 0.98 < FOM <= 1          -- run sets in FOM_DETECTORS, see below

(The LM875A/B/C > 1e-2 loss-monitor requirement of the earlier SPINE tech-note
definition is no longer applied, 2026-10-20.)

The FOM is the CORRECTED one (FOM_best) from the beam-quality ROOT files in data/
(FOM_FILES, written by bnb_fom_caf.py), joined to each bnb spill on the
exact (spill_time_sec, spill_time_nsec) stamp as in match_fom.py. Spills those
files do not cover keep the production FOM. See correct_fom().

POT is counted as the sum of TOR875 over surviving spills (spills with a
non-finite or non-positive TOR875 -- IFBeam query failures -- never count).
TOR875 is the correct POT counter for both detectors. NB the POT quoted in
data/good_beam_lists.md is the beam-quality files' `tor_used` (TOR860 for
ICARUS), about 0.5% lower; ignore it -- those numbers are not a reference for
the POT computed here.

Data quality (per run) uses the good-run lists in data/run_lists.md (a copy of
the conclusive prefiltered / tagged-bad / good lists for ICARUS Run 2, ICARUS
Run 4 and SBND Run 1; 2026-10-20 version, with ICARUS Run 4 and SBND tagged on
per-event DQ metrics and SBND covering the 1e20 and post-1e20 periods). A run is
good iff it is in the "Good" list -- for SBND the combined "Good, SBND Run 1
total" list. The reason a run is cut is its prefilter category or the DQ
metrics that tagged it.

Both cuts apply to DATA ONLY. The beam-quality cut needs spill information, so
it applies to on-beam data only; the good-run cut applies to on-beam and
off-beam data alike.

Some productions (EXTERNAL_BEAM, e.g. the ICARUS Run 4 ReCAF2026 sample) carry
no spill information at all: their bnb tables are empty and their hdr POT is 0.
For those the per-spill table and each event's spill come from the
beam-quality ROOT file instead (its bnbinfo and spill trees). See spill_table()
and external_event_spill().
"""
import functools
import importlib.util
import os
import re

import h5py
import numpy as np
import pandas as pd

# match_fom.py sits next to this file. Load it by path rather than putting this
# directory first on sys.path, where gump's makedf.py would shadow the package.
_spec = importlib.util.spec_from_file_location(
    "match_fom", os.path.join(os.path.dirname(os.path.abspath(__file__)), "match_fom.py"))
mf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mf)

RUN_LIST_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "run_lists.md")

# loaddf detector label -> section heading in run_lists.md
RUN_SETS = {"SBND": "SBND Run 1", "ICARUS Run2": "ICARUS Run 2", "ICARUS Run4": "ICARUS Run 4"}

# (prefiltered, tagged bad, good) run counts quoted in the summary table of
# run_lists.md -- the parser is checked against them (SBND: the Run 1 total).
_EXPECTED_COUNTS = {"ICARUS Run2": (500, 18, 229), "ICARUS Run4": (112, 12, 109), "SBND": (382, 9, 102)}

# Header-documented hand-curated exception to the ICARUS Run 2 good list
# (the other one, 9610, is voltage-tagged but kept, and so is in the Good list).
_EXTRA_REASONS = {"ICARUS Run2": {9435: ("Excluded", "Handscan run")}}

NOT_IN_LIST = ("Not in run list", "Run absent from the run list")

# ============================================================
# Run list parsing
# ============================================================
def _ints(block):
    return [int(x) for x in block.split()]

def parse_run_lists(fname=RUN_LIST_FILE):
    """Parse run_lists.md.

    Returns {detector: dict(good=set, reasons={run: (category, reason)})}, where
    `reasons` covers every run in the list that is NOT good: prefiltered runs
    map to ("Prefiltered", <prefilter category>) and tagged runs to
    ("Tagged bad", <DQ metrics>). Hand-curated Run 2 exceptions are included.
    """
    with open(fname) as f:
        text = f.read()

    # split into run-set sections on the level-2 headings
    sections = re.split(r"^## ", text, flags=re.M)
    out = {}
    for det, heading in RUN_SETS.items():
        sec = [s for s in sections if s.startswith(heading + "\n")]
        assert len(sec) == 1, "run set %r not found in %s" % (heading, fname)
        sec = sec[0]
        # level-3 and level-4 subsections. SBND has one Prefiltered block, a
        # Tagged bad / Good pair per period (1e20, post-1e20) and the combined
        # "Good, SBND Run 1 total" list, which is the one applied.
        sub = [(s.split("\n", 1)[0].strip(), s) for s in re.split(r"^#{3,4} ", sec, flags=re.M)[1:]]
        pre = [v for k, v in sub if k.startswith("Prefiltered")]
        tags = [v for k, v in sub if k.startswith("Tagged bad")]
        good = [v for k, v in sub if k.startswith("Good, ")] or [v for k, v in sub if k.startswith("Good")]
        assert len(pre) == 1 and len(tags) >= 1 and len(good) == 1, (det, [k for k, _ in sub])
        pre, good = pre[0], good[0]

        reasons = {}
        # prefilter categories: **Category** (N) followed by a code block
        for cat, n, block in re.findall(r"^\*\*(.+?)\*\* \((\d+)\)\s*\n+```\n(.*?)```", pre, flags=re.M | re.S):
            runs = _ints(block)
            assert len(runs) == int(n), (det, cat, len(runs), n)
            for r in runs:
                # a run may appear in several categories -- keep them all
                if r in reasons:
                    reasons[r] = ("Prefiltered", reasons[r][1] + "; " + cat)
                else:
                    reasons[r] = ("Prefiltered", cat)
        n_pre = len(reasons)

        # tagged: the "By run" tables, rows "| run | nmetric | metrics |"
        n_tag = 0
        for tag in tags:
            byrun = tag.split("**By metric**")[0]
            for r, _nm, metrics in re.findall(r"^\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(.+?)\s*\|\s*$", byrun, flags=re.M):
                reasons[int(r)] = ("Tagged bad", metrics)
                n_tag += 1

        goodruns = set(_ints(re.search(r"```\n(.*?)```", good, flags=re.S).group(1)))

        assert (n_pre, n_tag, len(goodruns)) == _EXPECTED_COUNTS[det], \
            (det, (n_pre, n_tag, len(goodruns)), _EXPECTED_COUNTS[det])

        for r, why in _EXTRA_REASONS.get(det, {}).items():
            reasons[r] = why
        # a tagged run kept by hand (ICARUS Run 2 9610) is good, not cut
        for r in goodruns:
            reasons.pop(r, None)

        out[det] = dict(good=goodruns, reasons=reasons)
    return out

# Run sets the figure-of-merit cut, 0.98 < FOM <= 1, is applied to. The cut is
# bounded on both sides so the failure flags -1/2/3/4/-999 and the ICARUS Run 4
# +100 assumed-width tier of the production FOM are rejected (spills not covered
# by FOM_FILES still carry the production value). With the corrected FOM the
# ICARUS -999s (34% of the Run 2 POT in sbn-rewgted-24) and the Run 4 +100 tier
# are gone, so the cut applies to every run set. A spill with no formable
# corrected FOM (NaN) fails it. (With the *_fx files Run 2 runs 9941-9945, which
# had no multiwire fits at all in the 2026-10-02 file, have a FOM again.)
FOM_DETECTORS = {"SBND", "ICARUS Run2", "ICARUS Run4"}

# Beam-quality ROOT files carrying the corrected FOM (FOM_best), per run set
# (the *_fx files of the 2026-10 recalculation behind data/good_beam_lists.md,
# shipped with sbn-rewgted-25 and linked into data/). The three SBND files are
# disjoint -- the dev file holds the FixedDev spills (runs 18255, 18259),
# part1/part2 the FullOnBeam ones; SBND RollingDev is not covered (keeps the
# production FOM). The Run 2 full file contains every spill of the 10% (unblind)
# sample, and the Run 4 full file covers the ReCAF2026 sample and the unblind one.
_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
FOM_FILES = {
    "SBND": [os.path.join(_DATA_DIR, f) for f in
             ("sbnd_run1_dev_fixed_fx.root", "sbnd_run1_part1_fixed_fx.root", "sbnd_run1_part2_fixed_fx.root")],
    "ICARUS Run2": [os.path.join(_DATA_DIR, "icarus_run2_full_fixed_fx_hptg.root")],
    "ICARUS Run4": [os.path.join(_DATA_DIR, "icarus_run4_full_fixed_fx.root")],
}

# On-beam productions that carry no spill information (empty bnb tables, NaN
# evt spill_* columns, hdr POT 0), keyed by .df basename: the beam-quality ROOT
# file to take the spills (bnbinfo tree) and each event's spill (spill tree)
# from instead. ICARUS_SpringRun4BNB_ReCAF2026.df is the full Run 4 dataset
# (2963936 triggers, 136 runs), exactly the triggered-event set of that file.
EXTERNAL_BEAM = {
    "ICARUS_SpringRun4BNB_ReCAF2026.df": os.path.join(_DATA_DIR, "icarus_run4_full_fixed_fx.root"),
}

def external_beam_file(fname):
    """The EXTERNAL_BEAM ROOT file for the .df `fname`, or None."""
    return EXTERNAL_BEAM.get(os.path.basename(fname))

RUN_LISTS = parse_run_lists()
GOOD_RUNS = {det: v["good"] for det, v in RUN_LISTS.items()}

def cut_reason(run, detector):
    """(category, reason) for a run the good-run cut removes; None for a good run."""
    if run in GOOD_RUNS[detector]:
        return None
    return RUN_LISTS[detector]["reasons"].get(run, NOT_IN_LIST)

# Identifies the run list (and the FOM run sets and files) in loaddf cache keys,
# so editing any of them busts the data caches.
with open(RUN_LIST_FILE, "rb") as _f:
    import hashlib
    _fsig = lambda p: (os.path.basename(p), os.path.getsize(p) if os.path.exists(p) else None)
    RUN_LIST_HASH = hashlib.sha256(_f.read() + repr((
        sorted((d, sorted(r)) for d, r in GOOD_RUNS.items()),
        sorted(FOM_DETECTORS),
        "no LM875 cut",
        sorted((d, [_fsig(p) for p in ps]) for d, ps in FOM_FILES.items()),
        sorted((k, _fsig(p)) for k, p in EXTERNAL_BEAM.items()))).encode()).hexdigest()[:12]

# ============================================================
# Cuts
# ============================================================
def data_quality_cut(runs, detector):
    """True for runs on the good-run list. `runs` is array-like of run numbers."""
    return np.isin(np.asarray(runs), list(GOOD_RUNS[detector]))

def beam_quality_cut(df, detector, prefix="spill_"):
    """Per-spill beam-quality cut (good_beam_lists.md), on columns <prefix>TOR860 etc.

    The FOM requirement is applied only for run sets in FOM_DETECTORS.

    Evaluated on the per-spill bnb table with prefix="" or on the evt frame
    (whose spill_* columns carry the spill matched to each event) with the
    default prefix. NaN (no matched spill) fails.
    """
    v = lambda n: df[prefix + n]
    intensity = (v("TOR860") > 1e11) & (v("TOR875") > 1e11)
    horn = (v("THCURR") >= 173) & (v("THCURR") <= 175)
    if detector in FOM_DETECTORS:
        fom = (v("FOM") > 0.98) & (v("FOM") <= 1)
        return intensity & horn & fom
    return intensity & horn

def valid_tor(bnb):
    """Spills whose TOR875 can be counted towards POT. TOR875 is the correct
    POT counter (not the TOR860-based `tor_used` behind the good_beam_lists.md
    POT numbers, which are to be ignored)."""
    return np.isfinite(bnb.TOR875) & (bnb.TOR875 > 0)

# ============================================================
# Corrected FOM
# ============================================================
@functools.lru_cache(maxsize=None)
def fom_lookup(detector):
    """(sorted unique spill stamp, {column: array}) from FOM_FILES[detector], or
    None (no correction) when the run set has no files or one is missing."""
    paths = FOM_FILES.get(detector, [])
    missing = [p for p in paths if not os.path.exists(p)]
    if not paths or missing:
        print("dataquality: no corrected FOM for %s (missing %s); keeping the production FOM"
              % (detector, missing or "FOM_FILES entry"))
        return None
    stamp, values, _ = mf.load_lookup(paths, "bnbinfo", "FOM_best", ["TOR875"])
    return stamp, values

def correct_fom(bnb, detector, log=print):
    """Overwrite bnb.FOM with the corrected FOM (FOM_best) of each matched spill.

    Joined on the exact (spill_time_sec, spill_time_nsec) stamp, as in
    match_fom.py. Unmatched spills keep the production FOM; a matched spill with
    no formable FOM gets NaN (and so fails the FOM cut). Returns the matched mask.

    The join is verified on TOR875 (and run, when present), which the
    beam-quality files carry unaltered. The production FOM is NOT a join check:
    the *_fx files' own `FOM` column no longer reproduces it bug-for-bug (e.g.
    ~84% on ICARUS Run 4 unblind, where the 10% file gave 99.999%), so its
    agreement is reported for information only.
    """
    lut = fom_lookup(detector)
    if lut is None or len(bnb) == 0:
        return np.zeros(len(bnb), dtype=bool)
    stamp_lut, values = lut
    stamp = mf.ns_stamp(bnb.spill_time_sec, bnb.spill_time_nsec)
    idx = mf.match(stamp, stamp_lut, 0)
    m = idx >= 0
    if m.any():
        old = bnb.FOM.to_numpy(dtype=np.float64)[m]
        ref = values["FOM"][idx[m]].astype(np.float64)
        both = np.isfinite(old) & np.isfinite(ref)
        repro = np.isclose(old[both], ref[both], atol=1e-4, rtol=0).mean() if both.any() else np.nan
        tor = bnb.TOR875.to_numpy(dtype=np.float64)[m]
        tref = values["TOR875"][idx[m]].astype(np.float64)
        tboth = np.isfinite(tor) & np.isfinite(tref)
        tor_ok = np.isclose(tor[tboth], tref[tboth], rtol=1e-5, atol=0).mean() if tboth.any() else np.nan
        msg = "corrected FOM: matched %d/%d spills (%.2f%%), TOR875 agreement %.4f%% (production FOM reproduced %.2f%%)" % (
            m.sum(), len(bnb), 100*m.mean(), 100*tor_ok, 100*repro)
        if "run" in bnb.columns:
            agree = (bnb.run.to_numpy()[m].astype(np.int64) == values["run"][idx[m]].astype(np.int64)).mean()
            msg += ", run agreement %.4f%%" % (100*agree)
            if agree < 1:
                msg += "  <-- WRONG JOIN"
        if tor_ok < 0.999:
            msg += "  <-- join or file mismatch"
        log(msg)
        fom = bnb.FOM.to_numpy(dtype=np.float64).copy()
        fom[m] = values["FOM_best"][idx[m]]
        bnb["FOM"] = fom
    else:
        log("corrected FOM: matched 0/%d spills; keeping the production FOM" % len(bnb))
    return m

def event_spill_fom(bnb, hdr, cols=("FOM", "TOR875")):
    """The (corrected) FOM and TOR875 (or the bnb columns `cols`) of each event's
    spill, indexed (__ntuple, entry), with `mapped` False for events that have no
    spill row.

    The evt frame's spill_* columns carry the LAST bnb row of the event's
    (run, subrun, event) group -- checked on sbn-rewgted-24 SBND FixedDev and
    ICARUS Run 4 to hold for every event. `bnb` is a spill_table() of the same split.
    """
    last = bnb.drop_duplicates(["run", "subrun", "event"], keep="last")
    last = last.set_index(["run", "subrun", "event"])[list(cols)].assign(mapped=True)
    key = pd.MultiIndex.from_frame(hdr[["run", "subrun", "evt"]].astype(np.int64))
    out = last.reindex(key)
    out.index = hdr.index
    out["mapped"] = out["mapped"].eq(True)
    return out

# ============================================================
# Per-spill tables and normalization
# ============================================================
def _keys(fname, prefix):
    with h5py.File(fname, "r") as f:
        return sorted([k for k in f.keys() if re.fullmatch(prefix + r"_\d+", k)],
                      key=lambda k: int(k.split("_")[-1]))

_BEAM_COLS = ["TOR860", "TOR875", "THCURR", "FOM_best", "run", "subrun", "spill_time_sec", "spill_time_nsec"]

@functools.lru_cache(maxsize=None)
def _external_spills(path):
    """Every spill of an EXTERNAL_BEAM file (its bnbinfo tree), FOM = FOM_best.

    Duplicate spill records (same exact time stamp) are kept once, preferring a
    finite FOM_best -- the same resolution as match_fom.load_lookup -- so a
    repeated record never counts its TOR875 twice.
    """
    import uproot
    t = uproot.open(path)["bnbinfo"]
    a = t.arrays(_BEAM_COLS, library="np")
    stamp = mf.ns_stamp(a["spill_time_sec"], a["spill_time_nsec"])
    order = np.lexsort((~np.isfinite(a["FOM_best"]), stamp))
    srt = stamp[order]
    keep = order[np.concatenate(([True], srt[1:] != srt[:-1]))]
    keep.sort()
    bnb = pd.DataFrame({k: a[k][keep] for k in _BEAM_COLS}).rename(columns={"FOM_best": "FOM"})
    bnb["run"] = bnb.run.astype(int)
    bnb["subrun"] = bnb.subrun.astype(int)
    print("dataquality: %s: %d spills (%d duplicate records dropped)"
          % (os.path.basename(path), len(bnb), len(stamp) - len(keep)))
    return bnb

@functools.lru_cache(maxsize=None)
def _external_event_spills(path):
    """Each triggered event's spill from an EXTERNAL_BEAM file (its spill tree),
    indexed (run, subrun, evt). The tree repeats some events with identical
    spill time and FOM_best (good_beam_lists.md); those are kept once."""
    import uproot
    t = uproot.open(path)["spill"]
    cols = ["TOR860", "TOR875", "THCURR", "FOM_best", "run", "subrun", "evt"]
    ev = pd.DataFrame(t.arrays(cols, library="np")).rename(columns={"FOM_best": "FOM"})
    for c in ("run", "subrun", "evt"):
        ev[c] = ev[c].astype(np.int64)
    n = len(ev)
    ev = ev.drop_duplicates(["run", "subrun", "evt"]).set_index(["run", "subrun", "evt"])
    print("dataquality: %s: %d triggered events (%d duplicate records dropped)"
          % (os.path.basename(path), len(ev), n - len(ev)))
    return ev

def external_event_spill(fname, hdr):
    """The beam readings (TOR860, TOR875, THCURR, corrected FOM) of each event's
    spill for an EXTERNAL_BEAM file, indexed like `hdr`, with `mapped` False
    for events absent from the spill tree (these count as not good)."""
    ev = _external_event_spills(external_beam_file(fname))
    key = pd.MultiIndex.from_frame(hdr[["run", "subrun", "evt"]].astype(np.int64))
    out = ev.reindex(key)
    out.index = hdr.index
    out["mapped"] = np.isfinite(out.TOR875.to_numpy(dtype=np.float64))
    return out

def _split_owner(fname, idf):
    """(run, subrun) pairs whose spills split `idf` owns: those in its header
    and in no earlier split's, so a subrun spread over splits counts once."""
    seen = set()
    for k in _keys(fname, "hdr"):
        i = int(k.split("_")[-1])
        rs = set(map(tuple, pd.read_hdf(fname, k)[["run", "subrun"]].drop_duplicates().astype(int).to_numpy()))
        if i == idf:
            return rs - seen
        seen |= rs
    return set()

def spill_table(fname, idf, detector):
    """The per-spill bnb table of one split, with run/subrun, the corrected FOM
    (correct_fom) and the cut flags.

    bnb rows are indexed (__ntuple, entry, spill) with (__ntuple, entry) the
    first_in_subrun header record that carries them. For an EXTERNAL_BEAM file
    the rows are instead the beam-quality file's spills of the (run, subrun)s
    the split owns, already carrying the corrected FOM.
    """
    bnb = pd.read_hdf(fname, "bnb_%i" % idf)
    ext = external_beam_file(fname)
    if ext is not None and len(bnb) == 0:
        allsp = _external_spills(ext)
        own = pd.MultiIndex.from_tuples(sorted(_split_owner(fname, idf)), names=["run", "subrun"])
        bnb = allsp[pd.MultiIndex.from_frame(allsp[["run", "subrun"]]).isin(own)].copy()
        print("[%s idf=%i] spills from %s: %d (%d runs)"
              % (os.path.basename(fname), idf, os.path.basename(ext), len(bnb), bnb.run.nunique()))
        bnb["valid"] = valid_tor(bnb)
        bnb["bq"] = beam_quality_cut(bnb, detector, prefix="") & bnb.valid
        return bnb
    hdr = pd.read_hdf(fname, "hdr_%i" % idf)
    rs = hdr[["run", "subrun"]]
    bnb = bnb.join(rs, on=["__ntuple", "entry"])
    assert not bnb.run.isna().any(), "bnb rows without a header record in %s split %i" % (fname, idf)
    bnb["run"] = bnb.run.astype(int)
    bnb["subrun"] = bnb.subrun.astype(int)
    correct_fom(bnb, detector, log=lambda m: print("[%s idf=%i] %s" % (os.path.basename(fname), idf, m)))
    bnb["valid"] = valid_tor(bnb)
    bnb["bq"] = beam_quality_cut(bnb, detector, prefix="") & bnb.valid
    return bnb

def split_onbeam_pot(fname, idf, detector, beam_quality=True, data_quality=True):
    """POT of one on-beam split after the requested cuts (sum of TOR875)."""
    bnb = spill_table(fname, idf, detector)
    m = bnb.valid.copy()
    if beam_quality:
        m &= bnb.bq
    if data_quality:
        m &= data_quality_cut(bnb.run, detector)
    return float(bnb.TOR875[m].astype(float).sum())

def per_run_spills(fname, detector):
    """Per-run spill counts and POT over all splits of ONE on-beam file.

    Columns: nspill (all), nspill_valid, nspill_bq, pot (valid TOR875), pot_bq,
    good (on the good-run list). NB different on-beam samples of a run set (e.g.
    SBND FullOnBeam and FixedDev) overlap, so never sum these across samples.
    """
    rows = []
    for k in _keys(fname, "bnb"):
        b = spill_table(fname, int(k.split("_")[-1]), detector)
        tor = b.TOR875.astype(float)
        rows.append(pd.DataFrame({
            "run": b.run, "nspill": 1, "nspill_valid": b.valid.astype(int), "nspill_bq": b.bq.astype(int),
            "pot": tor.where(b.valid, 0.), "pot_bq": tor.where(b.bq, 0.)}))
    per = pd.concat(rows).groupby("run").sum()
    per["good"] = data_quality_cut(per.index, detector)
    return per

def pot_stages(fname, detector):
    """Cumulative POT (and spills) initially, after beam quality, after data quality."""
    per = per_run_spills(fname, detector)
    hdrpot = sum(pd.read_hdf(fname, k).pot.sum() for k in _keys(fname, "hdr"))
    g = per[per.good]
    return dict(
        pot_hdr=float(hdrpot),
        pot_initial=float(per.pot.sum()), pot_bq=float(per.pot_bq.sum()), pot_dq=float(g.pot_bq.sum()),
        nspill_initial=int(per.nspill.sum()), nspill_bq=int(per.nspill_bq.sum()), nspill_dq=int(g.nspill_bq.sum()),
        nrun_initial=int(len(per)), nrun_dq=int(len(g)))

def _read_all(fname, prefix):
    return pd.concat([pd.read_hdf(fname, k) for k in _keys(fname, prefix)])

def _with_run(tab, hdr):
    return tab.join(hdr[["run"]], on=["__ntuple", "entry"]) if "run" not in tab.columns else tab

def data_norm(detector, onbeam, offbeam_files, log=print, beam_quality=True, data_quality=True,
              runs=None):
    """Gate counts, POT and off-beam weight for the on-beam file `onbeam`, after cuts.

    The same normalization mcdata_comparison_gumple.compute_norm always did, with
    the cuts folded in:
      * POT = sum TOR875 over spills passing beam quality in good runs.
      * SBND gates ON = number of such spills.
      * ICARUS gates ON = sum trig.gate_delta over good-run events, times the
        fraction of (good-run) spills passing beam quality -- gate_delta counts
        gates between triggers, so it cannot be split per spill directly.
      * gates OFF = the usual off-beam counts restricted to good runs.
      * NEVT_ON / NEVT_OFF (triggers) likewise, with the beam-quality spill fraction
        applied to NEVT_ON.

    runs, if given, restricts everything (on- and off-beam) to those run numbers,
    e.g. to normalize one time-ordered slice of a run period.
    """
    def in_runs(r):
        return np.ones(len(r), bool) if runs is None else np.isin(np.asarray(r), list(runs))

    hdr_on = _read_all(onbeam, "hdr")
    per = per_run_spills(onbeam, detector)
    good_on = data_quality_cut(hdr_on.run, detector) if data_quality else np.ones(len(hdr_on), bool)
    good_on &= in_runs(hdr_on.run)
    per = per[in_runs(per.index)]
    pr = per[per.good] if data_quality else per
    bq_frac = (pr.nspill_bq.sum() / pr.nspill.sum()) if beam_quality else 1.
    pot = float(pr.pot_bq.sum() if beam_quality else pr.pot.sum())

    offs = []
    for f in offbeam_files:
        h = _read_all(f, "hdr")
        good = data_quality_cut(h.run, detector) if data_quality else np.ones(len(h), bool)
        good &= in_runs(h.run)
        offs.append((f, h, good))

    if "ICARUS" in detector:
        trig_on = _with_run(_read_all(onbeam, "trig"), hdr_on)
        tgood = data_quality_cut(trig_on.run, detector) if data_quality else np.ones(len(trig_on), bool)
        tgood &= in_runs(trig_on.run)
        # (1-1/200.) / (1-1/40.): gate-loss corrections, as before
        gate_corr = (1.-1/200.) if detector == "ICARUS Run2" else (1.-1/40.)
        ngates_ON = float(trig_on.gate_delta[tgood].sum())*bq_frac*gate_corr
        ngates_OFF = 0.
        for f, h, good in offs:
            t = _with_run(_read_all(f, "trig"), h)
            tg = data_quality_cut(t.run, detector) if data_quality else np.ones(len(t), bool)
            tg &= in_runs(t.run)
            ngates_OFF += float(t.gate_delta[tg].sum())*(1-1/20.)
        off_w = ngates_ON / ngates_OFF
    else:
        ngates_ON = float(pr.nspill_bq.sum() if beam_quality else pr.nspill.sum())
        ngates_OFF = float(sum(h.noffbeambnb[good].sum() for _, h, good in offs))
        f_factor = 0.0754
        off_w = (1. - f_factor) * ngates_ON / ngates_OFF

    nevt_ON = float(good_on.sum())*bq_frac
    nevt_OFF = float(sum(good.sum() for _, _, good in offs))
    nevt = nevt_ON - nevt_OFF*off_w

    log("data cuts: beam_quality=%s data_quality=%s (FOM cut %s)"
        % (beam_quality, data_quality, "applied" if detector in FOM_DETECTORS else "not applied, see dataquality.FOM_DETECTORS"))
    log("  on-beam: %d runs -> %d good; spill BQ fraction (good runs) = %.4f"
        % (len(per), int(per.good.sum()), bq_frac))
    log("  on-beam POT: initial %.4e  after BQ %.4e  after BQ+DQ %.4e"
        % (per.pot.sum(), per.pot_bq.sum(), per[per.good].pot_bq.sum()))
    for f, h, good in offs:
        log("  off-beam %s: %d/%d triggers in good runs" % (os.path.basename(f), int(good.sum()), len(h)))
    log("ngates_ON = %r, ngates_OFF = %r, OFF_w = %r" % (ngates_ON, ngates_OFF, off_w))
    log("POT = %r" % pot)
    log("N GATES ON / 5e12 POT")
    log("%r" % (5e12*ngates_ON/pot))
    log("NEVT_ON = %r, POT = %r, NEVT_ON/1e15 POT = %r" % (nevt_ON, pot, nevt_ON / (pot / 1e15)))
    log("NEVT = %r, POT = %r, NEVT/1e15 POT = %r" % (nevt, pot, nevt / (pot / 1e15)))

    return dict(ngates_ON=ngates_ON, ngates_OFF=ngates_OFF, OFF_w=off_w, POT=pot,
                NEVT_ON=nevt_ON, NEVT_OFF=nevt_OFF, NEVT=nevt, BQ_FRAC=bq_frac)
