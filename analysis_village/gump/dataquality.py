"""Beam-quality and data-quality (good-run) cuts for the GUMPLE data samples.

Beam quality (per spill) follows Sec. 6.1 of the SBN SPINE numu-disappearance
technical note (itself following the ICARUS single-detector analysis):

  1. beam intensity: TOR860 > 1e11 and TOR875 > 1e11 protons,
                     LM875A, LM875B, LM875C > 1e-2
  2. horn current:   173 <= THCURR <= 175 kA
  3. figure of merit: 0.98 < FOM <= 1          -- NOT APPLIED, see TODO below

POT is counted as the sum of TOR875 over surviving spills (spills with a
non-finite or non-positive TOR875 -- IFBeam query failures -- never count).

Data quality (per run) uses the good-run lists in data/run_lists.md (a copy of
the conclusive prefiltered / tagged-bad / good lists for ICARUS Run 2, ICARUS
Run 4 and SBND Run 1). A run is good iff it is in the "Good" list; the reason a
run is cut is its prefilter category or the DQ metrics that tagged it. Prefilter
categories listed in IGNORED_PREFILTER (ICARUS "No End Time") are not applied:
runs removed only by them count as good. Runs in EXTRA_GOOD_RUNS (SBND 18255)
are added to the good list by hand.

Both cuts apply to DATA ONLY. The beam-quality cut needs spill information, so
it applies to on-beam data only; the good-run cut applies to on-beam and
off-beam data alike.
"""
import os
import re

import h5py
import numpy as np
import pandas as pd

RUN_LIST_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "run_lists.md")

# loaddf detector label -> section heading in run_lists.md
RUN_SETS = {"SBND": "SBND Run 1", "ICARUS Run2": "ICARUS Run 2", "ICARUS Run4": "ICARUS Run 4"}

# (prefiltered, tagged bad, good) run counts quoted in the summary table of
# run_lists.md -- the parser is checked against them.
_EXPECTED_COUNTS = {"ICARUS Run2": (500, 18, 229), "ICARUS Run4": (162, 9, 62), "SBND": (318, 7, 59)}

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
        sub = {s.split("\n", 1)[0].strip(): s for s in re.split(r"^### ", sec, flags=re.M)[1:]}
        pre = [v for k, v in sub.items() if k.startswith("Prefiltered")][0]
        tag = [v for k, v in sub.items() if k.startswith("Tagged bad")][0]
        good = [v for k, v in sub.items() if k.startswith("Good")][0]

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

        # tagged: the "By run" table, rows "| run | nmetric | metrics |"
        byrun = tag.split("**By metric**")[0]
        n_tag = 0
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

# Prefilter categories NOT applied as a cut, per run set. ICARUS "No End Time" runs
# (all in Run 4) are missing an end time in the run database -- a bookkeeping gap,
# not a detector problem: in sbn-rewgted-24 they carry as much beam as the good
# runs, at the same trigger and slice rates per POT. They are treated as good.
IGNORED_PREFILTER = {"ICARUS Run2": {"No End Time"}, "ICARUS Run4": {"No End Time"}}

def _apply_ignored(lists):
    rescued = {}
    for det, v in lists.items():
        ign = IGNORED_PREFILTER.get(det, set())
        rescued[det] = set()
        for r, (cat, why) in list(v["reasons"].items()):
            # rescue only runs whose every prefilter category is ignored
            if cat == "Prefiltered" and set(why.split("; ")) <= ign:
                rescued[det].add(r)
                del v["reasons"][r]
        v["good"] = v["good"] | rescued[det]
    return rescued

# Runs added to the good-run list by hand, per run set: {run: reason}. SBND 18255
# (inside the SBND Run 1 range 18250-18412) is absent from every category of the
# run list -- it was apparently never evaluated -- and looks like a normal physics
# run in sbn-rewgted-24 (19.6 triggers / 1e15 POT, 99% of its POT passes beam
# quality). It carries 74% of the SBND FixedDev POT.
EXTRA_GOOD_RUNS = {"SBND": {18255: "Absent from the run list; physics-like, added to the good runs by hand"}}

RUN_LISTS = parse_run_lists()
RESCUED_RUNS = _apply_ignored(RUN_LISTS)
for _det, _runs in EXTRA_GOOD_RUNS.items():
    RUN_LISTS[_det]["good"] = RUN_LISTS[_det]["good"] | set(_runs)
    for _r in _runs:
        RUN_LISTS[_det]["reasons"].pop(_r, None)
GOOD_RUNS = {det: v["good"] for det, v in RUN_LISTS.items()}

def cut_reason(run, detector):
    """(category, reason) for a run the good-run cut removes; None for a good run."""
    if run in GOOD_RUNS[detector]:
        return None
    return RUN_LISTS[detector]["reasons"].get(run, NOT_IN_LIST)

# Identifies the run list in loaddf cache keys, so editing it busts the data caches.
with open(RUN_LIST_FILE, "rb") as _f:
    import hashlib
    RUN_LIST_HASH = hashlib.sha256(_f.read() + repr((
        sorted((d, sorted(c)) for d, c in IGNORED_PREFILTER.items()),
        sorted((d, sorted(r)) for d, r in EXTRA_GOOD_RUNS.items()))).encode()).hexdigest()[:12]

# ============================================================
# Cuts
# ============================================================
def data_quality_cut(runs, detector):
    """True for runs on the good-run list. `runs` is array-like of run numbers."""
    return np.isin(np.asarray(runs), list(GOOD_RUNS[detector]))

def beam_quality_cut(df, prefix="spill_"):
    """Per-spill beam-quality cut (tech note Sec. 6.1.2), on columns <prefix>TOR860 etc.

    Evaluated on the per-spill bnb table with prefix="" or on the evt frame
    (whose spill_* columns carry the spill matched to each event) with the
    default prefix. NaN (no matched spill) fails.
    """
    v = lambda n: df[prefix + n]
    intensity = (v("TOR860") > 1e11) & (v("TOR875") > 1e11) & \
                (v("LM875A") > 1e-2) & (v("LM875B") > 1e-2) & (v("LM875C") > 1e-2)
    horn = (v("THCURR") >= 173) & (v("THCURR") <= 175)
    # TODO: apply the figure-of-merit cut, 0.98 < FOM <= 1 (bounded on both sides so
    # the failure flags -1/2/3/4/-999 and the ICARUS Run 4 +100 assumed-width tier
    # are rejected). Disabled for now: in sbn-rewgted-24 34% of the ICARUS Run 2 POT
    # carries FOM = -999 (IFBeam query failed), so the cut would keep only ~56% of
    # the Run 2 exposure (vs 90% in the tech note), and Run 4 would lose ~25% to the
    # +100 tier. Re-enable once the FOM values are understood:
    #   fom = (v("FOM") > 0.98) & (v("FOM") <= 1)
    #   return intensity & horn & fom
    return intensity & horn

def valid_tor(bnb):
    """Spills whose TOR875 can be counted towards POT."""
    return np.isfinite(bnb.TOR875) & (bnb.TOR875 > 0)

# ============================================================
# Per-spill tables and normalization
# ============================================================
def _keys(fname, prefix):
    with h5py.File(fname, "r") as f:
        return sorted([k for k in f.keys() if re.fullmatch(prefix + r"_\d+", k)],
                      key=lambda k: int(k.split("_")[-1]))

def spill_table(fname, idf):
    """The per-spill bnb table of one split, with run/subrun and the cut flags.

    bnb rows are indexed (__ntuple, entry, spill) with (__ntuple, entry) the
    first_in_subrun header record that carries them.
    """
    bnb = pd.read_hdf(fname, "bnb_%i" % idf)
    hdr = pd.read_hdf(fname, "hdr_%i" % idf)
    rs = hdr[["run", "subrun"]]
    bnb = bnb.join(rs, on=["__ntuple", "entry"])
    assert not bnb.run.isna().any(), "bnb rows without a header record in %s split %i" % (fname, idf)
    bnb["run"] = bnb.run.astype(int)
    bnb["valid"] = valid_tor(bnb)
    bnb["bq"] = beam_quality_cut(bnb, prefix="") & bnb.valid
    return bnb

def split_onbeam_pot(fname, idf, detector, beam_quality=True, data_quality=True):
    """POT of one on-beam split after the requested cuts (sum of TOR875)."""
    bnb = spill_table(fname, idf)
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
        b = spill_table(fname, int(k.split("_")[-1]))
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

def data_norm(detector, onbeam, offbeam_files, log=print, beam_quality=True, data_quality=True):
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
    """
    hdr_on = _read_all(onbeam, "hdr")
    per = per_run_spills(onbeam, detector)
    good_on = data_quality_cut(hdr_on.run, detector) if data_quality else np.ones(len(hdr_on), bool)
    pr = per[per.good] if data_quality else per
    bq_frac = (pr.nspill_bq.sum() / pr.nspill.sum()) if beam_quality else 1.
    pot = float(pr.pot_bq.sum() if beam_quality else pr.pot.sum())

    offs = []
    for f in offbeam_files:
        h = _read_all(f, "hdr")
        good = data_quality_cut(h.run, detector) if data_quality else np.ones(len(h), bool)
        offs.append((f, h, good))

    if "ICARUS" in detector:
        trig_on = _with_run(_read_all(onbeam, "trig"), hdr_on)
        tgood = data_quality_cut(trig_on.run, detector) if data_quality else np.ones(len(trig_on), bool)
        # (1-1/200.) / (1-1/40.): gate-loss corrections, as before
        gate_corr = (1.-1/200.) if detector == "ICARUS Run2" else (1.-1/40.)
        ngates_ON = float(trig_on.gate_delta[tgood].sum())*bq_frac*gate_corr
        ngates_OFF = 0.
        for f, h, good in offs:
            t = _with_run(_read_all(f, "trig"), h)
            tg = data_quality_cut(t.run, detector) if data_quality else np.ones(len(t), bool)
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

    log("data cuts: beam_quality=%s data_quality=%s (FOM cut disabled, see dataquality.py TODO)"
        % (beam_quality, data_quality))
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
