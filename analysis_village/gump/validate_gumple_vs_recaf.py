"""Validate the GUMPLE production (sbn-rewgted-19) against reCAF (sbn-rewgted-14).

Compares, per sample present in both productions, three tiers of quantities:

  INVARIANTS -- same underlying generated events, so these should agree closely:
    hdr POT sum, hdr row count, mcnu row count, mean true nu energy,
    genie interaction-mode fractions. Disagreement here signals a production
    problem, not a format change.

  BOOKKEEPING -- percent-level agreement expected:
    evt (slice) row count, slices passing the geometric slice-vertex FV cut
    (recomputed identically on both formats from slc_vtx_*).

  EXPECTED SHIFTS -- reported side by side, no pass/fail:
    slices passing nu_score > 0.6; slices passing the chi2 PID cut; mean
    mu_len / mu_E of PID-passing slices. The GUMPLE production redefines the
    muon candidate (longest passing track, was best chi2 ratio) and aggregates
    the proton-candidate chi2 over all candidates, so event-by-event shifts in
    these ARE expected.

Reads evt/mcnu/hdr HDF5 keys directly (no loaddf, no caching); streams the evt
chunks so nothing large stays resident.

    python validate_gumple_vs_recaf.py [--out summary.txt]
"""

import argparse
import os
import sys

import h5py
import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (_HERE, os.path.abspath(os.path.join(_HERE, "..", "..")),
           os.path.join(_HERE, "..", "gumple")):
    #if _p not in sys.path:
     sys.path.insert(0, _p)

import gumple_cuts as gc

OLD_DIR = "/Users/gputnam/Work/osc/sbn-rewgted-14/"
NEW_DIR = "/Users/gputnam/Work/osc/sbn-rewgted-20/"

# (label, [old files], [new files]). Shard layouts differ between productions,
# so every metric is aggregated over the sample's file list.
SAMPLES = [
    ("ICARUS Run2 MC CV",
     ["ICARUSRun2_SpringMCOverlay_rewgt.df"],
     ["ICARUSRun2_SpringMCOverlay_rewgt_%i.df" % i for i in range(3)]),
    ("ICARUS Run4 MC CV",
     ["ICARUSRun4_SpringMCOverlay_rewgt_%i.df" % i for i in range(2)],
     ["ICARUSRun4_SpringMCOverlay_rewgt_%i.df" % i for i in range(4)]),
    # NB: -14 also carries SBNDMCCV_Nom.df; -19 extends the CV to 13 shards.
    # Totals therefore NEED NOT match -- the per-event invariants (mean Enu,
    # mode fractions) still must.
    ("SBND MC CV",
     ["SBNDMCCV_%i.df" % i for i in range(4)],
     ["SBNDMCCV_%i.df" % i for i in range(13)]),
    ("ICARUS Run2 dirt",
     ["ICARUS_Spring_Overlay_Dirt.df"],
     ["ICARUSRun2_Spring_Overlay_Dirt.df"]),
    ("ICARUS Run4 dirt",
     ["ICARUSRun4_Spring_Overlay_Dirt.df"],
     ["ICARUSRun4_Spring_Overlay_Dirt.df"]),
    ("ICARUS Run2 onbeam",
     ["ICARUS_SpringRun2BNB_unblind.df"],
     ["ICARUS_SpringRun2BNB_unblind.df"]),
    ("ICARUS Run2 offbeam",
     ["ICARUS_SpringRun2BNBOff_unblind.df"],
     ["ICARUS_SpringRun2BNBOff_unblind.df"]),
    ("SBND onbeam (full)",
     ["SBND_SpringBNBData_FullOnBeam.df"],
     ["SBND_SpringBNBData_FullOnBeam.df"]),
    ("SBND offbeam",
     ["SBND_SpringBNBOffData.df"],
     ["SBND_SpringBNBOffData.df"]),
]

MODE_LIST = [0, 10, 1, 2, 3]  # QE, MEC, RES, SIS/DIS, COH


def _keys(path, base):
    with h5py.File(path, "r") as f:
        return sorted((k for k in f.keys() if k.startswith(base + "_")
                       and k[len(base) + 1:].isdigit()),
                      key=lambda k: int(k.rsplit("_", 1)[1]))


def _fv_mask(evt):
    """Geometric slice-vertex FV cut, applied identically to both formats.

    Uses gumple_cuts._fv_cut (vertex margins by default), which includes the
    ICARUS WW cable exclusion -- the same geometry for old and new, so the
    comparison stays apples-to-apples.
    """
    vtx = pd.DataFrame({
        "detector": evt["detector"], "Run": evt["Run"],
        "x": evt["slc_vtx_x"], "y": evt["slc_vtx_y"], "z": evt["slc_vtx_z"],
    }, index=evt.index)
    return gc._fv_cut(vtx)


def _pid_proxy(evt):
    """Chi2-only PID proxy that runs on BOTH formats.

    gumple_cuts.pid_cut needs the base-muon-mask columns (mu_trackScore,
    mu_dist_start, ...) which the old (-14) frames do not carry, so the
    validation uses the legacy fixed-threshold chi2 cut instead.
    """
    return (evt["mu_chi2_of_mu_cand"] < 30) & \
           (evt["prot_chi2_of_mu_cand"] > 80) & \
           (evt["mu_len"] > 25) & \
           (evt["mu_chi2_of_prot_cand"] > 0) & \
           (evt["prot_chi2_of_prot_cand"] < 90)


def sample_metrics(files):
    m = dict(pot=0.0, n_hdr=0, n_mcnu=0, sum_nuE=0.0, n_evt=0, n_fv=0,
             n_nuscore=0, n_pid=0, sum_mu_len=0.0, sum_mu_E=0.0,
             n_gump_sel=0, has_gump_sel=False,
             mode_counts=np.zeros(len(MODE_LIST) + 1))
    for path in files:
        for k in _keys(path, "hdr"):
            hdr = pd.read_hdf(path, k)
            m["pot"] += float(hdr["pot"].sum())
            m["n_hdr"] += len(hdr)
            del hdr
        for k in _keys(path, "mcnu"):
            mc = pd.read_hdf(path, k)
            m["n_mcnu"] += len(mc)
            # mcnu frames carry 2-level MultiIndex columns; selecting a name can
            # return a single-column DataFrame
            nuE = mc["nu_E"]
            nuE = nuE.iloc[:, 0] if isinstance(nuE, pd.DataFrame) else nuE
            m["sum_nuE"] += float(nuE.sum())
            gm = mc["genie_mode"]
            gm = gm.iloc[:, 0] if isinstance(gm, pd.DataFrame) else gm
            for i, mode in enumerate(MODE_LIST):
                m["mode_counts"][i] += int((gm == mode).sum())
            m["mode_counts"][-1] += int((~gm.isin(MODE_LIST)).sum())
            del mc
        for k in _keys(path, "evt"):
            evt = pd.read_hdf(path, k)
            m["n_evt"] += len(evt)
            fv = _fv_mask(evt)
            m["n_fv"] += int(fv.sum())
            m["n_nuscore"] += int((fv & (evt["nu_score"] > 0.6)).sum())
            pid = fv & _pid_proxy(evt)
            m["n_pid"] += int(pid.sum())
            m["sum_mu_len"] += float(evt.loc[pid, "mu_len"].sum())
            m["sum_mu_E"] += float(evt.loc[pid, "mu_E"].sum())
            if "gump_sel" in evt.columns:
                m["has_gump_sel"] = True
                m["n_gump_sel"] += int(evt["gump_sel"].sum())
            del evt
    return m


def fmt_row(name, old, new, kind, is_count=True, flags=None):
    if old in (None, 0) and new in (None, 0):
        rel = ""
    elif old:
        rel = " (%+.2f%%)" % (100.0 * (new - old) / old)
    else:
        rel = " (old=0)"
    line = "  %-34s %-12s old=%-16s new=%-16s%s" % (
        name, kind,
        ("%d" % old) if is_count else ("%.6g" % old),
        ("%d" % new) if is_count else ("%.6g" % new), rel)
    if flags is not None and old and kind == "INVARIANT":
        if abs(new - old) / max(abs(old), 1e-30) > 0.01:
            flags.append("INVARIANT mismatch: %s (old=%r new=%r)" % (name, old, new))
            line += "   <<< FLAG"
    return line


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=None, help="Also write the report to this file")
    ap.add_argument("--old-dir", default=OLD_DIR)
    ap.add_argument("--new-dir", default=NEW_DIR)
    args = ap.parse_args(argv)

    lines, flags = [], []

    def emit(s):
        print(s, flush=True)
        lines.append(s)

    emit("GUMPLE (-19) vs reCAF (-14) production validation")
    emit("old: %s" % args.old_dir)
    emit("new: %s" % args.new_dir)
    emit("")
    emit("Expected-shift metrics move with the new muon-candidate definition")
    emit("(longest passing track, was best chi2 ratio) and the proton-candidate")
    emit("chi2 aggregation -- differences there are not errors.")

    for label, old_files, new_files in SAMPLES:
        olds = [os.path.join(args.old_dir, f) for f in old_files]
        news = [os.path.join(args.new_dir, f) for f in new_files]
        missing = [p for p in olds + news if not os.path.exists(p)]
        emit("")
        emit("== %s (old: %d file(s), new: %d file(s))" % (label, len(olds), len(news)))
        if missing:
            emit("  SKIPPED -- missing: %s" % ", ".join(os.path.basename(p) for p in missing))
            continue
        try:
            mo = sample_metrics(olds)
            mn = sample_metrics(news)
        except Exception as e:
            emit("  ERROR reading sample: %r" % e)
            flags.append("read error on %s: %r" % (label, e))
            continue

        same_shards = old_files == new_files or (
            label not in ("SBND MC CV",))  # SBND CV totals expected to differ
        inv = "INVARIANT" if same_shards else "info"
        emit(fmt_row("hdr POT sum", mo["pot"], mn["pot"], inv, is_count=False, flags=flags))
        emit(fmt_row("hdr rows", mo["n_hdr"], mn["n_hdr"], inv, flags=flags))
        if mo["n_mcnu"] or mn["n_mcnu"]:
            emit(fmt_row("mcnu rows", mo["n_mcnu"], mn["n_mcnu"], inv, flags=flags))
            meano = mo["sum_nuE"] / max(mo["n_mcnu"], 1)
            meann = mn["sum_nuE"] / max(mn["n_mcnu"], 1)
            emit(fmt_row("mean true Enu [GeV]", meano, meann, "INVARIANT",
                         is_count=False, flags=flags))
            fo = mo["mode_counts"] / max(mo["mode_counts"].sum(), 1)
            fn = mn["mode_counts"] / max(mn["mode_counts"].sum(), 1)
            emit("  %-34s %-12s old=%s new=%s" % (
                "mode fractions QE/MEC/RES/DIS/COH/oth", "INVARIANT",
                "/".join("%.3f" % x for x in fo), "/".join("%.3f" % x for x in fn)))
            if np.abs(fo - fn).max() > 0.01:
                flags.append("mode-fraction mismatch on %s" % label)
        emit(fmt_row("evt (slice) rows", mo["n_evt"], mn["n_evt"], "bookkeeping"))
        emit(fmt_row("slices in vtx FV", mo["n_fv"], mn["n_fv"], "bookkeeping"))
        emit(fmt_row("FV & nu_score>0.6", mo["n_nuscore"], mn["n_nuscore"], "exp-shift"))
        emit(fmt_row("FV & PID cut", mo["n_pid"], mn["n_pid"], "exp-shift"))
        emit(fmt_row("mean mu_len (PID) [cm]",
                     mo["sum_mu_len"] / max(mo["n_pid"], 1),
                     mn["sum_mu_len"] / max(mn["n_pid"], 1), "exp-shift", is_count=False))
        emit(fmt_row("mean mu_E (PID) [GeV]",
                     mo["sum_mu_E"] / max(mo["n_pid"], 1),
                     mn["sum_mu_E"] / max(mn["n_pid"], 1), "exp-shift", is_count=False))
        if mn["has_gump_sel"]:
            emit("  %-34s %-12s new=%d (no old-format equivalent)"
                 % ("gump_sel (new only)", "info", mn["n_gump_sel"]))

    emit("")
    emit("=" * 70)
    if flags:
        emit("FLAGS (%d):" % len(flags))
        for f in flags:
            emit("  - %s" % f)
    else:
        emit("FLAGS: none -- invariants agree between the productions.")

    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w") as f:
            f.write("\n".join(lines) + "\n")
        print("wrote %s" % args.out, flush=True)

    return 0


if __name__ == "__main__":
    sys.exit(main())
