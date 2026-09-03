"""Serial runner for the GUMPLE (sbn-rewgted-20-calovarB) analyses.

Waits for the rsync of
/exp/sbnd/data/users/gputnam/GUMPLE/sbn-rewgted-20-calovarB ->
../sbn-rewgted-20-calovarB to finish (every remote file present locally with a
matching size, and no rsync partials left), then runs, serially, with the loaddf
cache writing into the -calovarB directory (loaddf's CACHE_WITH_INPUT default):

  1. validate_gumple_vs_recaf.py           (-14 vs -20 production validation)
  2. TrackSplittingCorrection_GUMPLE.py    (first: warms the shared-FV caches)
  3. mcdata_comparison_gumple.py -d all
  4. mcdata_constraint_gumple.py           (SBND 1p constraint, incl. dp bins)
  5. nb/MCDataComparisonPID-GUMPLE.ipynb   (headless, one run per detector
                                            x chi2 flavor -- see CHI2_FLAVORS)
  6. nb/SignalBoxSystematics-GUMPLE.ipynb  (headless; the longest job)
  7. offbeam_cosmicmc_gumple.py            (SBND off-beam data vs CORSIKA
                                            in-time cosmic MC; technote only)
  8. selection_plots_gumple.py             (CV cut breakdowns, efficiencies,
                                            near/far signal-box distributions)

--df-dir / --plotbase reach every step: the scripts by argument, the PID and
SignalBox notebooks through the GUMPLE_DF_DIR / GUMPLE_PLOTBASE environment
variables (the validate step gets --new-dir).

Every subprocess runs with cwd = analysis_village/gump (the cwd-relative
sys.path hacks in loaddf.py / kinematics.py need it) and logs to
<plotbase>/logs/. Failures do not stop the chain (the analyses are
independent); a summary table of exit codes is printed at the end.

The PID step runs once per (detector, chi2 flavor) pair. The calovarB production
stores the chi2 PID under several flavors -- which hits of the collection-plane
calorimetry the chi2 is computed from: the full hit set, or one of three
progressively harder trims -- and each flavor gets its own
<plotbase>/pid/<flavor>/ sub-directory. Flavors run innermost so a detector's
four runs share one warm loaddf cache. NB the flavors CHANGED from -calovar,
which had the three wire planes plus a best-plane pick.

    python run_gumple_serial.py [--plotbase DIR] [--skip-wait] [--only STEP ...]
                                [--chi2-flavor FLAVOR ...] [--df-dir DIR]
                                [--calo-model scan|nominal] [--no-angles DET ...]

--df-dir points the whole chain at a different production (the download gate is
sized from the remote directory of the -calovarB default, so pass --skip-wait
with it). --calo-model and --no-angles are pid-step-only knobs, passed to the
notebook as GUMPLE_CALO_MODEL / GUMPLE_ANGLE_BINS:

  --calo-model nominal   one tuned calorimetric budget per detector instead of
                         the 10/11-model scan (ICARUS 1% gain + dE/dx smearing +
                         dE/dx bias; SBND 2% gain + 2x dE/dx smearing + EMB)
  --calo-model smearscan that nominal budget with the dE/dx smearing term at 0x,
                         1x and 2x and every other calorimetric term held fixed,
                         on the angle-inclusive plots; the angular quintile plots
                         keep the single nominal budget
  --no-angles DET        drop the angular quintile stages for that detector,
                         leaving the inclusive comparison only. Repeatable.
  --cut-stages TAGS      comma-separated cut-stage filename tags to plot
                         (simplecosrej, simplecosrej_trkqual_nop,
                         simplecosrej_trkqual_p3), or "all"
"""

import argparse
import glob
import os
import subprocess
import sys
import time

GUMP_DIR = os.path.dirname(os.path.abspath(__file__))
PY = "/Users/gputnam-local/Work/fitter/env/bin/python"
JUPYTER = "/Users/gputnam-local/Work/fitter/env/bin/jupyter"

DF_DIR = "/Users/gputnam/Work/osc/sbn-rewgted-20-calovarB/"
REMOTE = "gputnam@sbndgpvm04.fnal.gov"
REMOTE_DIR = "/exp/sbnd/data/users/gputnam/GUMPLE/sbn-rewgted-20-calovarB"
DEFAULT_PLOTBASE = "/Users/gputnam/Work/osc/cafpyana/plots-gumple-2026-08-31-calovarB"

# chi2 PID flavors the PID notebook can plot, in the order they are run. Must
# match CHI2_FLAVORS in nb/MCDataComparisonPID-GUMPLE.ipynb.
#
# "nominal" is the single calculation stored by sbn-rewgted-21, which turned the
# alternates off in production; it is the angle-blended default PID (untrimmed
# collection plane above theta_x = 47 deg, p2trim2 below it, on the proton
# candidate). The other four are calovarB's: the untrimmed collection plane plus
# three progressively harder hit trims, each carrying the cafpyana calo/smear
# variation columns, so each can be run with the full systematic budget.
#
# "all" expands to the calovarB four -- "nominal" reads the same columns as "p2"
# under a different label, so running both would just duplicate the plots.
CHI2_FLAVORS = ["nominal", "p2", "p2trim", "p2trim2", "p2trim3"]
CALOVARB_FLAVORS = ["p2", "p2trim", "p2trim2", "p2trim3"]

# Fallback expected file set (name -> bytes), from the remote -calovarB listing
# on 2026-08-29 (47 files, 79.4 GB). Used only if the live ssh listing fails; the
# live listing wins so late additions on the remote side are waited for too.
EXPECTED_FILES = {
    "ICARUSRun2_SpringMCOverlay_rewgt_0.df":         3515245718,
    "ICARUSRun2_SpringMCOverlay_rewgt_1.df":         3509886290,
    "ICARUSRun2_SpringMCOverlay_rewgt_2.df":         622005505,
    "ICARUSRun2_Spring_Overlay_Dirt.df":             149640724,
    "ICARUSRun2_Spring_Overlay_SCE.df":              557217901,
    "ICARUSRun2_Spring_Overlay_WMXThXW.df":          511784229,
    "ICARUSRun2_Spring_Overlay_WMYZ.df":             546881541,
    "ICARUSRun4_SpringMCOverlay_rewgt_0.df":         4810748949,
    "ICARUSRun4_SpringMCOverlay_rewgt_1.df":         4801912033,
    "ICARUSRun4_SpringMCOverlay_rewgt_2.df":         4807228417,
    "ICARUSRun4_SpringMCOverlay_rewgt_3.df":         4412590963,
    "ICARUSRun4_Spring_Overlay_Dirt.df":             412792772,
    "ICARUSRun4_Spring_Overlay_SCE.df":              1254441567,
    "ICARUSRun4_Spring_Overlay_WMXThXW.df":          710808104,
    "ICARUSRun4_Spring_Overlay_WMYZ.df":             1251012068,
    "ICARUS_SpringRun2BNBOff_unblind.df":            199157637,
    "ICARUS_SpringRun2BNB_FullOnBeam.df":            1073208475,
    "ICARUS_SpringRun2BNB_unblind.df":               93081328,
    "ICARUS_SpringRun4BNBOff_unblind.df":            60861822,
    "ICARUS_SpringRun4BNB_unblind.df":               82419474,
    "SBNDAr25.df":                                   575370964,
    "SBNDIntimeMC.df":                               265086611,
    "SBNDMCCV_0.df":                                 2782417595,
    "SBNDMCCV_1.df":                                 2791788698,
    "SBNDMCCV_10.df":                                2785102364,
    "SBNDMCCV_11.df":                                2771211591,
    "SBNDMCCV_12.df":                                858746141,
    "SBNDMCCV_2.df":                                 2773203488,
    "SBNDMCCV_3.df":                                 2791773552,
    "SBNDMCCV_4.df":                                 2793369040,
    "SBNDMCCV_5.df":                                 2793792005,
    "SBNDMCCV_6.df":                                 2790264184,
    "SBNDMCCV_7.df":                                 2790388256,
    "SBNDMCCV_8.df":                                 2768104225,
    "SBNDMCCV_9.df":                                 2779260668,
    "SBND_SpringBNBData_FixedDev.df":                26572445,
    "SBND_SpringBNBData_FullOnBeam.df":              524959262,
    "SBND_SpringBNBData_RollingDev.df":              9035985,
    "SBND_SpringBNBOffData.df":                      199030544,
    "SBND_SpringLowEMC.df":                          494223966,
    "SBND_SpringMC_0xSCE.df":                        561134789,
    "SBND_SpringMC_2xSCE.df":                        584703115,
    "SBND_SpringMC_DENT.df":                         573177968,
    "SBND_SpringMC_Nom.df":                          610104664,
    "SBND_SpringMC_WMNom.df":                        2299001350,
    "SBND_SpringMC_WMXThetaXW.df":                   2541929540,
    "SBND_SpringMC_WMYZ.df":                         2492150024,
}


def remote_listing():
    """{name: size} from the remote directory, or None if ssh fails."""
    try:
        out = subprocess.run(
            ["ssh", "-K", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15", REMOTE,
             "ls -l %s" % REMOTE_DIR],
            capture_output=True, text=True, timeout=60)
    except Exception:
        return None
    if out.returncode != 0:
        return None
    files = {}
    for line in out.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 9 and parts[-1].endswith(".df"):
            files[parts[-1]] = int(parts[4])
    return files or None


def rsync_running():
    r = subprocess.run(["pgrep", "-f", "rsync.*sbn-rewgted-20-calovar"],
                       capture_output=True, text=True)
    return r.returncode == 0


def download_status(expected, df_dir):
    """(missing, mismatched, partials) for the local dir vs `expected`."""
    missing, mismatched = [], []
    for name, size in expected.items():
        p = os.path.join(df_dir, name)
        if not os.path.exists(p):
            missing.append(name)
        elif os.path.getsize(p) != size:
            mismatched.append(name)
    partials = [os.path.basename(p) for p in glob.glob(os.path.join(df_dir, ".*.df.*"))]
    return missing, mismatched, partials


def wait_for_download(df_dir, poll=60):
    expected = remote_listing()
    src = "remote listing"
    if expected is None:
        expected = EXPECTED_FILES
        src = "embedded fallback list"
    print("[gate] expecting %d files (%.1f GB) from the %s"
          % (len(expected), sum(expected.values()) / 1e9, src), flush=True)

    while True:
        missing, mismatched, partials = download_status(expected, df_dir)
        if not missing and not mismatched and not partials:
            print("[gate] download complete: all %d files present with matching sizes"
                  % len(expected), flush=True)
            return
        n_done = len(expected) - len(missing) - len(mismatched)
        print("[gate] %d/%d done, %d incomplete/missing, partials: %s"
              % (n_done, len(expected), len(missing) + len(mismatched),
                 ", ".join(partials) or "none"), flush=True)
        if not rsync_running():
            raise RuntimeError(
                "rsync is no longer running but the download is incomplete "
                "(missing: %s; size-mismatched: %s). Restart the rsync and re-run."
                % (missing[:5], mismatched[:5]))
        time.sleep(poll)


def run_step(name, cmd, logdir, env_extra=None):
    log = os.path.join(logdir, "%s.log" % name)
    env = dict(os.environ)
    if env_extra:
        env.update(env_extra)
    print("\n[run] %s\n      %s\n      log: %s" % (name, " ".join(cmd), log), flush=True)
    t0 = time.time()
    with open(log, "w") as lf:
        r = subprocess.run(cmd, cwd=GUMP_DIR, env=env, stdout=lf,
                           stderr=subprocess.STDOUT)
    dt = time.time() - t0
    status = "OK" if r.returncode == 0 else "FAILED (exit %d)" % r.returncode
    print("[run] %s: %s in %.1f min" % (name, status, dt / 60.), flush=True)
    return r.returncode


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plotbase", default=DEFAULT_PLOTBASE)
    ap.add_argument("--skip-wait", action="store_true",
                    help="Do not wait for the rsync (files already complete)")
    ap.add_argument("--poll", type=int, default=60, help="Download poll interval [s]")
    ap.add_argument("--only", action="append", default=None,
                    help="Run only these steps (validate, tracksplit, mcdata, "
                         "constraint, pid, signalbox, corsika, selection). "
                         "Repeatable.")
    ap.add_argument("--chi2-flavor", action="append", default=None,
                    help="chi2 PID flavor(s) for the pid step: %s, or 'all' "
                         "(= the calovarB four). Repeatable. Default: all."
                         % ", ".join(CHI2_FLAVORS))
    ap.add_argument("--df-dir", default=DF_DIR,
                    help="Input production directory (default: %(default)s). "
                         "Pass --skip-wait with a non-default one -- the "
                         "download gate is sized from the -calovarB remote.")
    ap.add_argument("--calo-model", default="scan",
                    choices=["scan", "nominal", "smearscan"],
                    help="pid step: 'scan' (default) runs every smearing model, "
                         "'nominal' runs the single tuned per-detector "
                         "calorimetric budget, 'smearscan' runs that budget at "
                         "0x/1x/2x dE/dx smearing on the angle-inclusive plots")
    ap.add_argument("--no-angles", action="append", default=None, metavar="DET",
                    help="pid step: drop the angular quintile stages for this "
                         "detector (e.g. --no-angles 'ICARUS Run4'). Repeatable.")
    ap.add_argument("--cut-stages", default="all", metavar="TAGS",
                    help="pid step: comma-separated cut-stage filename tags to "
                         "plot -- simplecosrej, simplecosrej_trkqual_nop, "
                         "simplecosrej_trkqual_p3 -- or 'all' (default)")
    args = ap.parse_args(argv)

    if not args.chi2_flavor:
        args.chi2_flavor = ["all"]
    flavors = []
    for f in args.chi2_flavor:
        if f == "all":
            flavors += [x for x in CALOVARB_FLAVORS if x not in flavors]
        elif f in CHI2_FLAVORS:
            if f not in flavors:
                flavors.append(f)
        else:
            ap.error("unknown chi2 flavor %r (choose from %s, or 'all')"
                     % (f, CHI2_FLAVORS))
    args.chi2_flavor = flavors

    df_dir = args.df_dir
    if not df_dir.endswith("/"):
        df_dir += "/"
    if not os.path.isdir(df_dir):
        ap.error("--df-dir %s does not exist" % df_dir)

    detectors = ["ICARUS Run2", "ICARUS Run4", "SBND"]
    no_angles = args.no_angles or []
    for d in no_angles:
        if d not in detectors:
            ap.error("unknown --no-angles detector %r (choose from %s)"
                     % (d, detectors))

    plotbase = args.plotbase
    logdir = os.path.join(plotbase, "logs")
    os.makedirs(logdir, exist_ok=True)

    print("[cfg] df-dir      %s" % df_dir)
    print("[cfg] plotbase    %s" % plotbase)
    print("[cfg] calo-model  %s" % args.calo_model)
    print("[cfg] flavors     %s" % ", ".join(args.chi2_flavor))
    print("[cfg] no-angles   %s" % (", ".join(no_angles) or "none"))
    print("[cfg] cut-stages  %s" % args.cut_stages, flush=True)

    if not args.skip_wait:
        wait_for_download(df_dir, poll=args.poll)

    steps = []

    steps.append(("validate", [
        PY, "validate_gumple_vs_recaf.py",
        "--new-dir", df_dir,
        "--out", os.path.join(logdir, "validate.txt")], None))

    steps.append(("tracksplit", [
        PY, "TrackSplittingCorrection_GUMPLE.py",
        "--df-dir", df_dir,
        "--plotdir", os.path.join(plotbase, "tracksplit") + "/"], None))

    steps.append(("mcdata", [
        PY, "mcdata_comparison_gumple.py", "-d", "all",
        "--df-dir", df_dir,
        "--plotdir", os.path.join(plotbase, "evtsel") + "/",
        "--nproc", "8"], None))

    steps.append(("constraint", [
        PY, "mcdata_constraint_gumple.py",
        "--df-dir", df_dir,
        "--plotdir", os.path.join(plotbase, "constraint") + "/"], None))

    # One PID run per (detector, chi2 flavor). Flavors innermost: the notebook's
    # CV MC load is cached by loaddf into the .df directory with a
    # flavor-independent key, so only a detector's first flavor pays the full load.
    for det in detectors:
        tag = det.replace(" ", "-")
        for flavor in args.chi2_flavor:
            steps.append(("pid-%s-%s" % (tag, flavor), [
                JUPYTER, "nbconvert", "--to", "notebook", "--execute",
                "--ExecutePreprocessor.timeout=-1",
                "--output", "MCDataComparisonPID-GUMPLE-%s-%s.ipynb" % (tag, flavor),
                "nb/MCDataComparisonPID-GUMPLE.ipynb"],
                {"GUMPLE_DETECTOR": det, "GUMPLE_CHI2_FLAVOR": flavor,
                 "GUMPLE_DF_DIR": df_dir, "GUMPLE_PLOTBASE": plotbase,
                 "GUMPLE_CALO_MODEL": args.calo_model,
                 "GUMPLE_CUT_STAGES": args.cut_stages,
                 "GUMPLE_ANGLE_BINS": "0" if det in no_angles else "1"}))

    steps.append(("signalbox", [
        JUPYTER, "nbconvert", "--to", "notebook", "--execute", "--inplace",
        "--ExecutePreprocessor.timeout=-1",
        "nb/SignalBoxSystematics-GUMPLE.ipynb"],
        {"GUMPLE_DF_DIR": df_dir, "GUMPLE_PLOTBASE": plotbase}))

    # Technote-only steps (2026-09-02): the CORSIKA / off-beam cosmic MC
    # comparison and the CV-only selection plots (cut breakdowns, efficiencies,
    # near/far signal-box distributions). Both are cheap next to the above.
    steps.append(("corsika", [
        PY, "offbeam_cosmicmc_gumple.py",
        "--df-dir", df_dir,
        "--plotdir", os.path.join(plotbase, "corsika") + "/"], None))

    steps.append(("selection", [
        PY, "selection_plots_gumple.py",
        "--df-dir", df_dir,
        "--plotdir", os.path.join(plotbase, "selection") + "/"], None))

    if args.only:
        keep = set(args.only)
        steps = [s for s in steps if any(s[0] == k or s[0].startswith(k + "-")
                                         for k in keep)]

    results = {}
    for name, cmd, env_extra in steps:
        results[name] = run_step(name, cmd, logdir, env_extra)

    print("\n" + "=" * 60)
    print("SERIAL RUN SUMMARY  (plots under %s)" % plotbase)
    for name, rc in results.items():
        print("  %-24s %s" % (name, "OK" if rc == 0 else "FAILED (exit %d)" % rc))
    print("=" * 60, flush=True)
    return 0 if all(rc == 0 for rc in results.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
