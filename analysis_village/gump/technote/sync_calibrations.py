#!/usr/bin/env python3
"""Compare the constants that the technote's calibration steps measure with the
values hard-coded in analysis_village/gump/loaddf.py, and optionally write the
measured values into loaddf.py.

Stdlib only. Calibration step -> output under --plotbase -> loaddf constant:

  flashfit    flashpe_fit/trigger_scale_summary.tsv   s_alltrig  -> MC_PE_SCALE
                                                       trig_unc   -> MC_PE_SCALE_UNC
  tracksplit  tracksplit/tracksplit_summary.csv       f_split of the Simple-stage
                                                       Run2+Run4 rows -> SPLIT_FRAC

Every MC figure depends on MC_PE_SCALE (the MC flash PE is scaled at load), and
the trigger / track-splitting systematics use the other two. Run the downstream
steps only once these agree.

A value is in sync when it equals the loaddf literal at the rounding loaddf
quotes: 3 decimals for the PE scale and its uncertainty, 4 for split fractions.

Exit status: 0 = in sync, 3 = at least one constant differs (after --apply: 0
once loaddf.py was rewritten), 2 = a summary file is missing.

Typical:
  sync_calibrations.py --plotbase ../../../plots-gumple-XX            # report
  sync_calibrations.py --plotbase ../../../plots-gumple-XX --apply    # rewrite loaddf.py
"""

import argparse
import ast
import csv
import datetime
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LOADDF = os.path.join(os.path.dirname(HERE), "loaddf.py")

# loaddf constant -> number of decimals it is quoted to
DECIMALS = {"MC_PE_SCALE": 3, "MC_PE_SCALE_UNC": 3, "SPLIT_FRAC": 4}


def measured_trigger(plotbase):
    path = os.path.join(plotbase, "flashpe_fit", "trigger_scale_summary.tsv")
    if not os.path.isfile(path):
        return path, None
    scale, unc = {}, {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            scale[row["detector"]] = float(row["s_alltrig"])
            unc[row["detector"]] = float(row["trig_unc"])
    return path, {"MC_PE_SCALE": scale, "MC_PE_SCALE_UNC": unc}


def measured_split(plotbase):
    path = os.path.join(plotbase, "tracksplit", "tracksplit_summary.csv")
    if not os.path.isfile(path):
        return path, None
    frac = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            if row["detector"] == "Run2+Run4" and row["stage"] == "simple":
                # summary names the planes "Z = 0"; loaddf uses "Z=0"
                frac[row["point"].replace(" = ", "=")] = float(row["f_split"])
    return path, {"SPLIT_FRAC": frac}


def loaddf_literals(text):
    """{name: (line index, dict)} for the single-line dict constants."""
    out = {}
    lines = text.split("\n")
    for name in DECIMALS:
        for i, line in enumerate(lines):
            m = re.match(r"^%s\s*=\s*(\{[^}]*\})" % name, line)
            if m:
                out[name] = (i, ast.literal_eval(m.group(1)))
                break
        else:
            raise SystemExit("no single-line `%s = {...}` in %s" % (name, LOADDF))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plotbase", required=True, help="plot directory the calibration steps wrote to")
    ap.add_argument("--loaddf", default=LOADDF, help="loaddf.py to compare / rewrite (default: %(default)s)")
    ap.add_argument("--apply", action="store_true", help="write the measured values into loaddf.py")
    args = ap.parse_args(argv)

    measured, missing = {}, []
    for fn in (measured_trigger, measured_split):
        path, vals = fn(args.plotbase)
        if vals is None:
            missing.append(path)
        else:
            measured.update(vals)
    for p in missing:
        print("MISSING calibration output: %s" % p)
    if missing:
        return 2

    text = open(args.loaddf).read()
    current = loaddf_literals(text)
    ndiff = 0
    print("%-16s %-14s %10s %10s %10s" % ("constant", "key", "loaddf", "measured", "rounded"))
    new = {}
    for name, nd in DECIMALS.items():
        cur = current[name][1]
        upd = dict(cur)
        for key in cur:
            if key not in measured[name]:
                print("%-16s %-14s %10s %10s   NOT MEASURED" % (name, key, cur[key], "-"))
                continue
            r = round(measured[name][key], nd)
            flag = "" if abs(r - cur[key]) < 0.5 * 10**-nd else "   <-- DIFFERS"
            ndiff += bool(flag)
            upd[key] = r
            print("%-16s %-14s %10.*f %10.*f %10.*f%s"
                  % (name, key, nd, cur[key], nd + 2, measured[name][key], nd, r, flag))
        new[name] = upd
    print("# %d constant(s) differ from loaddf.py" % ndiff)

    if ndiff and args.apply:
        lines = text.split("\n")
        stamp = "# sync_calibrations.py %s from %s" % (
            datetime.date.today().isoformat(), os.path.basename(os.path.normpath(args.plotbase)))
        for name, (i, _) in current.items():
            if new[name] != current[name][1]:
                body = ", ".join('"%s": %s' % (k, v) for k, v in new[name].items())
                lines[i] = "%s = {%s}  %s" % (name, body, stamp)
        open(args.loaddf, "w").write("\n".join(lines))
        print("# rewrote %s -- re-run every downstream step" % args.loaddf)
        return 0
    return 3 if ndiff else 0


if __name__ == "__main__":
    sys.exit(main())
