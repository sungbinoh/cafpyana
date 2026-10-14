"""POT bookkeeping for the data-only beam-quality and data-quality cuts.

For each run set (SBND Run 1, ICARUS Run 2, ICARUS Run 4) and each of its
on-beam data samples:

  * a bar chart of the POT initially (sum of valid TOR875), after the beam-quality
    cut, and after the data-quality (good-run) cut, applied cumulatively in that
    order -- <det>_<sample>_pot.{png,pdf}; the run set's headline sample (the
    full-statistics stream where one exists) is also written as <det>_pot.{png,pdf}
  * dq_tables.tex, one file holding, per run set, a LaTeX table of the runs in its
    on-beam data (DQ_TABLE_SAMPLES: SBND FullOnBeam and FixedDev in one table with a
    POT column per sample -- they overlap and are never summed -- ICARUS Run 2
    FullOnBeam, ICARUS Run 4 unblind) that the data-quality cut removes, with the
    reason for each; the ICARUS tables group adjacent same-reason runs into ranges
  * bq_tables.tex: two tables, the beam-quality criteria (with which are
    applied) and the cumulative spill / POT cutflow per sample (needs
    \\usepackage{booktabs})
  * pot_summary.csv / pot_summary.txt

The cuts are those of dataquality.py (NB the FOM cut is applied only for the
run sets in dataquality.FOM_DETECTORS -- SBND -- see the TODO there).

Usage:
    python pot_quality_summary.py --df-dir /path/to/sbn-rewgted-24/ --outdir OUT
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import dataquality as dq

plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), "dune.mplstyle"))

# run set -> [(sample label, on-beam file)], the headline sample first
SAMPLES = {
    "SBND": [("FullOnBeam", "SBND_SpringBNBData_FullOnBeam.df"),
             ("FixedDev", "SBND_SpringBNBData_FixedDev.df")],
    "ICARUS Run2": [("FullOnBeam", "ICARUS_SpringRun2BNB_FullOnBeam.df"),
                    ("unblind", "ICARUS_SpringRun2BNB_unblind.df")],
    "ICARUS Run4": [("unblind", "ICARUS_SpringRun4BNB_unblind.df")],
}

TITLES = {"SBND": "SBND Run 1", "ICARUS Run2": "ICARUS Run 2", "ICARUS Run4": "ICARUS Run 4"}

# run set -> (label, [on-beam files]) for the dq_tables "runs removed" table. With
# several files (SBND) the table has one row per removed run and a separate POT
# column and total per sample -- the samples OVERLAP, so they are never summed.
# ICARUS Run 2 shows only the full-statistics stream.
DQ_TABLE_SAMPLES = {
    "SBND": ("FullOnBeam+FixedDev", ["SBND_SpringBNBData_FullOnBeam.df", "SBND_SpringBNBData_FixedDev.df"]),
    "ICARUS Run2": ("FullOnBeam", ["ICARUS_SpringRun2BNB_FullOnBeam.df"]),
    "ICARUS Run4": ("unblind", ["ICARUS_SpringRun4BNB_unblind.df"]),
}

# run sets whose "runs removed" table groups adjacent rows with the same reason
# into a run range
GROUP_RANGES = {"ICARUS Run2", "ICARUS Run4"}


def tex_escape(s):
    return s.replace("&", "\\&").replace("%", "\\%").replace("_", "\\_").replace("#", "\\#")


def bar_chart(det, label, fname, st, outbase):
    vals = np.array([st["pot_initial"], st["pot_bq"], st["pot_dq"]])
    exp = int(np.floor(np.log10(vals[0])))
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    names = ["Initial", "After\nBeam Quality", "After\nData Quality"]
    bars = ax.bar(names, vals / 10**exp, color=["#7f7f7f", "#1f77b4", "#2ca02c"], width=0.6)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width()/2, b.get_height(), "%.3g\n(%.1f%%)" % (v / 10**exp, 100*v/vals[0]),
                ha="center", va="bottom", fontsize=11)
    ax.set_ylim(0, 1.25*vals[0] / 10**exp)
    ax.set_ylabel(r"POT [$\times 10^{%d}$]" % exp)
    ax.set_title("%s: %s" % (TITLES[det], os.path.basename(fname).replace(".df", "")), fontsize=12)
    bqtext = ("intensity + horn current + FOM" if det in dq.FOM_DETECTORS
              else "intensity + horn current\n(FOM cut not applied)")
    ax.text(0.98, 0.97, "Beam quality: " + bqtext,
            transform=ax.transAxes, ha="right", va="top", fontsize=9, color="dimgray")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig("%s.%s" % (outbase, ext))
    plt.close(fig)


def extra_good_note(det, pers, names):
    """Table note for runs added to the good-run list by hand."""
    extra = sorted(dq.EXTRA_GOOD_RUNS.get(det, {}))
    if not extra:
        return ""
    txt = []
    for r in extra:
        pots = ["%.2f$\\times 10^{17}$ POT in %s" % (p.pot_bq[r] / 1e17, tex_escape(n))
                for p, n in zip(pers, names) if r in p.index]
        txt.append("run %d%s" % (r, (" (%s, after the beam-quality cut)" % "; ".join(pots)) if pots else ""))
    return ("Note: %s %s absent from the run list and %s been added to the good runs by hand; "
            "%s not removed." % (", ".join(txt), "is" if len(extra) == 1 else "are",
                                 "has" if len(extra) == 1 else "have", "it is" if len(extra) == 1 else "they are"))


def _group_rows(bad, det):
    """[(run label, category, reason, [runs])] for the sorted removed runs `bad`, with
    adjacent same-reason runs merged into ranges for the run sets in GROUP_RANGES."""
    out = []
    for run in sorted(bad):
        cat, why = dq.cut_reason(run, det)
        if det in GROUP_RANGES and out and out[-1][1:3] == [cat, why]:
            out[-1][3].append(run)
        else:
            out.append([run, cat, why, [run]])
    return [("%d" % rs[0] if len(rs) == 1 else "%d--%d" % (rs[0], rs[-1]), cat, why, rs) for _, cat, why, rs in out]


def runs_cut_table(det, label, fnames, pers):
    """Runs in the sample(s) that the good-run cut removes (LaTeX lines).

    `pers` holds one per_run_spills frame per file in `fnames`; each sample keeps
    its own POT column and total (the samples overlap and are not summed)."""
    names = [os.path.basename(f).replace(".df", "").split("_")[-1] for f in fnames]
    bad = sorted(set().union(*[set(p.index[~p.good]) for p in pers]))
    npot = len(pers)
    colspec = "r l p{%s\\textwidth} %s" % ("0.5" if npot == 1 else "0.34", " ".join(["r"]*npot))
    if npot == 1:
        head = "Run & Category & Reason & POT [$10^{17}$] \\\\"
    else:
        head = (" & & & " + " & ".join(tex_escape(n) for n in names) + " \\\\\n"
                "Run & Category & Reason & " + " & ".join(["POT [$10^{17}$]"]*npot) + " \\\\")
    lines = [
        "% Runs in " + ", ".join(os.path.basename(f) for f in fnames) + " removed by the data-quality (good-run) cut."
        + (" The samples overlap: POT is given per sample and never summed." if npot > 1 else ""),
        "% Generated by analysis_village/gump/pot_quality_summary.py",
        "\\begin{longtable}{%s}" % colspec,
        "\\caption{%s (%s): runs removed by the data-quality cut. POT is after the beam-quality cut%s.}"
        % (TITLES[det], tex_escape(label), ", per sample (the samples overlap and are not summed)" if npot > 1 else "")
        + "\\label{tab:dq_runs_cut_%s_%s}\\\\" % (det.replace(" ", ""), label.replace("+", "")),
        "\\hline", head, "\\hline", "\\endfirsthead",
        "\\hline", head, "\\hline", "\\endhead",
    ]
    for runlab, cat, why, runs in _group_rows(bad, det):
        pots = []
        for p in pers:
            present = [r for r in runs if r in p.index and not p.good[r]]
            pots.append("%.2f" % (p.pot_bq[present].sum() / 1e17) if present else "--")
        lines.append("%s & %s & %s & %s \\\\" % (runlab, tex_escape(cat), tex_escape(why), " & ".join(pots)))
    lines.append("\\hline")
    for i, (n, p) in enumerate(zip(names, pers)):
        b_ = p[~p.good]
        tot = p.pot_bq.sum()
        cells = ["%.2f" % (b_.pot_bq.sum() / 1e17) if j == i else "" for j in range(npot)]
        lines.append("\\multicolumn{3}{l}{%s: %d of %d runs removed (%.1f\\%% of POT)} & %s \\\\"
                     % (tex_escape(n) if npot > 1 else "Total", len(b_), len(p),
                        100*b_.pot_bq.sum()/tot if tot > 0 else 0., " & ".join(cells)))
    lines.append("\\hline")
    # ICARUS tables carry no notes (removed at the user's request); SBND keeps the hand-added-run note
    for note in (extra_good_note(det, pers, names),):
        if note:
            lines.append("\\multicolumn{%d}{p{0.95\\textwidth}}{\\footnotesize %s} \\\\" % (3 + npot, note))
    lines.append("\\end{longtable}")
    return lines

def bq_cutflow(fname, det):
    """Cumulative spills / POT after each beam-quality criterion, over all splits.

    Rows: all spills, + beam intensity, + horn current, and + FOM (applied only
    for dataquality.FOM_DETECTORS; for the others it is shown for reference).
    POT is the sum of valid TOR875 (finite and > 0) over the surviving spills.
    """
    acc = {}
    for k in dq._keys(fname, "bnb"):
        b = dq.spill_table(fname, int(k.split("_")[-1]), det)
        tor = b.TOR875.astype(float).where(b.valid, 0.)
        inten = (b.TOR860 > 1e11) & (b.TOR875 > 1e11) & (b.LM875A > 1e-2) & (b.LM875B > 1e-2) & (b.LM875C > 1e-2)
        horn = (b.THCURR >= 173) & (b.THCURR <= 175)
        fom = (b.FOM > 0.98) & (b.FOM <= 1)
        for name, m in (("all", np.ones(len(b), bool)), ("intensity", inten), ("horn", inten & horn),
                        ("fom", inten & horn & fom)):
            n, pot = acc.get(name, (0, 0.))
            acc[name] = (n + int(np.sum(m)), pot + float(tor[m].sum()))
    return acc


def bq_table(cutflows):
    """LaTeX: the beam-quality criteria, and the cumulative spill / POT cutflow per sample."""
    # Two floats: the criteria (tab:beam_quality_cuts) and the cutflow
    # (tab:beam_quality_cutflow, cited from Samples.tex). The FOM footnote was
    # dropped on Overleaf in favour of the Samples.tex text, so it is not generated.
    L = ["% Beam-quality cut summary (SBN SPINE numu-disappearance technical note, Sec. 6.1).",
         "% Generated by analysis_village/gump/pot_quality_summary.py. Cuts implemented in dataquality.beam_quality_cut.",
         "",
         "\\begin{table}[htbp]", "\\centering",
         "\\caption{Beam-quality criteria, applied per spill to on-beam data only, cumulatively in the order listed. "
         "POT is the sum of TOR875 over the surviving spills (spills with a non-finite or non-positive TOR875 never count).}",
         "\\label{tab:beam_quality_cuts}",
         "\\begin{tabular}{l p{0.55\\textwidth} c}", "\\toprule",
         "Criterion & Requirement & Applied \\\\", "\\midrule",
         "Beam intensity & TOR860 $> 10^{11}$ and TOR875 $> 10^{11}$ protons; LM875A, LM875B, LM875C $> 10^{-2}$ & yes \\\\",
         "Horn current & $173 \\le$ THCURR $\\le 175$~kA & yes \\\\",
         "Figure of merit & $0.98 <$ FOM $\\le 1$ & %s only \\\\" % ", ".join(TITLES[d] for d in sorted(dq.FOM_DETECTORS)),
         "\\bottomrule", "\\end{tabular}",
         "\\end{table}", "",
         "\\begin{table}[htbp]", "\\centering",
         "\\caption{Beam-quality cutflow: spills and POT surviving each cumulative criterion, as a fraction of all "
         "stored spills and of all valid POT. The FOM requirement is applied to %s only; for the other samples "
         "its row is shown for reference.}" % ", ".join(TITLES[d] for d in sorted(dq.FOM_DETECTORS)),
         "\\label{tab:beam_quality_cutflow}",
         "\\begin{tabular}{l r r r r}", "\\toprule",
         "Selection & Spills & Spill frac. & POT & POT frac. \\\\"]
    for i, (det, sample, acc) in enumerate(cutflows):
        fomlab = "+ FOM" if det in dq.FOM_DETECTORS else "\\textit{+ FOM (not applied)}"
        names = [("all", "All"), ("intensity", "Beam intensity"), ("horn", "+ Horn current"), ("fom", fomlab)]
        n0, p0 = acc["all"]
        L += ["\\midrule", "\\multicolumn{5}{l}{\\textbf{%s}} \\\\" % tex_escape(sample)]
        for j, (k, lab) in enumerate(names):
            n, pot = acc[k]
            row = "\\quad %s & %s & %.3f & %s & %.3f \\\\" % (
                lab, "{:,}".format(n).replace(",", "\\,"), n / n0,
                "$%.3f \\times 10^{%d}$" % (pot / 10**int(np.floor(np.log10(pot))), int(np.floor(np.log10(pot)))) if pot > 0 else "0",
                pot / p0)
            L.append(row)
    L += ["\\bottomrule", "\\end{tabular}", "\\end{table}"]
    return L


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--df-dir", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    rows = []
    cutflows = []
    tex = ["% Data-quality (good-run) cut tables for SBND Run 1, ICARUS Run 2 and ICARUS Run 4.",
           "% Generated by analysis_village/gump/pot_quality_summary.py. Requires \\usepackage{longtable}.",
           "% Prefilter categories NOT applied as cuts: " +
           "; ".join("%s: %s" % (d, ", ".join(sorted(c))) for d, c in dq.IGNORED_PREFILTER.items()), "",
           # The longtables below do not float, but the beam-quality tables before them do: a
           # deferred float (e.g. the cutflow) would otherwise print AFTER the first longtable
           # while keeping its lower table number. Flush pending floats first -- the same test
           # as placeins' \FloatBarrier (not loaded in the note), so no page break is forced
           # when nothing is pending.
           "% Float barrier: flush any pending (deferred / bottom) floats before the longtables.",
           "\\par\\makeatletter",
           "\\begingroup\\let\\@elt\\relax\\edef\\@tempa{\\@botlist\\@deferlist\\@dbldeferlist}"
           "\\ifx\\@tempa\\@empty\\endgroup\\else\\endgroup\\clearpage\\fi",
           "\\makeatother", ""]
    for det, samples in SAMPLES.items():
        detname = det.replace(" ", "-")
        tlabel, tfiles = DQ_TABLE_SAMPLES[det]
        tfiles = [os.path.join(args.df_dir, f) for f in tfiles]
        tex += runs_cut_table(det, tlabel, tfiles, [dq.per_run_spills(f, det) for f in tfiles]) + ["", "\\clearpage", ""]
        for i, (label, f) in enumerate(samples):
            fname = os.path.join(args.df_dir, f)
            if not os.path.exists(fname):
                print("MISSING %s -- skipped" % fname)
                continue
            per = dq.per_run_spills(fname, det)
            st = dq.pot_stages(fname, det)
            bases = ["%s_%s" % (detname, label)] + (["%s" % detname] if i == 0 else [])
            for b in bases:
                bar_chart(det, label, fname, st, os.path.join(args.outdir, b + "_pot"))
            cutflows.append((det, "%s %s" % (TITLES[det], label), bq_cutflow(fname, det)))
            rows.append(dict(detector=det, sample=label, file=f, **st))
    with open(os.path.join(args.outdir, "dq_tables.tex"), "w") as out:
        out.write("\n".join(tex) + "\n")
    with open(os.path.join(args.outdir, "bq_tables.tex"), "w") as out:
        out.write("\n".join(bq_table(cutflows)) + "\n")

    tab = pd.DataFrame(rows)
    tab["frac_bq"] = tab.pot_bq / tab.pot_initial
    tab["frac_dq"] = tab.pot_dq / tab.pot_initial
    tab.to_csv(os.path.join(args.outdir, "pot_summary.csv"), index=False)
    with open(os.path.join(args.outdir, "pot_summary.txt"), "w") as out:
        out.write("Beam quality: TOR860,TOR875 > 1e11; LM875A/B/C > 1e-2; 173 <= THCURR <= 175 kA; "
                  "0.98 < FOM <= 1 for %s only (TODO in dataquality.py)\n" % ", ".join(sorted(dq.FOM_DETECTORS)))
        out.write("Data quality: good-run list data/run_lists.md; prefilter categories NOT applied: %s\n\n"
                  % "; ".join("%s: %s" % (d, ", ".join(sorted(c))) for d, c in dq.IGNORED_PREFILTER.items()))
        for r in rows:
            out.write("%-12s %-10s  initial %.4e  after BQ %.4e (%.1f%%)  after DQ %.4e (%.1f%%)  "
                      "[hdr.pot %.4e]  runs %d -> %d  spills %d -> %d -> %d\n"
                      % (r["detector"], r["sample"], r["pot_initial"], r["pot_bq"],
                         100*r["pot_bq"]/r["pot_initial"], r["pot_dq"], 100*r["pot_dq"]/r["pot_initial"],
                         r["pot_hdr"], r["nrun_initial"], r["nrun_dq"],
                         r["nspill_initial"], r["nspill_bq"], r["nspill_dq"]))
    print(open(os.path.join(args.outdir, "pot_summary.txt")).read())


if __name__ == "__main__":
    main()
