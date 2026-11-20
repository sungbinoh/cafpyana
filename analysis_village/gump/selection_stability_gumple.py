"""Stability of the GUMPLE selection rate: selected events / POT vs run number and date.

For each run set the on-beam data (after the data-only beam-quality and good-run
cuts of dataquality.py, as applied by loaddf.load_one) is split into groups of
consecutive runs, each holding enough POT that its expected number of selected
events reaches min(max(1/rel_err^2, N_total/max_groups), N_total/min_groups) -- at most rel_err
statistical error per point, and no more than max_groups points per plot (N_total summed over the
plot's run sets), unless that would leave fewer than min_groups points (small samples, e.g. ICARUS
Run 4 Mp), where the groups are made smaller instead. Each group's selected events / POT is plotted
divided by the run set's average, so the absolute rate stays blinded.

Off-beam subtraction (default; --no-offbeam-sub turns it off): the off-beam data
(same runs as the on-beam data, separate stream) are selected identically (DQ cut,
no BQ cut) and subtracted per GROUP, i.e. per date range, with that group's own
gate ratio w_g = ngates_ON_g / ngates_OFF_g:
    net_g = N_on_g - w_g N_off_g,   var_g = N_on_g + w_g^2 N_off_g.
The gates follow dataquality.data_norm, per run: SBND ON = BQ spills, OFF =
noffbeambnb, with w x (1 - 0.0754); ICARUS ON = trig.gate_delta x the run's BQ spill fraction
x gate-loss correction, OFF = trig.gate_delta x (1 - 1/20). Groups are still formed
from the on-beam counts (so they are the same as without subtraction).

Plots (vs date, the date taken from the BNB spill times; <name>_stability_date.{png,pdf}),
on the full on-beam data samples:
  * SBND                 -- SBND_SpringBNBData_FullOnBeam.df
  * ICARUS-Run2          -- ICARUS_SpringRun2BNB_FullOnBeam.df, both cryostats
  * ICARUS-Run4          -- ICARUS_SpringRun4BNB_ReCAF2026.df, both cryostats
  * ICARUS-West-Run2+4   -- West cryostat (slc_vtx_x > 0) of both runs, divided by
                            the common Run 2+4 average, with a broken x-axis

Selections: PID Np, i.e. FV & cosmic rejection & two-prong & PID & flash & n_pfp >= 2
(mcdata_comparison_gumple.pid_cut_np plus the flash cut applied at selection time), and
its exclusive split into 1p (n_pfp == 2) and Mp (n_pfp > 2), as in
signalbox_systematics_gumple.py. Every plot is made for each selection; Np keeps the
<name>_stability_date names, 1p / Mp are <name>_1p_stability_date / <name>_Mp_stability_date.
Each selection is grouped by its own counts (same rule), so Mp has fewer, larger groups.

Beam quality includes the FOM
cut for the run sets in dataquality.FOM_DETECTORS (corrected FOM, dataquality.correct_fom).

NB `FV` must stay a module-level function in `__main__`: loaddf's cache key uses
the preselection's qualname, and this keeps the key shared with the other GUMPLE
tools.

Usage:
    python selection_stability_gumple.py --df-dir /path/to/sbn-rewgted-24/ --outdir OUT
"""
import argparse
import datetime
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import warnings
warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

# repo root first, the gump dir LAST (analysis_village/gump/makedf.py shadows the makedf package)
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..", "..")))
sys.path.insert(0, os.path.join(_HERE, "..", "gumple"))
sys.path.append(_HERE)

import gumple_cuts as gc
import loaddf
import dataquality as dq

plt.style.use(os.path.join(_HERE, "dune.mplstyle"))

# run set -> on-beam files: the full on-beam data samples. Do NOT add the SBND
# FixedDev/RollingDev or ICARUS Run 2 *_unblind.df samples: they overlap the
# FullOnBeam data. ICARUS Run 4 uses the full-dataset ReCAF2026 file, whose spill
# information comes from dataquality.EXTERNAL_BEAM. (A run set may list several
# files, whose per-run POT and counts are summed; they must be disjoint.)
SAMPLES = {
    "SBND": ["SBND_SpringBNBData_FullOnBeam.df"],
    "ICARUS Run2": ["ICARUS_SpringRun2BNB_FullOnBeam.df"],
    "ICARUS Run4": ["ICARUS_SpringRun4BNB_ReCAF2026.df"],
}
# run set -> off-beam files (as in mcdata_comparison_gumple.get_files). ICARUS Run 4
# pairs the ReCAF2026 on-beam file with the ReCAF2026 off-beam file (same processing).
OFFBEAM_SAMPLES = {
    "SBND": ["SBND_SpringBNBOffData.df"],
    "ICARUS Run2": ["ICARUS_SpringRun2BNBOff_unblind.df"],
    "ICARUS Run4": ["ICARUS_SpringRun4BNBOff_ReCAF2026.df"],
}
# gate bookkeeping constants of dataquality.data_norm
SBND_F_FACTOR = 0.0754
ICARUS_GATE_CORR = {"ICARUS Run2": 1. - 1/200., "ICARUS Run4": 1. - 1/40.}
ICARUS_OFF_GATE_CORR = 1. - 1/20.
TITLES = {"SBND": "SBND Run 1", "ICARUS Run2": "ICARUS Run 2", "ICARUS Run4": "ICARUS Run 4"}
COLORS = {"SBND": "#1e3f54", "ICARUS Run2": "#d54c28", "ICARUS Run4": "#315031"}
# selection label -> (column suffix, n_pfp mask); Np = 1p + Mp
SELECTIONS = {"Np": ("", lambda n: n >= 2), "1p": ("_1p", lambda n: n == 2), "Mp": ("_mp", lambda n: n > 2)}


# ============================================================
# Selection
# ============================================================
# NB: keep at module scope, named FV (same body as mcdata_comparison_gumple.FV).
def FV(df):
    return gc.sanity_cut(df) & gc.slcfv_cut(df) & df.cut_contained & df.cut_cathode


def selection_np(df):
    return FV(df) & gc.cosmic_cut(df) & gc.trk_cut(df) & gc.pid_cut(df) & gc.flash_cut(df) & (df.n_pfp >= 2)


# ============================================================
# Per-run tables
# ============================================================
def selected_events(fname, detector, offbeam=False):
    """(run, evt, west, n_pfp) of every selected slice in the on- (or off-) beam file, after the data cuts."""
    out = []
    keys = []
    for k in dq._keys(fname, "hdr"):
        idf = int(k.split("_")[-1])
        df, _, _ = loaddf.load_one(fname, idf, load_truth=False, include_syst=False, detector=detector,
                                   preselection=FV, match_Enu=False, offbeampot=offbeam)
        hdr = pd.read_hdf(fname, "hdr_%i" % idf)
        keys.append(hdr[["run", "evt"]].assign(idf=idf))
        if len(df) == 0:  # gumple cuts crash on zero-row frames
            continue
        df = df[selection_np(df)]
        if len(df) == 0:
            continue
        ev = pd.MultiIndex.from_arrays([df.index.get_level_values("__ntuple"),
                                        df.index.get_level_values("entry")])
        re = hdr[["run", "evt"]].reindex(ev)
        assert not re.run.isna().any(), "selected slices without a header record in %s split %i" % (fname, idf)
        out.append(pd.DataFrame({"run": re.run.astype(int).to_numpy(), "evt": re.evt.astype(int).to_numpy(),
                                 "west": (df.slc_vtx_x > 0).to_numpy(), "n_pfp": df.n_pfp.astype(int).to_numpy()}))
    sel = pd.concat(out, ignore_index=True) if out else pd.DataFrame({"run": [], "evt": [], "west": [], "n_pfp": []})

    # cross-split dedup, as in loaddf.load: drop every event whose (run, evt) shows up in >1 split
    keys = pd.concat(keys)
    dup = keys[keys.duplicated(["run", "evt"], keep=False)]
    if len(dup):
        bad = pd.MultiIndex.from_frame(dup[["run", "evt"]]).unique()
        n0 = len(sel)
        sel = sel[~pd.MultiIndex.from_frame(sel[["run", "evt"]]).isin(bad)]
        print("[%s] cross-split dedup: %d duplicated (run, evt) keys, dropped %d selected slices"
              % (os.path.basename(fname), len(bad), n0 - len(sel)))
    return sel


def _gate_delta_per_run(fname):
    """sum trig.gate_delta per run (ICARUS gate counting), over all splits, without
    the duplicated (run, evt) events that loaddf drops (dq.dup_events)."""
    hdrs = {int(k.split("_")[-1]): pd.read_hdf(fname, k) for k in dq._keys(fname, "hdr")}
    dup = dq.dup_events(pd.concat(hdrs.values()))
    rows = []
    for idf, h in hdrs.items():
        h["_dup"], dup = dup[:len(h)], dup[len(h):]
        t = dq._with_run(pd.read_hdf(fname, "trig_%i" % idf), h)
        assert not t.run.isna().any(), "trig rows without a header record in %s split %i" % (fname, idf)
        t = t[~t._dup.astype(bool)]
        rows.append(pd.DataFrame({"run": t.run.astype(int).to_numpy(), "gd": t.gate_delta.astype(float).to_numpy()}))
    return pd.concat(rows).groupby("run").gd.sum()


def per_run_exposure(fname, detector):
    """Per good run: POT (TOR875 of beam-quality spills), the spill time span / POT-weighted mean time,
    and the on-beam gates of dataquality.data_norm (SBND: BQ spills; ICARUS: gate_delta x BQ spill fraction
    x gate-loss correction, the BQ fraction taken per run)."""
    rows = []
    for k in dq._keys(fname, "bnb"):
        b = dq.spill_table(fname, int(k.split("_")[-1]), detector)
        b = b[dq.data_quality_cut(b.run, detector)]
        bq = b.bq.to_numpy()
        t = b.spill_time_sec.astype(float) + 1e-9*b.spill_time_nsec.astype(float)
        tor = b.TOR875.astype(float).where(b.bq, 0.)
        rows.append(pd.DataFrame({"run": b.run.to_numpy(), "pot": tor.to_numpy(), "tpot": (t*tor).to_numpy(),
                                  "t": np.where(bq, t.to_numpy(), np.nan), "nspill": 1, "nspill_bq": bq.astype(int)}))
    r = pd.concat(rows).groupby("run").agg(pot=("pot", "sum"), tpot=("tpot", "sum"),
                                           t_start=("t", "min"), t_end=("t", "max"),
                                           nspill=("nspill", "sum"), nspill_bq=("nspill_bq", "sum"))
    r = r[r.nspill_bq > 0]  # runs with no beam-quality spill carry no POT
    if "ICARUS" in detector:
        gd = _gate_delta_per_run(fname).reindex(r.index).fillna(0.)
        r["ngates_on"] = gd * r.nspill_bq / r.nspill * ICARUS_GATE_CORR[detector]
    else:
        r["ngates_on"] = r.nspill_bq.astype(float)
    return r.drop(columns=["nspill", "nspill_bq"])


def per_run_offbeam_gates(fname, detector):
    """Off-beam gates per good run (dataquality.data_norm conventions; for SBND the
    (1 - f_factor) of the off-beam weight is folded in here)."""
    if "ICARUS" in detector:
        g = _gate_delta_per_run(fname) * ICARUS_OFF_GATE_CORR
    else:
        h = dq._read_all(fname, "hdr")
        h = h[~dq.dup_events(h)]
        g = h.groupby("run").noffbeambnb.sum().astype(float) / (1. - SBND_F_FACTOR)
    return g[dq.data_quality_cut(g.index, detector)].rename("ngates_off")


def combine_exposure(exposures):
    """Sum per-run exposure tables over files (a run can be split across files)."""
    r = pd.concat(exposures).groupby(level="run").agg(pot=("pot", "sum"), tpot=("tpot", "sum"),
                                                      t_start=("t_start", "min"), t_end=("t_end", "max"),
                                                      ngates_on=("ngates_on", "sum"))
    r["t_mean"] = r.tpot / r.pot
    return r.drop(columns="tpot")


def run_table(exposure, sel, sel_off=None, ngates_off=None):
    """Per-run table: exposure, on-beam counts n_sel / n_west and (if given) off-beam
    counts n_off / n_off_west and gates ngates_off, for each selection in SELECTIONS
    (column suffix "" for Np, "_1p", "_mp"). Off-beam runs with no on-beam POT
    are dropped (they cannot be paired with any on-beam data)."""
    tab = exposure.copy()
    cnt = lambda s: s.groupby("run").size().reindex(tab.index).fillna(0).astype(int)
    extra = set(sel.run) - set(tab.index)
    assert not extra, "selected events in runs with no counted POT: %s" % sorted(extra)
    for suf, mask in SELECTIONS.values():
        s = sel[mask(sel.n_pfp)]
        tab["n_sel" + suf] = cnt(s)
        tab["n_west" + suf] = cnt(s[s.west.astype(bool)])
        if sel_off is None:
            tab["n_off" + suf] = tab["n_off_west" + suf] = 0
        else:
            so = sel_off[mask(sel_off.n_pfp)]
            tab["n_off" + suf] = cnt(so)
            tab["n_off_west" + suf] = cnt(so[so.west.astype(bool)])
    for c in ("n_sel", "n_west", "n_off", "n_off_west"):
        assert (tab[c] == tab[c + "_1p"] + tab[c + "_mp"]).all(), "1p + Mp != Np in " + c
    if sel_off is None:
        tab["ngates_off"] = np.nan
    else:
        tab["ngates_off"] = ngates_off.reindex(tab.index)
        nogate = tab.index[~(tab.ngates_off > 0)]
        assert not len(nogate), "on-beam runs with no off-beam gates: %s" % list(nogate)
        orphan = sorted(set(ngates_off.index) - set(tab.index))
        if orphan:
            print("  off-beam runs without on-beam POT (dropped): %s (%d selected off-beam slices)"
                  % (orphan, int(sel_off.run.isin(orphan).sum())))
    return tab.sort_index()


# ============================================================
# Grouping
# ============================================================
def group_runs(tab, rate, nmin):
    """Consecutive runs into groups whose EXPECTED count (POT x rate) reaches nmin.

    A tail below threshold is merged into the last group. Returns one row per group.
    """
    gid = np.zeros(len(tab), dtype=int)
    g, acc = 0, 0.
    for i, pot in enumerate(tab.pot.to_numpy()):
        gid[i] = g
        acc += pot*rate
        if acc >= nmin:
            g, acc = g + 1, 0.
    if acc > 0 and acc < nmin and g > 0:
        gid[gid == g] = g - 1  # merge the short tail
    t = tab.assign(group=gid, run=tab.index, rpot=tab.index*tab.pot, tpot=tab.t_mean*tab.pot)
    grp = t.groupby("group").agg(run_min=("run", "min"), run_max=("run", "max"), nruns=("run", "size"),
                                 t_start=("t_start", "min"), t_end=("t_end", "max"),
                                 pot=("pot", "sum"), rpot=("rpot", "sum"), tpot=("tpot", "sum"),
                                 n_sel=("n_sel", "sum"), n_off=("n_off", "sum"),
                                 ngates_on=("ngates_on", "sum"), ngates_off=("ngates_off", "sum"))
    grp["run_mean"] = grp.rpot / grp.pot
    grp["t_mean"] = grp.tpot / grp.pot
    return grp.drop(columns=["rpot", "tpot"])


def offbeam_subtract(groups, subtract):
    """Per-group off-beam weight w, net count and its variance (w = 0 without subtraction)."""
    groups["w_off"] = groups.ngates_on / groups.ngates_off if subtract else 0.
    groups["n_net"] = groups.n_sel - groups.w_off*groups.n_off
    groups["var"] = groups.n_sel + groups.w_off**2*groups.n_off
    groups["offbeam_frac"] = groups.w_off*groups.n_off / groups.n_sel
    return groups


def normalize(groups, rate):
    """Add the blinded rate y = (net/POT)/rate, its statistical error, and the expected count."""
    groups["n_exp"] = groups.pot * rate
    groups["y"] = groups.n_net / groups.n_exp
    groups["yerr"] = np.sqrt(groups["var"]) / groups.n_exp
    return groups


def flatness(segments, subtract):
    """chi2 of the (net) counts against the common average rate (1 fitted normalization).

    Without subtraction the Pearson variance (expected count) is used, as before; with
    it the observed variance N_on + w^2 N_off (the off-beam term has no expectation here)."""
    g = pd.concat([s["groups"] for s in segments])
    chi2 = float((((g.n_net - g.n_exp)**2) / (g["var"] if subtract else g.n_exp)).sum())
    ndof = len(g) - 1
    p = float(stats.chi2.sf(chi2, ndof)) if ndof > 0 else np.nan
    return chi2, ndof, p


# ============================================================
# Plotting
# ============================================================
def _to_date(t):
    return mdates.date2num([datetime.datetime.fromtimestamp(x, datetime.timezone.utc) for x in np.atleast_1d(t)])


def _xvals(groups):
    x = _to_date(groups.t_mean.to_numpy())
    lo, hi = _to_date(groups.t_start.to_numpy()), _to_date(groups.t_end.to_numpy())
    return x, lo, hi


def _break_marks(ax_l, ax_r):
    d = 0.015
    kw = dict(color="k", clip_on=False, lw=1)
    for ax, xs in ((ax_l, 1), (ax_r, 0)):
        for y in (0, 1):
            ax.plot([xs - d, xs + d], [y - 2*d, y + 2*d], transform=ax.transAxes, **kw)


def _sub_tag(subtract):
    return "off-beam subtracted per group" if subtract else "on-beam, no off-beam subtraction"


def plot_stability(segments, title, outbase, subtract, selname="Np"):
    """Rate vs date, one panel per segment (x-axis broken between them), sharing the y-axis."""
    spans = []
    for s in segments:
        x, lo, hi = _xvals(s["groups"])
        spans.append((lo.min(), hi.max()))
    widths = np.array([hi - lo for lo, hi in spans], dtype=float)
    widths = np.maximum(widths, 0.25*widths.sum())

    fig, axes = plt.subplots(1, len(segments), sharey=True, figsize=(9, 5),
                             gridspec_kw=dict(width_ratios=widths, wspace=0.06), squeeze=False)
    axes = axes[0]
    for ax, s, (lo_s, hi_s), w in zip(axes, segments, spans, widths):
        g = s["groups"]
        x, lo, hi = _xvals(g)
        ax.errorbar(x, g.y, yerr=g.yerr, xerr=[x - lo, hi - x], fmt="o", ms=4, capsize=0,
                    color=s["color"], label=s["label"])
        pad = 0.04*w
        ax.set_xlim(lo_s - pad, hi_s + pad)
        ax.axhline(1., color="gray", ls="--", lw=1)
        loc = mdates.AutoDateLocator(maxticks=max(3, int(6*w/widths.sum())))
        ax.xaxis.set_major_locator(loc)
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))
        ax.xaxis.get_offset_text().set_fontsize(12)
        ax.tick_params(axis="x", labelsize=12)
        ax.tick_params(axis="y", labelsize=14)
    for ax_l, ax_r in zip(axes[:-1], axes[1:]):
        ax_l.spines["right"].set_visible(False)
        ax_r.spines["left"].set_visible(False)
        ax_r.tick_params(axis="y", which="both", left=False, labelleft=False)
        _break_marks(ax_l, ax_r)

    axes[0].set_ylabel("Selected events / POT\n(normalized to average)", fontsize=15)
    if len(axes) > 1:
        fig.supxlabel("Date (UTC)", fontsize=15, y=0.01)
    else:
        axes[0].set_xlabel("Date (UTC)", fontsize=15)

    chi2, ndof, p = flatness(segments, subtract)
    txt = "$\\chi^2$/ndf = %.1f/%d\np = %.2g" % (chi2, ndof, p)
    axes[0].text(0.03, 0.97, txt, transform=axes[0].transAxes, ha="left", va="top", fontsize=10)
    allg = pd.concat([s["groups"] for s in segments])
    ylo = min(0.8, (allg.y - allg.yerr).min() - 0.05)
    yhi = max(1.2, (allg.y + allg.yerr).max() + 0.12)
    axes[0].set_ylim(ylo, yhi)
    if len(segments) > 1:
        handles = [ax.get_legend_handles_labels() for ax in axes]
        axes[-1].legend([h for hs, _ in handles for h in hs], [l for _, ls in handles for l in ls],
                        loc="upper right", fontsize=10)
    fig.suptitle(title + ": GUMPLE PID %s selection" % selname, fontsize=12)
    fig.subplots_adjust(left=0.14, right=0.97, bottom=0.13, top=0.91)
    for ext in ("png", "pdf"):
        os.makedirs(os.path.join(os.path.dirname(outbase), ext), exist_ok=True)
        fig.savefig(os.path.join(os.path.dirname(outbase), ext, "%s.%s" % (os.path.basename(outbase), ext)))
    plt.close(fig)
    return chi2, ndof, p


# ============================================================
# Main
# ============================================================
def make_segments(tabs, nmin, max_groups, min_groups, col, offcol, subtract):
    """Group each run set's table (on-beam counts from column `col`, off-beam from `offcol`).

    Groups are formed from the on-beam rate; all groups are then divided by one common
    (net) rate."""
    rate_on = sum(t[col].sum() for t in tabs.values()) / sum(t.pot.sum() for t in tabs.values())
    # one threshold for the whole plot, so a plot has at most max_groups groups in total
    # ...and, for a small sample (e.g. ICARUS Run 4 Mp), at least min_groups groups at the cost of rel_err
    ntot = sum(t[col].sum() for t in tabs.values())
    nthr = min(max(nmin, ntot / max_groups), ntot / min_groups)
    segs = []
    for det, t in tabs.items():
        t = t.assign(n_sel=t[col], n_off=t[offcol])
        g = offbeam_subtract(group_runs(t, rate_on, nthr), subtract)
        segs.append(dict(det=det, label=TITLES[det], color=COLORS[det], groups=g))
    rate = sum(s["groups"].n_net.sum() for s in segs) / sum(s["groups"].pot.sum() for s in segs)
    for s in segs:
        normalize(s["groups"], rate)
    return segs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--df-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--rel-err", type=float, default=0.10,
                    help="max statistical error per group (on-beam); groups need >= 1/rel_err^2 expected events")
    ap.add_argument("--max-groups", type=int, default=15,
                    help="groups are made larger if a plot would otherwise have more than this many")
    ap.add_argument("--min-groups", type=int, default=5,
                    help="groups are made smaller (above rel_err) if a plot would otherwise have fewer than this many")
    ap.add_argument("--no-offbeam-sub", action="store_true",
                    help="do not subtract the off-beam data (the previous, on-beam-only behaviour)")
    args = ap.parse_args()
    subtract = not args.no_offbeam_sub
    os.makedirs(args.outdir, exist_ok=True)
    nmin = 1. / args.rel_err**2

    tabs, gate_check = {}, []
    for det, flist in SAMPLES.items():
        fnames = [os.path.join(args.df_dir, f) for f in flist]
        sel = pd.concat([selected_events(f, det) for f in fnames], ignore_index=True)
        # the files of a run set must be disjoint, or their POT would be double counted
        if len(fnames) > 1:
            keys = [pd.MultiIndex.from_frame(pd.concat([pd.read_hdf(f, k)[["run", "evt"]]
                                                        for k in dq._keys(f, "hdr")])) for f in fnames]
            for i in range(len(keys)):
                for j in range(i):
                    n = len(keys[i].intersection(keys[j]))
                    assert n == 0, "%s and %s share %d events" % (flist[i], flist[j], n)
        exp = combine_exposure([per_run_exposure(f, det) for f in fnames])
        if subtract:
            offs = [os.path.join(args.df_dir, f) for f in OFFBEAM_SAMPLES[det]]
            sel_off = pd.concat([selected_events(f, det, offbeam=True) for f in offs], ignore_index=True)
            ngates_off = pd.concat([per_run_offbeam_gates(f, det) for f in offs]).groupby(level=0).sum()
            tab = run_table(exp, sel, sel_off, ngates_off)
            # cross-check of the per-run gates against the global normalization
            if len(fnames) == 1:
                norm = dq.data_norm(det, fnames[0], offs, log=lambda m: None)
                w_glob = norm["OFF_w"]
                w_run = tab.ngates_on.sum() / tab.ngates_off.sum()
                gate_check.append("%-12s ngates_ON  per-run sum %.6e  data_norm %.6e  (ratio %.5f)\n"
                                  "%-12s ngates_OFF per-run sum %.6e  data_norm %.6e  (ratio %.5f)  [SBND: /(1-f)]\n"
                                  "%-12s OFF_w      per-run sums %.6g  data_norm %.6g  (ratio %.5f)"
                                  % (det, tab.ngates_on.sum(), norm["ngates_ON"], tab.ngates_on.sum()/norm["ngates_ON"],
                                     "", tab.ngates_off.sum(), norm["ngates_OFF"], tab.ngates_off.sum()/norm["ngates_OFF"],
                                     "", w_run, w_glob, w_run/w_glob))
                print(gate_check[-1])
        else:
            tab = run_table(exp, sel)
        tabs[det] = tab
        print("%s: %d good runs with POT, POT = %.4e" % (det, len(tab), tab.pot.sum()))

    seg = lambda dets, col, offcol: make_segments({d: tabs[d] for d in dets}, nmin, args.max_groups, args.min_groups, col, offcol, subtract)
    plots = []
    for selname, (suf, _) in SELECTIONS.items():
        tag = "" if selname == "Np" else "_" + selname
        plots += [(selname, "SBND" + tag, "SBND Run 1", seg(["SBND"], "n_sel" + suf, "n_off" + suf)),
                  (selname, "ICARUS-Run2" + tag, "ICARUS Run 2", seg(["ICARUS Run2"], "n_sel" + suf, "n_off" + suf)),
                  (selname, "ICARUS-Run4" + tag, "ICARUS Run 4", seg(["ICARUS Run4"], "n_sel" + suf, "n_off" + suf)),
                  (selname, "ICARUS-West-Run2+4" + tag, "ICARUS West Cryostat",
                   seg(["ICARUS Run2", "ICARUS Run4"], "n_west" + suf, "n_off_west" + suf))]

    rows, summary = [], []
    for selname, name, title, segs in plots:
        chi2, ndof, p = plot_stability(segs, title, os.path.join(args.outdir, "%s_stability_date" % name), subtract,
                                       selname)
        summary.append("%-3s %-23s groups %3d  chi2/ndf = %7.1f/%-3d  p = %.3g  min expected/group = %.0f"
                       % (selname, name, sum(len(s["groups"]) for s in segs), chi2, ndof, p,
                          min(s["groups"].n_exp.min() for s in segs)))
        if subtract:
            allg = pd.concat([s["groups"] for s in segs])
            summary[-1] += "  off-beam frac %.3f (groups %.3f-%.3f)" % (
                (allg.w_off*allg.n_off).sum() / allg.n_sel.sum(), allg.offbeam_frac.min(), allg.offbeam_frac.max())
        pot_tot = sum(s["groups"].pot.sum() for s in segs)
        for s in segs:
            g = s["groups"]
            for i, r in g.iterrows():
                rows.append(dict(selection=selname, plot=name, run_set=s["det"], group=i, run_min=r.run_min, run_max=r.run_max,
                                 nruns=r.nruns,
                                 t_start=datetime.datetime.fromtimestamp(r.t_start, datetime.timezone.utc).isoformat(),
                                 t_end=datetime.datetime.fromtimestamp(r.t_end, datetime.timezone.utc).isoformat(),
                                 pot_frac=r.pot/pot_tot, rel_rate=r.y, rel_rate_err=r.yerr,
                                 **(dict(offbeam_frac=r.offbeam_frac) if subtract else {})))
    # NB: no absolute counts or rates are written (blinded) -- only POT fractions, relative rates and
    # off-beam fractions.
    pd.DataFrame(rows).to_csv(os.path.join(args.outdir, "stability_groups.csv"), index=False)
    with open(os.path.join(args.outdir, "stability_summary.txt"), "w") as out:
        out.write("GUMPLE PID Np / 1p / Mp selected events / POT, normalized to the average (%s)\n" % _sub_tag(subtract))
        out.write("Data cuts: dataquality.py beam quality (FOM cut %s) + good-run lists data/run_lists.md\n" % (
                      "applied to " + ", ".join(sorted(dq.FOM_DETECTORS)) if dq.FOM_DETECTORS else "NOT applied"))
        out.write("Samples (full on-beam data): %s\n"
                  % "; ".join("%s=%s" % (d, " + ".join(fl)) for d, fl in SAMPLES.items()))
        if subtract:
            out.write("Off-beam samples: %s\n"
                      % "; ".join("%s=%s" % (d, " + ".join(fl)) for d, fl in OFFBEAM_SAMPLES.items()))
            out.write("Off-beam subtracted per group with the group's gate ratio; "
                      "chi2 uses the observed variance N_on + w^2 N_off.\n")
        out.write("Grouping: consecutive runs, >= min(max(%.0f, N_plot/%d), N_plot/%d) expected on-beam events per group "
                  "(<= %.1f%% stat. err. on-beam unless the plot would have < %d groups)\n"
                  "chi2 is against a flat line at the (common) average, 1 fitted normalization.\n\n"
                  % (nmin, args.max_groups, args.min_groups, 100*args.rel_err, args.min_groups))
        out.write("\n".join(summary) + "\n")
        if gate_check:
            out.write("\nGate bookkeeping (per-run sums vs dataquality.data_norm):\n" + "\n".join(gate_check) + "\n")
    print(open(os.path.join(args.outdir, "stability_summary.txt")).read())


if __name__ == "__main__":
    main()
