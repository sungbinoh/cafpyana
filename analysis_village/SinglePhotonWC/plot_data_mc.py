#!/usr/bin/env python3
"""Data / MC comparisons for the single-photon selection (cuts.py).

For each variable in VARIABLES and each selection stage in STAGES: on-beam
data (points, Poisson errors) over the prediction, stacked by truth category
(MC, scaled to the on-beam POT) plus off-beam data (scaled to the on-beam
spills), with the prediction's statistical uncertainty as a hatched band and a
data / prediction ratio panel. No under/overflow bins; the number of entries
outside each range is printed. Linear axes.

Normalisation, read from the dataframes (as analysis_village/gump/loaddf.py):
  MC        on-beam POT / MC POT. MC POT: histpotdf. On-beam POT: the sum of
            rec.hdr.pot in the hdr table (histpotdf if there is no hdr table)
  off-beam  on-beam POT / off-beam equivalent POT, the latter the number of
            off-beam BNB gates (sum of rec.hdr.noffbeambnb) x
            N_GATES_ON_PER_5E12POT x 5e12
MC events already contain cosmics in the neutrino spill ("Cosmic /
unmatched"); off-beam data models the spills without a neutrino.

--area-norm: shape comparison; the prediction (MC + off-beam) is scaled to the
on-beam data total in each plot (blinded stages: to unit area). Works without
the POT / gate information: if hdr.pot or hdr.noffbeambnb is missing, the
off-beam share of the prediction is fitted to the on-beam data at preselection
(binned Poisson likelihood over FIT_VARIABLES, see fit_offbeam_fraction), or set with
--offbeam-frac. Those fit variables are then not an independent check.

Blinding: off by default; stages listed in --blind (e.g. --blind final) are
drawn without on-beam data.

Inputs: MC .df from configs/sp_mc.py (wc + mcnu + wcevt + wgt), data .df from
configs/sp_data.py (wc + wcevt + hdr). The MicroBooNE BDT XMLs are read from
UBOONE_BDT_DIR.

  python analysis_village/SinglePhotonWC/plot_data_mc.py \\
      --mc <MC dir or .df> --onbeam <on-beam dir or .df> --offbeam <off-beam dir or .df> \\
      -o data_mc_plots/ [--area-norm]
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from analysis_village.SinglePhotonWC import cuts, dfio  # noqa: E402
from analysis_village.SinglePhotonWC.cutflow import load as load_mc  # noqa: E402
from analysis_village.SinglePhotonWC.plot_selection import GROUPS  # noqa: E402

plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), "dune.mplstyle"))

OFFBEAM_COLOR = "#a3a29d"
UBOONE_BDT_DIR = "/home/nathanielerowe/claude-devcontainer/photons/bdt/for_claude"
N_GATES_ON_PER_5E12POT = 1.05104   # SBND BNB, as in gump/loaddf.py

_full = cuts.CHAINS["presel+uboone"]
_names = list(_full)
STAGES = {   # name: (cut chain, title)
    "presel": (cuts.PRESELECTION, "preselection"),
    "bdt": ({k: _full[k] for k in _names[:_names.index("ub_other") + 1]}, "after the uB other BDT"),
    "final": (_full, "full selection"),
}

def _sp(col):
    """Single-photon-tagger column, NaN (not drawn) where the tagger did not run."""
    return lambda d: d[col].where(d.shw_sp_filled == 1)


def _mip(col):
    """make_wc_candidates track dQ/dx (T_rec_charge q / nq) in MIP units."""
    from analysis_village.SinglePhotonWC.make_wc_candidates import DQDX_MIP
    return lambda d: d[col] / DQDX_MIP


VARIABLES = [   # name, column or function of the candidate frame, bins, axis label; NaN entries are not drawn
    ("nu_x", "nu_x", np.linspace(-200, 200, 41), "reco vertex x [cm]"),
    ("nu_y", "nu_y", np.linspace(-200, 200, 41), "reco vertex y [cm]"),
    ("nu_z", "nu_z", np.linspace(0, 500, 51), "reco vertex z [cm]"),
    ("shw_energy", "shw_sp_energy", np.arange(0, 2001, 100), "counted shower reco energy [MeV]"),
    ("shw_a68", "shw1_a68", np.arange(0, 40.1, 2), "counted shower 68% opening angle [deg]"),
    ("shw_lobe_sep", "shw1_lobe_sep", np.arange(0, 60.1, 3), "counted shower two-lobe angle [deg]"),
    ("shw_costheta", "shw1_dir_z", np.linspace(-1, 1, 21), "counted shower cos(θ) to the beam"),
    ("shw_back_dist", "shw1_back_dist", np.arange(0, 501, 25), "counted shower traced back to the AV wall [cm]"),
    # dQ/dx of the counted shower's stem (SP tagger, first 20 fit points from the vertex end, MIP units)
    ("shw_stem_dqdx", _sp("shw_sp_medium_dQ_dx"), np.arange(0, 5.01, 0.2), "shower stem median dQ/dx [MIP]"),
    ("shw_stem_dqdx_0", _sp("shw_sp_vec_dQ_dx_0"), np.arange(0, 6.01, 0.25), "shower stem dQ/dx, 1st point [MIP]"),
    ("shw_stem_dqdx_1", _sp("shw_sp_vec_dQ_dx_1"), np.arange(0, 6.01, 0.25), "shower stem dQ/dx, 2nd point [MIP]"),
    ("shw_stem_dqdx_2", _sp("shw_sp_vec_dQ_dx_2"), np.arange(0, 6.01, 0.25), "shower stem dQ/dx, 3rd point [MIP]"),
    ("shw_stem_dqdx_max", _sp("shw_sp_max_dQ_dx_sample"), np.arange(0, 6.01, 0.25), "shower stem max dQ/dx [MIP]"),
    ("shw_stem_dedx", _sp("shw_sp_vec_median_dedx"), np.arange(0, 10.01, 0.4), "shower stem median dE/dx [MeV/cm]"),
    # dQ/dx of tracks (make_wc_candidates._tracks)
    ("trk1_dqdx", _mip("trk1_dqdx"), np.arange(0, 6.01, 0.2), "longest track median dQ/dx [MIP]"),
    ("mu1_dqdx", _mip("mu1_dqdx"), np.arange(0, 6.01, 0.2), "longest μ-labelled track median dQ/dx [MIP]"),
    ("mu1_dqdx_start", _mip("mu1_dqdx_start"), np.arange(0, 6.01, 0.2), "longest μ-labelled track dQ/dx, first 4 cm [MIP]"),
    ("mu1_end_ratio", "mu1_end_ratio", np.arange(0, 5.01, 0.2), "longest μ-labelled track: last 3 cm / median dQ/dx"),
    ("n_sp_showers", "shw_sp_n_20br1_showers", np.arange(-0.5, 6.5), "SP-tagger showers > 20 MeV at the vertex"),
    ("n_em", "n_em", np.arange(-0.5, 8.5), "EM particles in the candidate"),
    ("E_em_2", lambda d: d.E_em_all_2.fillna(0) if "E_em_all_2" in d else d.E_em_2.fillna(0),
     np.arange(0, 201, 10), "second EM shower energy (none: 0) [MeV]"),
    ("trk1_len", lambda d: d.trk1_len.fillna(0), np.arange(0, 201, 10), "longest track (none: 0) [cm]"),
    ("reco_Enu", "kine_reco_Enu", np.arange(0, 3001, 150), "reco neutrino energy [MeV]"),
    ("numu_score", "numu_score", np.linspace(-3, 3, 31), "Wire-Cell νμ CC score"),
    ("nue_score", "nue_score", np.arange(-15.5, 17.6, 1.0), "Wire-Cell νe score"),
    ("ub_numu", "ub_numu", np.arange(-3.6, 4.01, 0.4), "uB νμ BDT"),
    ("ub_other", "ub_other", np.arange(-4.4, 3.01, 0.4), "uB other BDT"),
    ("ub_ncpi0", "ub_ncpi0", np.arange(-3.0, 3.51, 0.5), "uB NC π⁰ BDT"),
    ("flash_pe", "flash_pe", np.linspace(0, 20000, 41), "matched flash PE"),
    ("flash_time", "flash_time_us", np.linspace(0, 2.5, 26), "matched flash time [µs]"),
]


def load_data(inputs, nproc=None):
    """On- or off-beam data: (candidates with the uB scores, POT from the hdr
    table, off-beam BNB gates)."""
    from analysis_village.SinglePhotonWC import microboone_bdt
    frames, totals = dfio.load(dfio.expand_inputs(inputs), keys=["wc", "wcevt", "hdr"], nproc=nproc)
    w = frames["wc"]
    w["wc_valid"] = 1   # a candidate implies Wire-Cell output
    w = w.join(microboone_bdt.score(w, UBOONE_BDT_DIR, nproc=nproc))
    hdr = frames["hdr"]
    if len(hdr):
        pot, gates = hdr.pot.sum(), hdr.noffbeambnb.sum()
    else:
        print("Warning: no hdr table (made without configs/sp_data.py?); POT from histpotdf, no off-beam gates.")
        pot, gates = totals.pot.sum(), 0
    return w, pot, gates


# variables whose preselection shapes fix the off-beam share when the gate count is unknown:
# cosmics enter from the top and leave little light in the beam window
FIT_VARIABLES = ["nu_y", "flash_pe", "numu_score"]


def fit_offbeam_fraction(mc, off, on):
    """Off-beam share f of the prediction minimising the summed Poisson
    deviance of the unit-area shapes (1 - f) MC + f off-beam vs on-beam over FIT_VARIABLES."""
    spec = {v[0]: v for v in VARIABLES}
    fs = np.linspace(0, 1, 201); chi2 = np.zeros(len(fs))
    for name in FIT_VARIABLES:
        _, var, bins, _ = spec[name]
        hm, ho, hd = (np.histogram(values(x, var), bins)[0].astype(float) for x in (mc, off, on))
        hm /= max(len(mc), 1); ho /= max(len(off), 1); n = len(on)   # per candidate, out-of-range ones included
        for i, f in enumerate(fs):
            p = np.maximum(n * ((1 - f) * hm + f * ho), 1e-9)
            chi2[i] += 2 * np.sum(p - hd + hd * np.log(np.where(hd > 0, hd / p, 1)))   # Poisson deviance
    return fs[np.argmin(chi2)]


def values(df, var):
    return np.asarray(var(df) if callable(var) else df[var], dtype=float)


def plot_one(name, var, bins, label, stage, title, mc, off, on, sc_mc, sc_off, out, blind, area_norm=False):
    x_mc = values(mc, var); x_off = values(off, var)
    fin = np.isfinite(x_mc); mc, x_mc = mc[fin], x_mc[fin]; x_off = x_off[np.isfinite(x_off)]
    stacks, weights, labels, colors = [], [], [], []
    for lab, cats, col in GROUPS:
        m = mc.cat.isin(cats).to_numpy()
        stacks.append(x_mc[m]); weights.append(np.full(m.sum(), sc_mc)); colors.append(col)
        labels.append("%s (%.1f)" % (lab, m.sum() * sc_mc))
    stacks.append(x_off); weights.append(np.full(len(x_off), sc_off)); colors.append(OFFBEAM_COLOR)
    labels.append("Off-beam data (%.1f)" % (len(x_off) * sc_off))
    pred = sum(np.histogram(x, bins, weights=w_)[0] for x, w_ in zip(stacks, weights))
    var_pred = (np.histogram(x_mc, bins)[0] * sc_mc ** 2 + np.histogram(x_off, bins)[0] * sc_off ** 2)
    if not blind:
        x_on = values(on, var); x_on = x_on[np.isfinite(x_on)]; n_on = np.histogram(x_on, bins)[0]
    if area_norm:   # scale the prediction to the data total (blinded: to unit area)
        k = (n_on.sum() if not blind else 1.0) / max(pred.sum(), 1e-12)
        weights = [w_ * k for w_ in weights]; pred = pred * k; var_pred = var_pred * k ** 2
        labels = [l.rsplit(" (", 1)[0] for l in labels]
    err = np.sqrt(var_pred)
    n_out = int(((x_mc < bins[0]) | (x_mc >= bins[-1])).sum() + ((x_off < bins[0]) | (x_off >= bins[-1])).sum())

    fig, (a, r) = plt.subplots(2, 1, figsize=(14, 8), sharex=True, gridspec_kw=dict(height_ratios=[3, 1], hspace=0.05))
    note = ""
    a.hist(stacks, bins=bins, weights=weights, stacked=True, color=colors, label=labels, edgecolor="white", linewidth=0.5)
    a.fill_between(bins, np.r_[pred - err, (pred - err)[-1]], np.r_[pred + err, (pred + err)[-1]], step="post",
                   facecolor="none", edgecolor="#3d3c39", hatch="////", linewidth=0, label="pred. stat. uncertainty")
    c = 0.5 * (bins[1:] + bins[:-1]); w_bin = np.diff(bins) / 2
    if blind:
        note = "   BLINDED: on-beam data not shown"
        r.text(0.5, 0.5, "blinded", transform=r.transAxes, ha="center", va="center", fontsize=12, color="#b00020")
    else:
        n_out += int(((x_on < bins[0]) | (x_on >= bins[-1])).sum())
        a.errorbar(c, n_on, yerr=np.sqrt(n_on), xerr=w_bin, fmt="o", color="k", ms=5, label="On-beam data (%d)" % n_on.sum())
        ok = pred > 0
        chi2 = np.sum(((n_on - pred) ** 2 / np.maximum(n_on + var_pred, 1e-9))[ok])
        note = "   χ²/ndf (stat.) = %.1f / %d, data/pred = %.2f" % (chi2, ok.sum(), n_on.sum() / max(pred.sum(), 1e-9))
        r.errorbar(c[ok], n_on[ok] / pred[ok], yerr=np.sqrt(n_on[ok]) / pred[ok], xerr=w_bin[ok], fmt="o", color="k", ms=4)
    rel = np.where(pred > 0, err / np.maximum(pred, 1e-9), 0)
    r.fill_between(bins, np.r_[1 - rel, 1 - rel[-1]], np.r_[1 + rel, 1 + rel[-1]], step="post", color="#c9c8c3", alpha=0.6, linewidth=0)
    r.axhline(1, color="k", lw=1); r.set_ylim(0, 2); r.set_ylabel("data / pred.")
    r.set_xlabel(label)
    a.set_ylabel("area-normalised" if area_norm and blind else
                 "candidates / bin" if area_norm else "candidates / bin")
    a.set_ylim(bottom=0)
    a.legend(fontsize=10, loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)
    a.set_title("%s%s%s" % (title, (", prediction area-normalised to data" if not blind else ", area-normalised") if area_norm else "", note), loc="left", fontsize=12,
                color="#b00020" if blind else "k")
    fig.savefig(os.path.join(out, "%s_%s%s.png" % (stage, name, "_areanorm" if area_norm else "")), dpi=120, bbox_inches="tight")
    plt.close(fig)
    return n_out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mc", nargs="+", required=True, help="MC .df files, directories or globs")
    ap.add_argument("--onbeam", nargs="+", required=True, help="on-beam data .df files, directories or globs")
    ap.add_argument("--offbeam", nargs="+", required=True, help="off-beam data .df files, directories or globs")
    ap.add_argument("-o", "--outdir", default="data_mc_plots")
    ap.add_argument("--area-norm", action="store_true", help="scale the prediction to the data total in each plot")
    ap.add_argument("--offbeam-frac", type=float, help="--area-norm without POT / gate info: off-beam share of the "
                    "prediction at preselection (default: fitted to the on-beam data)")
    ap.add_argument("-j", "--nproc", type=int, help="worker processes for reading the files and the BDT scoring (default: up to 8)")
    ap.add_argument("--stages", nargs="+", default=["presel", "bdt", "final"], choices=list(STAGES))
    ap.add_argument("--blind", nargs="*", default=[], choices=list(STAGES),
                    help="stages drawn without on-beam data (default: none)")
    ap.add_argument("--variables", nargs="+", help="subset of VARIABLES names (default: all)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    mc, _, pot_mc = load_mc(args.mc, UBOONE_BDT_DIR, nproc=args.nproc)
    on, pot_on, _ = load_data(args.onbeam, args.nproc)
    off, _, gates_off = load_data(args.offbeam, args.nproc)
    have_norm = pot_on > 0 and gates_off > 0
    if have_norm:
        pot_off = gates_off * N_GATES_ON_PER_5E12POT * 5e12
        sc_mc, sc_off = pot_on / pot_mc, pot_on / pot_off
        norm_line = "MC POT %.4g, on-beam POT %.4g, off-beam %d gates = %.4g POT-equivalent" % (pot_mc, pot_on, gates_off, pot_off)
    elif args.area_norm:
        pm, po, pd_ = (df[cuts.passes(df, cuts.PRESELECTION)] for df in (mc, off, on))
        f = args.offbeam_frac if args.offbeam_frac is not None else fit_offbeam_fraction(pm, po, pd_)
        sc_mc, sc_off = (1 - f) * len(pd_) / max(len(pm), 1), f * len(pd_) / max(len(po), 1)
        norm_line = ("no %s: off-beam share at preselection %s f = %.3f (MC weight %.4g, off-beam weight %.4g, "
                     "matched to the on-beam preselection total)" % (
                         "on-beam POT" if not pot_on > 0 else "off-beam gate count",
                         "set" if args.offbeam_frac is not None else "fitted over %s:" % "/".join(FIT_VARIABLES), f, sc_mc, sc_off))
        print("Warning: " + norm_line)
    else:
        sys.exit("No on-beam POT (hdr.pot / histpotdf) or no off-beam gates (hdr.noffbeambnb); make the data dataframes "
                 "with configs/sp_data.py, or use --area-norm for a shape-only comparison.")
    variables = [v for v in VARIABLES if not args.variables or v[0] in args.variables]

    lines = [norm_line,
             "weights: MC %.4g, off-beam %.4g%s" % (sc_mc, sc_off, "; plots area-normalised" if args.area_norm else ""), ""]
    for stage in args.stages:
        chain, title = STAGES[stage]
        m, o, d = (df[cuts.passes(df, chain)] for df in (mc, off, on))
        blind = stage in args.blind
        pred = len(m) * sc_mc + len(o) * sc_off
        lines.append("%-7s MC %8.1f (signal %.1f)  off-beam %8.1f  prediction %8.1f  on-beam %s" % (
            stage, len(m) * sc_mc, cuts.is_signal_cat(m.cat).sum() * sc_mc, len(o) * sc_off, pred,
            "blinded" if blind else "%d  (data/pred %.3f)" % (len(d), len(d) / max(pred, 1e-9))))
        for name, var, bins, label in variables:
            if not callable(var) and var not in mc:
                print("skip %s: column %s not in the MC dataframe" % (name, var))
                continue
            n_out = plot_one(name, var, bins, label, stage, title, m, o, d, sc_mc, sc_off, args.outdir, blind, args.area_norm)
            if n_out:
                print("%s %s: %d entries outside the range, not drawn" % (stage, name, n_out))
    open(os.path.join(args.outdir, "normalisation.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("wrote the plots to %s" % args.outdir)


if __name__ == "__main__":
    main()
