"""SBND off-beam data vs. CORSIKA in-time cosmic MC (GUMPLE / sbn-rewgted-2x).

Script port of nb/OffBeamCosmicMCComparison.ipynb (which read the reCAF
sbn-rewgted-10/-11 productions through the deleted gump_cuts module). The
comparison behind the CORSIKA-modeling systematic of the GUMP technote:

  * SBND_SpringBNBOffData.df  -- off-beam data, one gate per off-beam trigger
                                 record (hdr.noffbeambnb)
  * SBNDIntimeMC.df           -- CORSIKA in-time cosmic MC, one gate per
                                 generated event (hdr.ngenevt, once per subrun)

MC is gate-normalized to the data (MC_w = ngates_OFF / ngates_MC) and the two
are compared at the mcdata_comparison_gumple.py selection stages -- Contained,
Cosmic Rej., Two Prong Cut, PID 1p, PID Np -- for the muon-candidate length,
cos(theta), the calorimetric neutrino energy and delta p. Each figure comes in
two flavors: plain (MC-stat error only) and "_wnorm", which folds in a
fully-correlated normalization uncertainty. By default (--norm-unc nominal)
that is the analysis' CORSIKA normalization uncertainty,
mcdata_comparison_gumple.SBND_COSMIC_NORM (10.7%); --norm-unc measured uses
|MC/data - 1| after the full 1p selection of the production being plotted
(always printed, so the nominal can be re-checked after every production);
a number gives it explicitly.

Outputs (in --plotdir):
  pdf/SBND_<stage>_<var>[_wnorm].pdf, png/...   stage in contained, cosmicrej,
                                                twoprongcut, pid1p, pidnp
  SBND_offbeam_cosmicmc_counts.csv              per-stage counts + chi2 table
  SBND_offbeam_cosmicmc_summary.txt

    python offbeam_cosmicmc_gumple.py --df-dir ../sbn-rewgted-21/ --plotdir <dir>
                                      [--norm-unc nominal|measured|<float>]

Run with cwd = analysis_village/gump (kinematics.py's cwd-relative sys.path).
"""

import argparse
import os
import sys
import warnings

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(_HERE, "..", "gumple"))
import gumple_cuts as gc  # noqa: E402

# The analysis' nominal CORSIKA normalization uncertainty (mcdata_comparison_gumple
# .SBND_COSMIC_NORM = 0.107); read from the module so the two cannot drift apart.
from mcdata_comparison_gumple import SBND_COSMIC_NORM as NOMINAL_NORM_UNC  # noqa: E402

DETECTOR = "SBND"
DEFAULT_DF_DIR = "/Users/gputnam/Work/osc/sbn-rewgted-21/"
DEFAULT_PLOTDIR = "/Users/gputnam/Work/osc/cafpyana/plots-gumple-2026-09-02-rewgted21/corsika/"

OFFBEAM_FILE = "SBND_SpringBNBOffData.df"
COSMIC_MC_FILE = "SBNDIntimeMC.df"

# Best-fit MC flash-PE scale for SBND, as applied in loaddf.load_one
MC_PE_SCALE = 0.642

FONTSIZE = 14
MC_COLOR = "#95af8b"

PLOTVARS = ["mu_len", "mu_costh", "nu_E_calo", "del_p"]
PLOTBINS = [np.linspace(0, 500, 21), np.linspace(-1, 1, 21), np.linspace(0, 2, 21),
            np.linspace(0, 1, 11)]
PLOTLABELS = ["Muon Candidate Length [cm]", "Muon Candidate $\\cos\\theta$",
              "Reco. Neutrino Energy [GeV]", "$\\delta p$ [GeV/c]"]


# ----------------------------------------------------------------------------
# I/O
# ----------------------------------------------------------------------------
def read_dfs(file, key):
    with h5py.File(file, "r") as f:
        keys = [k for k in f.keys() if k.startswith(key)]
    return pd.concat([pd.read_hdf(file, k) for k in keys])


def FV(df):
    """The mcdata_comparison_gumple.py 'Contained' preselection."""
    return gc.sanity_cut(df) & gc.slcfv_cut(df) & df.cut_contained & df.cut_cathode


def load_evt(fname, pe_scale):
    with h5py.File(fname, "r") as f:
        keys = [k for k in f.keys() if k.startswith("evt")]
    dfs = []
    for k in keys:
        df = pd.read_hdf(fname, k)
        if "Run" not in df.columns:
            df["Run"] = 1
        if "detector" not in df.columns:
            df["detector"] = DETECTOR
        df["flash_maxpe"] = (df["flash_maxpe"] * pe_scale).fillna(0.).astype(float)
        df = df[FV(df)]
        dfs.append(df)
    return pd.concat(dfs)


# ----------------------------------------------------------------------------
# Selection stages -- mirror mcdata_comparison_gumple.CUTS / CUTNAMES
# ----------------------------------------------------------------------------
def contained(d):
    return FV(d)


def cosmic_rej(d):
    return contained(d) & gc.cosmic_cut(d)


def twoprong_cut(d):
    return cosmic_rej(d) & gc.trk_cut(d)


def pid_cut(d):
    return twoprong_cut(d) & gc.pid_cut(d)


def pid_cut_1p(d):
    return pid_cut(d) & (d.n_pfp == 2)


def pid_cut_np(d):
    return pid_cut(d) & (d.n_pfp > 2)


CUTS = [contained, cosmic_rej, twoprong_cut, pid_cut_1p, pid_cut_np]
CUTNAMES = ["Contained", "Cosmic Rej.", "Two Prong Cut", "PID 1p", "PID Np"]
NORM_CUT = "PID 1p"   # stage that defines the normalization uncertainty


def cut_stem(cname):
    return cname.replace(" ", "").replace(".", "").lower()


# ----------------------------------------------------------------------------
# Histogramming / chi2
# ----------------------------------------------------------------------------
def f_chi2(NMC, Ndata, cov):
    which_bin = NMC > 0  # ignore singular entries
    NMC = NMC[which_bin]
    Ndata = Ndata[which_bin]
    cov = cov[which_bin, :][:, which_bin]
    delta = NMC - Ndata
    try:
        cov_inv = np.linalg.inv(cov)
    except np.linalg.LinAlgError:
        return -1, which_bin.sum()
    return delta @ cov_inv @ delta, which_bin.sum()


def make_plot_data(MCdf, OFFdf, MC_w, var, bins, cut, xlabel, title, norm_unc=None):
    NMC_raw, bins = np.histogram(MCdf.loc[MCdf[cut], var], bins=bins)
    NMC = NMC_raw * MC_w
    NMCstaterr = np.sqrt(NMC_raw) * MC_w

    NData, _ = np.histogram(OFFdf.loc[OFFdf[cut], var], bins=bins)
    NData = NData.astype(float)
    NDataErr = np.sqrt(NData)

    cov = np.diag(NMCstaterr**2)  # MC-stat
    if norm_unc is not None:
        # fully-correlated normalization uncertainty
        cov = cov + norm_unc**2 * np.outer(NMC, NMC)
    NMCerr = np.sqrt(np.diag(cov))

    cov_w_stat = cov + np.diag(NDataErr**2)
    chi2, ndof = f_chi2(NMC, NData, cov_w_stat)

    return {"title": title, "xlabel": xlabel, "bins": bins,
            "NMC": NMC, "NMCerr": NMCerr, "NData": NData, "NDataErr": NDataErr,
            "chi2": chi2, "ndof": ndof}


def ratio_plot(plotdata, ngates_OFF):
    fig, (ax0, ax1) = plt.subplots(2, 1, height_ratios=[3, 1], sharex=True)
    bins = plotdata["bins"]
    centers = (bins[:-1] + bins[1:]) / 2

    NMC = plotdata["NMC"]
    err = plotdata["NMCerr"]
    ax0.hist(centers, bins=bins, weights=NMC, color=MC_COLOR, label="InTime Cosmic MC")

    NData = plotdata["NData"]
    NDataErr = plotdata["NDataErr"]
    line = ax0.errorbar(centers, NData, NDataErr, color="black", linestyle="none", marker=".")

    # step="post" holds each y from x[i] to x[i+1]: step over the full edge
    # array with the last value repeated so the last bin is shaded too
    _hi, _lo = NMC + err, NMC - err
    ax0.fill_between(bins, np.append(_hi, _hi[-1]), np.append(_lo, _lo[-1]),
                     facecolor="none", hatch="//", edgecolor="gray", linewidth=0.0, step="post")

    with np.errstate(divide="ignore", invalid="ignore"):
        ax1.errorbar(centers, NData / NMC, NDataErr / NMC, color="black", linestyle="none", marker=".")
        _hi, _lo = 1 + err / NMC, 1 - err / NMC
        ax1.fill_between(bins, np.append(_hi, _hi[-1]), np.append(_lo, _lo[-1]),
                         facecolor="none", hatch="//", edgecolor="gray", linewidth=0.0, step="post")
    ax1.set_ylim([0.5, 1.5])
    ax1.axhline(1, color="red", linestyle="--")

    for ax in (ax0, ax1):
        ax.tick_params(axis="both", which="both", direction="in", length=6, width=1.5,
                       labelsize=FONTSIZE, top=True, right=True)
    for spine in ax0.spines.values():
        spine.set_linewidth(1.5)

    ax1.set_xlabel(plotdata["xlabel"], fontsize=FONTSIZE, fontweight="bold")
    ax0.set_ylabel("Events / %.2f$\\times 10^{7}$ Gates" % (ngates_OFF / 1e7),
                   fontsize=FONTSIZE, fontweight="bold")
    ax1.set_ylabel("Data / MC", fontsize=FONTSIZE - 2)

    title = plotdata["title"]
    ax0.set_title(f"$\\bf{{{DETECTOR}}}$ {title}", fontsize=FONTSIZE + 2)

    ld = ax0.legend([line], ["OFF-Beam Data"], frameon=False, loc="upper left", fontsize=10)
    ax0_lo, ax0_hi = ax0.get_ylim()
    ax0.set_ylim([ax0_lo, ax0_hi * 1.7])
    ax0.legend(fontsize=12, loc="upper right")
    ax0.add_artist(ld)

    ax0.text(0.05, 0.8, "$\\chi^2$: %.1f / %i" % (plotdata["chi2"], plotdata["ndof"]),
             transform=ax0.transAxes, fontsize=FONTSIZE - 2)
    return fig, ax0, ax1


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--df-dir", default=DEFAULT_DF_DIR)
    p.add_argument("--plotdir", default=DEFAULT_PLOTDIR)
    p.add_argument("--norm-unc", default="nominal",
                   help="normalization uncertainty for the _wnorm figures: 'nominal' "
                        "(default; mcdata_comparison_gumple.SBND_COSMIC_NORM = %.3f), "
                        "'measured' (|MC/data - 1| after PID 1p on this production), "
                        "or a number" % NOMINAL_NORM_UNC)
    p.add_argument("--no-save", action="store_true")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    df_dir = args.df_dir if args.df_dir.endswith("/") else args.df_dir + "/"
    plotdir = args.plotdir
    if not args.no_save:
        os.makedirs(plotdir, exist_ok=True)
        os.makedirs(os.path.join(plotdir, "png"), exist_ok=True)
        os.makedirs(os.path.join(plotdir, "pdf"), exist_ok=True)

    lines = []

    def log(msg):
        print(msg, flush=True)
        lines.append(str(msg))

    OFFBEAM = df_dir + OFFBEAM_FILE
    COSMIC_MC = df_dir + COSMIC_MC_FILE
    log("off-beam data:  %s" % OFFBEAM)
    log("cosmic MC:      %s" % COSMIC_MC)

    # ---- gate normalization ----
    hdr_OFF = read_dfs(OFFBEAM, "hdr")
    hdr_MC = read_dfs(COSMIC_MC, "hdr")
    ngates_OFF = hdr_OFF.noffbeambnb.sum()
    ngates_MC = hdr_MC.loc[hdr_MC.first_in_subrun == 1, "ngenevt"].sum()
    MC_w = ngates_OFF / ngates_MC
    log("N gates OFF-beam: %.4g" % ngates_OFF)
    log("N gates cosmic MC: %.4g" % ngates_MC)
    log("MC weight (gate ratio): %.4f" % MC_w)
    log("Events per gate -- OFF-beam: %.4f%%, cosmic MC: %.4f%%" %
        (100. * len(hdr_OFF) / ngates_OFF, 100. * len(hdr_MC) / ngates_MC))

    # ---- events ----
    OFFdf = load_evt(OFFBEAM, 1.0)
    MCdf = load_evt(COSMIC_MC, MC_PE_SCALE)
    log("OFF-beam candidates after FV: %i" % len(OFFdf))
    log("Cosmic MC candidates after FV: %i" % len(MCdf))

    for c, cname in zip(CUTS, CUTNAMES):
        for frame in (MCdf, OFFdf):
            frame[cname] = c(frame) & gc.flash_cut(frame)

    # the track direction is a unit vector, so dir_z is already cos(theta)
    for frame in (MCdf, OFFdf):
        frame["mu_costh"] = frame.mu_dir_z

    # ---- normalization uncertainty ----
    n_data_sel = OFFdf[NORM_CUT].sum()
    n_mc_sel = MCdf[NORM_CUT].sum() * MC_w
    MEASURED_NORM_UNC = n_mc_sel / n_data_sel - 1.
    log("Selected (%s) -- data: %i, MC prediction: %.1f" % (NORM_CUT, n_data_sel, n_mc_sel))
    log("Measured normalization offset (MC/data - 1): %.4f   "
        "(compare mcdata_comparison_gumple.SBND_COSMIC_NORM = %.4f)"
        % (MEASURED_NORM_UNC, NOMINAL_NORM_UNC))
    if args.norm_unc == "nominal":
        NORM_UNC = NOMINAL_NORM_UNC
    elif args.norm_unc == "measured":
        NORM_UNC = abs(MEASURED_NORM_UNC)
    else:
        NORM_UNC = abs(float(args.norm_unc))
    log("Normalization uncertainty applied in the _wnorm figures: %.4f (--norm-unc %s)"
        % (NORM_UNC, args.norm_unc))

    # ---- plots ----
    plt.rcParams["figure.max_open_warning"] = 0
    all_plotdata = {}
    nfig = 0
    for cname in CUTNAMES:
        all_plotdata[cname] = {}
        for v, b, l in zip(PLOTVARS, PLOTBINS, PLOTLABELS):
            for suffix, nunc in [("", None), ("_wnorm", NORM_UNC)]:
                pdat = make_plot_data(MCdf, OFFdf, MC_w, v, b, cname, l, cname, norm_unc=nunc)
                all_plotdata[cname][v + suffix] = pdat
                fig, _, _ = ratio_plot(pdat, ngates_OFF)
                if not args.no_save:
                    stem = "%s_%s_%s%s" % (DETECTOR, cut_stem(cname), v, suffix)
                    fig.savefig(os.path.join(plotdir, "pdf", stem + ".pdf"), bbox_inches="tight")
                    fig.savefig(os.path.join(plotdir, "png", stem + ".png"), bbox_inches="tight")
                    nfig += 1
                plt.close(fig)

    # ---- counts / chi2 table ----
    rows = []
    for cname in CUTNAMES:
        n_data = float(OFFdf[cname].sum())
        n_mc_raw = float(MCdf[cname].sum())
        n_mc_pred = n_mc_raw * MC_w
        row = {"cut": cname, "n_offbeam": n_data, "n_mc_raw": n_mc_raw, "n_mc_pred": n_mc_pred,
               "data_over_mc": n_data / n_mc_pred if n_mc_pred > 0 else np.nan}
        for v in PLOTVARS:
            row["chi2_%s" % v] = all_plotdata[cname][v]["chi2"]
            row["ndof_%s" % v] = all_plotdata[cname][v]["ndof"]
            row["chi2_wnorm_%s" % v] = all_plotdata[cname][v + "_wnorm"]["chi2"]
        rows.append(row)
    countdf = pd.DataFrame(rows)
    log("")
    log(countdf.to_string(index=False))
    if not args.no_save:
        countdf.to_csv(os.path.join(plotdir, "%s_offbeam_cosmicmc_counts.csv" % DETECTOR), index=False)
        with open(os.path.join(plotdir, "%s_offbeam_cosmicmc_summary.txt" % DETECTOR), "w") as f:
            f.write("\n".join(lines) + "\n")
        log("\nSaved %d figures (pdf+png) to %s" % (nfig, plotdir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
