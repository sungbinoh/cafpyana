#!/usr/bin/env python3
"""One panel per variable of the per-candidate Wire-Cell dataframe, for
learning what each column holds before building a selection.

Input: the run_df_maker output of configs/wc_candidates.py (key "wc", one row
per Wire-Cell neutrino candidate, see make_wc_candidates.py), or CAFs, on
which make_wc_candidate_df is run here. Writes <outdir>/wc_variables.pdf and
a PNG per page, the variables grouped by topic (GROUPS). Each panel has the
column's meaning and units (DESC, from the Wire-Cell source:
HaiwangYu/wire-cell-toolkit branch sbnd-wc-flatcaf), N and the range.

  * columns with few integer values are drawn as counts per value
  * continuous columns are histogrammed between the 0.5 and 99.5 percentiles;
    the rows outside are given as under/overflow
  * sentinel values (SENTINELS, e.g. nue_score = -15 "not evaluated") and
    NaN (e.g. no proton -> no E_proton_1) are counted in the panel and left
    out of the histogram

Nothing is selected: every candidate is used.

  python run_df_maker.py -c analysis_village/NCDeltaWC/configs/wc_candidates.py -o wc_candidates -i <cafs>
  python analysis_village/NCDeltaWC/plot_wc_variables.py -i wc_candidates.df -o plots/
"""

import argparse
import glob
import os
import sys
import textwrap

import numpy as np
import pandas as pd
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

COLOR = "#2a78d6"
INK2 = "#52514e"
plt.rcParams.update({
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": "#e6e5e1", "grid.linewidth": 0.6,
    "axes.axisbelow": True, "font.size": 9,
})

KEY = "wc"
NCOL, NROW = 4, 3
MAX_DISCRETE = 25   # integer columns with at most this many values: counts per value

# Values that mean "not set"; counted, not histogrammed
SENTINELS = {"nue_score": [-15], "has_vertex": [-1], "tmatch_idx": [-1], "main_cluster_lm": [-1],
             "nu_n_flashes": [-1], "nu_n_flashes_in_window": [-1], "nu_n_bundles": [-1]}

# Columns only meaningful on some rows: drawn for those rows only
PIO = lambda df: df.kine_pio_flag > 0
CONDITIONS = {c: ("kine_pio_flag > 0", PIO) for c in [
    "kine_pio_mass", "kine_pio_vtx_dis", "kine_pio_angle", "kine_pio_energy_1", "kine_pio_energy_2",
    "kine_pio_theta_1", "kine_pio_theta_2", "kine_pio_phi_1", "kine_pio_phi_2",
    "kine_pio_dis_1", "kine_pio_dis_2"]}

DESC = {}
GROUPS = []   # (title, [columns]) in page order


def group(title, rows):
    cols = []
    for names, text in rows:
        for c in names.split():
            DESC[c] = text
            cols.append(c)
    GROUPS.append((title, cols))


group("Event", [
    ("nu_index", "Candidate index within the event (row of T_tagger / T_kine)"),
    ("run subrun evt", "Run / subrun / event number of the CAF event"),
    ("n_nu", "Candidates in the event (this one included)"),
    ("nu_n_flashes", "Trun: optical flashes in the event (both TPCs, whole readout)"),
    ("nu_n_flashes_in_window", "Trun: flashes in the beam window [0.2, 2.2) us"),
    ("nu_n_bundles", "Trun: in-window flash bundles examined (merged flash groups count once)"),
])
group("Candidate selection: flash, main cluster, bundle", [
    ("matched_flash_gid", "Global id of the matched flash: flash index + 1000000 x TPC"),
    ("flash_time_us", "Matched flash time [us]"),
    ("flash_pe", "Matched flash light [PE]"),
    ("flash_tpc", "TPC of the matched flash"),
    ("flash_both_tpc", "1 = the matched flash was also seen in the other TPC (its flash group has a flash of each TPC within 0.05 us)"),
    ("cluster_id", "Main cluster of the candidate, after the vertex finder (vertex, particles, segments belong to it)"),
    ("main_cluster_npoints", "3D imaging points in the main cluster"),
    ("main_cluster_length_cm", "Main cluster length [cm]"),
    ("main_cluster_fc", "Main cluster fully contained (TaggerCheck FC; differs from match_isFC)"),
    ("main_cluster_tgm", "Main cluster tagged through-going muon"),
    ("main_cluster_stm", "Main cluster tagged stopping muon"),
    ("main_cluster_lm", "Q/L lm_flag of the main cluster: -1 unset, 1 low energy, 2 light mismatch"),
    ("bundle_n_main", "In-window main clusters in the candidate's flash bundle"),
])
group("Vertex and tagger flags", [
    ("nu_x nu_y nu_z", "Reconstructed neutrino vertex [cm]"),
    ("kine_nu_x_corr kine_nu_y_corr kine_nu_z_corr", "Corrected vertex position [cm]"),
    ("has_vertex", "1 = vertex found, 0 = none (taggers do not run), -1 = not filled"),
    ("match_isFC", "1 = main cluster fully contained in the PR fiducial volume"),
    ("neutrino_type", "Tagger bitmask: 2 cosmic, 4 numu CC, 8 NC (numu tagger: not CC), 32 nue; 0 = no vertex"),
    ("cosmict_flag", "Overall cosmic verdict: OR of the 10 cosmic sub-taggers (1 = cosmic)"),
    ("numu_cc_flag", "1 = numu CC tagger fired (any of its 3 criteria)"),
    ("photon_flag", "1 = single-photon tagger fired"),
    ("shw_sp_filled", "1 = single-photon tagger variables (shw_sp_*) were computed"),
    ("numu_score", "numu CC BDT score (log-odds; higher = more numu CC)"),
    ("nue_score", "nue CC BDT score (log-odds); -15 = not evaluated (no good vertex shower > 60-80 MeV)"),
])
group("Energy and pi0 (T_kine)", [
    ("kine_reco_Enu", "Reconstructed neutrino energy: all particle energies + add_energy [MeV]"),
    ("kine_reco_add_energy", "Rest masses and binding energy added to kine_reco_Enu [MeV]"),
    ("kine_energy_excluded", "Energy of activity left out of the candidate [MeV]"),
    ("kine_n_excluded", "Objects left out of the candidate"),
    ("kine_main_vertex_id", "Id of the main vertex in the PR graph"),
    ("kine_pio_flag", "pi0 found: 0 no, 1 with vertex, 2 without vertex"),
    ("kine_pio_mass", "pi0 candidate diphoton mass [MeV] (often a pairing of small showers)"),
    ("kine_pio_vtx_dis", "pi0 candidate distance to the neutrino vertex [cm]"),
    ("kine_pio_angle", "Opening angle of the two pi0 photons [deg]"),
    ("kine_pio_energy_1 kine_pio_energy_2", "Energy of pi0 photon 1 / 2 [MeV]"),
    ("kine_pio_theta_1 kine_pio_theta_2", "Polar angle of pi0 photon 1 / 2 [deg]"),
    ("kine_pio_phi_1 kine_pio_phi_2", "Azimuthal angle of pi0 photon 1 / 2 [deg]"),
    ("kine_pio_dis_1 kine_pio_dis_2", "Conversion distance of pi0 photon 1 / 2 [cm]"),
])
group("Particles (T_kine list, connected particles: kine_energy_included == 1)", [
    ("n_particles", "Connected particles"),
    ("n_shower", "Connected shower objects (any PID label)"),
    ("n_em", "Connected particles labelled 11 (e / EM shower; photons included)"),
    ("n_muon", "Connected particles labelled 13"),
    ("n_pion", "Connected particles labelled 211"),
    ("n_proton", "Connected particles labelled 2212"),
    ("n_nopid", "Connected particles without a PID label (0)"),
    ("n_unclear", "Leftover showers with an unclear link (included == 3; still in kine_reco_Enu)"),
    ("E_unclear", "Summed energy of the unclear-link showers [MeV]"),
    ("E_em_1", "Energy of the most energetic connected EM shower (label 11) [MeV]; NaN = no EM shower"),
    ("E_em_2", "Energy of the second most energetic connected EM shower (label 11) [MeV]; NaN = fewer than 2"),
    ("E_muon_1", "Energy of the leading muon-labelled particle [MeV]; NaN = none"),
    ("E_pion_1", "Energy of the leading pion-labelled particle [MeV]; NaN = none"),
    ("E_proton_1", "Energy of the leading proton-labelled particle [MeV]; NaN = none"),
])
group("Segments and charge points", [
    ("n_segments", "T_segment rows of the candidate (fitted or not)"),
    ("n_segments_in_enu", "Segments whose energy is in kine_reco_Enu"),
    ("trk_len_max", "Longest track-like segment [cm]; NaN = no track segment"),
    ("trk_len_sum", "Summed length of track-like segments [cm]"),
    ("shw_len_max", "Longest shower-like segment [cm]; NaN = no shower segment"),
    ("n_charge_points", "T_rec_charge fitted 3D points of the candidate"),
])
group("MC truth match (T_match level 3, T_pf)", [
    ("tmatch_idx", "Matched true neutrino (row of T_truth_nu = rec.mc.nu); -1 = none"),
    ("tmatch_iou_n", "Point IoU with the matched true interaction; partners ranked by this"),
    ("tmatch_iou_q", "Charge IoU with the matched true interaction"),
    ("pf_n_nodes", "Nodes in the candidate's T_pf particle-flow tree"),
    ("pf_n_primary", "Primary nodes (depth 0)"),
    ("pf_n_pseudo", "Inferred pseudo-particles (photon of a gap shower, neutron of a detached proton)"),
])
group("Matched true neutrino (rec.mc.nu, GENIE)", [
    ("true_pdg", "Neutrino PDG: +-12 nue, +-14 numu"),
    ("true_iscc", "1 = CC"),
    ("true_isnc", "1 = NC"),
    ("true_genie_mode", "GENIE mode (0 QE, 1 RES, 2 DIS, 3 COH, 10 MEC)"),
    ("true_resnum", "GENIE resonance; 0 = Delta(1232) (only meaningful for RES)"),
    ("true_E", "True neutrino energy [GeV]"),
    ("true_w", "Hadronic invariant mass W [GeV]"),
    ("true_vtx_x true_vtx_y true_vtx_z", "True vertex [cm]"),
    ("true_npi0", "Primary pi0s"),
    ("true_ngamma", "Primary photons, any energy (incl. ~3-9 MeV nuclear de-excitation)"),
    ("true_ngamma_20MeV true_ngamma_100MeV", "Primary photons above 20 / 100 MeV"),
    ("true_leading_gamma_E", "Most energetic primary photon energy [GeV]"),
    ("true_np", "Primary protons, any energy"),
    ("true_np_20MeV true_np_35MeV true_np_50MeV", "Primary protons above 20 / 35 / 50 MeV kinetic energy"),
    ("true_nn", "Primary neutrons"),
    ("true_npi_30MeV", "Primary charged pions above 30 MeV kinetic energy"),
    ("true_nmu_27MeV", "Primary muons above 27 MeV kinetic energy"),
])
group("Matched true photons, pi0 and proton (make_truth.py)", [
    ("true_g1_E true_g2_E", "Leading / 2nd primary photon (> 20 MeV) energy [GeV]"),
    ("true_g1_costh true_g2_costh", "Primary photon direction, cos(angle to beam)"),
    ("true_g1_gap true_g2_gap", "Distance from creation to conversion [cm]; NaN unless made and ending in the TPC"),
    ("true_g1_conv true_g2_conv", "1 = trajectory ends in pair conversion (0: photoabsorption)"),
    ("true_g1_lost true_g2_lost", "1 = made in the TPC but leaves before it ends"),
    ("true_p1_KE", "Leading primary proton kinetic energy [GeV]"),
    ("true_p1_len", "Leading proton in-TPC track length [cm]"),
    ("true_p1_costh", "Leading proton cos(angle to beam)"),
    ("true_p1_cos_gp", "cos(angle) between leading photon and leading proton"),
    ("true_pi0_E", "Leading primary pi0 total energy [GeV]"),
    ("true_pi0_g1_E true_pi0_g2_E", "pi0 decay photons, leading / sub-leading energy [GeV]"),
    ("true_pi0_g1_costh true_pi0_g2_costh", "pi0 decay photon cos(angle to beam)"),
    ("true_pi0_g1_gap true_pi0_g2_gap", "pi0 decay photon conversion distance [cm]; NaN unless made and ending in the TPC"),
    ("true_pi0_g1_conv true_pi0_g2_conv", "pi0 decay photon ends in pair conversion"),
    ("true_pi0_g1_lost true_pi0_g2_lost", "pi0 decay photon leaves the TPC before it ends"),
    ("true_pi0_open", "pi0 decay photon opening angle [deg]"),
    ("true_pi0_asym", "pi0 energy asymmetry (E1 - E2) / (E1 + E2)"),
])


def expand(patterns):
    files = []
    for p in patterns:
        files += sorted(glob.glob(p)) if any(ch in p for ch in "*?[") else [p]
    return files


def read(files):
    parts = []
    for i, p in enumerate(files):
        if p.endswith(".df"):
            with pd.HDFStore(p, mode="r") as st:
                keys = set(st.keys())
                n = int(st["split"].n_split.iloc[0]) if "/split" in keys else 1
                parts += [st["%s_%d" % (KEY, j)] for j in range(n) if "/%s_%d" % (KEY, j) in keys]
        else:
            # imported here: makedf.makedf needs CVMFS, and .df inputs do not
            from analysis_village.NCDeltaWC.make_wc_candidates import make_wc_candidate_df
            parts.append(make_wc_candidate_df(uproot.open(p)))
    return pd.concat(parts).reset_index() if parts else pd.DataFrame()


def panel(ax, df, col):
    s = pd.to_numeric(df[col], errors="coerce")
    notes = []
    if col in CONDITIONS:
        text, cond = CONDITIONS[col]
        s = s[cond(df)]
        notes.append("only %s" % text)
    n_all = len(s)
    notes.insert(0, "N=%d" % n_all)
    if s.isna().any():
        notes.append("NaN: %d (%.0f%%)" % (s.isna().sum(), 100 * s.isna().mean()))
    s = s.dropna()
    for v in SENTINELS.get(col, []):
        k = int((s == v).sum())
        if k:
            notes.append("%g: %d (%.0f%%)" % (v, k, 100 * k / n_all))
    x = s[~s.isin(SENTINELS.get(col, []))]

    if not len(x):
        ax.text(0.5, 0.5, "no values", transform=ax.transAxes, ha="center", va="center", color=INK2)
    elif (x == np.round(x)).all() and x.nunique() <= MAX_DISCRETE:
        vc = x.astype(int).value_counts().sort_index()
        ax.bar(vc.index.astype(str), vc.values, color=COLOR, width=0.7)
        ax.tick_params(axis="x", labelsize=6, rotation=90 if len(vc) > 10 else 0)
        ax.grid(axis="x", visible=False)
        if vc.max() > 50 * max(vc.min(), 1):
            ax.set_yscale("log")
    else:
        lo, hi = np.percentile(x, [0.5, 99.5])
        if lo == hi:
            lo, hi = lo - 0.5, hi + 0.5
        h = x[(x >= lo) & (x <= hi)]
        if (x == np.round(x)).all() and hi - lo < 200:
            bins = np.arange(np.floor(lo) - 0.5, np.ceil(hi) + 1.5)
        else:
            bins = np.linspace(lo, hi, 51)
        counts, _, _ = ax.hist(h, bins=bins, histtype="stepfilled", color=COLOR, alpha=0.25, lw=0)
        ax.hist(h, bins=bins, histtype="step", color=COLOR, lw=1.5)
        if counts.max() > 100 * max(np.median(counts[counts > 0]), 1):
            ax.set_yscale("log")
        uf, of = int((x < lo).sum()), int((x > hi).sum())
        if uf or of:
            notes.append("under %d / over %d" % (uf, of))
        ax.tick_params(axis="x", labelsize=7)
    if len(x):
        notes.append("range [%.4g, %.4g]" % (x.min(), x.max()))

    ax.set_title(col, fontsize=9, fontweight="bold", loc="left")
    ax.tick_params(axis="y", labelsize=7)
    ax.text(0.0, -0.22, "\n".join(textwrap.wrap(DESC.get(col, "(no description)"), 62)),
            transform=ax.transAxes, fontsize=6.5, va="top", color="#0b0b0b")
    ax.text(0.99, 0.97, "\n".join(notes), transform=ax.transAxes, fontsize=6, ha="right", va="top",
            color=INK2, bbox=dict(fc="white", ec="none", alpha=0.8, pad=1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--inputs", nargs="+", required=True, help=".df files or CAFs, or glob patterns")
    ap.add_argument("-o", "--outdir", default="wc_variable_plots")
    args = ap.parse_args()

    df = read(expand(args.inputs))
    if not len(df):
        sys.exit("No candidates read.")
    os.makedirs(args.outdir, exist_ok=True)

    skip = {"__ntuple", "entry", "file"}
    known = {c for _, cols in GROUPS for c in cols}
    groups = [(t, [c for c in cols if c in df]) for t, cols in GROUPS]
    groups.append(("Other", [c for c in df.columns if c not in known and c not in skip]))

    per = NCOL * NROW
    page = 0
    path = os.path.join(args.outdir, "wc_variables.pdf")
    with PdfPages(path) as pdf:
        for title, cols in groups:
            chunks = [cols[i:i + per] for i in range(0, len(cols), per)]
            for ic, chunk in enumerate(chunks):
                page += 1
                fig, axes = plt.subplots(NROW, NCOL, figsize=(14, 9.5))
                fig.subplots_adjust(hspace=0.75, wspace=0.28, top=0.91, bottom=0.08, left=0.05, right=0.98)
                fig.suptitle("%s%s   (one row per candidate, %d candidates)"
                             % (title, " [%d/%d]" % (ic + 1, len(chunks)) if len(chunks) > 1 else "", len(df)),
                             fontsize=11, fontweight="bold")
                for ax, c in zip(axes.ravel(), chunk):
                    panel(ax, df, c)
                for ax in axes.ravel()[len(chunk):]:
                    ax.axis("off")
                pdf.savefig(fig)
                fig.savefig(os.path.join(args.outdir, "var_%02d.png" % page))
                plt.close(fig)
            if cols:
                print("%-70s %3d variables" % (title, len(cols)))
    if groups[-1][1]:
        print("No description for:", ", ".join(groups[-1][1]))
    print("Wrote %s (%d pages)" % (path, page))


if __name__ == "__main__":
    main()
