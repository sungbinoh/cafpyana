#!/usr/bin/env python
"""LEGACY: compare CAFANA-MAPLE output (SB_GUMP*.root) with a cafpyana-maple .df.

This script validated the original chi2-in-candidate MAPLE port against the
CAFANA macro (see DIFFERENCES.md, 2026-08-11).  The evt-df builder has since
moved to chi2-free candidate finding with post-hoc chi2 cuts (maple_sel.py)
and GUMP column names, which intentionally breaks bit-compatibility with
CAFANA and removed the maple_evt_cafanapid.py config this script requires.
It can only be run against dataframes produced before that rework.

Usage:
  python compare_cafana_cafpyana.py <SB_GUMP.root> <maple.df> [--rtol 1e-4] [--atol 1e-4]
"""
import argparse
import sys

import h5py
import numpy as np
import pandas as pd
import uproot


def load_df(fname, table):
    with h5py.File(fname, "r") as f:
        keys = [k for k in f.keys() if k.rsplit("_", 1)[0] == table]
    if not keys:
        return None
    dfs = [pd.read_hdf(fname, k) for k in sorted(keys, key=lambda k: int(k.rsplit("_", 1)[1]))]
    return pd.concat(dfs)


def evt_key_frame(evt, hdr):
    """Attach run/subrun/evt from hdr to the evt table (index __ntuple, entry, slc)."""
    h = hdr[["run", "subrun", "evt"]].copy()
    h.index = h.index.droplevel(list(range(2, h.index.nlevels))) if h.index.nlevels > 2 else h.index
    out = evt.join(h, on=list(h.index.names))
    return out


def compare_var(name, a, b, rtol, atol):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    both_nan = np.isnan(a) & np.isnan(b)
    close = np.isclose(a, b, rtol=rtol, atol=atol) | both_nan
    n_bad = int((~close).sum())
    with np.errstate(invalid="ignore"):
        maxdiff = float(np.nanmax(np.abs(a - b))) if len(a) else 0.0
    status = "OK " if n_bad == 0 else "BAD"
    print("  %s %-28s n=%4d  n_mismatch=%4d  max|diff|=%.3e" % (status, name, len(a), n_bad, maxdiff))
    return n_bad, close


RECO_VARS = [
    # (cafana name, cafpyana column)
    ("recoE", "recoE"),
    ("Reco_class", "Reco_class"),
    ("Muon_length", "Muon_length"),
    ("Proton_length_leading", "Proton_length_leading"),
    ("Proton_kinetic_leading", "Proton_kinetic_leading"),
    ("Proton_chi2mu", "Proton_chi2mu"),
    ("Proton_chi2pro", "Proton_chi2pro"),
    ("Muon_chi2mu", "Muon_chi2mu"),
    ("Muon_chi2pro", "Muon_chi2pro"),
    ("Vertex_x", "vtx_x"),
    ("Vertex_y", "vtx_y"),
    ("Vertex_z", "vtx_z"),
    ("Muon_endx", "Muon_endx"),
    ("Muon_endy", "Muon_endy"),
    ("Muon_endz", "Muon_endz"),
    ("Proton_endx", "Proton_endx"),
    ("Proton_endy", "Proton_endy"),
    ("Proton_endz", "Proton_endz"),
    ("Transverse_angle", "Transverse_angle"),
    ("T3D_angle_mup", "T3D_angle_mup"),
    ("deltaPt", "deltaPt"),
    ("Transverse_mom_reco_mu", "Transverse_mom_reco_mu"),
    ("Transverse_mom_reco_pro", "Transverse_mom_reco_pro"),
    ("Number_protons", "Number_protons"),
    ("Barycenter_delta", "Barycenter_delta"),
    ("Muon_trackScore", "Muon_trackScore"),
    ("Proton_trackScore", "Proton_trackScore"),
    ("Nu_score", "nu_score"),
    ("CryoSel", "CryoSel"),
    ("E_nu_true", "E_nu_true"),
    ("FlashPE", "FlashPE"),
]

CUTCOLS = ["cut_sanity", "cut_fv", "cut_crtveto", "cut_cryo", "cut_contained",
           "cut_muon", "cut_np", "cut_0pi", "cut_0shwother", "maxcut",
           "n_proton", "n_pion", "n_shower", "n_other",
           "Muon_chi2mu", "Muon_chi2pro"]


def rasters_for_nu(slc_group):
    """Replicate the CAFANA Eff_raster loop over the slices matched to one nu."""
    maxcut = 1
    angle = -9.0
    nuscore = -9.0
    found_best = False
    best_vtx_x = 0.0
    for _, row in slc_group.iterrows():
        if row.maxcut > maxcut:
            maxcut = row.maxcut
            if row.cut_muon:
                nuscore = row.nu_score
                best_vtx_x = row.vtx_x
                found_best = True
                if not np.isnan(row.T3D_angle_mup):
                    angle = row.T3D_angle_mup
    return maxcut, angle, nuscore, found_best, best_vtx_x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cafana")
    ap.add_argument("df")
    ap.add_argument("--rtol", type=float, default=1e-4)
    ap.add_argument("--atol", type=float, default=1e-4)
    args = ap.parse_args()

    fr = uproot.open(args.cafana)
    creco = fr["events/selectedReco"].arrays(library="pd")
    cnu = fr["events/selectedNu"].arrays(library="pd")

    evt = load_df(args.df, "evt")
    mcnu = load_df(args.df, "mcnu")
    hdr = load_df(args.df, "hdr")

    evt = evt_key_frame(evt, hdr)
    if str(evt.pid_mode.iloc[0]) != "cafana":
        print("WARNING: df built with pid_mode=%s; validation expects 'cafana'" % evt.pid_mode.iloc[0])

    nbad_total = 0

    # =================================================================
    # selectedReco: event/slice set comparison
    # =================================================================
    print("=" * 70)
    print("selectedReco comparison")
    print("=" * 70)
    csel = creco.copy()
    csel["key"] = list(zip(csel.Run.astype(int), csel.Subrun.astype(int),
                           csel.Evt.astype(int), csel.Slice_index.astype(int)))
    psel = evt[evt.maple_sel].copy()
    psel["key"] = list(zip(psel.run.astype(int), psel.subrun.astype(int),
                           psel.evt.astype(int), psel.slice_index.astype(int)))

    ckeys = set(csel.key)
    pkeys = set(psel.key)
    print("CAFANA selected slices:   %d" % len(ckeys))
    print("cafpyana selected slices: %d" % len(pkeys))
    only_c = ckeys - pkeys
    only_p = pkeys - ckeys
    print("selected by CAFANA only:   %d" % len(only_c))
    print("selected by cafpyana only: %d" % len(only_p))
    nbad_total += len(only_c) + len(only_p)

    evt["key"] = list(zip(evt.run.astype(int), evt.subrun.astype(int),
                          evt.evt.astype(int), evt.slice_index.astype(int)))
    for label, keys in [("CAFANA-only", sorted(only_c)), ("cafpyana-only", sorted(only_p))]:
        for k in keys:
            print("--- %s slice %s ---" % (label, str(k)))
            rows = evt[evt.key == k]
            if len(rows):
                print(rows[CUTCOLS].to_string())
            else:
                print("  (slice not present in cafpyana evt table)")

    # per-variable comparison on the intersection
    common = sorted(ckeys & pkeys)
    if common:
        cm = csel.set_index("key").loc[common]
        pm = psel.set_index("key").loc[common]
        print("\nper-variable comparison on %d matched slices:" % len(common))
        for cname, pname in RECO_VARS:
            if cname not in cm.columns:
                print("  SKIP %s (not in CAFANA tree)" % cname)
                continue
            nb, close = compare_var(cname, cm[cname], pm[pname], args.rtol, args.atol)
            nbad_total += nb
            if nb:
                bad = np.where(~close)[0][:5]
                for i in bad:
                    print("      e.g. key=%s cafana=%.6g cafpyana=%.6g" %
                          (common[i], cm[cname].iloc[i], pm[pname].iloc[i]))

    # =================================================================
    # selectedNu: truth-level comparison
    # =================================================================
    print("=" * 70)
    print("selectedNu comparison")
    print("=" * 70)
    h = hdr[["run", "subrun", "evt"]].copy()
    pn = mcnu[mcnu.is_1muNp_maple].join(h, on=list(h.index.names))

    # order nus within an event by mc.nu index on both sides
    cnu = cnu.copy()
    cnu["ord"] = cnu.groupby(["Run", "Subrun", "Evt"]).cumcount()
    pn = pn.sort_index()
    pn["ord"] = pn.groupby(["run", "subrun", "evt"]).cumcount()

    ckeys = set(zip(cnu.Run.astype(int), cnu.Subrun.astype(int), cnu.Evt.astype(int), cnu["ord"]))
    pkeys = set(zip(pn.run.astype(int), pn.subrun.astype(int), pn.evt.astype(int), pn["ord"]))
    print("CAFANA selectedNu entries:   %d" % len(ckeys))
    print("cafpyana 1muNp mcnu entries: %d" % len(pkeys))
    print("CAFANA only:   %d %s" % (len(ckeys - pkeys), sorted(ckeys - pkeys)[:10]))
    print("cafpyana only: %d %s" % (len(pkeys - ckeys), sorted(pkeys - ckeys)[:10]))
    nbad_total += len(ckeys - pkeys) + len(pkeys - ckeys)

    common = sorted(ckeys & pkeys)
    if common:
        cm = cnu.set_index(["Run", "Subrun", "Evt", "ord"]).loc[common]
        pm = pn.reset_index().set_index(["run", "subrun", "evt", "ord"]).loc[common]

        # Pass_cut and rasters from the evt table
        evt_eff = evt[~(evt.tmatch_eff < 0.5)]
        keycols = ["run", "subrun", "evt", "tmatch_idx"]
        grouped = dict(list(evt_eff.groupby(keycols)))

        pass_cut, angle, nuscore, cryosel = [], [], [], []
        for (r, s, e, o), row in zip(common, pm.itertuples()):
            g = grouped.get((r, s, e, float(row.ind)), None)
            if g is None or not len(g):
                pass_cut.append(1)
                angle.append(-9.0)
                nuscore.append(-9.0)
                cryosel.append(-1)
                continue
            g = g.sort_values("slice_index")
            mc, an, ns, found, bvx = rasters_for_nu(g)
            pass_cut.append(mc)
            angle.append(an)
            nuscore.append(ns)
            if not found:
                cryosel.append(-1)
            else:
                cl = g.cryo_light.iloc[0]
                if cl < 0:
                    cryosel.append(-1)
                elif bvx > 0 and cl in (1, 2):
                    cryosel.append(1)
                elif bvx < 0 and cl in (0, 2):
                    cryosel.append(0)
                else:
                    cryosel.append(-1)

        print("\nper-variable comparison on %d matched truth interactions:" % len(common))
        nb, _ = compare_var("true_Enu", cm.true_Enu, pm.nu_E, args.rtol, args.atol); nbad_total += nb
        nb, _ = compare_var("True_visible_Enu", cm.True_visible_Enu, pm.true_visible_Enu, args.rtol, args.atol); nbad_total += nb
        nb, _ = compare_var("Genie_mode", cm.Genie_mode, pm.genie_mode, args.rtol, args.atol); nbad_total += nb
        nb, _ = compare_var("True_protons", cm.True_protons, np.ones(len(common)), args.rtol, args.atol); nbad_total += nb
        nb, close = compare_var("Pass_cut", cm.Pass_cut, pass_cut, args.rtol, args.atol); nbad_total += nb
        if nb:
            for i in np.where(~close)[0][:10]:
                print("      e.g. key=%s cafana=%g cafpyana=%g" % (common[i], cm.Pass_cut.iloc[i], pass_cut[i]))
        nb, _ = compare_var("Eff_raster_angle", cm.Eff_raster_angle, angle, args.rtol, args.atol); nbad_total += nb
        nb, _ = compare_var("Eff_raster_nuscore", cm.Eff_raster_nuscore, nuscore, args.rtol, args.atol); nbad_total += nb
        nb, _ = compare_var("Eff_cryosel", cm.Eff_cryosel, cryosel, args.rtol, args.atol); nbad_total += nb

    print("=" * 70)
    print("TOTAL mismatches: %d" % nbad_total)
    print("=" * 70)
    return 0 if nbad_total == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
