"""SBND Wire-Cell numu CC BDT, release numu_round3af (sbnd-wirecell-bdt,
data/numu_round3af/README_release.md), outside the Wire-Cell chain. It replaces
the shipped-weight numu_score of the dataframes.

Two levels:
  small models numu_1, numu_2, cosmict_10: one score per object (the T_tagger
    vectors numu_cc_1_*, numu_cc_2_*, cosmict_10_*), pooled per candidate
    (max, max, min), float32(tanh(scale x margin + offset)), a fixed value
    for a candidate without objects;
  final model: 74 per-candidate T_tagger scalars (all in the wc table) with
    the three small-model values; numu_score = 0.868589 x margin.
Simulation correction (part of the model, MC only, correction.json): the seven
dQ/dx inputs x 1.071 (numu_cc_1_medium_dQ_dx per segment, cosmict_{2,3,7}_
dQ_dx_{end,front} where cosmict_N_filled == 1; filled values > 0 only); a
numu_1 segment keeps numu_cc_flag_1 only if its corrected dQ/dx is still
below numu_cc_1_dQ_dx_cut, and numu_cc_flag is switched off where the
candidate lost every flag_1 segment, has no numu_cc_2 flag and
numu_cc_flag_3 != 1 (training/numu/make_corrected_tables.py --update-flag1).
Scoring follows training/numu/score_numu_release.py.

make_wc_numuobj_df (dataframe maker, key "numuobj"): the small models' object
inputs, one row per (entry, nu_index, l1_model, obj); l1_model 0 numu_1,
1 numu_2, 2 cosmict_10; the model's inputs and its flag, NaN elsewhere.
score(): the candidates' round3af values from the wc and numuobj tables
(needs xgboost; the cut is left to the analysis).
"""

import json
import os

import numpy as np
import pandas as pd

RELEASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "numu_round3af")
L1 = {   # model: (code, object lists, flag) as sbnd-wirecell-bdt spec/l1_models.py
    "numu_1": (0, ["numu_cc_1_particle_type", "numu_cc_1_length", "numu_cc_1_medium_dQ_dx", "numu_cc_1_dQ_dx_cut",
                   "numu_cc_1_direct_length", "numu_cc_1_n_daughter_tracks", "numu_cc_1_n_daughter_all"], "numu_cc_flag_1"),
    "numu_2": (1, ["numu_cc_2_length", "numu_cc_2_total_length", "numu_cc_2_n_daughter_tracks", "numu_cc_2_n_daughter_all"],
               "numu_cc_flag_2"),
    "cosmict_10": (2, ["cosmict_10_vtx_z", "cosmict_10_flag_shower", "cosmict_10_flag_dir_weak", "cosmict_10_angle_beam",
                       "cosmict_10_length"], "cosmict_flag_10"),
}
OBJ_LEVELS = ["entry", "nu_index", "l1_model", "obj"]
OBJ_COLUMNS = list(dict.fromkeys(c for _, lists, _ in L1.values() for c in lists)) + ["flag"]
COSMIC_DQDX = ["cosmict_%d_dQ_dx_%s" % (n, e) for n in (2, 3, 7) for e in ("end", "front")]
SCORE_SCALE = 0.868589   # numu_score = margin x 2 / ln(10), the shipped scale


def _empty_objs():
    return pd.DataFrame(columns=OBJ_COLUMNS, index=pd.MultiIndex.from_arrays([[]] * 4, names=OBJ_LEVELS), dtype=np.float32)


def make_wc_numuobj_df(f):
    """Never raises: a CAF it cannot read gives an empty table (with a warning)."""
    try:
        return _numuobjs(f)
    except Exception as e:
        print("make_wc_numuobj_df: skipped (%s: %s)" % (type(e).__name__, e))
        return _empty_objs()


def _numuobjs(f):
    if "recTreeWireCell" not in f:
        return _empty_objs()
    t = f["recTreeWireCell"]
    keys = set(t.keys())
    br = ["wc.T_tagger.nu_index"]
    for _, lists, flag in L1.values():
        for n in lists + [flag]:
            b = "wc.T_tagger." + n
            if b not in keys:
                raise KeyError(b)
            br += [b, b + "..idx", b + "..length"]
    a = t.arrays(br, library="np")
    rows, idx = [], []
    for i, nus in enumerate(a["wc.T_tagger.nu_index"]):
        for c, nu in enumerate(nus):
            for model, (code, lists, flag) in L1.items():
                b0 = "wc.T_tagger." + lists[0]
                if c >= len(a[b0 + "..idx"][i]):
                    continue
                n = int(a[b0 + "..length"][i][c])
                for j in range(n):
                    r = dict.fromkeys(OBJ_COLUMNS, np.nan)
                    for name in lists + [flag]:
                        b = "wc.T_tagger." + name
                        if c < len(a[b + "..idx"][i]) and j < int(a[b + "..length"][i][c]):
                            r["flag" if name == flag else name] = a[b][i][int(a[b + "..idx"][i][c]) + j]
                    rows.append(r)
                    idx.append((i, int(nu), code, j))
    if not rows:
        return _empty_objs()
    return pd.DataFrame(rows, index=pd.MultiIndex.from_tuples(idx, names=OBJ_LEVELS)).astype(np.float32)


def load_release(path=RELEASE):
    import xgboost as xgb
    man = json.load(open(os.path.join(path, "manifest.json")))
    corr = json.load(open(os.path.join(path, "correction.json")))

    def booster(fn):
        b = xgb.Booster()
        b.load_model(os.path.join(path, fn))
        return b
    l1 = {m: dict(spec, booster=booster(spec["file"])) for m, spec in man["level1"]["models"].items()}
    return dict(manifest=man, correction=corr, l1=l1, final=booster(man["models"]["d6_f0"]["file"]))


def _scale(q, corr):
    if corr["mode"] == "uniform":
        return np.full_like(q, corr["k"])
    if corr["mode"] == "linear":
        return corr["a"] + corr["b"] * np.clip(q, corr["q_min"], corr["q_max"])
    raise ValueError(corr["mode"])


def _correct(q, mask, corr):
    q = np.asarray(q, np.float64).copy()
    ok = mask & np.isfinite(q) & (q > 0)
    q[ok] = q[ok] * _scale(q[ok], corr)
    return q


def score(cands, objs, ismc, release=None, correct=True):
    """round3af values per candidate (index as cands): numu_score_3af and the small-model values
    numu_1_score_3af, numu_2_score_3af, cosmict_10_score_3af; numu_cc_flag_3af (the flag after the correction).
    cands: the wc table; objs: the numuobj table (same outer index levels + l1_model, obj); ismc: apply the
    simulation correction (correct=False skips it)."""
    import xgboost as xgb
    rel = release or load_release()
    man, corr = rel["manifest"], rel["correction"]
    inputs = man["inputs"]
    key = list(cands.index.names)
    C = cands[list(dict.fromkeys(inputs + ["numu_cc_flag_3"] + ["cosmict_%d_filled" % n for n in (2, 3, 7)]))].astype(np.float64)
    O = objs.reset_index(["l1_model", "obj"])
    sim = ismc and correct
    if sim:
        for c in COSMIC_DQDX:
            C[c] = _correct(C[c].to_numpy(), C["cosmict_%s_filled" % c.split("_")[1]].to_numpy() == 1, corr)
        o1 = (O.l1_model == L1["numu_1"][0]).to_numpy()
        q = O["numu_cc_1_medium_dQ_dx"].to_numpy(np.float64)
        q[o1] = _correct(q[o1], np.ones(o1.sum(), bool), corr)
        O["numu_cc_1_medium_dQ_dx"] = q
        # the numu_1 muon flag follows the corrected dQ/dx; numu_cc_flag recomputed where flag_1 was lost
        fl = O["flag"].to_numpy(np.float64)
        cut = O["numu_cc_1_dQ_dx_cut"].to_numpy(np.float64)
        newfl = np.where(o1 & (fl == 1) & ~(q < cut), 0.0, fl)
        seg = pd.DataFrame({"old": o1 & (fl == 1), "new": o1 & (newfl == 1)}, index=O.index).groupby(level=key).any()
        lost = seg.index[seg.old & ~seg.new]
        o2 = O[(O.l1_model == L1["numu_2"][0]) & (O.flag == 1)].index.unique()
        hit = C.index.isin(lost) & ~C.index.isin(o2) & (C.numu_cc_flag_3.to_numpy() != 1)
        C.loc[hit, "numu_cc_flag"] = 0.0
        O["flag"] = newfl
    for c in man["nan_to_0"]:
        C[c] = C[c].fillna(0.0)
    out = pd.DataFrame(index=cands.index)
    for m, spec in rel["l1"].items():
        o = O[O.l1_model == L1[m][0]]
        x = o[spec["inputs"]].copy()
        for c, (kind, v) in spec["protect"].items():
            val = x[c].to_numpy(np.float64)
            x[c] = np.where(np.isinf(val) if kind == "inf" else np.isnan(val), v, val)
        if len(o):
            p = spec["booster"].predict(xgb.DMatrix(x.to_numpy(np.float32), feature_names=spec["inputs"]), output_margin=True)
            mg = pd.Series(p, index=o.index).groupby(level=key).agg(spec["pool"]).reindex(C.index).to_numpy(np.float64)
        else:
            mg = np.full(len(C), np.nan)
        v = np.where(np.isnan(mg), np.float32(spec["empty_value"]),
                     np.tanh(spec["scale"] * mg + spec["offset"]).astype(np.float32)).astype(np.float32)
        C[m + "_score"] = v
        out[m + "_score_3af"] = v
    X = xgb.DMatrix(C[inputs].to_numpy(np.float32), feature_names=inputs)
    out["numu_score_3af"] = SCORE_SCALE * rel["final"].predict(X, output_margin=True)
    out["numu_cc_flag_3af"] = C["numu_cc_flag"].to_numpy()
    return out
