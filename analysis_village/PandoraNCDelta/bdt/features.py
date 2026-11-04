"""Reconstructed features for the NCDelta BDT.

Truth columns are excluded here rather than in the training script, so there
is one place to check for leakage (CLAUDE.md rule 1: never cut on truth).
"""

import numpy as np

TRUTH_PREFIXES = ("true_", "shw_true", "trk_true", "tmatch_")

# Reconstructed but not physics: bookkeeping, or an index into something else.
EXCLUDE = {
    "shw_pfp_idx", "trk_pfp_idx", "Run", "ismc", "detector",
    "is_clear_cosmic",          # already applied in the pre-selection
    "nu_score_rank",            # rank within the event, not a slice property
    # PID labels: the argmax of score columns already in the list, so a tree
    # would only be splitting on the numeric PDG code.
    "shw_razzle_pdg", "shw_razzled_pdg", "trk_razzled_pdg",
}


def feature_columns(df):
    """Numeric reconstructed columns of `df`, truth and bookkeeping removed."""
    out = []
    for c in df.columns:
        if c in EXCLUDE or any(c.startswith(p) or p in c for p in TRUTH_PREFIXES):
            continue
        if df[c].dtype.kind not in "fiu":
            continue
        if df[c].nunique(dropna=True) < 2:      # constant, nothing to learn
            continue
        out.append(c)
    return out


def derived(df, cuts, shower_energy):
    """The handful of quantities the cut-based selection computes on the fly."""
    bp = shower_energy.shw_best_plane(df)
    nh = np.take_along_axis(
        df[["shw_nhits_I0", "shw_nhits_I1", "shw_nhits_I2"]].to_numpy(dtype=float),
        bp[:, None], axis=1).ravel()
    e = shower_energy.shw_energy_corrected(df).to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return {
            "shw_energy_corrected": e,
            "shw_nhits_bestplane": nh,
            "hits_per_GeV": nh / np.where(e > 0, e, np.nan),
            "extrap_det_path": cuts.extrap_det_path(df).to_numpy(dtype=float),
        }


def sanitize(X):
    """+-inf -> NaN, which XGBoost treats as missing and learns around.

    The dataframes keep the CAF's unphysical values on purpose, so they stay
    visible (opt0_measPE reaches 9e20 with 32 infinities, shw_dedx_best 1e13).
    Trees split on rank, so a finite outlier only lands at one end of a split
    and is tolerable; an infinity is a hard error in DMatrix. Sentinel values
    (-5, -999, -10000) are left alone -- "no measurement" is a real state and
    the model may legitimately use it.
    """
    n = int(np.isinf(X.to_numpy(dtype=float)).sum())
    if n:
        print("   %d infinite feature values -> NaN" % n)
    return X.replace([np.inf, -np.inf], np.nan)
