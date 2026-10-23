import os

import numpy as np
import pandas as pd

########## VERSION HISTORY #############
#   v0: Initial Launch of parameters
#
########################################
llpid_version = "v0" 
llpid_tables_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "analysis_village", "ll_pid", "tables")


#### Default
#### == every number the method uses; copy it, change what you want and pass it as params=
LLPID_PARAMS = {
    # --- which hits count
    "rr_min": 1.0,
    "dedx_min": 0.5,            # MeV/cm, below this the hit is noise
    "dedx_max": 16.,            # MeV/cm, above this a delta ray...
    "dedx_max_free_rr": 0.,     # ...unless within this many cm of either track end
    "pitch_max": None,            # cm, None disables
    "skip_first": 1,            # usable hits dropped at the low-rr end
    "skip_last": 1,             # and at the high-rr end
    "min_hits": 9,              # fewer usable hits (all planes): no dL fit

    # --- the dL (TLE) scan
    "scan_coarse": 5.,          # cm
    "scan_fine": 1.,            # cm, +-fine_steps*scan_coarse around each track's coarse minimum
    "fine_steps": 4,
    "scan_max": 600.,           # cm
    "flat_abs": 10.,            # -2lnL a well edge must rise above its floor, else dL = 0. FIXME: tuned on per-hit lambda
    "allow_flip": True,         # also fit the track reversed
    "flip_margin": 25.,         # -2lnL the reversed fit must win by. FIXME: tuned on per-hit lambda
    "orient_shared": False,     # one flip decision per track instead of per hypothesis
    "prior_file": "dl_prior_kde_truthlen.npz",


    # --- the per-hit dE/dx pdf
    "tables": "mc_bysp@mc_pdf_binned_n2_end10.npz",
    "tables_bare": "mc_bysp@mc_bare_pdf_binned_n2_end10.npz",
}

# variations
# calorimetry variations go through dedxname instead: dedx from chi2pid.dedx(..., new_calo_params=chi2pid.CALO_VARIATIONS[...])
LLPID_VARIATIONS = {
    "cv": LLPID_PARAMS,
}

MASS = {"mu": 105.6583755, "pi": 139.57039, "p": 938.27208816} # MeV
HYPS = ("mu", "pi", "p")
PAR_TO_HYP = {"muon": "mu", "pion": "pi", "proton": "p"}


# Bethe-Bloch / Landau constants, liquid argon
K_COEF = 0.307075   # MeV cm^2 / mol
Z_AR, A_AR = 18.0, 39.9484
RHO_AR = 1.38434       # g/cm^3
I_AR = 0.000188     # MeV
M_E = 0.511         # MeV
EULER = 0.5772156649015329


def likelihood(hitdfs, dedxname="dedx", rrcorr=True, params=None, hyps=HYPS):
    params = LLPID_PARAMS if params is None else params
    frames = [hitdfs] if isinstance(hitdfs, pd.DataFrame) else list(hitdfs)

    tracks = _track_index(frames)
    cols = _columns(rrcorr, hyps)
    if len(tracks) == 0:
        return pd.DataFrame({c: pd.Series(dtype=float) for c in cols}).set_axis(tracks)

    planes = [_prepare_plane(df, tracks, _plane_of(df, i, len(frames)), dedxname, params)
              for i, df in enumerate(frames)]
    n_hits = sum(h["n_usable"] for h in planes)

    out = {}

    tables_bare = _load_tables(params["tables_bare"])
    ell0 = {hyp: _ell_at_zero(planes, tables_bare, hyp) for hyp in hyps}

    if rrcorr:
        tables = _load_tables(params["tables"])
        priors, dl_min = _load_priors(params["prior_file"])
        fits = _fit_dl(planes, tables, priors, dl_min, n_hits, params, hyps)

    for hyp in hyps:
        out["ell0_" + hyp] = ell0[hyp]
    for hyp in hyps:
        out["nhits_" + hyp] = n_hits

    if rrcorr:
        for hyp in hyps:
            out["ell_" + hyp] = fits[hyp]["cost"]
        for name in ("dl", "flip", "fitok", "flattened", "onedge"):
            for hyp in hyps:
                out[name + "_" + hyp] = fits[hyp][name]


    out["version"] = llpid_version
    return pd.DataFrame(out, index=tracks)[cols]


def rr_corr(hitdf, lldf, par, rr_floor=0.3):
    """Per-hit residual range after the dL fit, for one plane's hits and one hypothesis."""
    hyp = PAR_TO_HYP.get(par, par)
    trk = hitdf.index.droplevel(-1)
    dl = np.nan_to_num(lldf["dl_" + hyp].reindex(trk).to_numpy(float), nan=0.0)
    flip = lldf["flip_" + hyp].reindex(trk).to_numpy()
    flip = np.where(pd.isna(flip), False, flip).astype(bool)

    rr = hitdf.rr.to_numpy(float)
    rr_max = hitdf.rr.groupby(level=list(range(hitdf.index.nlevels-1))).transform("max").to_numpy(float)
    # a flipped track measures its range from the other end
    return pd.Series(np.maximum(np.where(flip, rr_max - rr, rr) + dl, rr_floor), index=hitdf.index, name="rr_corr_" + hyp)

def _columns(rrcorr, hyps):
    cols = ["ell0_" + h for h in hyps]
    cols += ["nhits_" + h for h in hyps]
    if rrcorr:
        cols += ["ell_" + h for h in hyps]
        cols += [n + "_" + h for n in ("dl", "flip", "fitok", "flattened", "onedge") for h in hyps]

    cols += ["version"]
    return cols



def _pdf_params(KE, mass, pitch):
    """Returns the parameters read in by the hit PDF: xi_per_cm, kappa, beta2, mpv_shift"""
    with np.errstate(invalid="ignore", divide="ignore"):
        gamma = KE / mass + 1.0
        beta2 = 1.0 - 1.0 / (gamma * gamma)
        beta = np.sqrt(beta2)
        bg2 = beta2 * gamma * gamma

        w_max = 2.0 * M_E * bg2 / (1.0 + 2.0 * M_E * gamma / mass + (M_E / mass) ** 2)

        xi = 0.5 * RHO_AR * pitch * Z_AR * K_COEF / (A_AR * beta2)
        kappa = np.clip(xi / w_max, 1e-12, None)
        xi_per_cm = xi / pitch

        # density effect (Sternheimer, argon)
        y = np.log10(np.maximum(beta * gamma, 1e-30))
        delta = 2.0 * np.log(10.0) * y - 5.2146 + 0.19559 * np.clip(3.0 - y, 0, None) ** 3
        delta = np.where(y < 0.2, 0.0, delta)

        mean_dedx = (RHO_AR * K_COEF * (Z_AR / A_AR) / beta2) * (
            0.5 * np.log(2.0 * M_E * bg2 * w_max / I_AR ** 2) - beta2 - delta / 2.0)

        mpv_shift = (1.0 - EULER + beta2 + np.log(kappa)) * xi_per_cm + mean_dedx

    return xi_per_cm, kappa, beta2, mpv_shift



def _hit_pdf(dedx, pitch, KE, mass, table):
    "Returns the hit likelihood pdf(dE/dx | KE, mass, pitch) for one plane's hits under one hypothesis."
    xi_per_cm, kappa, beta2, mpv_shift = _pdf_params(KE, mass, pitch)
    lam = np.clip((dedx - mpv_shift) / xi_per_cm, -50.0, 100.0)

    shape = np.broadcast(lam, kappa, beta2).shape
    lam = np.broadcast_to(lam, shape).ravel()
    kap = np.broadcast_to(kappa, shape).ravel()
    bet = np.broadcast_to(beta2, shape).ravel()
    jac = np.broadcast_to(xi_per_cm, shape).ravel()

    return np.clip(_lookup_binned(lam, kap, bet, table) / jac, 1e-12, None).reshape(shape)


def _lookup_binned(lam, kappa, beta2, table):
    """pdf of the bin each hit falls in; values outside the edges take the first/last bin"""
    idx = tuple(np.clip(np.searchsorted(e, x, side="right") - 1, 0, len(e) - 2)
                for e, x in zip(table["edges"], (lam, kappa, beta2)))
    return table["vav_pdf"][idx]



##############################
# RESIDUAL RANGE -> KINETIC ENERGY
##############################
# the only place the mass hypothesis enters

_KE_CACHE = {}

def _ke_table(hyp, ke_max=4000.0, step=0.05):
    if hyp in _KE_CACHE:
        return _KE_CACHE[hyp]
    mass = MASS[hyp]
    ke = np.arange(step, ke_max, step)
    gamma = ke / mass + 1.0
    beta2 = 1.0 - 1.0 / (gamma * gamma)
    beta = np.sqrt(beta2)
    bg2 = beta2 * gamma * gamma
    w_max = 2.0 * M_E * bg2 / (1.0 + 2.0 * M_E * gamma / mass + (M_E / mass) ** 2)
    y = np.log10(beta * gamma)
    delta = 2.0 * np.log(10.0) * y - 5.2146 + 0.19559 * np.clip(3.0 - y, 0, None) ** 3
    delta = np.where(y < 0.2, 0.0, delta)
    dedx = (RHO_AR * K_COEF * (Z_AR / A_AR) / beta2) * (
        0.5 * np.log(2.0 * M_E * bg2 * w_max / I_AR ** 2) - beta2 - delta / 2.0)
    _KE_CACHE[hyp] = {"range": np.cumsum(step / dedx), "ke": ke}
    return _KE_CACHE[hyp]

def _ke_from_rr(table, rr):
    return np.interp(np.clip(rr, 1e-3, None), table["range"], table["ke"])

def _rr_from_ke(table, ke):
    return np.interp(ke, table["ke"], table["range"])

##############################
# HITS -> TRACKS
##############################
# each plane's hits are sorted by (track, rr); every "per track" quantity is a reduceat
# over that plane's non-empty tracks, scattered into the union of tracks of all planes,
# so a track missing a plane just gets nothing from it

def _track_index(frames):
    tracks = None
    for df in frames:
        if len(df) == 0:
            continue
        t = df.index.droplevel(-1).unique()
        tracks = t if tracks is None else tracks.union(t)
    if tracks is None:
        names = frames[0].index.names[:-1] if frames else []
        return pd.MultiIndex.from_tuples([], names=names)
    return tracks.sort_values()

def _plane_of(df, position, n_frames):
    if "plane" in df.columns and len(df):
        pl = df["plane"].to_numpy()
        if pl.dtype.kind == "f":
            pl = pl[np.isfinite(pl)]
        if pl.size:
            return int(np.bincount(pl.astype(np.intp)).argmax())
    if n_frames > 1:
        return position
    raise ValueError("a single hit frame needs a 'plane' column (make_trkhitdf adds one)")

def _prepare_plane(df, tracks, plane, dedxname, params):
    rr = df.rr.to_numpy(float)
    keep = rr > params["rr_min"]
    tid = tracks.get_indexer(df.index.droplevel(-1))[keep]
    rr = rr[keep]
    dedx = df[dedxname].to_numpy(float)[keep]
    pitch = df.pitch.to_numpy(float)[keep]

    order = np.lexsort((rr, tid))
    counts = np.bincount(tid, minlength=len(tracks))
    first = np.r_[0, np.cumsum(counts)[:-1]]
    hits = {
        "tid": tid[order], "rr": rr[order], "dedx": dedx[order], "pitch": pitch[order],
        "plane": plane, "n_tracks": len(tracks), "first": first,
        "nonempty": np.flatnonzero(counts), "starts": first[counts > 0],
    }
    hits["rr_max"] = _per_track_max(hits, hits["rr"])[hits["tid"]]
    hits["usable"] = _usable_hits(hits, params)
    hits["n_usable"] = _per_track_sum(hits, hits["usable"])
    return hits

def _per_track_sum(hits, values):
    out = np.zeros(hits["n_tracks"])
    if hits["starts"].size:
        out[hits["nonempty"]] = np.add.reduceat(np.asarray(values, float), hits["starts"])
    return out

def _per_track_max(hits, values):
    out = np.full(hits["n_tracks"], np.nan)
    if hits["starts"].size:
        out[hits["nonempty"]] = np.maximum.reduceat(np.asarray(values, float), hits["starts"])
    return out

def _rank_within_track(hits, mask):
    c = np.cumsum(mask) - mask
    return c - c[hits["first"][hits["tid"]]]

def _usable_hits(hits, params):
    rr, dedx, pitch = hits["rr"], hits["dedx"], hits["pitch"]
    finite = np.isfinite(dedx) & np.isfinite(pitch) & np.isfinite(rr)

    free = params["dedx_max_free_rr"]
    near_an_end = (rr < free) | (rr > hits["rr_max"] - free)
    in_window = (dedx > params["dedx_min"]) & ((dedx < params["dedx_max"]) | near_an_end)
    keep = finite & in_window

    pmax = params.get("pitch_max")
    if pmax is not None and np.isfinite(pmax):
        keep &= pitch <= pmax

    rank = _rank_within_track(hits, keep)
    n_keep = _per_track_sum(hits, keep)[hits["tid"]]
    keep &= (rank >= params["skip_first"]) & (rank < n_keep - params["skip_last"])
    return keep

def _track_ell(hits, tables, hyp, offset=0.0, flipped=None):

    """ Track likelihood with possibilities for residual range correction:
    Returns (usable hits, summed -2lnL) per track on one plane, under one hypothesis."""
    rr = hits["rr"]
    if flipped is not None:
        rr = np.maximum(np.where(flipped[hits["tid"]], hits["rr_max"] - rr, rr), 0.0)
    if np.ndim(offset) > 0:
        offset = np.asarray(offset)[hits["tid"]]

    KE = _ke_from_rr(_ke_table(hyp), rr + offset)
    pdf = _hit_pdf(hits["dedx"], hits["pitch"], KE, MASS[hyp], _table_for(tables, hits["plane"], hyp))

    ln_pdf = np.log(np.maximum(pdf, 1e-300))
    n = hits["n_usable"]
    raw = -2.0 * _per_track_sum(hits, np.where(hits["usable"], ln_pdf, 0.0))
    return n, np.where(n > 0, raw, np.nan)

def _ell_at_zero(planes, tables, hyp):
    """Track likelihood without residual range correction"""
    total, n_total = 0.0, 0
    for hits in planes:
        n, raw = _track_ell(hits, tables, hyp)
        total = total + np.nan_to_num(raw)
        n_total = n_total + n
    return np.where(n_total > 0, total, np.nan)

##############################
# THE dL (TLE) FIT
##############################
# dL = range the particle still had past its last hit.  Per hypothesis and per
# direction: a coarse scan, a fine scan around its minimum, then
#   flattened  the well's edges never rise flat_abs above its floor -> dE/dx does not
#              determine dL, keep dL = 0
#   on edge    the minimum sits at scan_max -> runaway, keep dL = 0

def _fit_dl(planes, tables, priors, dl_min, n_hits, params, hyps):
    """Choose the optimal residual range correction for each track, under each hypothesis and direction"""
    n_tracks = n_hits.size
    enough = n_hits >= params["min_hits"]
    directions = (False, True) if params["allow_flip"] else (False,)

    step = float(params["scan_coarse"])
    reach = float(params["fine_steps"]) * step
    fine = np.arange(-reach, reach + params["scan_fine"], params["scan_fine"])
    flat_abs = float(params["flat_abs"])
    cols = np.arange(n_tracks)

    best = {hyp: {} for hyp in hyps}
    for hyp in hyps:
        coarse = np.arange(dl_min[hyp], params["scan_max"] + step, step)
        for flipped in directions:
            fine_ell, fine_dl, zero_ell = _scan(planes, tables, hyp, coarse, fine, priors[hyp], n_tracks, flipped)

            edge_lo, edge_hi = _well_edges(fine_ell)
            ell_min = np.nanmin(np.where(np.isnan(fine_ell), np.inf, fine_ell), axis=0)
            flattened = ((edge_lo - ell_min) < flat_abs) | ((edge_hi - ell_min) < flat_abs)

            total = np.where((fine_dl < dl_min[hyp]) | (fine_dl > params["scan_max"]), np.inf, fine_ell)
            k = np.argmin(total, axis=0)
            best[hyp][flipped] = {"dl": np.where(flattened, 0.0, fine_dl[k, cols]),
                                  "cost": np.where(flattened, zero_ell, total[k, cols]),
                                  "zero_ell": zero_ell, "flattened": flattened}

    fits = {}
    for hyp in hyps:
        flipped = _orientation(best[hyp], params) if params["allow_flip"] else np.zeros(n_tracks, bool)
        fits[hyp] = _apply_orientation(best[hyp], flipped, enough, params)

    if params["allow_flip"] and params["orient_shared"]:
        fwd = np.nanmin([best[hyp][False]["cost"] for hyp in hyps], axis=0)
        rev = np.nanmin([best[hyp][True]["cost"] for hyp in hyps], axis=0)
        shared = ((fwd - rev) > params["flip_margin"]) & np.isfinite(fwd) & np.isfinite(rev)
        for hyp in hyps:
            fits[hyp] = _apply_orientation(best[hyp], shared, enough, params)
    return fits

def _scan(planes, tables, hyp, coarse, fine, prior, n_tracks, flipped):
    """Scans through residual range offsets and returns the track likelihood distribution for each track as a function of residual range corrected"""
    flip = np.full(n_tracks, True) if flipped else None

    def surface(grid, base):
        ell = np.full((len(grid), n_tracks), np.nan)
        for g, step in enumerate(grid):
            off = np.full(n_tracks, step, float) if base is None else base + step
            tot = np.zeros(n_tracks)
            cnt = np.zeros(n_tracks)
            for hits in planes:
                n_p, raw_p = _track_ell(hits, tables, hyp, off, flip)
                ok = n_p > 0
                tot[ok] += raw_p[ok]
                cnt[ok] += n_p[ok]
            good = cnt > 0
            ell[g, good] = tot[good]
        return ell

    # the uncorrected ell: dL = 0, no prior term (zero there by construction)
    zero_ell = surface(np.array([0.0]), None)[0]

    # Use a two step scan to allow for fine scan without taking too much time; coarse scan finds approximate minimum, and fine scan searches around the minimum. 
    coarse_ell = surface(coarse, None)
    coarse_ell = coarse_ell + _prior_charge(prior, np.broadcast_to(coarse[:, None], coarse_ell.shape))
    centre = coarse[np.argmin(np.where(np.isnan(coarse_ell), np.inf, coarse_ell), axis=0)]

    fine_dl = centre[None, :] + fine[:, None]
    fine_ell = surface(fine, centre) + _prior_charge(prior, fine_dl)
    return fine_ell, fine_dl, zero_ell

def _well_edges(ell, n_edge=5):
    G, T = ell.shape
    y = np.where(np.isnan(ell), np.inf, ell)
    edge_lo = np.full(T, np.nan)
    edge_hi = np.full(T, np.nan)
    for t in range(T):
        idx = np.flatnonzero(np.isfinite(y[:, t]))
        if idx.size < 2 * n_edge:
            continue
        edge_lo[t] = y[idx[:n_edge], t].mean()
        edge_hi[t] = y[idx[-n_edge:], t].mean()
    return edge_lo, edge_hi

def _orientation(best, params):
    """Check what orientation is preferred for each track, based on the best -2lnL(rr+dL) for each orientation.  
    Returns a boolean array: True = flipped, False = forward."""
    if True not in best:
        return np.zeros(np.shape(best[False]["cost"]), bool)
    f, r = best[False]["cost"], best[True]["cost"]
    return ((f - r) > params["flip_margin"]) & np.isfinite(f) & np.isfinite(r)

def _apply_orientation(best, flipped, enough, params):
    """Given the best choice direction for the track, returns the final offset, likelihood, and flags for each of the tracks.
    If poor fits, returns default values (dL=0, cost=zero_ell, fitok=False)."""
    if True in best:
        def choose(name):
            return np.where(flipped, best[True][name], best[False][name])
    else:
        def choose(name):
            return best[False][name]

    dl, flattened = choose("dl"), choose("flattened")
    on_edge = dl >= params["scan_max"] - params["scan_coarse"]
    keep = enough & ~flattened & ~on_edge & np.isfinite(dl)
    return {"dl": np.where(keep, dl, 0.0),
            "cost": np.where(keep, choose("cost"), choose("zero_ell")),
            "fitok": keep, "flattened": flattened, "flip": flipped, "onedge": on_edge}

##############################
# THE dL PRIOR
##############################
# -2 ln P(dL) from the KDE table made in table_maker.ipynb, applied to both directions.
def _prior_charge(prior, dl):
    x, nlp = prior
    return np.interp(dl, x, nlp)



##############################
# LIKELIHOOD TABLES
##############################
# Files in llpid_tables/ next to this module (committed, so grid jobs get them from git archive).
# params["tables"]:
#   "mc@<file>"        MC-measured pdf(lambda | kappa, beta2) per plane, species pooled
#   "mc_bysp@<file>"   the same per plane AND per hypothesis
#   "mc_<hyp>@<file>"  one species' template for every hypothesis
# Loaded on first use and cached per process, not at import like chi2pid's templates: they are ~40-90 MB.

_TABLE_CACHE = {}
_PRIOR_CACHE = {}

def _table_path(fname):
    path = fname if os.path.isabs(fname) else os.path.join(llpid_tables_dir, fname)
    if not os.path.exists(path):
        raise FileNotFoundError(f"likelihood PID table {path} not found")
    return path


def _load_tables(spec):
    if spec in _TABLE_CACHE:
        return _TABLE_CACHE[spec]

    which, _, fname = spec.partition("@")
    which = which[2:].lstrip("_") or "all"
    z = np.load(_table_path(fname))
    edges = (z["lam_edges"], z["kap_edges"], z["beta2_edges"])
    def template(tag):
        return {"edges": edges, "vav_pdf": np.ascontiguousarray(z["pdf_" + tag])}

    planes = [int(p) for p in z["planes"]]
    if which == "bysp":
        tables = {pl: {hyp: template(f"{pl}_{hyp}") for hyp in HYPS} for pl in planes}
    else:
        tables = {pl: template(f"{pl}_{which}") for pl in planes}

    _TABLE_CACHE[spec] = tables
    return tables


def _table_for(tables, plane, hyp):
    t = tables
    if "vav_pdf" not in t:
        if plane not in t:
            raise KeyError(f"no likelihood table for plane {plane}; have {sorted(t)}")
        t = t[plane]
    if "vav_pdf" not in t:
        t = t[hyp]
    return t

def _load_priors(fname):
    if fname in _PRIOR_CACHE:
        return _PRIOR_CACHE[fname]
    z = np.load(_table_path(fname))
    priors = {hyp: (z["x"], z[hyp]) for hyp in HYPS}
    dl_min = {hyp: float(z["x"][0]) for hyp in HYPS}
    _PRIOR_CACHE[fname] = (priors, dl_min)
    return priors, dl_min
