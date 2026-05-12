#!/usr/bin/env python3
"""
Read SBND voxel flux histogram ROOT files, merge jobs with POT weights from
``/run/beamOn`` in ``production*.in``, and write ``full_voxel_flux_potmacro.pkl``.
"""

from __future__ import annotations

import argparse
import pickle
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import uproot
from tqdm.auto import tqdm

# ---------------------------------------------------------------------------
# Job direcotories
# ---------------------------------------------------------------------------
DEFAULT_BASE_DIR = Path("/pnfs/sbnd/scratch/users/munjung/beammc_20cm")
DEFAULT_PRODUCTION = "production_BooNE_50m_I174000A"
DEFAULT_CACHE_DIR = Path("/exp/sbnd/data/users/munjung/flux/SBND_dk2nu")
DEFAULT_VOXEL_SIZE = 20
DEFAULT_N_JOBS = 100

# ---------------------------------------------------------------------------
# Nu Flavors
# ---------------------------------------------------------------------------
FLAVORS = ("nue", "nuebar", "numu", "numubar")
FLAVOR_HIST = {f: f"h50{i + 1}" for i, f in enumerate(FLAVORS)}

_BEAMON_RE = re.compile(r"/run/beamOn\s+(\d+)")


# ---------------------------------------------------------------------------
# Binnings
# ---------------------------------------------------------------------------
def edges_z(voxel_size: int = DEFAULT_VOXEL_SIZE) -> np.ndarray:
    return np.arange(11000, 11500 + voxel_size, voxel_size)

def edges_y(voxel_size: int = DEFAULT_VOXEL_SIZE) -> np.ndarray:
    # jobs include (200, 220) -- exclude this
    return np.arange(-200, 200 + voxel_size, voxel_size)

def edges_x(voxel_size: int = DEFAULT_VOXEL_SIZE) -> np.ndarray:
    return np.arange(-200, 200 + voxel_size, voxel_size)

# ---------------------------------------------------------------------------
# File name configs
# ---------------------------------------------------------------------------

def hist_filename(
    z_lo: int,
    y_lo: int,
    job_id: int,
    production: str,
    voxel_size: int,
) -> str:
    z_hi = z_lo + voxel_size
    y_hi = y_lo + voxel_size
    return (
        f"hist_sbnd_NuBeam_z{z_lo}-{z_hi}_y{y_lo}-{y_hi}"
        f"_{production}_{job_id}_voxelsize{voxel_size}cm.root"
    )


_timestamp_dir_cache: dict[int, Path | None] = {}

def timestamp_dir_for_z(z_lo: int, base_dir: Path, production: str) -> Path | None:
    if z_lo in _timestamp_dir_cache:
        return _timestamp_dir_cache[z_lo]
    parent = base_dir / str(z_lo) / production
    if not parent.exists():
        _timestamp_dir_cache[z_lo] = None
        return None
    candidates = [d for d in parent.iterdir() if d.is_dir()]
    if not candidates:
        _timestamp_dir_cache[z_lo] = None
        return None
    chosen = sorted(candidates, key=lambda p: p.name)[-1]
    _timestamp_dir_cache[z_lo] = chosen
    return chosen


def hist_path(
    z_lo: int,
    y_lo: int,
    job_id: int,
    *,
    base_dir: Path,
    production: str,
    voxel_size: int,
) -> Path | None:
    ts = timestamp_dir_for_z(z_lo, base_dir, production)
    if ts is None:
        return None
    return ts / str(job_id) / hist_filename(z_lo, y_lo, job_id, production, voxel_size)


def total_pot_from_production_macro(job_dir: Path) -> float:
    if not job_dir.is_dir():
        return 0.0
    for macro in sorted(job_dir.glob("production*.in")):
        try:
            text = macro.read_text(errors="replace")
        except OSError:
            continue
        m = _BEAMON_RE.search(text)
        if m:
            return float(m.group(1))
    return 0.0


def check_files(
    *,
    base_dir: Path,
    production: str,
    voxel_size: int,
    z_lo_arr: np.ndarray,
    y_lo_arr: np.ndarray,
    n_jobs: int,
    ny: int,
    verbose: bool = True,
):
    """Verify (z, y, job) hist file coverage."""
    missing: dict[tuple[int, int], list[int]] = defaultdict(list)
    missing_z: list[int] = []
    n_expected = 0
    n_present = 0
    for z_lo in tqdm(z_lo_arr, desc="checking z"):
        ts = timestamp_dir_for_z(z_lo, base_dir, production)
        if ts is None:
            missing_z.append(int(z_lo))
            n_expected += ny * n_jobs
            continue
        for y_lo in y_lo_arr:
            for job_id in range(n_jobs):
                n_expected += 1
                fname = ts / str(job_id) / hist_filename(
                    int(z_lo), int(y_lo), job_id, production, voxel_size
                )
                if fname.exists():
                    n_present += 1
                else:
                    missing[(int(z_lo), int(y_lo))].append(job_id)
    if verbose:
        print(
            f"Expected {n_expected} hist files, found {n_present} "
            f"({n_expected - n_present} missing)"
        )
        if missing_z:
            print(f"Missing z directories entirely: {missing_z}")
        if missing:
            print(f"Missing files in {len(missing)} (z, y) pairs:")
            for (z_lo, y_lo), jobs in list(missing.items())[:20]:
                tail = "..." if len(jobs) > 10 else ""
                print(
                    f"  z={z_lo} y={y_lo} -> missing {len(jobs)} job(s): "
                    f"{jobs[:10]}{tail}"
                )
    return missing, missing_z


def cache_path_zy(z_lo: int, y_lo: int, cache_dir: Path) -> Path:
    return cache_dir / f"yz_z{int(z_lo)}_y{int(y_lo)}_potmacro.pkl"


def load_zy_bin(
    z_lo: int,
    y_lo: int,
    *,
    base_dir: Path,
    production: str,
    voxel_size: int,
    cache_dir: Path,
    nx: int,
    n_jobs: int,
    force: bool = False,
    verbose: bool = False,
):
    cache = cache_path_zy(z_lo, y_lo, cache_dir)
    if cache.exists() and not force:
        with open(cache, "rb") as f:
            blob = pickle.load(f)
        if isinstance(blob, tuple) and len(blob) == 4:
            return blob
        cache.unlink(missing_ok=True)

    weighted = None
    bin_edges = None
    total_pot = 0.0
    n_merged = 0

    for job_id in range(n_jobs):
        fpath = hist_path(
            z_lo, y_lo, job_id, base_dir=base_dir, production=production, voxel_size=voxel_size
        )
        if fpath is None or not fpath.exists():
            continue
        pot = total_pot_from_production_macro(fpath.parent)
        if pot <= 0.0:
            if verbose:
                print(f"No POT from production*.in next to {fpath}")
            continue

        try:
            with uproot.open(fpath) as f:
                job_ok = False
                for ix in range(nx):
                    for fi, flav in enumerate(FLAVORS):
                        key = (
                            f"vox_ix{ix:04d}_iy0000_iz0000/"
                            f"{FLAVOR_HIST[flav]}_vox{ix:06d}"
                        )
                        try:
                            h = f[key]
                        except (KeyError, uproot.exceptions.KeyInFileError):
                            continue
                        vals = np.asarray(h.values(), dtype=np.float64)
                        if weighted is None:
                            bin_edges = h.axis().edges()
                            weighted = np.zeros(
                                (nx, len(FLAVORS), len(vals)), dtype=np.float64
                            )
                        weighted[ix, fi] += vals * pot
                        job_ok = True
        except Exception as exc:
            print(f"  warning: could not open {fpath}: {exc}")
            continue

        if job_ok:
            total_pot += pot
            n_merged += 1

    if weighted is None or total_pot <= 0.0:
        return None, None, 0, 0.0

    contents = weighted / total_pot
    cache_dir.mkdir(parents=True, exist_ok=True)
    with open(cache, "wb") as f:
        pickle.dump((contents, bin_edges, n_merged, total_pot), f)
    return contents, bin_edges, n_merged, total_pot


def build_full_voxel_array(
    *,
    base_dir: Path,
    production: str,
    voxel_size: int,
    cache_dir: Path,
    z_edges: np.ndarray,
    y_edges: np.ndarray,
    x_edges: np.ndarray,
    n_jobs: int,
    force: bool = False,
    verbose: bool = False,
) -> dict:
    full_cache = cache_dir / "full_voxel_flux_potmacro.pkl"
    if full_cache.exists() and not force:
        with open(full_cache, "rb") as f:
            cached = pickle.load(f)
        if isinstance(cached, dict) and cached.get("pot_per_zy") is not None:
            return cached
        full_cache.unlink(missing_ok=True)

    nz = len(z_edges) - 1
    ny = len(y_edges) - 1
    nx = len(x_edges) - 1
    z_lo_arr = z_edges[:-1]
    y_lo_arr = y_edges[:-1]

    flux = None
    n_jobs_per_zy = np.zeros((nz, ny), dtype=np.int32)
    pot_per_zy = np.zeros((nz, ny), dtype=np.float64)
    bin_edges = None

    for iz, z_lo in enumerate(tqdm(z_lo_arr, desc="z slices")):
        for iy, y_lo in enumerate(y_lo_arr):
            contents, edges, n_merged, total_pot = load_zy_bin(
                int(z_lo),
                int(y_lo),
                base_dir=base_dir,
                production=production,
                voxel_size=voxel_size,
                cache_dir=cache_dir,
                nx=nx,
                n_jobs=n_jobs,
                force=force,
                verbose=verbose,
            )
            if contents is None:
                continue
            if flux is None:
                bin_edges = edges
                flux = np.zeros(
                    (nz, ny, nx, len(FLAVORS), len(bin_edges) - 1), dtype=np.float64
                )
            flux[iz, iy] = contents
            n_jobs_per_zy[iz, iy] = n_merged
            pot_per_zy[iz, iy] = total_pot

    if flux is None or bin_edges is None:
        raise RuntimeError("No voxel histograms were loaded; check --base-dir and inputs.")

    voxel_area_m2 = (voxel_size * 1e-2) ** 2
    out = {
        "flux": flux,
        "bin_edges": bin_edges,
        "n_jobs": n_jobs_per_zy,
        "pot_per_zy": pot_per_zy,
        "x_edges": x_edges,
        "y_edges": y_edges,
        "z_edges": z_edges,
        "flavors": list(FLAVORS),
        "voxel_size": int(voxel_size),
        "voxel_area_m2": float(voxel_area_m2),
    }
    cache_dir.mkdir(parents=True, exist_ok=True)
    with open(full_cache, "wb") as f:
        pickle.dump(out, f)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--base-dir",
        type=Path,
        default=DEFAULT_BASE_DIR,
        help="Top beammc_20cm directory (contains <z_lo>/...).",
    )
    ap.add_argument(
        "--production",
        default=DEFAULT_PRODUCTION,
        help="Subdirectory name under each z slice.",
    )
    ap.add_argument(
        "--cache-dir",
        type=Path,
        default=DEFAULT_CACHE_DIR,
        help="Where per-(z,y) pickles and the full bundle are written.",
    )
    ap.add_argument(
        "--voxel-size",
        type=int,
        default=DEFAULT_VOXEL_SIZE,
        help="Voxel size in cm (must match hist filenames).",
    )
    ap.add_argument(
        "--n-jobs",
        type=int,
        default=DEFAULT_N_JOBS,
        help="Number of job subdirectories (0 .. n-1) to merge per (z,y).",
    )
    ap.add_argument(
        "--force",
        action="store_true",
        help="Ignore on-disk caches and re-read all ROOT files.",
    )
    ap.add_argument(
        "--no-check",
        action="store_true",
        help="Skip the initial hist file inventory pass.",
    )
    ap.add_argument(
        "--verbose",
        action="store_true",
        help="Print extra diagnostics (e.g. missing POT next to a hist file).",
    )
    args = ap.parse_args()

    vs = args.voxel_size
    z_edges = edges_z(vs)
    y_edges = edges_y(vs)
    x_edges = edges_x(vs)
    z_lo_arr = z_edges[:-1]
    y_lo_arr = y_edges[:-1]
    nz, ny, nx = len(z_edges) - 1, len(y_edges) - 1, len(x_edges) - 1

    _timestamp_dir_cache.clear()

    if not args.no_check:
        check_files(
            base_dir=args.base_dir,
            production=args.production,
            voxel_size=vs,
            z_lo_arr=z_lo_arr,
            y_lo_arr=y_lo_arr,
            n_jobs=args.n_jobs,
            ny=ny,
            verbose=True,
        )

    data = build_full_voxel_array(
        base_dir=args.base_dir,
        production=args.production,
        voxel_size=vs,
        cache_dir=args.cache_dir,
        z_edges=z_edges,
        y_edges=y_edges,
        x_edges=x_edges,
        n_jobs=args.n_jobs,
        force=args.force,
        verbose=args.verbose,
    )

    flux = data["flux"]
    bin_edges = data["bin_edges"]
    if flux is None or bin_edges is None:
        raise SystemExit("Build produced no flux data (check inputs and logs).")

    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    n_jobs_arr = data["n_jobs"]
    pot_per_zy = data.get("pot_per_zy")
    if pot_per_zy is None:
        pot_per_zy = np.zeros((nz, ny), dtype=np.float64)

    out_path = args.cache_dir / "full_voxel_flux_potmacro.pkl"
    print("Wrote", out_path)
    print("flux shape (NZ, NY, NX, n_flav, n_E):", flux.shape)
    print("energy bins:", len(bin_centers), bin_edges[0], "...", bin_edges[-1], "GeV")
    if (n_jobs_arr > 0).any():
        print(
            "jobs merged per (z,y): min",
            int(n_jobs_arr[n_jobs_arr > 0].min()),
            "max",
            int(n_jobs_arr.max()),
        )
    if (pot_per_zy > 0).any():
        print(
            "POT sum per (z,y) slice: min",
            float(pot_per_zy[pot_per_zy > 0].min()),
            "max",
            float(pot_per_zy.max()),
        )


if __name__ == "__main__":
    main()
