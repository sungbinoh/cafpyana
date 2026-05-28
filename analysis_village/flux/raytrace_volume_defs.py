from __future__ import annotations

import numpy as np

# ==== Volume definitions ====
# nominal TPC F Vwith cathode cut: 10 < |x| < 190 cm, full |y| < 190, 10 < z < 450
FV_SPLIT_BOXES: list[dict] = [
    dict(x_range=(-190.0, -10.0), y_range=(-190.0, 190.0), z_range=(10.0, 450.0)),
    dict(x_range=(10.0, 190.0), y_range=(-190.0, 190.0), z_range=(10.0, 450.0)),
]

# + truncate highYZ -- split into 4 boxes
FV_SPLIT_TRUNCY_BOXES: list[dict] = [
    dict(x_range=(-190.0, -10.0), y_range=(-190.0, 190.0), z_range=(10.0, 250.0)),
    dict(x_range=(10.0, 190.0), y_range=(-190.0, 190.0), z_range=(10.0, 250.0)),
    dict(x_range=(-190.0, -10.0), y_range=(-190.0, 100.0), z_range=(250.0, 450.0)),
    dict(x_range=(10.0, 190.0), y_range=(-190.0, 100.0), z_range=(250.0, 450.0)),
]

RAYTRACE_VOLUME_DEFS: list[tuple[str, list[dict]]] = [
    ("FV_split", FV_SPLIT_BOXES),
    ("FV_split_truncY", FV_SPLIT_TRUNCY_BOXES),
]

RAYTRACE_VOLUME_LABEL: dict[str, str] = {
    "FV_split": r"$10<|x|<190,\;|y|<190,\;10<z<450$",
    "FV_split_truncY": r"$10<|x|<190$: $|y|<190$ ($10<z<250$), $-190<y<100$ ($250<z<450$)",
}

VOLUMES: dict[str, list[dict]] = dict(RAYTRACE_VOLUME_DEFS)
VOLUME_LABEL: dict[str, str] = dict(RAYTRACE_VOLUME_LABEL)

FV_SPLIT_VOLUME_CM3: float = 2 * 180.0 * 380.0 * 440.0  # 6.0192e7
FV_SPLIT_TRUNCY_VOLUME_CM3: float = (
    2 * 180.0 * 380.0 * 240.0 + 2 * 180.0 * 290.0 * 200.0
)  # 5.3712e7

# ==== Neutrino energy histogram configs ====
E_BIN_WIDTH_GEV: float = 0.05
E_MAX_GEV: float = 10.0
E_BINS: np.ndarray = np.arange(0.0, E_MAX_GEV + E_BIN_WIDTH_GEV, E_BIN_WIDTH_GEV)
BIN_CENTERS: np.ndarray = 0.5 * (E_BINS[:-1] + E_BINS[1:])
N_E_BINS: int = len(BIN_CENTERS)

# # ==== slab averages that match ray-traced FV volumes [cm] ====
# SLAB_Z_RANGE_FV_CM: tuple[float, float] = (10.0, 450.0)
# SLAB_COMPARE_Z_RANGE_CM: dict[str, tuple[float, float]] = {
#     "FV_split": SLAB_Z_RANGE_FV_CM,
#     "FV_split_truncY": SLAB_Z_RANGE_FV_CM,
# }


def _z_plane_in_box(z_cm: float, z_range: tuple[float, float]) -> bool:
    z_lo, z_hi = z_range
    return z_lo <= z_cm <= z_hi


def mask_xy_in_volume_at_z(
    x_cm: np.ndarray,
    y_cm: np.ndarray,
    z_cm: float,
    boxes: list[dict],
) -> np.ndarray:
    """Union of (x, y) footprints for boxes whose z_range contains ``z_cm``."""
    m = np.zeros_like(x_cm, dtype=bool)
    for box in boxes:
        if not _z_plane_in_box(z_cm, box["z_range"]):
            continue
        x_lo, x_hi = box["x_range"]
        y_lo, y_hi = box["y_range"]
        m |= (
            (x_cm >= x_lo)
            & (x_cm < x_hi)
            & (y_cm >= y_lo)
            & (y_cm < y_hi)
        )
    return m


def plane_area_cm2_at_z(z_cm: float, boxes: list[dict]) -> float:
    """Cross-sectional area [cm^2] of the volume union at plane z = z_cm."""
    area = 0.0
    for box in boxes:
        if not _z_plane_in_box(z_cm, box["z_range"]):
            continue
        dx = box["x_range"][1] - box["x_range"][0]
        dy = box["y_range"][1] - box["y_range"][0]
        area += dx * dy
    return area


def plane_area_m2_at_z(z_cm: float, boxes: list[dict]) -> float:
    return plane_area_cm2_at_z(z_cm, boxes) * 1.0e-4
