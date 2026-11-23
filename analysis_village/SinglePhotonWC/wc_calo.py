"""Calorimetry for Wire-Cell track points (T_rec_charge): dQ/dx -> dE/dx,
expected dE/dx vs residual range for mu / pi / p, and a per-track chi2.

dQ/dx -> dE/dx: ellipsoidal modified box recombination with the SBND
constants of makedf/calo.py. The flat CAF carries no gain, so Wire-Cell's
charge units are tied to electrons by a single scale, DQDX_SCALE (electrons
per cm per unit of q / nq): the median measured / expected dE/dx of true
stopping muons (2 < rr <= 26 cm) is 1 (Nov16 local CAFs, 103 muons; true
protons then come out 12% high, per-point spread ~20% for muons, ~27% for
protons).

Recombination: dqdx_to_dedx / dedx_to_dqdx take the box parameters
(default: the makedf/calo.py ones the scale was fitted with). RECOMB holds
the SBND MC and data sets of makedf/chi2pid.SBND_CALO_PARAMS, and
swap_recombination() re-derives MC dQ/dx as if the charge had recombined as
the data parameters say: measured dQ/dx -> dE/dx (MC set) -> dQ/dx (data set).

Residual range: Wire-Cell's T_rec_charge.rr restarts on every segment and
points the wrong way for ~1/4 of segments, so residual_range() uses the path
length along the particle's own points to the end farther from the vertex
(the true end of a stopping particle for ~80% of them).

Expected dE/dx(rr): kinetic energy at residual range rr from the CSDA range
(mean Bethe-Bloch with the Sternheimer density effect for argon), then the
Landau-Vavilov most-probable dE/dx (PDG eq. 34.12) for a step of PITCH_CM,
Wire-Cell's typical fit-point spacing. No electron-lifetime or yz correction.
"""

import numpy as np

# --- Recombination (makedf/calo.py, SBND)
MODA = 0.906
B90 = 0.203
R_EMB = 1.25
WION = 23.6e-6        # MeV per electron
RHO = 1.38434         # g/cm^3
EFIELD = 0.5          # kV/cm

# --- Argon and particles
K_BETHE = 0.307075    # MeV mol^-1 cm^2
Z_A = 18.0 / 39.948
I_EV = 188.0e-6       # MeV
ME = 0.51099895       # MeV
MASS = {"mu": 105.6584, "pi": 139.5704, "p": 938.2721}   # MeV
# Sternheimer density-effect parameters for liquid argon
STERN = dict(C=5.2146, x0=0.2, x1=3.0, a=0.19559, k=3.0)

PITCH_CM = 0.6        # cm, Wire-Cell fit-point spacing
RR_MAX_CM = 26.0      # cm, residual range used in the chi2 (as LArSoft chi2 PID)
DEDX_SIGMA_FRAC = 0.20   # fractional per-point resolution assumed in the chi2
CHI2_POINT_MAX = 25.0    # per-point chi2 cap (delta rays, mis-assigned charge)
DQDX_SCALE = 13.8     # electrons/cm per Wire-Cell q/nq unit; true stopping muons (Nov16)
# Pandora calorimetry gains, ADC per electron for planes 0-2 (makedf/chi2pid.SBND_CALO_PARAMS)
PANDORA_GAIN = {"mc": (0.0203521, 0.0202351, 0.0200727), "data": (0.02, 0.02, 0.02)}
RECOMB = {   # (alpha, beta_90, R), makedf/chi2pid.SBND_CALO_PARAMS
    "mc": (0.904, 0.204, 1.25),
    "data": (0.9366111025888879, 0.1835451204193374, 1.567685536351266),
}


def _beta(phi, b90=B90, r=R_EMB):
    return b90 / np.sqrt(np.sin(phi) ** 2 + np.cos(phi) ** 2 / r ** 2) / (RHO * EFIELD)


def dqdx_to_dedx(dqdx_wc, phi=np.pi / 2, scale=None, recomb=(MODA, B90, R_EMB)):
    """Wire-Cell dQ/dx (q / nq) -> dE/dx [MeV/cm]. phi: track angle to the
    drift (x) axis [rad]; recomb: (alpha, beta_90, R)."""
    s = DQDX_SCALE if scale is None else scale
    a, b90, r = recomb
    b = _beta(phi, b90, r)
    return (np.exp(np.asarray(dqdx_wc, float) * s * WION * b) - a) / b


def dedx_to_dqdx(dedx, phi=np.pi / 2, scale=None, recomb=(MODA, B90, R_EMB)):
    """dE/dx [MeV/cm] -> Wire-Cell dQ/dx (q / nq), the inverse of dqdx_to_dedx."""
    s = DQDX_SCALE if scale is None else scale
    a, b90, r = recomb
    b = _beta(phi, b90, r)
    with np.errstate(invalid="ignore"):
        return np.log(a + b * np.asarray(dedx, float)) / (b * WION * s)


def swap_recombination(dqdx_wc, phi, src="mc", dst="data", scale=None):
    """dQ/dx re-derived with the dst recombination: the measured charge's
    dE/dx under RECOMB[src], back to charge under RECOMB[dst]. scale:
    electrons per cm per unit of dqdx (default DQDX_SCALE, Wire-Cell q / nq;
    1 / gain for Pandora's ADC/cm)."""
    return dedx_to_dqdx(dqdx_to_dedx(dqdx_wc, phi, scale, recomb=RECOMB[src]), phi, scale, recomb=RECOMB[dst])


def _density(bg):
    x = np.log10(bg)
    c = STERN
    d = np.where(x >= c["x1"], 2 * np.log(10) * x - c["C"],
                 np.where(x >= c["x0"], 2 * np.log(10) * x - c["C"] + c["a"] * (c["x1"] - x) ** c["k"], 0.0))
    return d


def _bg(ke, mass):
    g = 1 + ke / mass
    return np.sqrt(g ** 2 - 1), g


def mean_dedx(ke, mass):
    """Mean Bethe-Bloch dE/dx [MeV/cm] at kinetic energy ke [MeV]."""
    bg, g = _bg(ke, mass)
    b2 = bg ** 2 / g ** 2
    wmax = 2 * ME * bg ** 2 / (1 + 2 * g * ME / mass + (ME / mass) ** 2)
    return K_BETHE * Z_A / b2 * (0.5 * np.log(2 * ME * bg ** 2 * wmax / I_EV ** 2) - b2 - _density(bg) / 2) * RHO


def mpv_dedx(ke, mass, dx=PITCH_CM):
    """Landau-Vavilov most-probable dE/dx [MeV/cm] for a step dx [cm]."""
    bg, g = _bg(ke, mass)
    b2 = bg ** 2 / g ** 2
    xi = K_BETHE / 2 * Z_A * RHO * dx / b2
    return xi / dx * (np.log(2 * ME * bg ** 2 / I_EV) + np.log(xi / I_EV) + 0.200 - b2 - _density(bg))


def _range_table(mass):
    ke = np.logspace(-1.5, 4, 4000)                     # MeV
    r = np.r_[0, np.cumsum(np.diff(ke) / mean_dedx(0.5 * (ke[1:] + ke[:-1]), mass))]
    return r, ke


_TABLES = {p: _range_table(m) for p, m in MASS.items()}


def expected_dedx(rr, particle):
    """Most-probable dE/dx [MeV/cm] at residual range rr [cm] for 'mu', 'pi' or 'p'."""
    r, ke = _TABLES[particle]
    return mpv_dedx(np.interp(np.asarray(rr, float), r, ke), MASS[particle])


def chi2(dedx, rr, particle):
    """Mean chi2 per point of measured dE/dx against the expectation (each
    point capped at CHI2_POINT_MAX), using points with 0 < rr <= RR_MAX_CM;
    NaN with fewer than 3 such points."""
    dedx, rr = np.asarray(dedx, float), np.asarray(rr, float)
    m = (rr > 0) & (rr <= RR_MAX_CM) & np.isfinite(dedx)
    if m.sum() < 3:
        return np.nan
    e = expected_dedx(rr[m], particle)
    return float(np.mean(np.minimum((dedx[m] - e) ** 2 / (DEDX_SIGMA_FRAC * e) ** 2, CHI2_POINT_MAX)))


def residual_range(xyz, vtx):
    """(order, rr, phi): the particle's points ordered from the vertex end
    along their principal axis, the path length from each to the far end
    (residual range, cm), and the axis angle to the drift (x) direction."""
    c = xyz.mean(0)
    u = np.linalg.svd(xyz - c, full_matrices=False)[2][0]
    t = (xyz - c) @ u
    if abs(t.min() - (vtx - c) @ u) > abs(t.max() - (vtx - c) @ u):
        t, u = -t, -u
    o = np.argsort(t)
    path = np.r_[0, np.cumsum(np.linalg.norm(np.diff(xyz[o], axis=0), axis=1))]
    return o, path[-1] - path, np.arccos(min(1.0, abs(u[0])))
