"""Pure-numpy quick-look petrophysics.

Every function takes array-likes (numpy arrays or pandas Series), propagates
NaN, and returns a numpy array of the same shape. Porosities are fractions
(0-1) unless the name says ``_pu``. Resistivities are ohm-m.
"""

from __future__ import annotations

from typing import Dict

import numpy as np

#: Grain (matrix) densities, g/cc. Textbook values; the original notebook
#: carried SS=2.45, LS=2.71, DL=2.85, Salt=2.30, of which only LS is standard.
MATRIX_DENSITY: Dict[str, float] = {
    "sandstone": 2.65,
    "limestone": 2.71,
    "dolomite": 2.87,
    "salt": 2.03,
    "anhydrite": 2.98,
}


def _arr(x) -> np.ndarray:
    return np.asarray(x, dtype=float)


def pu_to_frac(x) -> np.ndarray:
    """Porosity units (percent) -> fraction."""
    return _arr(x) / 100.0


def frac_to_pu(x) -> np.ndarray:
    """Fraction -> porosity units (percent)."""
    return _arr(x) * 100.0


def gamma_ray_index(gr, gr_clean: float, gr_dirty: float) -> np.ndarray:
    """IGR = (GR - GR_clean) / (GR_dirty - GR_clean), clipped to [0, 1]."""
    if gr_dirty <= gr_clean:
        raise ValueError("gr_dirty must exceed gr_clean")
    igr = (_arr(gr) - gr_clean) / (gr_dirty - gr_clean)
    return np.clip(igr, 0.0, 1.0)


def vshale_linear(gr, gr_clean: float, gr_dirty: float) -> np.ndarray:
    """Linear shale volume: Vsh = IGR."""
    return gamma_ray_index(gr, gr_clean, gr_dirty)


def vshale_larionov(gr, gr_clean: float, gr_dirty: float, older: bool = False) -> np.ndarray:
    """Larionov (1969) shale volume from the gamma-ray index.

    ``older=False`` uses the Tertiary-rock form ``0.083 * (2**(3.7 IGR) - 1)``;
    ``older=True`` the older-rock form ``0.33 * (2**(2 IGR) - 1)``.
    """
    igr = gamma_ray_index(gr, gr_clean, gr_dirty)
    if older:
        vsh = 0.33 * (2.0 ** (2.0 * igr) - 1.0)
    else:
        vsh = 0.083 * (2.0 ** (3.7 * igr) - 1.0)
    return np.clip(vsh, 0.0, 1.0)


def density_porosity(rhob, rho_matrix: float = 2.65, rho_fluid: float = 1.0) -> np.ndarray:
    """phi_D = (rho_ma - RHOB) / (rho_ma - rho_f), as a fraction (not clipped).

    ``rho_matrix`` may be a number or a key of :data:`MATRIX_DENSITY`.
    Negative values (RHOB above the matrix density) are returned as-is so the
    caller can see them; clip if you need a physical porosity.
    """
    if isinstance(rho_matrix, str):
        rho_matrix = MATRIX_DENSITY[rho_matrix]
    if rho_matrix <= rho_fluid:
        raise ValueError("rho_matrix must exceed rho_fluid")
    return (rho_matrix - _arr(rhob)) / (rho_matrix - rho_fluid)


def neutron_density_crossover(nphi, dphi, threshold: float = 0.0):
    """Gas flag from neutron-density separation.

    Both inputs are fractions. Returns ``(separation, gas_flag)`` where
    ``separation = dphi - nphi`` and ``gas_flag`` is True where the density
    porosity exceeds the neutron porosity by more than ``threshold`` (the
    classic gas crossover). NaN separation never flags.
    """
    sep = _arr(dphi) - _arr(nphi)
    with np.errstate(invalid="ignore"):
        flag = sep > threshold
    return sep, flag


def archie_sw(rt, phi, rw: float, a: float = 1.0, m: float = 2.0, n: float = 2.0, clip: bool = True) -> np.ndarray:
    """Archie water saturation: Sw = ((a Rw) / (phi^m Rt))^(1/n).

    ``phi`` is a fraction. Where ``phi <= 0`` or ``rt <= 0`` the result is NaN.
    With ``clip=True`` (default) values above 1 are clipped to 1, which is the
    usual presentation; pass ``clip=False`` to see the raw ratio.
    """
    if rw <= 0:
        raise ValueError("rw must be positive")
    rt = _arr(rt)
    phi = _arr(phi)
    with np.errstate(invalid="ignore", divide="ignore"):
        sw = ((a * rw) / (phi**m * rt)) ** (1.0 / n)
    sw = np.where((phi > 0) & (rt > 0), sw, np.nan)
    if clip:
        sw = np.minimum(sw, 1.0)
    return sw


def archie_sw_lines(phi, rw: float, sw_values=(1.0, 0.5, 0.25, 0.1), a: float = 1.0, m: float = 2.0, n: float = 2.0):
    """Rt as a function of phi for fixed Sw — the iso-saturation lines of a Pickett plot.

    Returns ``{sw: rt_array}``.
    """
    phi = _arr(phi)
    with np.errstate(divide="ignore", invalid="ignore"):
        return {sw: (a * rw) / (phi**m * sw**n) for sw in sw_values}
