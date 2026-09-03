"""Pure-numpy quick-look petrophysics.

Every function takes array-likes (numpy arrays or pandas Series), propagates
NaN, and returns a numpy array of the same shape. Porosities are fractions
(0-1) unless the name says ``_pu``. Resistivities are ohm-m.
"""

from __future__ import annotations

from typing import Dict, Optional

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


def matrix_density(matrix) -> float:
    """Grain density in g/cc from a :data:`MATRIX_DENSITY` key or a number.

    Raises ``ValueError`` naming the valid keys for an unknown name.
    """
    if isinstance(matrix, str):
        key = matrix.strip().lower()
        if key not in MATRIX_DENSITY:
            raise ValueError(f"unknown matrix {matrix!r}; choose one of {sorted(MATRIX_DENSITY)} or give a density in g/cc")
        return MATRIX_DENSITY[key]
    return float(matrix)


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
    rho_matrix = matrix_density(rho_matrix)
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


#: First-order neutron lithology correction, fraction, relative to a limestone-
#: scaled compensated-neutron reading: phi_true ~= phi_lime + offset[matrix].
#: Chart-derived round numbers (Schlumberger CNL, ~20 pu); good to about +-2 pu.
NEUTRON_MATRIX_OFFSET: Dict[str, float] = {"limestone": 0.0, "sandstone": 0.04, "dolomite": -0.06}


def neutron_lithology_correction(phin, from_matrix: str = "limestone", to_matrix: str = "limestone") -> np.ndarray:
    """Re-scale a neutron porosity (fraction) logged on ``from_matrix`` to ``to_matrix``.

    Linear chart approximation via :data:`NEUTRON_MATRIX_OFFSET`; raises
    ``ValueError`` for a matrix outside sandstone / limestone / dolomite —
    there is no sensible neutron scale for salt or anhydrite, and a numeric
    grain density says nothing about the tool's lithology response.
    """
    f, t = str(from_matrix).strip().lower(), str(to_matrix).strip().lower()
    for k in (f, t):
        if k not in NEUTRON_MATRIX_OFFSET:
            raise ValueError(f"no neutron lithology correction for {k!r}; supported: {sorted(NEUTRON_MATRIX_OFFSET)}")
    return _arr(phin) + (NEUTRON_MATRIX_OFFSET[t] - NEUTRON_MATRIX_OFFSET[f])


def rwa(rt, phi, a: float = 1.0, m: float = 2.0) -> np.ndarray:
    """Apparent water resistivity Rwa = Rt * phi^m / a (ohm-m). Equals Rw where Sw = 1."""
    rt = _arr(rt)
    phi = _arr(phi)
    with np.errstate(invalid="ignore"):
        out = rt * phi**m / a
    return np.where((phi > 0) & (rt > 0), out, np.nan)


def pick_rw_from_rwa(
    rt,
    phi,
    vsh=None,
    depth=None,
    q: float = 5.0,
    vsh_cut: float = 0.15,
    phi_min: float = 0.06,
    a: float = 1.0,
    m: float = 2.0,
    min_points: int = 20,
) -> Dict[str, object]:
    """Pick Rw as a low percentile of Rwa over clean, porous samples.

    The wettest clean rock has Rwa ~ Rw; everything hydrocarbon-bearing has
    Rwa > Rw. Taking the ``q``-th percentile rather than the minimum keeps a
    single bad sample from setting the pick. Unlike :func:`fit_water_line`
    this does not need a straight Pickett line, so it is the more robust
    choice in mixed-lithology sections — at the cost of assuming ``m``.

    Returns ``{"rw", "m", "a", "q", "n_points", "interval"}`` where
    ``interval`` is ``(top, base)`` of the samples at or below the pick when
    ``depth`` is given (else None). Raises ``ValueError`` with fewer than
    ``min_points`` usable samples.
    """
    rt = _arr(rt)
    phi = _arr(phi)
    ok = np.isfinite(rt) & np.isfinite(phi) & (rt > 0) & (phi >= phi_min)
    if vsh is not None:
        v = _arr(vsh)
        ok &= np.isfinite(v) & (v < vsh_cut)
    n = int(ok.sum())
    if n < min_points:
        raise ValueError(f"only {n} clean porous samples (need {min_points})")
    r = rwa(rt[ok], phi[ok], a=a, m=m)
    pick = float(np.nanpercentile(r, q))
    interval = None
    if depth is not None:
        d = _arr(depth)[ok][r <= pick]
        if d.size:
            interval = (float(np.nanmin(d)), float(np.nanmax(d)))
    return {"rw": pick, "m": float(m), "a": float(a), "q": float(q), "n_points": n, "interval": interval}


def pick_rsh(rt, vsh, vsh_min: float = 0.8, q: float = 50.0, min_points: int = 10) -> float:
    """Shale resistivity: the ``q``-th percentile of Rt where Vsh >= ``vsh_min``."""
    rt = _arr(rt)
    v = _arr(vsh)
    ok = np.isfinite(rt) & np.isfinite(v) & (rt > 0) & (v >= vsh_min)
    if int(ok.sum()) < min_points:
        raise ValueError(f"only {int(ok.sum())} shale samples with Vsh >= {vsh_min} (need {min_points})")
    return float(np.percentile(rt[ok], q))


def sw_simandoux(rt, phi, vsh, rw: float, rsh: float, a: float = 1.0, m: float = 2.0, n: float = 2.0, clip: bool = True) -> np.ndarray:
    """Modified Simandoux water saturation.

    Solves ``1/Rt = phi^m Sw^n / (a Rw) + Vsh Sw / Rsh`` for Sw. Closed form
    for ``n = 2``; Newton iterations from the Archie estimate otherwise. With
    ``Vsh = 0`` it reduces exactly to Archie.
    """
    if rw <= 0 or rsh <= 0:
        raise ValueError("rw and rsh must be positive")
    rt, phi, vsh = _arr(rt), _arr(phi), np.clip(_arr(vsh), 0.0, 1.0)
    valid = (phi > 0) & (rt > 0) & np.isfinite(vsh)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = phi**m / (a * rw)  # coefficient of Sw^n
        d = vsh / rsh          # coefficient of Sw
        if n == 2:
            sw = (np.sqrt(d**2 + 4.0 * c / rt) - d) / (2.0 * c)
        else:
            sw = archie_sw(rt, phi, rw, a=a, m=m, n=n, clip=False)
            sw = np.where(np.isfinite(sw), sw, 1.0)
            for _ in range(30):
                f = c * sw**n + d * sw - 1.0 / rt
                fp = n * c * sw ** (n - 1) + d
                sw = np.clip(sw - f / fp, 1e-6, 10.0)
    sw = np.where(valid, sw, np.nan)
    return np.minimum(sw, 1.0) if clip else sw


def sw_indonesia(rt, phi, vsh, rw: float, rsh: float, a: float = 1.0, m: float = 2.0, n: float = 2.0, clip: bool = True) -> np.ndarray:
    """Indonesia (Poupon-Leveaux) water saturation.

    ``1/sqrt(Rt) = [Vsh^(1 - Vsh/2) / sqrt(Rsh) + phi^(m/2) / sqrt(a Rw)] Sw^(n/2)``.
    With ``Vsh = 0`` it reduces exactly to Archie.
    """
    if rw <= 0 or rsh <= 0:
        raise ValueError("rw and rsh must be positive")
    rt, phi, vsh = _arr(rt), _arr(phi), np.clip(_arr(vsh), 0.0, 1.0)
    valid = (phi > 0) & (rt > 0) & np.isfinite(vsh)
    with np.errstate(invalid="ignore", divide="ignore"):
        term = vsh ** (1.0 - vsh / 2.0) / np.sqrt(rsh) + phi ** (m / 2.0) / np.sqrt(a * rw)
        sw = ((1.0 / np.sqrt(rt)) / term) ** (2.0 / n)
    sw = np.where(valid, sw, np.nan)
    return np.minimum(sw, 1.0) if clip else sw


SW_MODELS = ("archie", "simandoux", "indonesia")


def water_saturation(model: str, rt, phi, vsh, rw: float, rsh: Optional[float] = None, a: float = 1.0, m: float = 2.0, n: float = 2.0, clip: bool = True) -> np.ndarray:
    """Dispatch on ``model`` (one of :data:`SW_MODELS`). Shaly models need ``rsh``."""
    key = str(model).strip().lower()
    if key == "archie":
        return archie_sw(rt, phi, rw, a=a, m=m, n=n, clip=clip)
    if key not in SW_MODELS:
        raise ValueError(f"unknown Sw model {model!r}; choose from {SW_MODELS}")
    if rsh is None or not np.isfinite(rsh):
        raise ValueError(f"{key} needs rsh (shale resistivity); pick one with pick_rsh()")
    fn = sw_simandoux if key == "simandoux" else sw_indonesia
    return fn(rt, phi, vsh, rw, rsh, a=a, m=m, n=n, clip=clip)


def fit_water_line(
    rt,
    phi,
    q: float = 5.0,
    phi_min: float = 0.06,
    phi_max: float = 0.35,
    bin_width: float = 0.05,
    min_per_bin: int = 15,
    a: float = 1.0,
) -> Dict[str, object]:
    """Fit the Sw = 1 line of a Pickett plot from the low-Rt envelope.

    Points are binned in ``log10(phi)`` (``bin_width`` decades); in each bin
    the ``q``-th percentile of ``log10(Rt)`` is taken as the water line. A
    straight line through those envelope points has slope ``-m`` and
    intercept ``log10(a * Rw)``.

    Feed it *clean* points only (low Vsh). Keep ``phi_min`` at a few percent:
    below that, shale conductivity and matrix-density error flatten the
    envelope and drag ``m`` toward 1.

    Returns ``{"m", "rw", "a", "n_points", "envelope"}`` where ``envelope`` is
    an ``(n_bins, 3)`` array of ``[log10 phi, log10 Rt, count]``.
    Raises ``ValueError`` with fewer than three usable bins.
    """
    rt = _arr(rt)
    phi = _arr(phi)
    ok = np.isfinite(rt) & np.isfinite(phi) & (rt > 0) & (phi >= phi_min) & (phi <= phi_max)
    lp, lr = np.log10(phi[ok]), np.log10(rt[ok])
    edges = np.arange(np.log10(phi_min), np.log10(phi_max) + bin_width, bin_width)
    idx = np.digitize(lp, edges)
    env = []
    for b in np.unique(idx):
        sel = idx == b
        if sel.sum() < min_per_bin:
            continue
        env.append((lp[sel].mean(), np.percentile(lr[sel], q), int(sel.sum())))
    if len(env) < 3:
        raise ValueError(f"only {len(env)} porosity bins with >= {min_per_bin} points; need 3")
    env_arr = np.array(env)
    slope, intercept = np.polyfit(env_arr[:, 0], env_arr[:, 1], 1)
    return {
        "m": float(-slope),
        "rw": float(10**intercept / a),
        "a": float(a),
        "n_points": int(ok.sum()),
        "envelope": env_arr,
    }
