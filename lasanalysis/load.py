"""Read a LAS file and clean the values lasio cannot know are missing.

lasio masks the declared ``NULL`` value (``-999.25`` in the KGS files) to NaN.
It cannot mask sentinels the service company never declared, and the Kansas
files carry at least two:

* ``100000.0`` in the induction resistivity curves (RILD / RILM) where the
  tool went off scale;
* ``-999.0`` (not ``-999.25``) in RHOB, which then drives DPOR to tens of
  thousands of porosity units on the same rows.

``read_las`` applies three unit-driven rules and records what it masked so a
notebook or test can see exactly how much data was discarded.
"""

from __future__ import annotations

import os
from typing import Dict, Iterable, Optional

import lasio
import numpy as np
import pandas as pd

#: Resistivity values at or above this are treated as off-scale (ohm-m).
RESISTIVITY_CEILING = 1e5

#: Sentinels seen in the data that are not the declared NULL.
EXTRA_NULLS = (-999.0,)

#: Plausible range for porosity curves in porosity units (percent).
POROSITY_RANGE_PU = (-50.0, 100.0)

#: Plausible range for bulk-density curves (g/cc). RHOB in the KGS file also
#: carries exact 0.0 on a few rows.
DENSITY_RANGE = (1.0, 3.5)

#: Declared NULL in the KGS exports (used for CSV frames that have no header).
DEFAULT_NULL = -999.25

#: Standard curve name -> mnemonics that service companies use for it.
#: Matching is case-insensitive; first hit wins.
ALIASES: Dict[str, list] = {
    "DEPT": ["DEPT", "DEPTH", "MD"],
    "GR": ["GR", "GRC", "SGR", "GAM"],
    "RT": ["RILD", "ILD", "RT", "LLD", "RDEP", "AT90"],
    "RM": ["RILM", "ILM", "RMED", "AT30"],
    "RXO": ["RLL3", "LL3", "SFL", "RXO", "MSFL"],
    "SP": ["SP"],
    "RHOB": ["RHOB", "DEN", "ZDEN", "RHOZ"],
    "NPHI": ["CNPOR", "NPHI", "NPOR", "TNPH", "CNC"],
    "DPHI": ["DPOR", "DPHI", "PHID"],
    "DT": ["DT", "DTC", "AC", "DTCO"],
    "CALI": ["DCAL", "CALI", "CAL", "HCAL"],
}


def _is_resistivity(unit: str) -> bool:
    return "OHM" in unit.upper()


def _is_porosity_pu(unit: str) -> bool:
    return unit.strip().upper() in ("PU", "%")


def _is_density(unit: str) -> bool:
    return unit.strip().upper().replace(" ", "") in ("G/CC", "G/CM3", "GM/CC", "KG/M3")


def _clean_array(
    data: np.ndarray,
    kind: Optional[str],
    resistivity_ceiling: float,
    extra_nulls: tuple,
    porosity_range_pu: Optional[tuple],
    density_range: Optional[tuple],
) -> tuple:
    """Apply the masking rules to one array. Returns ``(cleaned, counts)``.

    ``kind`` is ``"resistivity"``, ``"porosity_pu"``, ``"density"`` or None.
    """
    counts: Dict[str, int] = {}
    mask = np.zeros(data.shape, dtype=bool)
    for v in extra_nulls:
        mask |= np.isclose(data, v, rtol=0, atol=1e-6)
    if mask.any():
        counts["extra_null"] = int(mask.sum())

    with np.errstate(invalid="ignore"):
        if kind == "resistivity":
            m = (data >= resistivity_ceiling) & ~mask
            if m.any():
                counts["off_scale"] = int(m.sum())
            mask |= m
        elif kind == "porosity_pu" and porosity_range_pu is not None:
            lo, hi = porosity_range_pu
            m = ((data < lo) | (data > hi)) & ~mask
            if m.any():
                counts["porosity"] = int(m.sum())
            mask |= m
        elif kind == "density" and density_range is not None:
            lo, hi = density_range
            m = ((data < lo) | (data > hi)) & ~mask
            if m.any():
                counts["density"] = int(m.sum())
            mask |= m

    if mask.any():
        data = data.copy()
        data[mask] = np.nan
    return data, counts


def _kind_from_unit(unit: str) -> Optional[str]:
    if _is_resistivity(unit):
        return "resistivity"
    if _is_porosity_pu(unit):
        return "porosity_pu"
    if _is_density(unit):
        return "density"
    return None


#: Standard name -> rule kind, for frames that carry no units (CSV exports).
_KIND_FROM_STANDARD = {
    "RT": "resistivity",
    "RM": "resistivity",
    "RXO": "resistivity",
    "NPHI": "porosity_pu",
    "DPHI": "porosity_pu",
    "RHOB": "density",
}


def clean_curves(
    las: lasio.LASFile,
    resistivity_ceiling: float = RESISTIVITY_CEILING,
    extra_nulls: Iterable[float] = EXTRA_NULLS,
    porosity_range_pu: Optional[tuple] = POROSITY_RANGE_PU,
    density_range: Optional[tuple] = DENSITY_RANGE,
) -> Dict[str, Dict[str, int]]:
    """Mask undeclared sentinels in place. Returns ``{mnemonic: {rule: n_masked}}``.

    Rules (all leave the depth curve untouched; the rule is chosen from the
    curve's unit):

    ``extra_null``   value equals one of ``extra_nulls`` (any curve)
    ``off_scale``    value >= ``resistivity_ceiling`` (curves with an OHM unit)
    ``porosity``     value outside ``porosity_range_pu`` (curves in PU / %)
    ``density``      value outside ``density_range`` (curves in g/cc)
    """
    extra_nulls = tuple(float(v) for v in extra_nulls)
    report: Dict[str, Dict[str, int]] = {}
    for i, curve in enumerate(las.curves):
        if i == 0:  # index / depth
            continue
        data, counts = _clean_array(
            np.asarray(curve.data, dtype=float),
            _kind_from_unit(curve.unit),
            resistivity_ceiling,
            extra_nulls,
            porosity_range_pu,
            density_range,
        )
        if counts:
            curve.data = data
            report[curve.mnemonic] = counts
    return report


def clean_frame(
    df: pd.DataFrame,
    null: Optional[float] = DEFAULT_NULL,
    resistivity_ceiling: float = RESISTIVITY_CEILING,
    extra_nulls: Iterable[float] = EXTRA_NULLS,
    porosity_range_pu: Optional[tuple] = POROSITY_RANGE_PU,
    density_range: Optional[tuple] = DENSITY_RANGE,
) -> tuple:
    """Same rules as :func:`clean_curves` for a unit-less DataFrame.

    The rule for each column is chosen by matching its name against
    :data:`ALIASES` (so ``RILD`` gets the resistivity ceiling, ``CNPOR`` the
    porosity range, ...). ``null`` is masked everywhere first. Returns
    ``(cleaned_df, report)``; the input is not modified.
    """
    nulls = tuple(float(v) for v in extra_nulls)
    if null is not None:
        nulls = (float(null),) + nulls
    out = df.copy()
    report: Dict[str, Dict[str, int]] = {}
    for col in out.columns:
        if not pd.api.types.is_numeric_dtype(out[col]):
            continue
        kind = None
        for standard, k in _KIND_FROM_STANDARD.items():
            if find_curve([col], standard) is not None:
                kind = k
                break
        data, counts = _clean_array(
            out[col].to_numpy(dtype=float), kind, resistivity_ceiling, nulls, porosity_range_pu, density_range
        )
        if counts:
            out[col] = data
            report[str(col)] = counts
    return out, report


def read_log_csv(path: "os.PathLike | str", depth_col: str = "Depth", clean: bool = True, **clean_kwargs) -> pd.DataFrame:
    """Read a comma-separated log export (one row per depth) into a frame indexed by depth.

    Applies :func:`clean_frame` by default; the report is attached as
    ``df.attrs["mask_report"]``.
    """
    df = pd.read_csv(path)
    if depth_col not in df.columns:
        raise KeyError(f"{depth_col!r} not in columns {list(df.columns)}")
    df = df.set_index(depth_col)
    df.index.name = "DEPT"
    report: Dict[str, Dict[str, int]] = {}
    if clean:
        df, report = clean_frame(df, **clean_kwargs)
    df.attrs["mask_report"] = report
    return df


def read_las(path: "os.PathLike | str", clean: bool = True, **clean_kwargs) -> lasio.LASFile:
    """Read a LAS file with lasio and (by default) mask undeclared sentinels.

    The mask report is attached as ``las.mask_report`` (empty dict when
    ``clean=False`` or nothing was masked).
    """
    las = lasio.read(os.fspath(path))
    las.mask_report = clean_curves(las, **clean_kwargs) if clean else {}
    return las


def curves(las: lasio.LASFile) -> pd.DataFrame:
    """Curve data as a DataFrame indexed by depth, with upper-case column names."""
    df = las.df()
    df.columns = [str(c).upper() for c in df.columns]
    df.index.name = "DEPT"
    return df


def find_curve(columns: Iterable[str], standard: str) -> Optional[str]:
    """Return the column that carries ``standard`` (e.g. ``"RT"``), or None."""
    if standard not in ALIASES:
        raise KeyError(f"unknown standard curve {standard!r}; known: {sorted(ALIASES)}")
    upper = {str(c).upper(): c for c in columns}
    for alias in ALIASES[standard]:
        if alias.upper() in upper:
            return upper[alias.upper()]
    return None


def standardize(df: pd.DataFrame, aliases: Optional[Dict[str, list]] = None) -> pd.DataFrame:
    """Rename recognised curves to their standard names; leave the rest alone.

    Returns a new DataFrame. Columns are upper-cased first so that ``RxoRt`` and
    ``RXORT`` are the same curve.
    """
    aliases = ALIASES if aliases is None else aliases
    out = df.copy()
    out.columns = [str(c).upper() for c in out.columns]
    rename = {}
    for standard, names in aliases.items():
        for alias in names:
            if alias.upper() in out.columns and standard not in rename.values():
                rename[alias.upper()] = standard
                break
    return out.rename(columns=rename)
