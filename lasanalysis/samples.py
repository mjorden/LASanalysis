"""Sample data model: lab results (SRA / Rock-Eval, XRD, core) joined to logs.

Geochemistry and mineralogy arrive as *samples* — tens per well, at a depth
or over an interval, from a lab, by a method — not as continuous curves. This
module keeps them in one long table so any analyte from any lab can be
carried, validated and joined the same way::

    api  depth_top  depth_base  sample_type  lab  method  analyte  value  unit  source

Identity is the **API number** (``api10`` = first ten digits), because lab
reports know API numbers and not KGS KIDs. Depths are as reported by the lab;
:func:`depth_shift` estimates a per-well core-to-log shift when a property
measured on the samples (core GR, core bulk density) is also logged.

Readers take a *wide* file (one row per sample, one column per analyte) —
see ``docs/templates/`` — and return the long form.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

SAMPLE_COLUMNS = ["api", "api10", "depth_top", "depth_base", "sample_type", "lab", "method", "analyte", "value", "unit", "source"]
SAMPLE_TYPES = ("core", "cuttings", "swc", "outcrop", "unknown")

#: Analyte -> (unit, low, high). None = unbounded on that side.
SRA_ANALYTES: Dict[str, tuple] = {
    "TOC": ("wt%", 0.0, 100.0),
    "S1": ("mg HC/g", 0.0, None),
    "S2": ("mg HC/g", 0.0, None),
    "S3": ("mg CO2/g", 0.0, None),
    "Tmax": ("degC", 300.0, 650.0),
    "HI": ("mg HC/g TOC", 0.0, 1200.0),
    "OI": ("mg CO2/g TOC", 0.0, 600.0),
    "PI": ("frac", 0.0, 1.0),
    "Ro": ("%", 0.0, 6.0),
}
XRD_MINERALS = ("quartz", "k_feldspar", "plagioclase", "calcite", "dolomite", "ankerite", "siderite", "pyrite", "apatite",
                "anhydrite", "halite", "illite", "smectite", "mixed_layer", "kaolinite", "chlorite", "total_clay", "other")
XRD_ANALYTES: Dict[str, tuple] = {m: ("wt%", 0.0, 100.0) for m in XRD_MINERALS}
#: Core-measured properties that can be compared with a log curve for depth shifting.
CORE_ANALYTES: Dict[str, tuple] = {"core_GR": ("GAPI", 0.0, 3000.0), "core_RHOB": ("g/cc", 1.0, 3.5), "core_phi": ("frac", 0.0, 1.0), "core_perm": ("mD", 0.0, None)}
ANALYTES: Dict[str, tuple] = {**SRA_ANALYTES, **XRD_ANALYTES, **CORE_ANALYTES}

#: Column names accepted for the sample metadata in a wide file (case-insensitive).
_META_ALIASES = {
    "api": ("api", "api_number", "api_no", "uwi", "well_api"),
    "depth_top": ("depth_top", "top", "depth", "top_depth", "from"),
    "depth_base": ("depth_base", "base", "bottom", "bottom_depth", "to"),
    "sample_type": ("sample_type", "type", "sample"),
    "lab": ("lab", "laboratory"),
    "method": ("method", "instrument", "technique"),
    "source": ("source", "report", "report_id", "reference"),
}
#: Analyte column aliases in a wide file (case-insensitive; the canonical name always matches).
_ANALYTE_ALIASES = {
    "TOC": ("toc", "toc_wt", "toc_wt%", "total_organic_carbon"),
    "Tmax": ("tmax", "t_max"),
    "Ro": ("ro", "vr", "vitrinite", "ro_pct", "vro"),
    "k_feldspar": ("kspar", "k-feldspar", "kfeldspar", "orthoclase"),
    "plagioclase": ("plag", "albite"),
    "mixed_layer": ("i/s", "illite_smectite", "mixed-layer", "ml"),
    "total_clay": ("clay", "clays", "total_clays"),
    "core_GR": ("core_gr", "coregr", "gr_core"),
    "core_RHOB": ("core_rhob", "rhob_core", "grain_density_bulk", "bulk_density"),
    "core_phi": ("core_phi", "phi_core", "core_porosity", "porosity"),
    "core_perm": ("core_perm", "perm", "permeability", "kair"),
}


def normalize_api(api) -> Tuple[str, Optional[str]]:
    """``("15-195-23011", "1519523011")`` from any API spelling; ``("", None)`` when empty.

    Keeps the string as given (stripped) and derives ``api10``, the ten-digit
    state-county-well number every source agrees on. Fewer than ten digits is
    not an API number and yields ``api10 = None``.
    """
    s = "" if api is None or (isinstance(api, float) and np.isnan(api)) else str(api).strip()
    digits = re.sub(r"\D", "", s)
    return s, (digits[:10] if len(digits) >= 10 else None)


def _find_col(columns: Iterable[str], names: Sequence[str]) -> Optional[str]:
    low = {str(c).strip().lower(): c for c in columns}
    for n in names:
        if n.lower() in low:
            return low[n.lower()]
    return None


def _analyte_for_column(col: str) -> Optional[str]:
    c = str(col).strip()
    for name in ANALYTES:
        if c.lower() == name.lower():
            return name
    for name, aliases in _ANALYTE_ALIASES.items():
        if c.lower() in aliases:
            return name
    return None


def read_samples(path_or_df, *, lab: str = "", method: str = "", source: str = "", sample_type: str = "",
                 analytes: Optional[Iterable[str]] = None) -> pd.DataFrame:
    """Read a wide sample file (CSV / Excel / DataFrame) into the long form.

    One row per sample; metadata columns are matched case-insensitively
    (``api``, ``depth_top`` or ``depth``, ``depth_base`` optional,
    ``sample_type``, ``lab``, ``method``, ``source``); every other column
    whose name is a known analyte (or alias) becomes rows. Keyword defaults
    fill metadata the file does not carry. Unknown columns are ignored —
    pass ``analytes`` to restrict to a subset. Raises ``ValueError`` when no
    API or depth column can be found, or no analyte column at all.
    """
    if isinstance(path_or_df, pd.DataFrame):
        wide = path_or_df.copy()
    else:
        p = Path(path_or_df)
        wide = pd.read_excel(p) if p.suffix.lower() in (".xlsx", ".xls") else pd.read_csv(p)
    cols = list(wide.columns)
    api_col = _find_col(cols, _META_ALIASES["api"])
    top_col = _find_col(cols, _META_ALIASES["depth_top"])
    if api_col is None or top_col is None:
        raise ValueError(f"need an API column and a depth column; found columns {cols}")
    base_col = _find_col(cols, _META_ALIASES["depth_base"])
    meta_cols = {k: _find_col(cols, v) for k, v in _META_ALIASES.items()}
    wanted = None if analytes is None else {a for a in analytes}
    analyte_cols = {}
    for c in cols:
        if c in meta_cols.values():
            continue
        a = _analyte_for_column(c)
        if a and (wanted is None or a in wanted):
            analyte_cols[c] = a
    if not analyte_cols:
        raise ValueError(f"no analyte columns recognised in {cols}; known analytes: {sorted(ANALYTES)}")

    rows = []
    for _, r in wide.iterrows():
        api_str, api10 = normalize_api(r[api_col])
        top = pd.to_numeric(r[top_col], errors="coerce")
        base = pd.to_numeric(r[base_col], errors="coerce") if base_col else np.nan
        if pd.isna(base):
            base = top
        meta = {
            "api": api_str, "api10": api10, "depth_top": float(top), "depth_base": float(base),
            "sample_type": str(r[meta_cols["sample_type"]]).strip().lower() if meta_cols["sample_type"] and not pd.isna(r[meta_cols["sample_type"]]) else (sample_type or "unknown"),
            "lab": str(r[meta_cols["lab"]]).strip() if meta_cols["lab"] and not pd.isna(r[meta_cols["lab"]]) else lab,
            "method": str(r[meta_cols["method"]]).strip() if meta_cols["method"] and not pd.isna(r[meta_cols["method"]]) else method,
            "source": str(r[meta_cols["source"]]).strip() if meta_cols["source"] and not pd.isna(r[meta_cols["source"]]) else source,
        }
        for c, a in analyte_cols.items():
            v = pd.to_numeric(r[c], errors="coerce")
            if pd.isna(v):
                continue
            rows.append({**meta, "analyte": a, "value": float(v), "unit": ANALYTES[a][0]})
    out = pd.DataFrame(rows, columns=SAMPLE_COLUMNS)
    out["sample_type"] = out["sample_type"].where(out["sample_type"].isin(SAMPLE_TYPES), "unknown")
    return out


def validate_samples(samples: pd.DataFrame, xrd_tolerance: float = 5.0, strict: bool = False) -> List[str]:
    """Check a long sample table. Returns a list of problems (empty when clean).

    Checks: required columns; every row has ``api10`` and a finite
    ``depth_top <= depth_base``; values inside the analyte's physical range
    (:data:`ANALYTES`); per-sample XRD totals within ``xrd_tolerance`` wt% of
    100; duplicate (api10, depth, analyte, lab) rows. ``strict=True`` raises
    ``ValueError`` with the list instead of returning it.
    """
    problems: List[str] = []
    missing = [c for c in SAMPLE_COLUMNS if c not in samples.columns]
    if missing:
        problems.append(f"missing columns {missing}")
        if strict:
            raise ValueError("; ".join(problems))
        return problems
    if samples.empty:
        return problems
    bad_api = samples["api10"].isna() | (samples["api10"] == "")
    if bad_api.any():
        problems.append(f"{int(bad_api.sum())} rows without a usable API number (apis: {sorted(samples.loc[bad_api, 'api'].astype(str).unique())[:5]})")
    top, base = samples["depth_top"], samples["depth_base"]
    bad_depth = ~np.isfinite(top) | ~np.isfinite(base) | (top > base)
    if bad_depth.any():
        problems.append(f"{int(bad_depth.sum())} rows with missing or inverted depths")
    unknown = ~samples["analyte"].isin(ANALYTES)
    if unknown.any():
        problems.append(f"unknown analytes {sorted(samples.loc[unknown, 'analyte'].unique())}")
    for a, (unit, lo, hi) in ANALYTES.items():
        v = samples.loc[samples["analyte"] == a, "value"]
        if v.empty:
            continue
        out = pd.Series(False, index=v.index)
        if lo is not None:
            out |= v < lo
        if hi is not None:
            out |= v > hi
        if out.any():
            problems.append(f"{a}: {int(out.sum())} values outside [{lo}, {hi}] {unit} (e.g. {v[out].iloc[0]:g})")
    xrd = samples[samples["analyte"].isin(XRD_MINERALS) & (samples["analyte"] != "total_clay")]
    if not xrd.empty:
        totals = xrd.groupby(["api10", "depth_top", "depth_base", "lab"])["value"].sum()
        off = totals[(totals - 100.0).abs() > xrd_tolerance]
        if not off.empty:
            problems.append(f"XRD totals off 100 by more than {xrd_tolerance} wt% for {len(off)} sample(s) (e.g. {off.iloc[0]:.1f} at {off.index[0][1]:g} ft)")
    dup = samples.duplicated(subset=["api10", "depth_top", "depth_base", "analyte", "lab"], keep=False)
    if dup.any():
        problems.append(f"{int(dup.sum())} duplicate rows (same well, depth, analyte, lab)")
    if strict and problems:
        raise ValueError("; ".join(problems))
    return problems


def to_wide(samples: pd.DataFrame, analytes: Optional[Iterable[str]] = None) -> pd.DataFrame:
    """One row per sample (api10, depth_top, depth_base, sample_type, lab, method, source), one column per analyte."""
    keys = ["api10", "api", "depth_top", "depth_base", "sample_type", "lab", "method", "source"]
    df = samples if analytes is None else samples[samples["analyte"].isin(list(analytes))]
    if df.empty:
        return pd.DataFrame(columns=keys)
    wide = df.pivot_table(index=keys, columns="analyte", values="value", aggfunc="mean").reset_index()
    wide.columns.name = None
    wide["depth_mid"] = (wide["depth_top"] + wide["depth_base"]) / 2
    return wide.sort_values(["api10", "depth_top"]).reset_index(drop=True)


def _log_at(log: pd.DataFrame, curve: str, depths: np.ndarray, tolerance: float) -> np.ndarray:
    """Log value nearest each depth (NaN beyond ``tolerance``)."""
    idx = log.index.to_numpy(dtype=float)
    vals = log[curve].to_numpy(dtype=float)
    order = np.argsort(idx)
    idx, vals = idx[order], vals[order]
    pos = np.searchsorted(idx, depths)
    pos = np.clip(pos, 1, len(idx) - 1)
    left, right = pos - 1, pos
    choose = np.where(np.abs(idx[left] - depths) <= np.abs(idx[right] - depths), left, right)
    out = vals[choose]
    out[np.abs(idx[choose] - depths) > tolerance] = np.nan
    return out


def join_to_log(samples_wide: pd.DataFrame, log: pd.DataFrame, curves: Optional[Iterable[str]] = None,
                tolerance: float = 0.5, shift: float = 0.0) -> pd.DataFrame:
    """Attach log curve values to samples (one well's log, one well's samples).

    Point samples (``depth_top == depth_base``) take the nearest log sample
    within ``tolerance``; interval samples take the mean of the log over the
    interval. ``shift`` is added to the sample depths first (see
    :func:`depth_shift`). Returns a copy with one ``log_<curve>`` column per
    curve (default: every numeric log column).
    """
    out = samples_wide.copy()
    curves = [c for c in (curves or log.columns) if c in log.columns and pd.api.types.is_numeric_dtype(log[c])]
    top = out["depth_top"].to_numpy(float) + shift
    base = out["depth_base"].to_numpy(float) + shift
    mid = (top + base) / 2
    idx = log.index.to_numpy(float)
    for c in curves:
        vals = _log_at(log, c, mid, tolerance)
        for i in np.where(base - top > 0)[0]:
            m = (idx >= top[i]) & (idx <= base[i])
            if m.any():
                v = log[c].to_numpy(float)[m]
                vals[i] = np.nanmean(v) if np.isfinite(v).any() else np.nan
        out[f"log_{c}"] = vals
    return out


def samples_on_grid(samples_wide: pd.DataFrame, log_index: pd.Index, analyte: str, shift: float = 0.0) -> pd.Series:
    """Place one analyte on the log's depth grid: the value over each sample's interval, NaN elsewhere.

    Point samples occupy the single nearest grid depth. Useful for plotting a
    sample track next to curves or for regressions against log values.
    """
    idx = log_index.to_numpy(float)
    out = np.full(idx.shape, np.nan)
    df = samples_wide.dropna(subset=[analyte]) if analyte in samples_wide else samples_wide.iloc[0:0]
    for _, r in df.iterrows():
        top, base = r["depth_top"] + shift, r["depth_base"] + shift
        if base > top:
            m = (idx >= top) & (idx <= base)
        else:
            m = np.zeros(idx.shape, bool)
            m[int(np.argmin(np.abs(idx - top)))] = True
        out[m] = r[analyte]
    return pd.Series(out, index=log_index, name=analyte)


def depth_shift(samples_wide: pd.DataFrame, log: pd.DataFrame, sample_curve: str, log_curve: str,
                max_shift: float = 20.0, step: float = 0.25, tolerance: float = 0.5, min_points: int = 5) -> Dict[str, object]:
    """Estimate the shift that best aligns a property measured on the samples with its log.

    Slides the samples by ``-max_shift .. +max_shift`` in ``step`` increments,
    reads the log at each shifted depth (:func:`join_to_log`), and keeps the
    shift with the highest Pearson correlation between sample and log values
    (RMSE breaks ties). A positive shift means the samples must move *deeper*
    to match the log. Returns ``{"shift", "r", "rmse", "n", "curve"}``; raises
    ``ValueError`` when fewer than ``min_points`` samples have the property.
    """
    df = samples_wide.dropna(subset=[sample_curve]) if sample_curve in samples_wide else samples_wide.iloc[0:0]
    if len(df) < min_points:
        raise ValueError(f"only {len(df)} samples carry {sample_curve!r} (need {min_points})")
    y = df[sample_curve].to_numpy(float)
    best = None
    for s in np.arange(-max_shift, max_shift + step / 2, step):
        x = join_to_log(df, log, [log_curve], tolerance=tolerance, shift=float(s))[f"log_{log_curve}"].to_numpy(float)
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() < min_points or np.std(x[ok]) == 0 or np.std(y[ok]) == 0:
            continue
        r = float(np.corrcoef(x[ok], y[ok])[0, 1])
        rmse = float(np.sqrt(np.mean((x[ok] - y[ok]) ** 2)))
        cand = (r, -rmse, float(s), int(ok.sum()))
        if best is None or cand[:2] > best[:2]:
            best = cand
    if best is None:
        raise ValueError(f"no shift in +-{max_shift} gave {min_points} overlapping points")
    r, neg_rmse, s, n = best
    return {"shift": s, "r": r, "rmse": -neg_rmse, "n": n, "curve": f"{sample_curve} vs {log_curve}"}


def apply_shift(samples: pd.DataFrame, shift: float) -> pd.DataFrame:
    """Return a copy with ``shift`` added to ``depth_top`` / ``depth_base`` (long or wide table)."""
    out = samples.copy()
    for c in ("depth_top", "depth_base", "depth_mid"):
        if c in out:
            out[c] = out[c] + shift
    return out


#: Marker colours for sample analytes on plots and in the viewer.
ANALYTE_COLORS = {"TOC": "#1f1f1f", "S2": "#8c564b", "Tmax": "#d62728", "HI": "#9467bd", "Ro": "#e377c2",
                  "quartz": "#ff7f0e", "calcite": "#1f77b4", "dolomite": "#17becf", "total_clay": "#2ca02c",
                  "core_GR": "#2ca02c", "core_RHOB": "#d62728", "core_phi": "#1f77b4", "core_perm": "#7f7f7f"}


def sample_tracks(samples, analytes: Sequence[str] = ("TOC",), shift: float = 0.0, title: Optional[str] = None,
                  one_track: bool = True) -> List[dict]:
    """Track specs (for :func:`~lasanalysis.plot.plot_tracks` / the viewer) that show samples as markers.

    ``samples`` may be the long or the wide table. Analytes with no values
    are skipped. With ``one_track`` every analyte shares one track (fine for
    TOC alone; pass ``one_track=False`` for a track per analyte when units
    differ). Each track carries ``points`` — ``[{label, depth, value, color, unit}]``
    — and an ``xlim`` from the data.
    """
    wide = samples if "analyte" not in samples.columns else to_wide(samples, analytes)
    tracks: List[dict] = []
    for a in analytes:
        pts = sample_points(wide, a, shift)
        if not pts["value"]:
            continue
        unit = ANALYTES.get(a, ("", None, None))[0]
        point = {"label": f"{a} [{unit}]" if unit else a, "depth": pts["depth"], "value": pts["value"], "color": ANALYTE_COLORS.get(a, "#000000"), "unit": unit}
        hi = max(pts["value"]) * 1.15 or 1.0
        if one_track and tracks:
            tracks[0]["points"].append(point)
            tracks[0]["xlim"] = (0, max(tracks[0]["xlim"][1], hi))
        else:
            tracks.append({"curves": [], "points": [point], "xlim": (0, hi), "title": title or (f"{a} (samples)" if not one_track else "samples")})
    return tracks


def sample_points(samples_wide: pd.DataFrame, analyte: str, shift: float = 0.0) -> Dict[str, list]:
    """``{"depth": [...], "value": [...]}`` at sample mid-depths — the shape the plot / viewer track specs take."""
    if analyte not in samples_wide:
        return {"depth": [], "value": []}
    df = samples_wide.dropna(subset=[analyte])
    mid = ((df["depth_top"] + df["depth_base"]) / 2 + shift).to_numpy(float)
    return {"depth": [float(v) for v in mid], "value": [float(v) for v in df[analyte].to_numpy(float)]}
