"""Run the quick-look workflow over many wells.

``run_well`` does one LAS file: load, standardise mnemonics, compute what the
available curves allow (Vsh, porosities, Sw), save a track plot, and return a
one-row summary. ``run_search`` chains a KGS search, ``fetch_las`` per hit and
``run_well`` per file, and writes ``summary.csv``. Both take the search and
fetch functions as arguments so they can be tested without a network.

Command line::

    python -m lasanalysis.multiwell --township 13 --range 22 --ew W --section 35 --out output/T13S_R22W_35
    python -m lasanalysis.multiwell --lease PEARSON --out output/pearson
    python -m lasanalysis.multiwell --las data/1046139243.las --out output/pearson --depth 3400 4200
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Callable, Dict, List, Optional

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import kgs
from .load import curves, read_las, standardize
from .petro import archie_sw, density_porosity, matrix_density, pu_to_frac, vshale_larionov
from .plot import plot_tracks

#: Petrophysical parameters. Rw / m are the Pearson #1-35 picks (see notebook
#: section 5); override per area.
DEFAULT_PARAMS: Dict[str, float] = {
    "gr_clean": 20.0,
    "gr_dirty": 110.0,
    "matrix": "limestone",
    "rho_fluid": 1.0,
    "rw": 0.03,
    "a": 1.0,
    "m": 2.0,
    "n": 2.0,
    # pay flags: porous, clean, hydrocarbon-bearing
    "phi_cut": 0.08,
    "vsh_cut": 0.30,
    "sw_cut": 0.50,
}


def parse_params(items: List[str]) -> Dict[str, object]:
    """``["rw=0.04", "matrix=dolomite"]`` -> validated overrides for :data:`DEFAULT_PARAMS`.

    Raises ``ValueError`` for an unknown key, a non-numeric value, or a matrix
    name that is not in :data:`~lasanalysis.petro.MATRIX_DENSITY`.
    """
    params: Dict[str, object] = {}
    for kv in items:
        k, sep, v = kv.partition("=")
        if not sep or k not in DEFAULT_PARAMS:
            raise ValueError(f"unknown param {k!r}; choose from {sorted(DEFAULT_PARAMS)}")
        if k == "matrix":
            try:
                params[k] = float(v)
            except ValueError:
                matrix_density(v)  # raises with the valid names
                params[k] = v.strip().lower()
        else:
            try:
                params[k] = float(v)
            except ValueError:
                raise ValueError(f"param {k} must be a number, got {v!r}") from None
    return params


def analyze(df: pd.DataFrame, params: Optional[Dict[str, float]] = None) -> pd.DataFrame:
    """Add VSH, PHID, PHIN, PHIND, SW and PAY columns where the inputs exist.

    ``df`` must already be :func:`~lasanalysis.load.standardize`-d. Missing
    inputs leave the derived column absent rather than raising.
    """
    p = {**DEFAULT_PARAMS, **(params or {})}
    out = df.copy()
    if "GR" in out:
        out["VSH"] = vshale_larionov(out["GR"], p["gr_clean"], p["gr_dirty"], older=True)
    if "RHOB" in out:
        out["PHID"] = density_porosity(out["RHOB"], p["matrix"], p["rho_fluid"])
    if "NPHI" in out:
        out["PHIN"] = pu_to_frac(out["NPHI"])
    if "PHID" in out and "PHIN" in out:
        out["PHIND"] = (out["PHID"] + out["PHIN"]) / 2
    phi = out["PHIND"] if "PHIND" in out else out["PHID"] if "PHID" in out else None
    if phi is not None and "RT" in out:
        out["SW"] = archie_sw(out["RT"], phi, rw=p["rw"], a=p["a"], m=p["m"], n=p["n"])
    if phi is not None and "SW" in out and "VSH" in out:
        with np.errstate(invalid="ignore"):
            out["PAY"] = (phi > p["phi_cut"]) & (out["VSH"] < p["vsh_cut"]) & (out["SW"] < p["sw_cut"])
        out["PAY"] = out["PAY"].fillna(False).astype(bool)
    return out


def default_tracks(columns) -> List[dict]:
    """Track spec using only the curves present (standard names)."""
    cols = set(columns)
    tracks: List[dict] = []
    if "GR" in cols:
        tracks.append({"curves": ["GR"], "xlim": (0, 175), "fill": "left", "title": "GR"})
    res = [c for c in ("RT", "RM", "RXO") if c in cols]
    if res:
        tracks.append({"curves": res, "xlim": (0.5, 2000), "log": True, "title": "Res [ohm-m]"})
    if "RHOB" in cols and "NPHI" in cols:
        tracks.append({"curves": ["RHOB", "NPHI"], "xlim": [(1.95, 2.95), (45, -15)], "twin": True, "title": "RHOB / NPHI"})
    elif "RHOB" in cols:
        tracks.append({"curves": ["RHOB"], "xlim": (1.95, 2.95)})
    if "DT" in cols:
        tracks.append({"curves": ["DT"], "xlim": (140, 40), "title": "DT"})
    if "VSH" in cols:
        tracks.append({"curves": ["VSH"], "xlim": (0, 1), "title": "Vsh"})
    phis = [c for c in ("PHID", "PHIN") if c in cols]
    if phis:
        tracks.append({"curves": phis, "xlim": (0.45, -0.15), "title": "phi [frac]"})
    if "SW" in cols:
        tracks.append({"curves": ["SW"], "xlim": (1, 0), "fill": "right", "title": "Sw"})
    return tracks


def _well_name(las) -> str:
    try:
        return str(las.well["WELL"].value).strip()
    except (KeyError, AttributeError):
        return ""


def _step(index) -> float:
    d = np.diff(np.asarray(index, dtype=float))
    d = d[np.isfinite(d) & (d > 0)]
    return float(np.median(d)) if d.size else float("nan")


def summarize(df: pd.DataFrame, params: Optional[Dict[str, float]] = None, meta: Optional[dict] = None) -> dict:
    """One-row summary of an :func:`analyze`-d frame (whole frame; slice first for an interval)."""
    p = {**DEFAULT_PARAMS, **(params or {})}
    step = _step(df.index)
    row: Dict[str, object] = dict(meta or {})
    row.update(
        {
            "depth_top": float(df.index.min()) if len(df) else np.nan,
            "depth_base": float(df.index.max()) if len(df) else np.nan,
            "step": step,
            "n_samples": int(len(df)),
            "curves": ",".join(c for c in ("GR", "SP", "RT", "RM", "RXO", "RHOB", "NPHI", "DPHI", "DT", "CALI") if c in df),
        }
    )
    for col, name in (("VSH", "vsh_mean"), ("PHIND", "phind_mean"), ("PHID", "phid_mean"), ("SW", "sw_mean")):
        row[name] = float(df[col].mean()) if col in df and df[col].notna().any() else np.nan
    if "PAY" in df:
        n_pay = int(df["PAY"].sum())
        row["pay_ft"] = n_pay * step if np.isfinite(step) else np.nan
        pay = df.loc[df["PAY"]]
        row["pay_top"] = float(pay.index.min()) if n_pay else np.nan
        row["pay_base"] = float(pay.index.max()) if n_pay else np.nan
        row["pay_sw_mean"] = float(pay["SW"].mean()) if n_pay else np.nan
    else:
        row.update({"pay_ft": np.nan, "pay_top": np.nan, "pay_base": np.nan, "pay_sw_mean": np.nan})
    row["params"] = json.dumps({k: p[k] for k in ("matrix", "rw", "a", "m", "n", "gr_clean", "gr_dirty")})
    return row


def run_well(
    las_path,
    out_dir,
    params: Optional[Dict[str, float]] = None,
    depth_range: Optional[tuple] = None,
    meta: Optional[dict] = None,
    dpi: int = 110,
    plot: bool = True,
    html: bool = False,
) -> dict:
    """Load one LAS, analyse it, save ``<stem>.png`` (and ``.csv`` of derived curves), return the summary row.

    ``html=True`` also writes ``<stem>.html``, the interactive viewer (see :mod:`lasanalysis.viewer`).
    """
    las_path = Path(las_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    las = read_las(las_path)
    df = analyze(standardize(curves(las)), params)
    if depth_range is not None:
        df = df.loc[depth_range[0] : depth_range[1]]
    meta = {"file": las_path.name, "well_name": _well_name(las), **(meta or {})}
    row = summarize(df, params, meta)
    row["mask_report"] = json.dumps(las.mask_report)
    stem = las_path.stem
    df.to_csv(out_dir / f"{stem}_derived.csv", float_format="%.4f")
    row["png"] = ""
    tracks = default_tracks(df.columns)
    if plot and tracks and len(df):
        fig = plot_tracks(df, tracks, depth_range=depth_range, dpi=dpi, title=f"{meta['well_name'] or stem} ({stem})")
        png = out_dir / f"{stem}.png"
        fig.savefig(png)
        plt.close(fig)
        row["png"] = png.name
    row["html"] = ""
    if html:
        from .viewer import _well_meta, write_viewer

        vmeta = {**_well_meta(las), "file": las_path.name, **{k: str(v) for k, v in (meta or {}).items() if v not in (None, "")}}
        out = write_viewer(df.drop(columns=[c for c in ("VSH", "PHID", "PHIN", "PHIND", "SW", "PAY") if c in df]),
                           out_dir / f"{stem}.html", params, depth_range, meta=vmeta)
        row["html"] = out.name
    return row


def run_search(
    search_kwargs: dict,
    out_dir,
    cache_dir="data/cache",
    params: Optional[Dict[str, float]] = None,
    depth_range: Optional[tuple] = None,
    search: Callable[..., List[dict]] = kgs.search_wells,
    fetch: Callable[..., Path] = kgs.fetch_las,
    log: Callable[[str], None] = print,
    **run_kwargs,
) -> pd.DataFrame:
    """Search KGS, fetch every LAS hit, run each, write ``summary.csv``.

    A well that fails to download or parse gets a row with an ``error`` column
    instead of aborting the batch.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    hits = search(**search_kwargs)
    log(f"{len(hits)} LAS files match {search_kwargs}")
    rows = []
    for h in hits:
        meta = {k: h.get(k) for k in ("kid", "well", "api", "operator", "location", "las_url")}
        try:
            path = fetch(h["kid"], cache_dir, url=h.get("las_url"))
            log(f"  {h['kid']} {h.get('well', '')}: {path.name}")
            row = run_well(path, out_dir, params, depth_range, meta, **run_kwargs)
            row["error"] = ""
        except Exception as e:  # noqa: BLE001 - keep the batch going, record why
            log(f"  {h['kid']} {h.get('well', '')}: FAILED {type(e).__name__}: {e}")
            row = {**meta, "error": f"{type(e).__name__}: {e}"}
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary.to_csv(out_dir / "summary.csv", index=False)
    return summary


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_argument_group("source (pick KGS search filters or local files)")
    src.add_argument("--township", type=int)
    src.add_argument("--range", type=int, dest="range_")
    src.add_argument("--ew", choices=["E", "W"])
    src.add_argument("--section", type=int)
    src.add_argument("--lease", default="")
    src.add_argument("--operator", default="")
    src.add_argument("--county", default="")
    src.add_argument("--api", default="")
    src.add_argument("--las", nargs="*", help="local LAS files instead of a KGS search")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--cache", default="data/cache", help="where fetched LAS files go")
    ap.add_argument("--depth", nargs=2, type=float, metavar=("TOP", "BASE"))
    ap.add_argument("--param", action="append", default=[], metavar="KEY=VALUE", help=f"override any of {sorted(DEFAULT_PARAMS)}")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--html", action="store_true", help="also write an interactive viewer per well")
    args = ap.parse_args(argv)

    try:
        params = parse_params(args.param)
    except ValueError as e:
        ap.error(str(e))
    depth = tuple(args.depth) if args.depth else None
    matplotlib.use("Agg")

    if args.las:
        rows = [run_well(p, args.out, params, depth, plot=not args.no_plot, html=args.html) for p in args.las]
        summary = pd.DataFrame(rows)
        summary.to_csv(Path(args.out) / "summary.csv", index=False)
    else:
        search_kwargs = {
            k: v
            for k, v in dict(
                township=args.township, range_=args.range_, ew=args.ew, section=args.section,
                lease=args.lease, operator=args.operator, county=args.county, api=args.api,
            ).items()
            if v not in (None, "")
        }
        if not search_kwargs:
            ap.error("give KGS search filters or --las files")
        summary = run_search(search_kwargs, args.out, args.cache, params, depth, plot=not args.no_plot, html=args.html)

    cols = [c for c in ("kid", "well", "api", "depth_top", "depth_base", "phind_mean", "sw_mean", "pay_ft", "error") if c in summary]
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(summary[cols].to_string(index=False))
    print(f"wrote {Path(args.out) / 'summary.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
