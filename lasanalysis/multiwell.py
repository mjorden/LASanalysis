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
import warnings

from .petro import (
    NEUTRON_MATRIX_OFFSET,
    archie_sw,
    density_porosity,
    fit_water_line,
    matrix_density,
    neutron_lithology_correction,
    pick_rsh,
    pick_rw_from_rwa,
    pu_to_frac,
    vshale_larionov,
    water_saturation,
)
from .plot import plot_tracks

#: Petrophysical parameters. Rw / m are the Pearson #1-35 picks (see notebook
#: section 5); override per area.
DEFAULT_PARAMS: Dict[str, object] = {
    "gr_clean": 20.0,
    "gr_dirty": 110.0,
    "matrix": "limestone",
    "neutron_matrix": "limestone",  # lithology the neutron curve was scaled on (#26)
    "rho_fluid": 1.0,
    "rw": 0.03,
    "a": 1.0,
    "m": 2.0,
    "n": 2.0,
    "sw_model": "archie",           # archie | simandoux | indonesia (#28)
    "rsh": float("nan"),            # shale resistivity for the shaly models; NaN = pick from the log
    # pay flags: porous, clean, hydrocarbon-bearing
    "phi_cut": 0.08,
    "vsh_cut": 0.30,
    "sw_cut": 0.50,
}

#: Parameters that take a name rather than a number, with their validator.
_STRING_PARAMS = {
    "matrix": lambda v: (matrix_density(v), v.strip().lower())[1],
    "neutron_matrix": lambda v: (neutron_lithology_correction(0.0, v, v), v.strip().lower())[1],
    "sw_model": lambda v: (water_saturation(v, 1.0, 0.1, 0.0, 0.05, rsh=1.0), v.strip().lower())[1],
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
        if k in _STRING_PARAMS:
            try:
                params[k] = float(v) if k == "matrix" else _STRING_PARAMS[k](v)
            except ValueError:
                params[k] = _STRING_PARAMS[k](v)  # raises with the valid names
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
    out.attrs = dict(df.attrs)
    if "GR" in out:
        out["VSH"] = vshale_larionov(out["GR"], p["gr_clean"], p["gr_dirty"], older=True)
    if "RHOB" in out:
        out["PHID"] = density_porosity(out["RHOB"], p["matrix"], p["rho_fluid"])
    if "NPHI" in out:
        phin = pu_to_frac(out["NPHI"])
        # #26: the neutron curve is scaled on the logging matrix; bring it onto the
        # density matrix before averaging. Only meaningful for the three chart
        # lithologies — otherwise leave it and say so.
        corrected = False
        if isinstance(p["matrix"], str) and p["matrix"] in NEUTRON_MATRIX_OFFSET:
            phin = neutron_lithology_correction(phin, p["neutron_matrix"], p["matrix"])
            corrected = True
        elif str(p["neutron_matrix"]).lower() != "limestone" or not isinstance(p["matrix"], str) or p["matrix"] not in NEUTRON_MATRIX_OFFSET:
            warnings.warn(
                f"neutron porosity left on the {p['neutron_matrix']} scale: no lithology correction to matrix {p['matrix']!r}",
                stacklevel=2,
            )
        out["PHIN"] = phin
        out.attrs["phin_corrected"] = corrected
    if "PHID" in out and "PHIN" in out:
        out["PHIND"] = (out["PHID"] + out["PHIN"]) / 2
    phi = out["PHIND"] if "PHIND" in out else out["PHID"] if "PHID" in out else None
    if phi is not None and "RT" in out:
        model = str(p["sw_model"]).lower()
        rsh = p["rsh"]
        if model != "archie" and "VSH" in out:
            if rsh is None or not np.isfinite(rsh):
                try:
                    rsh = pick_rsh(out["RT"], out["VSH"])
                except ValueError as e:
                    warnings.warn(f"{model}: {e}; falling back to Archie", stacklevel=2)
                    model = "archie"
            out.attrs["rsh"] = rsh
        elif model != "archie":
            warnings.warn(f"{model} needs a GR/Vsh curve; falling back to Archie", stacklevel=2)
            model = "archie"
        out.attrs["sw_model"] = model
        out["SW"] = water_saturation(model, out["RT"], phi, out["VSH"] if "VSH" in out else 0.0, p["rw"], rsh=rsh, a=p["a"], m=p["m"], n=p["n"])
        if model != "archie":
            out["SW_ARCHIE"] = archie_sw(out["RT"], phi, rw=p["rw"], a=p["a"], m=p["m"], n=p["n"])
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
    # #27: two independent Rw picks so a disagreement is visible per well.
    row.update({"rw_envelope": np.nan, "m_envelope": np.nan, "rw_rwa": np.nan, "rwa_interval": ""})
    phi = df["PHIND"] if "PHIND" in df else df["PHID"] if "PHID" in df else None
    if phi is not None and "RT" in df and "VSH" in df:
        clean = df[df["VSH"] < 0.15]
        try:
            f = fit_water_line(clean["RT"], phi[clean.index], phi_min=0.06)
            row["rw_envelope"], row["m_envelope"] = round(f["rw"], 4), round(f["m"], 2)
        except ValueError:
            pass
        try:
            r = pick_rw_from_rwa(df["RT"], phi, df["VSH"], depth=df.index, m=p["m"], a=p["a"])
            row["rw_rwa"] = round(r["rw"], 4)
            row["rwa_interval"] = "" if r["interval"] is None else f"{r['interval'][0]:.1f}-{r['interval'][1]:.1f}"
        except ValueError:
            pass
    row["sw_model"] = df.attrs.get("sw_model", str(p["sw_model"]))
    row["rsh"] = df.attrs.get("rsh", np.nan)
    row["phin_corrected"] = bool(df.attrs.get("phin_corrected", False))
    row["params"] = json.dumps({k: (None if isinstance(p[k], float) and np.isnan(p[k]) else p[k]) for k in ("matrix", "neutron_matrix", "rw", "a", "m", "n", "sw_model", "rsh", "gr_clean", "gr_dirty")})
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
        attrs = dict(df.attrs)
        df = df.loc[depth_range[0] : depth_range[1]]
        df.attrs = attrs
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
        out = write_viewer(df.drop(columns=[c for c in ("VSH", "PHID", "PHIN", "PHIND", "SW", "SW_ARCHIE", "PAY") if c in df]),
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
    coords: bool = True,
    well_info: Optional[Callable[..., List[dict]]] = None,
    map_png: Optional[str] = "wells.png",
    **run_kwargs,
) -> pd.DataFrame:
    """Search KGS, fetch every LAS hit, run each, write ``summary.csv``.

    With ``coords=True`` (default) each well's KGS page is read for NAD83
    ``lat`` / ``lon`` and header fields (#31), which land in the summary, and a
    ``wells.png`` location map is written when at least one well has
    coordinates. A well that fails to download or parse gets a row with an
    ``error`` column instead of aborting the batch.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    hits = search(**search_kwargs)
    log(f"{len(hits)} LAS files match {search_kwargs}")
    if coords:
        (well_info or kgs.add_well_info)(hits, cache_dir=cache_dir, log=log)
    rows = []
    for h in hits:
        meta = {k: h.get(k) for k in ("kid", "well", "api", "operator", "location", "las_url", "well_kid", *kgs.WELL_ROW_FIELDS) if k in h}
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
    if map_png and "lat" in summary and summary["lat"].notna().any():
        plot_wells(summary, out_dir / map_png)
        log(f"wrote {out_dir / map_png}")
    return summary


def plot_wells(summary: pd.DataFrame, out_png, label: str = "well", color_by: Optional[str] = "pay_ft", dpi: int = 120):
    """Simple location map of a batch: lon/lat scatter, labelled, optionally coloured by a summary column."""
    matplotlib.use("Agg")
    df = summary.dropna(subset=["lat", "lon"])
    if df.empty:
        raise ValueError("no wells with coordinates")
    fig, ax = plt.subplots(figsize=(7, 6), dpi=dpi)
    c = df[color_by] if color_by and color_by in df and df[color_by].notna().any() else None
    sc = ax.scatter(df["lon"], df["lat"], c=c, cmap="viridis", s=45, edgecolor="k", linewidth=0.5, zorder=3)
    if c is not None:
        fig.colorbar(sc, ax=ax, label=color_by, shrink=0.8)
    for _, r in df.iterrows():
        ax.annotate(str(r.get(label, "") or r.get("kid", "")), (r["lon"], r["lat"]), xytext=(4, 4), textcoords="offset points", fontsize=7)
    ax.set_xlabel("Longitude (NAD83)")
    ax.set_ylabel("Latitude (NAD83)")
    ax.set_aspect(1.0 / np.cos(np.deg2rad(df["lat"].mean())))
    ax.grid(True, linewidth=0.3, alpha=0.6)
    ax.set_title(f"{len(df)} wells")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)
    return Path(out_png)


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
    src.add_argument("--index", nargs="?", const="data/cache/ks_las_files.zip", metavar="ZIP",
                     help="search KGS's offline index (ks_las_files.zip) instead of scraping; downloads it to the given path if missing")
    src.add_argument("--within", nargs=3, type=float, metavar=("LAT", "LON", "KM"), help="with --index: wells within KM of LAT,LON")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--cache", default="data/cache", help="where fetched LAS files go")
    ap.add_argument("--depth", nargs=2, type=float, metavar=("TOP", "BASE"))
    ap.add_argument("--param", action="append", default=[], metavar="KEY=VALUE", help=f"override any of {sorted(DEFAULT_PARAMS)}")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--html", action="store_true", help="also write an interactive viewer per well")
    ap.add_argument("--no-coords", action="store_true", help="skip the per-well KGS page lookup (lat/lon, TD, status) and the wells.png map")
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
        if args.within:
            search_kwargs["within"] = tuple(args.within)
        if not search_kwargs:
            ap.error("give KGS search filters or --las files")
        search = kgs.search_wells
        if args.index is not None:
            index_df = kgs.load_index(kgs.fetch_index(args.index))
            search_kwargs.pop("county", None)  # the index has no county column

            def search(**kw):  # noqa: F811 - offline variant of kgs.search_wells
                return kgs.search_index(index_df, **kw)

        elif args.within:
            ap.error("--within needs --index")
        summary = run_search(search_kwargs, args.out, args.cache, params, depth, plot=not args.no_plot, html=args.html, coords=not args.no_coords, search=search)

    cols = [c for c in ("kid", "well", "api", "lat", "lon", "depth_top", "depth_base", "phind_mean", "sw_mean", "pay_ft", "rw_envelope", "rw_rwa", "error") if c in summary]
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(summary[cols].to_string(index=False))
    print(f"wrote {Path(args.out) / 'summary.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
