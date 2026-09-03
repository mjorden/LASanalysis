"""Log plots built on the matplotlib object API only (no pyplot state).

``plot_tracks`` takes a declarative list of tracks and draws them side by side
on one shared, inverted depth axis. The depth axis is set exactly once, so
there is no ``invert_yaxis`` juggling and twin tracks cannot flip it back.
"""

from __future__ import annotations

from typing import Iterable, List, Optional, Sequence, Union

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .petro import archie_sw_lines, matrix_density

_DEFAULT_COLORS = ("tab:green", "tab:blue", "tab:red", "tab:orange", "tab:purple", "tab:brown")

TrackSpec = dict
"""One track. Keys (all optional except ``curves``):

curves   list of column names to draw
xlim     (lo, hi) shared by all curves, or a list of (lo, hi) per curve
log      True for a log x-axis
colors   list of matplotlib colours, one per curve
twin     True to give each curve after the first its own x-axis (twiny)
title    track title; defaults to the curve names joined with " / "
fill     "left" or "right": shade between the first curve and that edge
"""


def _as_df(data) -> pd.DataFrame:
    if isinstance(data, pd.DataFrame):
        return data
    if hasattr(data, "df"):  # lasio.LASFile
        from .load import curves

        return curves(data)
    raise TypeError("data must be a pandas DataFrame or a lasio.LASFile")


def _per_curve(value, n: int, name: str) -> list:
    if value is None:
        return [None] * n
    if isinstance(value, (list, tuple)) and value and isinstance(value[0], (list, tuple)):
        if len(value) != n:
            raise ValueError(f"{name}: expected {n} entries, got {len(value)}")
        return list(value)
    return [value] * n


def plot_tracks(
    data,
    tracks: Sequence[TrackSpec],
    depth_range: Optional[tuple] = None,
    figsize: Optional[tuple] = None,
    dpi: int = 150,
    title: Optional[str] = None,
    track_width: float = 2.5,
    height: float = 10.0,
) -> matplotlib.figure.Figure:
    """Draw ``tracks`` side by side against depth.

    ``data`` is a DataFrame indexed by depth or a ``lasio.LASFile``.
    ``depth_range=(top, bottom)`` in depth units (top < bottom); default is the
    full log. Depth increases downward on every track.
    """
    df = _as_df(data)
    depth = df.index.to_numpy(dtype=float)
    n = len(tracks)
    if n == 0:
        raise ValueError("no tracks given")
    if figsize is None:
        figsize = (track_width * n, height)

    fig, axes = plt.subplots(1, n, figsize=figsize, dpi=dpi, sharey=True, squeeze=False)
    axes = axes[0]

    if depth_range is None:
        top, bottom = float(np.nanmin(depth)), float(np.nanmax(depth))
    else:
        top, bottom = depth_range
        if top >= bottom:
            raise ValueError("depth_range must be (top, bottom) with top < bottom")
    # The one and only place the depth axis is set. sharey propagates it.
    axes[0].set_ylim(bottom, top)
    axes[0].set_ylabel(f"Depth [{df.index.name or 'DEPT'}]")

    for ax, spec in zip(axes, tracks):
        curves: List[str] = list(spec["curves"])
        missing = [c for c in curves if c not in df.columns]
        if missing:
            raise KeyError(f"curves not in data: {missing}")
        k = len(curves)
        xlims = _per_curve(spec.get("xlim"), k, "xlim")
        colors = spec.get("colors") or _DEFAULT_COLORS[:k]
        twin = bool(spec.get("twin", False))
        log = bool(spec.get("log", False))

        ax.set_title(spec.get("title", " / ".join(curves)))
        ax.grid(True, which="both", linewidth=0.3, alpha=0.6)
        if log:
            ax.set_xscale("log")

        first_ax = ax
        for i, (name, color) in enumerate(zip(curves, colors)):
            target = ax
            if twin and i > 0:
                target = first_ax.twiny()
                if log:
                    target.set_xscale("log")
                target.spines["top"].set_position(("axes", 1.0 + 0.08 * (i - 1)))
                target.tick_params(axis="x", colors=color)
            target.plot(df[name].to_numpy(dtype=float), depth, color=color, linewidth=0.6, label=name)
            if xlims[i] is not None:
                target.set_xlim(*xlims[i])
            if twin and i > 0:
                target.set_xlabel(name, color=color)

        fill = spec.get("fill")
        if fill in ("left", "right"):
            x = df[curves[0]].to_numpy(dtype=float)
            lo, hi = ax.get_xlim()
            edge = lo if fill == "left" else hi
            ax.fill_betweenx(depth, x, edge, where=~np.isnan(x), color=colors[0], alpha=0.15, linewidth=0)
        if not twin and k > 1:
            ax.legend(fontsize=7, loc="upper right")

    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig


def crossplot_neutron_density(
    nphi_pu,
    rhob,
    color_by=None,
    rho_fluid: float = 1.0,
    matrices: Iterable[str] = ("sandstone", "limestone", "dolomite"),
    ax: Optional[matplotlib.axes.Axes] = None,
    **scatter_kw,
) -> matplotlib.axes.Axes:
    """Neutron porosity (PU) vs bulk density with approximate matrix lines.

    The matrix lines join (0 pu, rho_ma) to (100 pu, rho_fluid); they ignore the
    neutron tool's lithology response and are a guide, not a chart.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6), dpi=120)
    scatter_kw.setdefault("s", 3)
    scatter_kw.setdefault("alpha", 0.5)
    sc = ax.scatter(np.asarray(nphi_pu, float), np.asarray(rhob, float), c=color_by, **scatter_kw)
    if color_by is not None:
        ax.figure.colorbar(sc, ax=ax, label=getattr(color_by, "name", None) or "")
    for name in matrices:
        rho_ma = matrix_density(name)
        ax.plot([0, 100], [rho_ma, rho_fluid], linewidth=0.8, label=f"{name} ({rho_ma:.2f})")
    ax.set_xlim(-5, 60)
    ax.set_ylim(3.0, 1.8)
    ax.set_xlabel("Neutron porosity [pu]")
    ax.set_ylabel("Bulk density [g/cc]")
    ax.grid(True, linewidth=0.3, alpha=0.6)
    ax.legend(fontsize=8)
    return ax


def pickett_plot(
    rt,
    phi,
    rw: float,
    a: float = 1.0,
    m: float = 2.0,
    n: float = 2.0,
    sw_lines: Sequence[float] = (1.0, 0.5, 0.25, 0.1),
    color_by=None,
    ax: Optional[matplotlib.axes.Axes] = None,
    **scatter_kw,
) -> matplotlib.axes.Axes:
    """log Rt vs log phi with Archie iso-Sw lines for the given Rw, a, m, n.

    Points falling on the Sw=1 line are water-bearing; a straight line through
    them has slope -m and intercept a*Rw at phi=1.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6), dpi=120)
    rt = np.asarray(rt, float)
    phi = np.asarray(phi, float)
    ok = (rt > 0) & (phi > 0)
    scatter_kw.setdefault("s", 3)
    scatter_kw.setdefault("alpha", 0.5)
    c = None if color_by is None else np.asarray(color_by)[ok]
    sc = ax.scatter(phi[ok], rt[ok], c=c, **scatter_kw)
    if color_by is not None:
        ax.figure.colorbar(sc, ax=ax, label=getattr(color_by, "name", None) or "")
    phi_line = np.logspace(-2, 0, 50)
    for sw, rt_line in archie_sw_lines(phi_line, rw, sw_lines, a=a, m=m, n=n).items():
        ax.plot(phi_line, rt_line, linewidth=0.8, label=f"Sw = {sw:g}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.01, 1.0)
    ax.set_xlabel("Porosity [frac]")
    ax.set_ylabel("Rt [ohm-m]")
    ax.set_title(f"Pickett plot (Rw={rw:g}, a={a:g}, m={m:g}, n={n:g})")
    ax.grid(True, which="both", linewidth=0.3, alpha=0.6)
    ax.legend(fontsize=8)
    return ax
