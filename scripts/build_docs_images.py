"""Render the figures embedded in README.md and docs/ into docs/images/.

    python scripts/build_docs_images.py

Every image is a real output of the package on the two wells in data/, so the
documentation cannot drift from what the code draws. The viewer screenshot
needs Playwright + Chromium (``python -m playwright install chromium``) and
network access for Plotly.js; it is skipped with a message when unavailable.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lasanalysis import (
    crossplot_neutron_density,
    curves,
    fit_water_line,
    kgs,
    pick_rw_from_rwa,
    pickett_plot,
    plot_tracks,
    read_las,
    read_log_csv,
    standardize,
)
from lasanalysis.multiwell import analyze, plot_wells
from lasanalysis.viewer import write_viewer

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "images"
PEARSON = ROOT / "data" / "1046139243.las"
PBW = ROOT / "data" / "1045399712.csv"
DPI = 110


def pearson_frame() -> pd.DataFrame:
    df = analyze(standardize(curves(read_las(PEARSON))), {"rw": 0.03, "m": 2.0})
    return df


def tracks(df: pd.DataFrame, path: Path, title: str, depth_range: tuple) -> None:
    spec = [
        {"curves": ["GR"], "xlim": (0, 175), "fill": "left", "title": "GR [GAPI]"},
        {"curves": ["RT", "RM", "RXO"], "xlim": (1, 1000), "log": True, "title": "Resistivity [ohm-m]"},
        {"curves": ["RHOB", "NPHI"], "xlim": [(1.95, 2.95), (45, -15)], "twin": True, "title": "RHOB / NPHI"},
        {"curves": ["VSH"], "xlim": (0, 1), "title": "Vsh"},
        {"curves": ["PHID", "PHIN"], "xlim": (0.45, -0.15), "title": "phi [frac]"},
        {"curves": ["SW"], "xlim": (1, 0), "fill": "right", "title": "Sw (Archie)"},
    ]
    fig = plot_tracks(df, spec, depth_range=depth_range, dpi=DPI, height=8.5, track_width=2.1, title=title)
    fig.savefig(path)
    plt.close(fig)


def pickett(df: pd.DataFrame, path: Path) -> None:
    sel = df.loc[3400:4200]
    clean = sel[sel["VSH"] < 0.15]
    fit = fit_water_line(clean["RT"], clean["PHIND"], phi_min=0.06)
    rwa = pick_rw_from_rwa(sel["RT"], sel["PHIND"], sel["VSH"], depth=sel.index, m=2.0)
    ax = pickett_plot(clean["RT"], clean["PHIND"], rw=round(fit["rw"], 3), m=round(fit["m"], 1),
                      color_by=clean.index.to_series().rename("depth [ft]"), s=6)
    env = fit["envelope"]
    ax.scatter(10 ** env[:, 0], 10 ** env[:, 1], marker="x", color="k", s=45, zorder=5, label="5th-pct envelope")
    ax.set_xlim(0.03, 0.4)
    ax.set_ylim(0.5, 500)
    ax.set_title(f"Pearson #1-35, 3400–4200 ft, Vsh < 0.15\nenvelope: m = {fit['m']:.2f}, a·Rw = {fit['rw']:.3f}   |   Rwa pick: Rw = {rwa['rw']:.3f} ({rwa['n_points']} samples)", fontsize=10)
    ax.legend(fontsize=8)
    ax.figure.set_size_inches(7, 6)
    ax.figure.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(ax.figure)


def crossplot(df: pd.DataFrame, path: Path) -> None:
    sel = df.loc[3400:4200]
    ax = crossplot_neutron_density(sel["NPHI"], sel["RHOB"], color_by=sel["GR"].rename("GR [GAPI]"), s=6)
    ax.set_title("Pearson #1-35, 3400–4200 ft: neutron–density, coloured by GR")
    ax.figure.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(ax.figure)


def wells_map(path: Path) -> None:
    fx = ROOT / "tests" / "fixtures"
    rows = []
    for page, kid, well, pay in (("displaywell_1046105344.html", 1046139243, "Pearson Family #1-35", 7.5),
                                 ("displaywell_1045079321.html", 1045399712, "PBW #1-32", 120.5)):
        w = kgs.parse_well_page((fx / page).read_text(encoding="utf-8"))
        rows.append({"kid": kid, "well": well, "lat": w["lat"], "lon": w["lon"], "pay_ft": pay})
    plot_wells(pd.DataFrame(rows), path, dpi=DPI)


def viewer_screenshot(path: Path) -> bool:
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        print("playwright not installed; skipping viewer screenshot")
        return False
    html = write_viewer(PEARSON, OUT / "_viewer_tmp.html", params={"rw": 0.03, "m": 2.0}, depth_range=(3400, 4200))
    try:
        with sync_playwright() as p:
            try:
                b = p.chromium.launch()
            except Exception as e:  # noqa: BLE001
                print(f"chromium not available ({e}); skipping viewer screenshot")
                return False
            page = b.new_page(viewport={"width": 1400, "height": 900}, device_scale_factor=1)
            page.goto(html.resolve().as_uri())
            page.wait_for_function("() => window.Plotly && document.getElementById('logs').data && document.getElementById('pickett').data", timeout=60000)
            page.wait_for_timeout(1000)
            page.screenshot(path=str(path))
            b.close()
    finally:
        html.unlink(missing_ok=True)
    return True


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    df = pearson_frame()
    tracks(df, OUT / "tracks_pearson.png", "Pearson Family #1-35 — quick-look (limestone matrix, Rw 0.03, m 2)", (3400, 4200))
    pbw = analyze(standardize(read_log_csv(PBW)), {"rw": 0.06, "m": 2.0})
    tracks(pbw, OUT / "tracks_pbw.png", "PBW #1-32 — quick-look (limestone matrix, Rw 0.06, m 2)", (3000, 3850))
    pickett(df, OUT / "pickett_pearson.png")
    crossplot(df, OUT / "crossplot_pearson.png")
    wells_map(OUT / "wells_map.png")
    viewer_screenshot(OUT / "viewer_pearson.png")
    for p in sorted(OUT.glob("*.png")):
        print(f"{p.relative_to(ROOT)}  {p.stat().st_size // 1024} KB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
