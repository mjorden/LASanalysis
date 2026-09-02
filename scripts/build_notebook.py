"""Regenerate KansasLAS.ipynb from source (outputs stripped).

Run:  python scripts/build_notebook.py
"""

from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]

cells = []


def md(text):
    cells.append(nbf.v4.new_markdown_cell(text.strip()))


def code(text):
    cells.append(nbf.v4.new_code_cell(text.strip()))


md(
    """
# Quick-look petrophysics: Pearson Family #1-35, Trego County, Kansas

Kansas Geological Survey LAS file KID **1046139243** (API 15-195-23011), logged
2016-11-21 by Casedhole Solutions for Downing Nelson Oil Company. Curves: GR,
DIL resistivity (deep / medium / shallow), SP, density (RHOB, DPOR), compensated
neutron (CNPOR), sonic (DT, SPOR), calipers and borehole volumes.

Everything here is driven by the `lasanalysis` package in this repo:
`read_las` cleans the sentinels the service company never declared,
`plot_tracks` draws the log, and `petro` holds the (pure numpy) calculations.
"""
)

code(
    """
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lasanalysis import (
    read_las, curves, standardize, plot_tracks,
    vshale_linear, vshale_larionov, density_porosity, neutron_density_crossover,
    archie_sw, crossplot_neutron_density, pickett_plot, MATRIX_DENSITY,
)
from lasanalysis.petro import pu_to_frac, frac_to_pu

%matplotlib inline
"""
)

md(
    """
## 1. Load and clean

lasio masks the declared `NULL` (`-999.25`). It cannot know about the
`100000.0` off-scale values in RILD/RILM, the `-999.0` and `0.0` in RHOB, or the
DPOR values in the tens of thousands of pu that those RHOB rows produce.
`read_las` masks all of them and reports what it touched.
"""
)

code(
    """
las = read_las("data/1046139243.las")
for curve in las.curves:
    print(f"{curve.mnemonic:8s} [{curve.unit:9s}] {curve.descr}")
print()
print("masked:", las.mask_report)
"""
)

code(
    """
df = standardize(curves(las))          # RILD -> RT, CNPOR -> NPHI, ... originals otherwise
df.describe().T[["count", "min", "50%", "max"]]
"""
)

md(
    """
## 2. Raw log

One shared depth axis, set once. Twin tracks get their own x-axis; the neutron
porosity axis is reversed (45 to -15 pu) so the density and neutron curves
overlay in the conventional way and gas crossover shows as density porosity
to the right of neutron.
"""
)

code(
    """
DEPTH_RANGE = (3400, 4200)   # ft, the interval of interest in this well

tracks = [
    {"curves": ["GR"],           "xlim": (0, 175),   "fill": "left", "title": "GR [GAPI]"},
    {"curves": ["RT", "RM", "RXO"], "xlim": (1, 1000), "log": True, "title": "Resistivity [ohm-m]"},
    {"curves": ["SP"],           "xlim": (-200, 140), "title": "SP [mV]"},
    {"curves": ["DT"],           "xlim": (140, 40),   "title": "DT [us/ft]"},
    {"curves": ["RHOB", "NPHI"], "xlim": [(1.95, 2.95), (45, -15)], "twin": True, "title": "RHOB / NPHI"},
    {"curves": ["DPHI", "NPHI"], "xlim": (45, -15), "title": "DPOR / CNPOR [pu]"},
]
fig = plot_tracks(df, tracks, depth_range=DEPTH_RANGE, title="Pearson Family #1-35")
"""
)

md(
    """
## 3. Shale volume from gamma ray

Clean and dirty GR picks are the ones the original notebook staged but never
used (20 / 110 GAPI). Larionov's older-rock form is appropriate for Kansas
Paleozoic section; the linear index is shown for comparison.
"""
)

code(
    """
GR_CLEAN, GR_DIRTY = 20, 110

df["VSH_LIN"] = vshale_linear(df["GR"], GR_CLEAN, GR_DIRTY)
df["VSH_LAR"] = vshale_larionov(df["GR"], GR_CLEAN, GR_DIRTY, older=True)
df[["VSH_LIN", "VSH_LAR"]].loc[DEPTH_RANGE[0]:DEPTH_RANGE[1]].describe().T
"""
)

md(
    """
## 4. Density porosity and neutron-density crossover

The matrix density is a parameter. Limestone (2.71 g/cc) is the working
assumption for this section; the other keys of `MATRIX_DENSITY` are one edit
away. Porosities below are **fractions**; the LAS porosity curves are in pu.
"""
)

code(
    """
MATRIX = "limestone"
RHO_FLUID = 1.0

df["PHID"] = density_porosity(df["RHOB"], MATRIX, RHO_FLUID)
df["PHIN"] = pu_to_frac(df["NPHI"])
sep, gas = neutron_density_crossover(df["PHIN"], df["PHID"], threshold=0.03)
df["ND_SEP"] = sep
df["GAS_FLAG"] = gas
print(f"matrix={MATRIX} ({MATRIX_DENSITY[MATRIX]} g/cc); "
      f"{int(gas[(df.index >= DEPTH_RANGE[0]) & (df.index <= DEPTH_RANGE[1])].sum())} samples flag crossover in the interval")
"""
)

code(
    """
ax = crossplot_neutron_density(df["NPHI"], df["RHOB"], color_by=df["GR"], rho_fluid=RHO_FLUID)
ax.set_title("Neutron-density crossplot, coloured by GR")
"""
)

md(
    """
## 5. Archie water saturation

`Rw` is not known for this well; the Pickett plot lets you pick it: water-bearing
points fall on the Sw = 1 line, whose intercept at phi = 1 is `a * Rw` and whose
slope is `-m`. Adjust `RW`, `M`, `N` until the wet limestone plots on the line,
then the Sw track is meaningful. The values below are a placeholder.
"""
)

code(
    """
RW, A, M, N = 0.05, 1.0, 2.0, 2.0

sel = df.loc[DEPTH_RANGE[0]:DEPTH_RANGE[1]]
ax = pickett_plot(sel["RT"], sel["PHID"], rw=RW, a=A, m=M, n=N, color_by=sel["GR"])
"""
)

code(
    """
df["SW"] = archie_sw(df["RT"], df["PHID"], rw=RW, a=A, m=M, n=N)

fig = plot_tracks(
    df,
    [
        {"curves": ["GR"], "xlim": (0, 175), "fill": "left"},
        {"curves": ["VSH_LAR", "VSH_LIN"], "xlim": (0, 1), "title": "Vsh"},
        {"curves": ["RT"], "xlim": (1, 1000), "log": True},
        {"curves": ["PHID", "PHIN"], "xlim": (0.45, -0.15), "title": "phi [frac]"},
        {"curves": ["SW"], "xlim": (1, 0), "fill": "right", "title": "Sw (Archie)"},
    ],
    depth_range=DEPTH_RANGE,
    title=f"Pearson #1-35 quick-look (matrix={MATRIX}, Rw={RW})",
)
"""
)

md(
    """
## 6. Second well from the same operator

`data/1045399712.csv` is PBW #1-32 (API 15-165-22116, Rush County, logged
2015-09-25). It arrived as a whitespace-delimited export with classic-Mac line
endings and has been converted to a real CSV; it has the same sentinels, and
`read_log_csv` applies the same rules by curve name since a CSV has no units.
"""
)

code(
    """
from lasanalysis import read_log_csv

pbw = standardize(read_log_csv("data/1045399712.csv"))
print("masked:", pbw.attrs["mask_report"])
fig = plot_tracks(pbw, tracks[:5], depth_range=(3000, 3850), title="PBW #1-32")
"""
)

md(
    """
## 7. More wells from KGS

`lasanalysis.kgs` fetches LAS files by KID and searches the KGS index by
township / range / section. See the README for a walk-through; network access
is deliberately not exercised in this notebook so it can run in CI.
"""
)

nb = nbf.v4.new_notebook()
nb["cells"] = cells
nb["metadata"] = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python"},
}
out = ROOT / "KansasLAS.ipynb"
nbf.write(nb, out)
print("wrote", out, len(cells), "cells")
