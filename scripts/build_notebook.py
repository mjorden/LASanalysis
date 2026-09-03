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
    archie_sw, fit_water_line, pick_rw_from_rwa, crossplot_neutron_density, pickett_plot, MATRIX_DENSITY,
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
## 5. Picking Rw and m from the Pickett plot

No water analysis is available for this well, so `Rw` and the cementation
exponent `m` come from the log itself. On a Pickett plot (log Rt vs log phi)
water-bearing rock falls on a straight line with slope `-m` and intercept
`a * Rw` at phi = 1; everything hydrocarbon-bearing plots above it.
`fit_water_line` takes the low-Rt envelope of the clean points (5th percentile
of Rt in each porosity bin) and fits that line.

Two choices matter and are shown below:

* **porosity source** — density porosity alone vs the neutron-density average
  (`PHIND`), the usual choice in carbonates;
* **low-porosity cutoff** — below ~6 % porosity, shale conductivity and matrix
  error flatten the envelope and drag the apparent `m` toward 1.3.
"""
)

code(
    """
df["PHIND"] = (df["PHID"] + df["PHIN"]) / 2

sel = df.loc[DEPTH_RANGE[0]:DEPTH_RANGE[1]]
clean = sel[sel["VSH_LAR"] < 0.15]

rows = []
for phicol in ("PHID", "PHIND"):
    for phi_min in (0.03, 0.06, 0.08):
        f = fit_water_line(clean["RT"], clean[phicol], phi_min=phi_min)
        rows.append({"porosity": phicol, "phi_min": phi_min, "m": round(f["m"], 2), "a*Rw": round(f["rw"], 4), "n": f["n_points"]})
pd.DataFrame(rows)
"""
)

md(
    """
With `PHIND` and a 6 % cutoff the fit gives `m ≈ 2.0` and `a·Rw ≈ 0.03`, stable
across Vsh cutoffs (0.10–0.25 all give 1.93–1.97). The cleanest wet zones
independently agree: 3580–3650 ft has Vsh < 0.15, phi ≈ 16 % and
Rt ≈ 0.7–1.5 ohm-m, which at `m = 2` implies Rw ≈ 0.02–0.04. A 100,000 ppm
NaCl brine at ~110 °F (roughly 4000 ft in central Kansas) has Rw ≈ 0.03, so
the pick is physically reasonable. `n = 2` and `a = 1` are assumed, not fitted.
"""
)

code(
    """
FIT = fit_water_line(clean["RT"], clean["PHIND"], phi_min=0.06)
RW, A, M, N = round(FIT["rw"], 3), 1.0, round(FIT["m"], 1), 2.0
print(f"pick: Rw = {RW}, a = {A}, m = {M}, n = {N}   (fit: m={FIT['m']:.2f}, a*Rw={FIT['rw']:.4f} from {FIT['n_points']} points)")

ax = pickett_plot(clean["RT"], clean["PHIND"], rw=RW, a=A, m=M, n=N, color_by=clean.index.to_series().rename("depth"))
env = FIT["envelope"]
ax.scatter(10 ** env[:, 0], 10 ** env[:, 1], marker="x", color="k", s=40, zorder=5, label="5th-pct envelope")
ax.set_xlim(0.02, 0.5); ax.set_ylim(0.5, 1000); ax.legend(fontsize=8)
"""
)

md(
    """
A second, independent pick: the apparent water resistivity `Rwa = Rt·phi^m / a`
equals Rw in water-bearing rock and exceeds it elsewhere, so a low percentile
of Rwa over clean, porous samples is an Rw estimate that needs no straight
Pickett line (it assumes `m` instead). The two should agree here; on a mixed
section such as PBW #1-32 the envelope fit wanders and Rwa is the one to trust.
"""
)

code(
    """
r = pick_rw_from_rwa(sel["RT"], sel["PHIND"], sel["VSH_LAR"], depth=sel.index, m=M, a=A)
print(f"Rwa pick: Rw = {r['rw']:.3f} from {r['n_points']} clean samples; wettest interval {r['interval'][0]:.0f}-{r['interval'][1]:.0f} ft")
print(f"envelope: Rw = {RW}, m = {M}")
"""
)

md(
    """
## 6. Archie water saturation

Archie ignores shale conductivity, so it overstates Sw in shaly intervals; the
pay flag later copes with a Vsh cutoff. `water_saturation(model, ...)` also
offers modified Simandoux and Indonesia (Poupon–Leveaux), which need a shale
resistivity (`pick_rsh` takes the median Rt of the shales). The viewer has the
same switch.
"""
)

code(
    """
df["SW"] = archie_sw(df["RT"], df["PHIND"], rw=RW, a=A, m=M, n=N)

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
    title=f"Pearson #1-35 quick-look (matrix={MATRIX}, Rw={RW}, m={M})",
)
"""
)

md(
    """
## 7. Second well from the same operator

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
## 8. Many wells at once

`lasanalysis.multiwell` runs this whole workflow over a batch: search the KGS
index, fetch each LAS, analyse it with the parameters picked above, save a
track plot, and write `summary.csv` (porosity, Sw, pay feet per well). From a
shell:

```
python -m lasanalysis.multiwell --township 13 --range 22 --ew W --section 35 --out output/T13S_R22W_35
python -m lasanalysis.multiwell --las data/1046139243.las --out output/pearson --depth 3400 4200
```

Network access is deliberately not exercised in this notebook so it can run
in CI; the cell below runs the local-file path.
"""
)

code(
    """
from lasanalysis.multiwell import run_well

row = run_well("data/1046139243.las", "output/pearson", params={"rw": RW, "m": M}, depth_range=DEPTH_RANGE, plot=False)
{k: row[k] for k in ("well_name", "depth_top", "depth_base", "phind_mean", "sw_mean", "pay_ft", "pay_top", "pay_base")}
"""
)

md(
    """
## 9. Interactive viewer

`lasanalysis.viewer` writes a self-contained HTML page: the same tracks on one
zoomable depth axis, with Rw / a / m / n, the matrix density, the GR picks and
the pay cutoffs as live sliders — Vsh, porosity, Sw, pay shading and a Pickett
panel recompute in the browser. Open the file below in a browser tab (it is
too dynamic to embed usefully here).
"""
)

code(
    """
from lasanalysis.viewer import write_viewer

write_viewer("data/1046139243.las", "output/pearson.html", params={"rw": RW, "m": M}, depth_range=DEPTH_RANGE)
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
