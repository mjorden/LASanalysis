# LASanalysis

Quick-look petrophysics for LAS well logs from the Kansas Geological Survey:
load a log, mask the sentinels the service company never declared, draw a
track plot on one shared depth axis, and compute Vshale, density porosity,
neutron–density crossover and Archie Sw. Includes a small client for the KGS
LAS-file index so the same workflow runs on any Kansas well.

## Setup

Python 3.10+ (verified on 3.12).

```bash
pip install -r requirements.txt -e .
jupyter notebook KansasLAS.ipynb
```

Run the tests with `pytest`. The notebook is regenerated from
`scripts/build_notebook.py` and committed with outputs stripped
(`pre-commit install` sets up `nbstripout`). CI executes the notebook end to
end on every push.

## What is in the box

| Path | Purpose |
|---|---|
| `lasanalysis/load.py` | `read_las` (lasio + sentinel masking, reports what it masked), `read_log_csv`, mnemonic aliases (`standardize`, `find_curve`) |
| `lasanalysis/petro.py` | Pure-numpy: `vshale_linear`, `vshale_larionov`, `density_porosity`, `neutron_density_crossover`, `archie_sw`, `MATRIX_DENSITY` |
| `lasanalysis/plot.py` | `plot_tracks` (declarative tracks, one inverted depth axis, twin x-axes), `crossplot_neutron_density`, `pickett_plot` |
| `lasanalysis/kgs.py` | `search_wells` (township / range / section / lease / operator), `fetch_las` (download a LAS by KID, cached) |
| `KansasLAS.ipynb` | The walk-through on the Pearson #1-35 well |
| `data/` | Two wells from KGS (see below) |

### Loading a log

```python
from lasanalysis import read_las, curves, standardize

las = read_las("data/1046139243.las")
print(las.mask_report)
# {'DPOR': {'porosity': 31}, 'RHOB': {'extra_null': 28, 'density': 3},
#  'RILD': {'off_scale': 268}, 'RILM': {'off_scale': 862}}
df = standardize(curves(las))   # RILD -> RT, CNPOR -> NPHI, DPOR -> DPHI, ...
```

`read_las` applies four rules, each chosen from the curve's unit:

| Rule | Trigger | Why |
|---|---|---|
| `extra_null` | value in `EXTRA_NULLS` (`-999.0`) | RHOB uses `-999.000`, not the declared `-999.25` |
| `off_scale` | resistivity (OHM unit) `>= 1e5` | RILD / RILM carry `100000.0` where the tool went off scale — 1,130 samples in the Pearson log |
| `porosity` | PU / % curve outside `[-50, 100]` | DPOR reaches 58,000 pu on the bad-RHOB rows |
| `density` | g/cc curve outside `[1.0, 3.5]` | RHOB has a few exact `0.0` |

Every threshold is a keyword argument; `clean=False` gives you lasio's
untouched output.

### Track plot

```python
from lasanalysis import plot_tracks

tracks = [
    {"curves": ["GR"],           "xlim": (0, 175), "fill": "left"},
    {"curves": ["RT", "RM"],     "xlim": (1, 1000), "log": True},
    {"curves": ["RHOB", "NPHI"], "xlim": [(1.95, 2.95), (45, -15)], "twin": True},
]
fig = plot_tracks(df, tracks, depth_range=(3400, 4200))
```

The depth axis is set exactly once (`set_ylim(bottom, top)` on the first
track; the rest share it). There is no `invert_yaxis` anywhere, so a twin
axis cannot flip it back.

### Petrophysics

All functions take arrays or Series, propagate NaN, and use **fractions** for
porosity. Convert the LAS porosity curves (pu) with `petro.pu_to_frac`.

```python
from lasanalysis import vshale_larionov, density_porosity, archie_sw

vsh  = vshale_larionov(df["GR"], gr_clean=20, gr_dirty=110, older=True)
phid = density_porosity(df["RHOB"], "limestone")        # or a number, g/cc
sw   = archie_sw(df["RT"], phid, rw=0.05, a=1, m=2, n=2)
```

Use `pickett_plot(rt, phi, rw, ...)` to pick `Rw` and `m`: water-bearing
points fall on the Sw = 1 line. The notebook's `RW = 0.05` is a placeholder,
not a calibrated value for this well.

### More wells from KGS

```python
from lasanalysis import kgs

rows = kgs.search_wells(township=13, range_=22, ew="W", section=35)
for r in rows:
    print(r["kid"], r["well"], r["api"], r["las_url"])

path = kgs.fetch_las(1046139243, dest_dir="data/cache")   # cached; returns the local path
```

`search_wells` wraps the KGS index search (township 1–35 S, range 1–43 W or
1–25 E, section 1–36, plus lease / operator / county / API filters) and
returns one dict per LAS file with its download URL. `fetch_las` resolves a
KID to its URL (pass `url=` from a search row to skip that; otherwise it reads
the header page for the log year and probes the per-year folders KGS files
under) and streams it to disk, refusing anything that does not start like a
LAS file. There is no documented API — this scrapes two KGS pages, so expect
it to need a fix when KGS changes their HTML. KGS also publishes a full index
of every well with an LAS file as
[`ks_las_files.zip`](https://www.kgs.ku.edu/PRS/Ora_Archive/ks_las_files.zip)
if you want to query offline.

## Data

Both files come from the [KGS LAS File Database](https://www.kgs.ku.edu/Magellan/Logs/)
and were filed with KGS by the Kansas Corporation Commission ("LAS File,
courtesy KCC"). Operator: Downing-Nelson Oil Co. Inc.; logging contractor:
Casedhole Solutions. They are redistributed here as-is; the code license
(MIT, see `LICENSE`) does not cover them — refer to the
[KGS data resources page](https://www.kgs.ku.edu/General/dataLib.html).

| File | KGS LAS KID | Well | API | County | Logged | Notes |
|---|---|---|---|---|---|---|
| `data/1046139243.las` | 1046139243 | Pearson Family #1-35, T13S R22W Sec 35 | 15-195-23011 | Trego | 2016-11-21 | LAS 2.0, 0.25 ft, 0–4360.75 ft, 20 curves |
| `data/1045399712.csv` | 1045399712 | PBW #1-32, T18S R17W Sec 32 | 15-165-22116 | Rush | 2015-09-25 | Converted from a whitespace-delimited, CR-terminated export; 0.5 ft, 0–3853 ft, 21 columns (adds `SSD`); `-999.25` nulls kept verbatim |

Curve mnemonics (both wells): `GR`, `SP`, DIL resistivity `RILD` / `RILM` /
`RLL3`, `RxoRt`, density `RHOB` / `RHOC` / `DPOR`, neutron `CNPOR`, sonic `DT` /
`ITT` / `SPOR`, microlog `MEL15` / `MEL20`, calipers `DCAL` / `MELCAL`, borehole
volumes `ABHV` / `TBHV`.

## History

This started in 2016–2017 as a single notebook plotting the Pearson log. The
2026 rework turned it into a package with tests and CI; see the closed issues
for what was wrong with the original (undeclared sentinels, a matplotlib
style name that no longer exists, 15,000 × 15,000-pixel figures committed
into the notebook, and petrophysics constants that were defined but never
used).
