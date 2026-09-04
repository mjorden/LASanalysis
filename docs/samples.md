# Lab samples: SRA / Rock-Eval, XRD and core data

Layer 1 of [RFC #32](https://github.com/mjorden/LASanalysis/issues/32).
`lasanalysis.samples` carries lab results — which arrive as tens of
*samples* per well, at a depth or over an interval — alongside the
continuous curves, without pretending they are curves.

## The table

One long table, one row per (sample, analyte):

| column | meaning |
|---|---|
| `api`, `api10` | API number as given, and its first ten digits — the join key. Lab reports know API numbers, not KGS KIDs |
| `depth_top`, `depth_base` | as reported by the lab; equal for a point sample |
| `sample_type` | `core`, `cuttings`, `swc`, `outcrop`, `unknown` |
| `lab`, `method` | who measured it and how (`Rock-Eval 6`, `SRA`, `LECO`, `XRD Rietveld`, …). Mixed methods must stay distinguishable |
| `analyte`, `value`, `unit` | one of `ANALYTES`: SRA (`TOC`, `S1`, `S2`, `S3`, `Tmax`, `HI`, `OI`, `PI`, `Ro`), XRD minerals (`quartz` … `total_clay`), core (`core_GR`, `core_RHOB`, `core_phi`, `core_perm`) |
| `source` | report id, so every number is traceable |

`read_samples(path_or_df, lab=…, method=…, source=…)` reads a *wide* CSV /
Excel file — one row per sample, one column per analyte — and returns the
long form. Column names are matched case-insensitively with aliases
(`Depth` or `depth_top`; `toc`, `kspar`, `GR_core`, …); unknown columns are
ignored. Templates: [`docs/templates/sra_template.csv`](templates/sra_template.csv),
[`docs/templates/xrd_template.csv`](templates/xrd_template.csv).

`validate_samples(df)` returns a list of problems: rows without a usable API,
missing or inverted depths, values outside each analyte's physical range
(TOC 0–100 wt%, Tmax 300–650 °C, …), XRD totals more than 5 wt% off 100, and
duplicate (well, depth, analyte, lab) rows. `strict=True` raises.

`to_wide(df)` pivots back to one row per sample with a `depth_mid` column.

## Depth alignment

Log depths are measured from the KB; core depths may be driller's depths,
and cuttings lag and cave. When a property measured on the samples is also
logged — core gamma (`core_GR` vs `GR`), core bulk density (`core_RHOB` vs
`RHOB`) — `depth_shift(samples_wide, log, "core_GR", "GR", max_shift=20)`
slides the samples through ±`max_shift` ft and keeps the shift with the best
correlation, reporting `r`, `rmse` and `n`. Positive means the samples move
deeper. `apply_shift` applies it; every join and plot below takes a `shift`
too. There is no shift for cuttings without a measured property — treat their
depths as ±one sample interval.

## Joins

- `join_to_log(samples_wide, log, curves, tolerance=0.5, shift=0)` — log
  values onto samples: nearest log sample for point samples, mean over the
  interval for interval samples. This is what a TOC-from-logs calibration
  regresses against.
- `samples_on_grid(samples_wide, log.index, "TOC")` — one analyte onto the
  log's depth grid (the value over each sample's interval, NaN elsewhere).
- `sample_tracks(samples, ("TOC",), shift=…)` — track specs with `points`
  that `plot_tracks` and the viewer draw as markers:

```python
from lasanalysis import read_las, curves, standardize, read_samples, sample_tracks, plot_tracks
from lasanalysis.viewer import write_viewer

df = standardize(curves(read_las("data/1046139243.las")))
sra = read_samples("my_sra_results.csv", lab="Acme", method="Rock-Eval 6", source="Report 42")
fig = plot_tracks(df, [{"curves": ["GR"], "xlim": (0, 175)}] + sample_tracks(sra, ("TOC",)), depth_range=(3400, 4200))
write_viewer(df, "out.html", samples=sra, sample_analytes=("TOC", "Tmax"))
```

## What is not here yet

The OFR 2000-64 parser (KGS's published Rock-Eval compilation, PDF tables),
the HI / OI / PI / Ro-equivalent derivations and kerogen-type plots
(Layer 2), and log calibration (Layer 3) — see the RFC.
