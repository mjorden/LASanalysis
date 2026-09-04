# Methods

What the package computes, in the order the pipeline runs it, with the
choices that were made and why. Function names refer to `lasanalysis.petro`
unless stated.

## Cleaning

lasio masks the declared `NULL` (`-999.25` in the KGS files). It cannot know
about the sentinels the service company never declared, and the Kansas files
carry several. `read_las` applies these rules and records what it touched in
`las.mask_report`:

| Rule | Trigger | Why it exists |
|---|---|---|
| `extra_null` | value in `EXTRA_NULLS` (`-999.0`) | RHOB uses `-999.000`, not the declared `-999.25` |
| `off_scale` | resistivity curve (OHM unit) `>= 1e5` | RILD / RILM carry `100000.0` where the tool went off scale — 1,130 samples in the Pearson log; the old notebook hid them with `xlim` |
| `porosity` | PU / % curve outside `[-50, 100]` | DPOR reaches 58,000 pu on the bad-RHOB rows |
| `density` | **bulk-density** curve outside `[1.0, 3.5]` g/cc (`[1000, 3500]` for kg/m³) | RHOB has a few exact `0.0` |

`rule_for(mnemonic, unit)` in `lasanalysis.load` is the single decision
point. The density rule is gated on the mnemonic as well as the unit: a
density *correction* curve such as `RHOC` / `DRHO` is in g/cc too and
legitimately sits near zero. The first version of this rule did not check
the name and silently wiped all 3,767 RHOC samples in the Pearson log
(issue #17) — the reason the mask report exists is so that this kind of
thing is visible. Units are whitespace-normalised. When a curve has no unit
(CSV exports) the mnemonic alone decides.

Pearson after cleaning:
`{'DPOR': {'porosity': 31}, 'RHOB': {'extra_null': 28, 'density': 3}, 'RILD': {'off_scale': 268}, 'RILM': {'off_scale': 862}}`.

## Standard names

`standardize(df)` renames curves via `ALIASES` (first match wins; unmatched
columns are kept):

| Standard | Aliases |
|---|---|
| `GR` | GR, GRC, SGR, GAM |
| `RT` / `RM` / `RXO` | RILD, ILD, RT, LLD, RDEP, AT90 / RILM, ILM, RMED, AT30 / RLL3, LL3, SFL, RXO, MSFL |
| `RHOB` | RHOB, DEN, ZDEN, RHOZ |
| `NPHI` / `DPHI` / `SPHI` | CNPOR, NPHI, NPOR, TNPH, CNC / DPOR, DPHI, PHID / SPOR, SPHI, PHIS |
| `DT`, `SP`, `CALI` | DT, DTC, AC, DTCO / SP / DCAL, CALI, CAL, HCAL |

## Shale volume

`vshale_linear` is the gamma-ray index `IGR = (GR − GR_clean) / (GR_dirty − GR_clean)`
clipped to [0, 1]. `vshale_larionov(older=True)` — `0.33·(2^(2·IGR) − 1)` —
is what the batch driver and viewer use, since the Kansas section is
Paleozoic; the Tertiary form (`0.083·(2^(3.7·IGR) − 1)`) is available. The
picks used throughout are `GR_clean = 20`, `GR_dirty = 110` GAPI; the
original 2017 notebook defined the same constants and never used them.

## Porosity

`density_porosity(rhob, matrix, rho_fluid)` = `(ρma − RHOB) / (ρma − ρf)`,
as a fraction, not clipped. `matrix` is a `MATRIX_DENSITY` key (sandstone
2.65, limestone 2.71, dolomite 2.87, salt 2.03, anhydrite 2.98) or a number.
Limestone is the working assumption for these wells.

**Neutron scale.** A compensated-neutron curve is reported on the lithology
the service company selected; `CNPOR` is limestone-scaled. Averaging it with
a density porosity computed on another matrix mixes two scales, so
`neutron_lithology_correction` shifts it first using first-order chart
offsets relative to limestone (sandstone +4 pu, dolomite −6 pu; good to
about ±2 pu). `neutron_matrix` is a parameter (default limestone). For
matrices with no neutron scale (salt, anhydrite, a bare density) the curve is
left as-is with a warning, and `summary.csv` records `phin_corrected = False`.
`PHIND` is the arithmetic mean of the two and is the porosity used for Sw.

## Rw and m

No water analysis exists for either well, so Rw and the cementation exponent
come from the log, two independent ways.

**Pickett envelope (`fit_water_line`).** Water-bearing rock plots on a line
of slope −m in log Rt – log φ space with intercept a·Rw at φ = 1. The
function bins clean points by log φ, takes the 5th percentile of Rt per bin
as the envelope, and fits a line through those points. Two choices matter:

| Well | Porosity | φ ≥ | m | a·Rw | n |
|---|---|---|---|---|---|
| Pearson | PHID | 0.03 | 1.30 | 0.081 | 810 |
| Pearson | PHID | 0.06 | 1.60 | 0.043 | 480 |
| Pearson | **PHIND** | **0.06** | **1.96** | **0.031** | 689 |
| Pearson | PHIND | 0.08 | 2.35 | 0.015 | 536 |
| PBW | PHIND | 0.06 (Vsh < 0.15) | 1.21 | 0.333 | 469 |
| PBW | PHIND | 0.06 (Vsh < 0.10) | 1.99 | 0.063 | 332 |
| PBW | PHIND | 0.08 (Vsh < 0.10) | 2.84 | 0.011 | 243 |

Below ~6 % porosity shale conductivity and matrix error flatten the
envelope and drag m toward 1.3, so the pick uses neutron-density porosity
with a 6 % cutoff. On Pearson that is stable across Vsh cutoffs (1.93–1.97);
on PBW it is not — the 3000–3850 ft section mixes lithologies (N-D porosity
runs to 39 %) and the fit swings from 0.4 to 2.8 with the cut.

**Rwa pick (`pick_rw_from_rwa`).** Apparent water resistivity
`Rwa = Rt·φ^m / a` equals Rw where Sw = 1 and exceeds it elsewhere, so a low
percentile (5th) of Rwa over clean, porous samples is an Rw estimate that
needs no straight line — at the cost of assuming m. It reports the depth
interval the supporting samples came from.

| Well | Envelope | Rwa pick (m = 2) | Clean wet zones | **Pick** |
|---|---|---|---|---|
| Pearson #1-35 | m 1.96, a·Rw 0.031 | Rw 0.027 (689 samples) | 3580–3650 ft: φ ≈ 16 %, Rt 0.7–1.5 → Rw 0.02–0.04 | **Rw 0.03, m 2.0** |
| PBW #1-32 | unstable | Rw ≈ 0.06 | 3250–3340 ft: φ 15–17 %, Rt 1.7–2.8 → Rw 0.05–0.066 | **Rw 0.06, m 2.0** |

`a = 1` and `n = 2` are assumed, not fitted. 0.03 Ω·m is what a ~100,000 ppm
NaCl brine looks like at ~110 °F (roughly 4,000 ft in central Kansas); PBW's
0.06 is consistent with being ~400 ft shallower. Neither is a substitute for
a produced-water analysis. `multiwell` writes `rw_envelope`, `m_envelope`,
`rw_rwa` and `rwa_interval` for every well so the two can disagree visibly.

## Water saturation

`water_saturation(model, rt, phi, vsh, rw, rsh, a, m, n)` dispatches on:

- **archie** — `Sw = ((a·Rw) / (φ^m·Rt))^(1/n)`. Default. Ignores shale
  conductivity, so it overstates Sw in shaly intervals.
- **simandoux** (modified) — solves `1/Rt = φ^m·Sw^n / (a·Rw) + Vsh·Sw / Rsh`;
  closed form at n = 2, Newton iterations from the Archie estimate otherwise.
- **indonesia** (Poupon–Leveaux) — `1/√Rt = [Vsh^(1−Vsh/2)/√Rsh + φ^(m/2)/√(a·Rw)]·Sw^(n/2)`.

Both shaly models reduce exactly to Archie at Vsh = 0 (tested to 1e-9) and
need a shale resistivity; `pick_rsh` takes the median Rt where Vsh ≥ 0.8
(Pearson: 4.24 Ω·m) when `rsh` is not given. Results are clipped to 1 by
default. When a shaly model is active the frame also carries `SW_ARCHIE` for
comparison. On Pearson, 3400–4200 ft, mean Sw is 0.75 (Archie), 0.64
(Simandoux) and 0.63 (Indonesia) with the auto-picked Rsh.

## Pay flag

`PAY = φ > phi_cut and Vsh < vsh_cut and Sw < sw_cut` with defaults 0.08 /
0.30 / 0.50, using PHIND (or PHID when there is no neutron). `pay_ft` is the
flagged sample count times the sample step. It is a screening flag over a
mixed section, not a net-pay count; the Vsh cutoff is also what keeps
Archie's shale bias out of it.

## `summary.csv`

One row per well from `multiwell`: identity (`kid`, `well`, `api`,
`well_kid`, `lat`, `lon`, `elevation`, `total_depth`, `status`,
`producing_formation`, dates); interval (`depth_top`, `depth_base`, `step`,
`n_samples`, `curves`); means of Vsh, PHIND, PHID, Sw; pay (`pay_ft`,
`pay_top`, `pay_base`, `pay_sw_mean`); the two Rw picks; `sw_model`, `rsh`,
`phin_corrected`; `params` (JSON of what was used); `mask_report`; the file
names of the PNG / HTML written; `error` when a well failed.

## Caveats

- Everything is single-well quick-look: no temperature correction of Rw, no
  invasion correction, no tool-specific environmental corrections.
- The neutron correction is a linear chart approximation.
- Larionov and the GR picks are held constant over the whole interval.
- The pay flag depends entirely on the three cutoffs and the Rw pick.
