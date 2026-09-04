# The KGS data sources

Everything `lasanalysis.kgs` knows about the Kansas Geological Survey's LAS
holdings, learned by probing in September 2026. There is no documented API
behind any of this; the client scrapes three pages and reads one zip file,
and will need a fix whenever KGS changes their HTML.

## Two kinds of KID

KGS assigns a KID to a *well* and a separate KID to each *LAS file* of that
well. The files in `data/` are named after LAS KIDs. The pages take
different ones:

| Page | Takes | Example |
|---|---|---|
| `qualified.well_page.DisplayWell?f_kid=` | well KID | 1046105344 (Pearson) |
| `las.lasd5.ViewLasHeader?f_kid=` | LAS KID | 1046139243 (Pearson's log) |

Passing a LAS KID to the well page returns an empty template, which
`well_info` recognises and refuses.

## Search

`https://chasm.kgs.ku.edu/ords/las.lasd5.SelectWells` accepts the form
fields `f_t` (township, 1–35, all south), `f_r` (range), `ew` (`E` ≤ 25 /
`W` ≤ 43), `f_s` (section 1–36), `f_l` (lease substring), `f_op` (operator
substring), `f_c` (county code), `f_api`, `f_st=15`. GET works. The result
is an HTML table; each LAS file row carries the well page link, the header
page link and a direct download link. `parse_search_html` reads it; the two
markers it checks for (`Select location of well`, `No wells found`) are the
guard against parsing something else.

## Download

LAS files live on Azure blob storage under a per-year folder:

```
https://kgsimages.blob.core.windows.net/web/web_1/WebDocs/WellLogs/kcc_logs_2016/1046139243.las
```

The year is the year KGS received the log (usually, not always, the log
year). Older logs are filed under township folders instead
(`…/WellLogs/01S02E/1020069094.las`). `resolve_las_url` reads the header
page's `DATE.` for a first guess and probes the per-year folders with HEAD
requests; it cannot find township-folder logs — the offline index can.

`fetch_las` streams to a temp file and promotes it only if the content starts
like a LAS file (`~V`). Any caller-supplied URL must be on the blob host with
the blob scheme (`check_las_url`); `las_url` values come out of KGS's HTML by
regex and must never become an arbitrary request target (issue #18).

## Well page

`DisplayWell` for a well KID gives: API, lease, well number, original and
current operator, field, location (T-R-S plus footages), **NAD27 and NAD83
latitude / longitude** with KGS's own statement of their source ("from GPS"
or "calculated from footages"), county, permit / spud / completion / plugging
dates, well type, status, total depth, elevation with datum (`2395 KB`),
producing formation. `parse_well_page` returns them; `well_info` caches the
result as JSON per well KID under `data/cache/wells/`. `multiwell` fetches
this for every well in a batch by default and draws `wells.png`.

## Offline index

`https://www.kgs.ku.edu/PRS/Ora_Archive/ks_las_files.zip` (~1.4 MB) holds
one row per LAS file KGS has — 29,089 on 2026-09-03 — with the columns
`KGS_ID` (well KID), `Latitude`, `Longitude`, `Location`, `Operator`,
`Lease`, `API`, `API_NUM_NODASH`, `Elevation`, `Elev_Ref`, `Depth_start`,
`Depth_stop`, `URL`. Three things to know:

1. The coordinates are **NAD27** (they match the well page's NAD27 values).
   `load_index` names them `lat_nad27` / `lon_nad27`.
2. The `URL` host, `www.kgs.ku.edu/b_1/WebDocs/WellLogs/…`, returned 404 for
   every row tried; the path tail (`<folder>/<las_kid>.las`) is exactly the
   blob store's. `load_index` keeps the tail as `las_path` and rebuilds
   `las_url` on the blob host, verified live on both year and township
   folders. This is how the index reaches logs the year probe cannot.
3. About 1,700 rows carry a folder-less `https://www.kgs.ku.edu//<kid>.las`.
   Those get `las_url = None` and fall back to `resolve_las_url`.

`search_index` filters the frame with the same arguments as `search_wells`
plus `within=(lat, lon, km)` (great-circle, on the NAD27 coordinates);
`fetch_index` caches the zip for 30 days. From the CLI: `multiwell --index
[ZIP] --within LAT LON KM`.

## Other KGS holdings relevant to the RFC

- Core library index: `https://www.kgs.ku.edu/PRS/Ora_Archive/coreLib.zip`
  (which wells have core).
- Open-file Report 2000-64 (Newell & Hatch): 493 Rock-Eval / TOC samples and
  63 vitrinite-reflectance values from 83 wells across the Salina, Forest
  City and Sedgwick basins and the CKU / Nemaha uplifts — PDF tables.
- Formation tops per well on the well page (not yet read by the client).

## Politeness

One request per page, a descriptive `User-Agent`, results cached (index for
30 days, well pages and LAS files indefinitely), and nothing downloaded that
the caller did not ask for. The test suite never touches KGS: it runs
against captured fixtures on a local stub server, with one HEAD-only live
check behind `LASANALYSIS_NETWORK=1`.
