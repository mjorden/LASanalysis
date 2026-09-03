"""Build the GitHub Pages site: one interactive viewer per well plus an index.

Run:  python scripts/build_site.py [site_dir]      (default: site/)
Deployed by .github/workflows/pages.yml on every push to master.
"""

from __future__ import annotations

import html
import sys
from pathlib import Path

from lasanalysis import read_log_csv, standardize
from lasanalysis.viewer import write_viewer

ROOT = Path(__file__).resolve().parents[1]
REPO_URL = "https://github.com/mjorden/LASanalysis"

WELLS = [
    {
        "page": "pearson.html",
        "source": ROOT / "data" / "1046139243.las",
        "kid": 1046139243,
        "well": "Pearson Family #1-35",
        "api": "15-195-23011",
        "location": "T13S R22W Sec 35, Trego County, KS",
        "logged": "2016-11-21",
        "depth_range": (3400, 4200),
        "meta": None,  # taken from the LAS header
        # Pickett envelope, Vsh < 0.15, N-D porosity >= 6 %: m = 1.96, a*Rw = 0.031;
        # clean wet zone 3580-3650 ft gives Rw 0.02-0.04 at m = 2.
        "params": {"rw": 0.03, "m": 2.0},
        "pick": "Rw 0.03, m 2.0 from the Pickett envelope; wet zone 3580–3650 ft agrees",
    },
    {
        "page": "pbw.html",
        "source": ROOT / "data" / "1045399712.csv",
        "kid": 1045399712,
        "well": "PBW #1-32",
        "api": "15-165-22116",
        "location": "T18S R17W Sec 32, Rush County, KS",
        "logged": "2015-09-25",
        "depth_range": (3000, 3850),
        "meta": {"well": "PBW #1-32", "uwi": "15-165-22116-00-00", "county": "RUSH", "state": "KANSAS",
                 "operator": "DOWNING NELSON OIL CO. INC.", "service": "CASEDHOLE SOLUTIONS", "date": "2015-09-25",
                 "file": "1045399712.csv"},
        # The envelope fit is unstable here (mixed lithology, N-D porosity to 39 %):
        # m ranges 0.4-2.8 with the cut. The cleanest cut (Vsh < 0.10, N-D porosity
        # >= 6 %) gives m = 1.99, a*Rw = 0.063, and the clean wet zones at
        # 3250-3340 ft (phi 15-17 %, Rt 1.7-2.8) give Rw 0.05-0.066 at m = 2.
        "params": {"rw": 0.06, "m": 2.0},
        "pick": "Rw 0.06, m 2.0 from wet zones 3250–3340 ft; the Pickett envelope is unstable in this mixed section",
    },
]

INDEX = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>LASanalysis — interactive well logs</title>
<style>
  body {{ margin:0; padding:32px 16px; background:#fafaf8; color:#222; font:15px/1.5 Georgia, "Times New Roman", serif; }}
  main {{ max-width:720px; margin:0 auto; }}
  h1 {{ font-weight:normal; font-size:26px; margin:0 0 4px; }}
  .sub {{ color:#666; margin-bottom:28px; }}
  a {{ color:#2a6f97; }}
  .well {{ display:block; padding:16px 18px; margin:12px 0; border:1px solid #ddd; background:#fff; text-decoration:none; color:inherit; }}
  .well:hover {{ border-color:#2a6f97; }}
  .well b {{ font-weight:normal; font-size:18px; color:#2a6f97; }}
  .well small {{ color:#666; display:block; margin-top:2px; }}
  .fine {{ color:#666; font-size:13px; margin-top:28px; }}
</style></head><body><main>
<h1>LASanalysis</h1>
<div class="sub">Quick-look petrophysics on Kansas Geological Survey well logs. Each page is a self-contained viewer:
shared-depth tracks, live Rw / m / n / matrix sliders, pay shading, Pickett plot.</div>
{cards}
<p class="fine">Source and method: <a href="{repo}">{repo}</a>. Data from the
<a href="https://www.kgs.ku.edu/Magellan/Logs/">KGS LAS File Database</a> (courtesy KCC), operator Downing-Nelson Oil Co. Inc.
Rw and m are picked per well from its own Pickett plot and wet zones (shown on each card); a = 1, n = 2 assumed. Not a reserves estimate.</p>
</main></body></html>
"""

CARD = """<a class="well" href="{page}"><b>{well}</b><small>API {api} · {location} · logged {logged} · KGS KID {kid} · viewer opens at {top}–{base} ft</small><small>{pick}</small></a>"""


def build(site_dir: Path) -> list:
    site_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for w in WELLS:
        src = w["source"]
        if src.suffix.lower() == ".csv":
            df = standardize(read_log_csv(src))
            out = write_viewer(df, site_dir / w["page"], params=w["params"], depth_range=w["depth_range"], meta=w["meta"], title=w["well"])
        else:
            out = write_viewer(src, site_dir / w["page"], params=w["params"], depth_range=w["depth_range"], title=w["well"])
        written.append(out)
    cards = "\n".join(
        CARD.format(page=w["page"], well=html.escape(w["well"]), api=w["api"], location=html.escape(w["location"]),
                    logged=w["logged"], kid=w["kid"], top=w["depth_range"][0], base=w["depth_range"][1], pick=html.escape(w["pick"]))
        for w in WELLS
    )
    index = site_dir / "index.html"
    index.write_text(INDEX.format(cards=cards, repo=REPO_URL), encoding="utf-8")
    written.append(index)
    (site_dir / ".nojekyll").write_text("", encoding="utf-8")
    return written


if __name__ == "__main__":
    target = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "site"
    for p in build(target):
        print(f"wrote {p} ({p.stat().st_size // 1024} KB)")
