"""Kansas Geological Survey LAS-file index client.

KGS indexes every LAS file it holds at https://www.kgs.ku.edu/Magellan/Logs/.
The search form posts to an Oracle ORDS endpoint (``las.lasd5.SelectWells``)
that returns an HTML table; each LAS file in it has a KGS *LAS KID* (the
number this repo's data files are named after), a header-preview page
(``las.lasd5.ViewLasHeader?f_kid=``) and a direct download link on Azure blob
storage filed under the year the log was received::

    https://kgsimages.blob.core.windows.net/web/web_1/WebDocs/WellLogs/kcc_logs_2016/1046139243.las

There is no documented API; this module scrapes the two pages above and is
therefore only as stable as KGS's HTML. Every function takes ``base_url`` /
``blob_url`` so the tests can point them at a local stub server.
"""

from __future__ import annotations

import datetime as _dt
import os
import re
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional
from urllib.parse import urlparse

import requests

CHASM_URL = "https://chasm.kgs.ku.edu/ords"
BLOB_URL = "https://kgsimages.blob.core.windows.net/web/web_1/WebDocs/WellLogs"
USER_AGENT = "lasanalysis (https://github.com/mjorden/LASanalysis)"

#: Kansas state FIPS code, the only value the KGS form accepts.
KANSAS_FIPS = 15

_TAG = re.compile(r"<[^>]+>")
_WS = re.compile(r"\s+")


class KGSError(RuntimeError):
    """The KGS service answered, but not with what we asked for."""


def _session(session: Optional[requests.Session]) -> requests.Session:
    if session is not None:
        return session
    s = requests.Session()
    s.headers["User-Agent"] = USER_AGENT
    return s


def _text(html: str) -> str:
    return _WS.sub(" ", _TAG.sub(" ", html)).strip()


# --------------------------------------------------------------------------- search


def _validate_location(township, range_, ew, section) -> None:
    if township is not None and not (1 <= int(township) <= 35):
        raise ValueError("township must be 1-35 (south)")
    if section is not None and not (1 <= int(section) <= 36):
        raise ValueError("section must be 1-36")
    if range_ is not None:
        if ew not in ("E", "W"):
            raise ValueError("ew must be 'E' or 'W' when range_ is given")
        hi = 25 if ew == "E" else 43
        if not (1 <= int(range_) <= hi):
            raise ValueError(f"range must be 1-{hi} for {ew}")


def search_params(
    township=None,
    range_=None,
    ew=None,
    section=None,
    lease: str = "",
    operator: str = "",
    county: str = "",
    api: str = "",
) -> Dict[str, str]:
    """Form fields for ``las.lasd5.SelectWells``. Validates the location fields."""
    _validate_location(township, range_, ew, section)
    return {
        "f_t": "" if township is None else str(int(township)),
        "f_r": "" if range_ is None else str(int(range_)),
        "ew": ew or "",
        "f_s": "" if section is None else str(int(section)),
        "f_l": lease,
        "f_op": operator,
        "f_st": str(KANSAS_FIPS),
        "f_c": county,
        "f_api": api,
    }


def parse_search_html(html: str) -> List[dict]:
    """Rows of a ``SelectWells`` result page, one dict per LAS file.

    Keys: ``kid`` (LAS KID, int), ``las_url``, ``header_url``, ``well_kid``,
    ``location``, ``operator``, ``well``, ``api``, ``spud``, ``plug``,
    ``depth_start``, ``depth_stop``, ``las_label``.
    A well with several LAS files yields several dicts sharing the well fields.
    """
    if "No wells found" in html:
        return []
    rows = re.findall(r"<tr[^>]*>(.*?)</tr>", html, flags=re.S | re.I)
    out: List[dict] = []
    well: Dict[str, object] = {}
    for row in rows:
        cells = re.findall(r"<td[^>]*>(.*?)</td>", row, flags=re.S | re.I)
        m = re.search(r'DisplayWell\?f_kid=(\d+)"[^>]*>(.*?)</a>', row, flags=re.S | re.I)
        if m and len(cells) >= 7:
            dates = re.split(r"<br\s*/?>", cells[4], flags=re.I)
            well = {
                "well_kid": int(m.group(1)),
                "location": _text(m.group(2)),
                "operator": _text(cells[1]),
                "well": _text(cells[2]),
                "api": _text(cells[3]),
                "spud": _text(dates[0]) if dates else "",
                "plug": _text(dates[1]) if len(dates) > 1 else "",
                "depth_start": _to_float(_text(cells[5])),
                "depth_stop": _to_float(_text(cells[6])),
            }
        for cell in cells:
            hm = re.search(r'ViewLasHeader\?f_kid=(\d+)', cell)
            if not hm:
                continue
            um = re.search(r'href="([^"]+\.las)"', cell, flags=re.I)
            label = _text(re.split(r"<br\s*/?>", cell, flags=re.I)[0])
            out.append(
                {
                    **well,
                    "kid": int(hm.group(1)),
                    "las_url": um.group(1) if um else None,
                    "header_url": f"{CHASM_URL}/las.lasd5.ViewLasHeader?f_kid={hm.group(1)}",
                    "las_label": label,
                }
            )
    return out


def _to_float(s: str) -> Optional[float]:
    try:
        return float(s)
    except ValueError:
        return None


def search_wells(
    township=None,
    range_=None,
    ew=None,
    section=None,
    lease: str = "",
    operator: str = "",
    county: str = "",
    api: str = "",
    *,
    base_url: str = CHASM_URL,
    session: Optional[requests.Session] = None,
    timeout: float = 60,
    with_coords: bool = False,
    cache_dir=None,
) -> List[dict]:
    """Search the KGS LAS index. Returns :func:`parse_search_html` rows.

    Township is 1-35 (all Kansas townships are south), range 1-43 W or 1-25 E,
    section 1-36. Text filters match a substring, case-insensitively, the way
    the web form does. All filters are optional but KGS refuses an empty query.
    ``with_coords=True`` also fetches each well's page (#31) and merges
    NAD83 ``lat`` / ``lon``, elevation, TD, status and dates into the rows.
    """
    params = search_params(township, range_, ew, section, lease, operator, county, api)
    if not any(v for k, v in params.items() if k != "f_st"):
        raise ValueError("give at least one search filter")
    s = _session(session)
    r = s.get(f"{base_url}/las.lasd5.SelectWells", params=params, timeout=timeout)
    r.raise_for_status()
    if "Select location of well" not in r.text and "No wells found" not in r.text:
        raise KGSError("unexpected SelectWells response (KGS page layout may have changed)")
    rows = parse_search_html(r.text)
    if with_coords:
        add_well_info(rows, base_url=base_url, session=s, timeout=timeout, cache_dir=cache_dir)
    return rows


# --------------------------------------------------------------------------- well page (#31)

_WELL_FIELDS = {
    # label on the page -> key in the returned dict
    "API": "api", "KID": "well_kid", "Lease": "lease", "Well": "well_no", "Original operator": "operator",
    "Current operator": "current_operator", "Field": "field", "Location": "location",
    "NAD27 Longitude": "lon_nad27", "NAD27 Latitude": "lat_nad27", "NAD83 Longitude": "lon", "NAD83 Latitude": "lat",
    "County": "county", "Permit Date": "permit_date", "Spud Date": "spud_date", "Completion Date": "completion_date",
    "Plugging Date": "plug_date", "Well Type": "well_type", "Status": "status", "Total Depth": "total_depth",
    "Elevation": "elevation", "Producing Formation": "producing_formation",
}
_NUMERIC_WELL_FIELDS = ("lon", "lat", "lon_nad27", "lat_nad27", "total_depth")


def parse_well_page(html: str) -> Dict[str, object]:
    """Fields of a ``qualified.well_page.DisplayWell`` page.

    Returns a dict with the keys in ``_WELL_FIELDS`` (missing ones absent),
    numbers parsed for the coordinate and depth fields, ``elevation`` as a
    number with ``elevation_datum`` (e.g. ``"KB"``) split off, and
    ``latlon_source`` (KGS states whether it came from GPS or footages).
    An empty template (unknown KID) yields an empty dict.
    """
    tokens = [t.strip() for t in re.split(r"<[^>]+>", html)]
    tokens = [t for t in tokens if t]
    out: Dict[str, object] = {}
    for i, tok in enumerate(tokens):
        label = tok.rstrip(":").strip()
        if label in _WELL_FIELDS and tok.endswith(":") and i + 1 < len(tokens):
            val = tokens[i + 1]
            if val.endswith(":"):
                continue  # empty field: the next token is another label ("IP Oil (bbl):", ...)
            out[_WELL_FIELDS[label]] = val
    m = re.search(r"Lat-long\s+([^<|]+?)\s*(?:<|$)", html)
    if m:
        out["latlon_source"] = m.group(1).strip()
    for k in _NUMERIC_WELL_FIELDS:
        if k in out:
            try:
                out[k] = float(str(out[k]).replace(",", ""))
            except ValueError:
                del out[k]
    if "elevation" in out:
        em = re.match(r"\s*([-\d.,]+)\s*([A-Za-z]*)", str(out["elevation"]))
        if em:
            out["elevation"] = float(em.group(1).replace(",", ""))
            out["elevation_datum"] = em.group(2).upper() or None
        else:
            del out["elevation"]
    if "well_kid" in out:
        try:
            out["well_kid"] = int(str(out["well_kid"]))
        except ValueError:
            pass
    if "location" in out and "lat" not in out and "lon" not in out and len(out) <= 1:
        return {}
    return out


def well_info(
    well_kid: int,
    *,
    base_url: str = CHASM_URL,
    session: Optional[requests.Session] = None,
    timeout: float = 60,
    cache_dir: "os.PathLike | str | None" = None,
) -> Dict[str, object]:
    """Header data for a *well* KID (not a LAS KID) from the KGS well page, cached as JSON.

    Raises ``KGSError`` when KGS returns its empty template for the KID.
    """
    import json

    well_kid = int(well_kid)
    cache = None if cache_dir is None else Path(cache_dir) / "wells" / f"{well_kid}.json"
    if cache is not None and cache.exists():
        return json.loads(cache.read_text(encoding="utf-8"))
    r = _session(session).get(f"{base_url}/qualified.well_page.DisplayWell", params={"f_kid": well_kid}, timeout=timeout)
    r.raise_for_status()
    info = parse_well_page(r.text)
    if not info.get("lat") and not info.get("api"):
        raise KGSError(f"KGS well page for KID {well_kid} is empty (is this a LAS KID rather than a well KID?)")
    info.setdefault("well_kid", well_kid)
    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(info, indent=1), encoding="utf-8")
    return info


#: Fields copied from :func:`well_info` onto a search row by :func:`add_well_info`.
WELL_ROW_FIELDS = ("lat", "lon", "latlon_source", "elevation", "elevation_datum", "total_depth", "well_type", "status",
                   "producing_formation", "spud_date", "completion_date", "county")


def add_well_info(rows: List[dict], *, base_url: str = CHASM_URL, session=None, timeout: float = 60,
                  cache_dir=None, log=None) -> List[dict]:
    """Merge coordinates and header fields into search rows in place (one request per distinct well_kid).

    A well page that fails leaves the row without those keys and, if ``log``
    is given, reports why — the batch never aborts on a missing page.
    """
    s = _session(session)
    seen: Dict[int, Optional[dict]] = {}
    for row in rows:
        wk = row.get("well_kid")
        if wk is None:
            continue
        if wk not in seen:
            try:
                seen[wk] = well_info(wk, base_url=base_url, session=s, timeout=timeout, cache_dir=cache_dir)
            except (requests.RequestException, KGSError) as e:
                seen[wk] = None
                if log:
                    log(f"  well {wk}: no header info ({type(e).__name__}: {e})")
        info = seen[wk]
        if info:
            for k in WELL_ROW_FIELDS:
                if k in info:
                    row[k] = info[k]
    return rows


# --------------------------------------------------------------------------- fetch


def las_header(kid: int, *, base_url: str = CHASM_URL, session=None, timeout: float = 60) -> str:
    """Plain text of the LAS header preview page for ``kid`` (through ``~A``)."""
    r = _session(session).get(f"{base_url}/las.lasd5.ViewLasHeader", params={"f_kid": int(kid)}, timeout=timeout)
    r.raise_for_status()
    text = _text(r.text)
    if "~Version" not in text:
        raise KGSError(f"KGS has no LAS header for KID {kid}")
    return text


def log_year_from_header(header_text: str) -> Optional[int]:
    """Four-digit year from the ``DATE.`` line of a LAS header, or None."""
    m = re.search(r"DATE\.\s*(.*?):\s*LOG DATE", header_text)
    if not m:
        m = re.search(r"DATE\.\s*([^:]*)", header_text)
    if not m:
        return None
    years = re.findall(r"(?<!\d)((?:19|20)\d{2})(?!\d)", m.group(1))
    if years:
        return int(years[-1])
    two = re.findall(r"\d{1,2}/\d{1,2}/(\d{2})(?!\d)", m.group(1))  # 11/21/16
    return 2000 + int(two[-1]) if two else None


def las_url(kid: int, year: int, *, blob_url: str = BLOB_URL) -> str:
    return f"{blob_url}/kcc_logs_{int(year)}/{int(kid)}.las"


def check_las_url(url: str, blob_url: str = BLOB_URL) -> str:
    """Refuse a download URL that does not live on the expected blob host.

    ``las_url`` values come out of KGS's HTML by regex; they must never become
    an arbitrary request target. Scheme and host must match ``blob_url``.
    """
    want, got = urlparse(blob_url), urlparse(str(url))
    if got.scheme != want.scheme or got.netloc != want.netloc:
        raise KGSError(f"refusing LAS URL {url!r}: expected {want.scheme}://{want.netloc}/...")
    return url


def resolve_las_url(
    kid: int,
    *,
    year: Optional[int] = None,
    years: Optional[Iterable[int]] = None,
    base_url: str = CHASM_URL,
    blob_url: str = BLOB_URL,
    session=None,
    timeout: float = 60,
) -> str:
    """Find the download URL for a LAS KID.

    KGS files each LAS under the year it was received, which is usually — not
    always — the log year. We read the header's ``DATE.`` to get a first guess
    and then ``HEAD`` the candidate URLs in order: guess, guess+1, guess-1,
    then every year from 2005 to next year.
    """
    s = _session(session)
    if year is None and years is None:
        try:
            year = log_year_from_header(las_header(kid, base_url=base_url, session=s, timeout=timeout))
        except (requests.RequestException, KGSError):
            year = None
    this_year = _dt.date.today().year
    if years is None:
        ordered: List[int] = []
        if year is not None:
            ordered += [year, year + 1, year - 1]
        ordered += range(this_year + 1, 2004, -1)
        years = list(dict.fromkeys(y for y in ordered if 1990 <= y <= this_year + 1))
    tried = []
    for y in years:
        url = las_url(kid, y, blob_url=blob_url)
        tried.append(url)
        try:
            r = s.head(url, timeout=timeout, allow_redirects=True)
        except requests.RequestException:
            continue  # transient failure on one candidate year; try the next
        if r.status_code == 200:
            return url
    raise KGSError(f"no LAS file found for KID {kid}; tried {len(tried)} URLs, first {tried[0]}")


def fetch_las(
    kid: int,
    dest_dir: "os.PathLike | str" = "data/cache",
    *,
    url: Optional[str] = None,
    overwrite: bool = False,
    base_url: str = CHASM_URL,
    blob_url: str = BLOB_URL,
    session=None,
    timeout: float = 120,
) -> Path:
    """Download the LAS file for ``kid`` to ``dest_dir/<kid>.las`` and return the path.

    Pass ``url`` (e.g. the ``las_url`` from a :func:`search_wells` row) to skip
    resolution. An existing file is reused unless ``overwrite``. The download
    goes to a temp file and is promoted only if it starts like a LAS file, so a
    200 HTML error page never lands in the cache.
    """
    dest_dir = Path(dest_dir)
    dest = dest_dir / f"{int(kid)}.las"
    if dest.exists() and not overwrite:
        return dest
    s = _session(session)
    if url is None:
        url = resolve_las_url(kid, base_url=base_url, blob_url=blob_url, session=s, timeout=timeout)
    else:
        check_las_url(url, blob_url)
    dest_dir.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=f"{kid}.", suffix=".part", dir=dest_dir)
    try:
        with os.fdopen(fd, "wb") as fh, s.get(url, stream=True, timeout=timeout) as r:
            r.raise_for_status()
            first = b""
            for chunk in r.iter_content(chunk_size=1 << 16):
                if not first:
                    first = chunk
                fh.write(chunk)
        if not _looks_like_las(first):
            raise KGSError(f"{url} did not return a LAS file (starts {first[:40]!r})")
        os.replace(tmp, dest)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    return dest


def _looks_like_las(head: bytes) -> bool:
    return head.lstrip().upper().startswith(b"~V")
