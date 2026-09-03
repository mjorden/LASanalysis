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
) -> List[dict]:
    """Search the KGS LAS index. Returns :func:`parse_search_html` rows.

    Township is 1-35 (all Kansas townships are south), range 1-43 W or 1-25 E,
    section 1-36. Text filters match a substring, case-insensitively, the way
    the web form does. All filters are optional but KGS refuses an empty query.
    """
    params = search_params(township, range_, ew, section, lease, operator, county, api)
    if not any(v for k, v in params.items() if k != "f_st"):
        raise ValueError("give at least one search filter")
    r = _session(session).get(f"{base_url}/las.lasd5.SelectWells", params=params, timeout=timeout)
    r.raise_for_status()
    if "Select location of well" not in r.text and "No wells found" not in r.text:
        raise KGSError("unexpected SelectWells response (KGS page layout may have changed)")
    return parse_search_html(r.text)


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
        r = s.head(url, timeout=timeout, allow_redirects=True)
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
