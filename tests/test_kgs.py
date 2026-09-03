"""kgs client tests against a local stub of the two KGS pages + the blob store.

Fixtures under tests/fixtures/ are real responses captured 2026-09-02.
Set LASANALYSIS_NETWORK=1 to also run the live HEAD check at the bottom.
"""

import os
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pytest
import requests

from lasanalysis import kgs

FIX = Path(__file__).parent / "fixtures"
SEARCH_HTML = (FIX / "selectwells_T13S_R22W_S35.html").read_text(encoding="utf-8")
EMPTY_HTML = (FIX / "selectwells_empty.html").read_text(encoding="utf-8")
HEADER_HTML = (FIX / "viewlasheader_1046139243.html").read_text(encoding="utf-8")
TINY_LAS = b"~Version Information\nVERS. 2.0\n~A DEPT GR\n0.0 10.0\n"
PEARSON_WELL_HTML = (FIX / "displaywell_1046105344.html").read_text(encoding="utf-8")
PBW_WELL_HTML = (FIX / "displaywell_1045079321.html").read_text(encoding="utf-8")
# what KGS returns for a KID it does not know (or a LAS KID passed by mistake): every field blank
EMPTY_WELL_HTML = "<html><body><table><tr><td>API:</td><td>KID:</td><td>Lease:</td><td>Well:</td></tr>" \
                  "<tr><td>NAD83 Longitude:</td><td>NAD83 Latitude:</td><td>County:</td></tr></table></body></html>"


class Stub(BaseHTTPRequestHandler):
    hits = []

    def log_message(self, *a):  # silence
        pass

    def _route(self):
        u = urlparse(self.path)
        q = {k: v[0] for k, v in parse_qs(u.query).items()}
        Stub.hits.append((self.command, u.path, q))
        if u.path.endswith("las.lasd5.SelectWells"):
            if q.get("f_t") == "35":
                return 200, "text/html", b"<html><body>Service temporarily unavailable</body></html>"
            return 200, "text/html", (SEARCH_HTML if q.get("f_t") == "13" else EMPTY_HTML).encode()
        if u.path == "/blob/kcc_logs_2017/4242.las":
            # simulate a transient network failure on this candidate year: drop the connection
            self.close_connection = True
            raise ConnectionAbortedError("simulated")
        if u.path == "/blob/kcc_logs_2016/4242.las":
            return 200, "text/plain", TINY_LAS
        if u.path.endswith("qualified.well_page.DisplayWell"):
            f = FIX / f"displaywell_{q.get('f_kid')}.html"
            if f.exists():
                return 200, "text/html", f.read_bytes()
            return 200, "text/html", EMPTY_WELL_HTML.encode()
        if u.path.endswith("las.lasd5.ViewLasHeader"):
            if q.get("f_kid") == "1046139243":
                return 200, "text/html", HEADER_HTML.encode()
            return 200, "text/html", b"<html>Kansas LAS files--No such file</html>"
        if u.path == "/blob/kcc_logs_2016/1046139243.las":
            return 200, "text/plain", TINY_LAS
        if u.path == "/blob/kcc_logs_2016/555.las":
            return 200, "text/html", b"<html>oops</html>"
        return 404, "text/plain", b"nope"

    def do_GET(self):
        code, ct, body = self._route()
        self.send_response(code)
        self.send_header("Content-Type", ct)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_HEAD(self):
        code, ct, body = self._route()
        self.send_response(code)
        self.send_header("Content-Type", ct)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()


@pytest.fixture(scope="module")
def stub():
    srv = ThreadingHTTPServer(("127.0.0.1", 0), Stub)
    t = threading.Thread(target=srv.serve_forever, daemon=True)
    t.start()
    base = f"http://127.0.0.1:{srv.server_port}"
    yield {"base_url": f"{base}/ords", "blob_url": f"{base}/blob"}
    srv.shutdown()


def test_parse_search_html():
    rows = kgs.parse_search_html(SEARCH_HTML)
    assert len(rows) == 2
    r = rows[0]
    assert r["kid"] == 1046139243
    assert r["well_kid"] == 1046105344
    assert r["well"] == "PEARSON FAMILY 1-35"
    assert r["api"] == "15-195-23011"
    assert r["operator"] == "Downing-Nelson Oil Co Inc"
    assert r["location"] == "T13S R22W, Sec. 35, SE NW SW NE"
    assert (r["spud"], r["plug"]) == ("14-NOV-2016", "21-NOV-2016")
    assert (r["depth_start"], r["depth_stop"]) == (0.0, 4360.75)
    assert r["las_url"].endswith("/kcc_logs_2016/1046139243.las")
    assert r["las_label"] == "LAS File 1, courtesy KCC"
    assert rows[1]["kid"] == 1046427082 and rows[1]["plug"] == "No Plug Date"
    assert rows[1]["las_url"].endswith("/kcc_logs_2017/1046427082.las")
    assert kgs.parse_search_html(EMPTY_HTML) == []


def test_parse_search_html_multiple_las_per_well():
    # A well with two LAS files spans two <tr>s; the second carries only the file cell.
    extra = (
        '<tr><td>LAS File 2, courtesy KCC<br><a href="x/ViewLasHeader?f_kid=777">h</a><br>'
        '<a href="http://b/kcc_logs_2018/777.las">d</a></td></tr>'
    )
    cut = SEARCH_HTML.rfind("</table>")  # the results table is the last one on the page
    html = SEARCH_HTML[:cut] + extra + SEARCH_HTML[cut:]
    rows = kgs.parse_search_html(html)
    assert [r["kid"] for r in rows] == [1046139243, 1046427082, 777]
    assert rows[2]["well"] == "FABRIZIUS-YOUNGER 'B' 1"  # inherits the last well
    assert rows[2]["las_label"] == "LAS File 2, courtesy KCC"


def test_search_params_validation():
    p = kgs.search_params(13, 22, "W", 35)
    assert (p["f_t"], p["f_r"], p["ew"], p["f_s"], p["f_st"]) == ("13", "22", "W", "35", "15")
    with pytest.raises(ValueError):
        kgs.search_params(township=36)
    with pytest.raises(ValueError):
        kgs.search_params(range_=22)  # no ew
    with pytest.raises(ValueError):
        kgs.search_params(range_=30, ew="E")  # east only goes to 25
    with pytest.raises(ValueError):
        kgs.search_params(section=0)


def test_search_wells(stub):
    base = stub["base_url"]
    rows = kgs.search_wells(13, 22, "W", 35, base_url=base)
    assert [r["kid"] for r in rows] == [1046139243, 1046427082]
    assert kgs.search_wells(1, 1, "W", 36, base_url=base) == []
    with pytest.raises(ValueError):
        kgs.search_wells(base_url=base)  # empty query


def test_search_wells_rejects_unexpected_page(stub):
    # #20: neither marker present -> fail closed rather than parse whatever came back
    with pytest.raises(kgs.KGSError, match="unexpected SelectWells response"):
        kgs.search_wells(35, 1, "W", 1, base_url=stub["base_url"])


def test_las_header_rejects_page_without_version_block():
    # #20: an error page must not be mistaken for a LAS header
    with pytest.raises(kgs.KGSError, match="no LAS header"):
        kgs.las_header(1, base_url="unused", session=_FakeSession("<html>Kansas LAS files--No such file</html>"))


def test_resolve_skips_a_candidate_year_that_errors(stub):
    # #21: a dropped connection on 2017 must not abort; 2016 is served next
    url = kgs.resolve_las_url(4242, years=[2017, 2016], **stub)
    assert url.endswith("/blob/kcc_logs_2016/4242.las")


def test_fetch_refuses_url_off_the_blob_host(stub, tmp_path):
    # #18: a scraped las_url must stay on the expected scheme+host
    evil = "http://169.254.169.254/latest/meta-data/x.las"
    with pytest.raises(kgs.KGSError, match="refusing LAS URL"):
        kgs.fetch_las(1046139243, tmp_path, url=evil, **stub)
    assert not (tmp_path / "1046139243.las").exists()
    with pytest.raises(kgs.KGSError):
        kgs.check_las_url("https://kgsimages.blob.core.windows.net.evil.example/x.las")
    with pytest.raises(kgs.KGSError):
        kgs.check_las_url("http://kgsimages.blob.core.windows.net/x.las")  # scheme downgrade
    assert kgs.check_las_url(kgs.las_url(1, 2016)) == kgs.las_url(1, 2016)


def test_log_year_from_header():
    text = kgs.las_header(1046139243, base_url="unused", session=_FakeSession(HEADER_HTML))
    assert "~Version" in text
    assert kgs.log_year_from_header(text) == 2016
    assert kgs.log_year_from_header("DATE. 11/21/16: LOG DATE") == 2016
    assert kgs.log_year_from_header("no date here") is None


class _FakeSession:
    def __init__(self, text):
        self.text = text

    def get(self, *a, **k):
        r = requests.Response()
        r.status_code = 200
        r._content = self.text.encode()
        return r


def test_resolve_and_fetch(stub, tmp_path):
    Stub.hits.clear()
    url = kgs.resolve_las_url(1046139243, **stub)
    assert url.endswith("/blob/kcc_logs_2016/1046139243.las")
    heads = [h for h in Stub.hits if h[0] == "HEAD"]
    assert len(heads) == 1  # header said 2016, first guess hit

    p = kgs.fetch_las(1046139243, tmp_path, **stub)
    assert p == tmp_path / "1046139243.las"
    assert p.read_bytes() == TINY_LAS
    assert not list(tmp_path.glob("*.part"))
    gets_before = len([h for h in Stub.hits if h[0] == "GET"])
    assert kgs.fetch_las(1046139243, tmp_path, **stub) == p  # cached
    assert len([h for h in Stub.hits if h[0] == "GET"]) == gets_before

    # explicit url skips resolution entirely
    Stub.hits.clear()
    kgs.fetch_las(1046139243, tmp_path, url=url, overwrite=True, **stub)
    assert all(h[0] == "GET" and h[1].endswith(".las") for h in Stub.hits)


def test_fetch_rejects_non_las_body(stub, tmp_path):
    with pytest.raises(kgs.KGSError, match="did not return a LAS"):
        kgs.fetch_las(555, tmp_path, url=f"{stub['blob_url']}/kcc_logs_2016/555.las", **stub)
    assert not (tmp_path / "555.las").exists()
    assert not list(tmp_path.glob("*.part"))


def test_resolve_unknown_kid_gives_up_loudly(stub):
    with pytest.raises(kgs.KGSError, match="no LAS file found"):
        kgs.resolve_las_url(999, years=[2016, 2017], **stub)


def test_parse_well_page():
    w = kgs.parse_well_page(PEARSON_WELL_HTML)
    assert w["well_kid"] == 1046105344 and w["api"].startswith("15-195-23011")
    assert (w["lease"], w["well_no"], w["operator"], w["field"], w["county"]) == ("PEARSON FAMILY", "1-35", "Downing-Nelson Oil Co Inc", "Wildcat", "Trego")
    assert w["lat"] == pytest.approx(38.880121) and w["lon"] == pytest.approx(-99.7321228)
    assert w["lat_nad27"] == pytest.approx(38.8801099)
    assert w["latlon_source"].startswith("calculated from footages")
    assert (w["elevation"], w["elevation_datum"], w["total_depth"]) == (2395.0, "KB", 4350.0)
    assert (w["well_type"], w["status"], w["spud_date"], w["completion_date"]) == ("D&A", "Plugged and Abandoned", "Nov-14-2016", "Nov-21-2016")
    assert "producing_formation" not in w  # blank on a D&A well: must not swallow the next label
    p = kgs.parse_well_page(PBW_WELL_HTML)
    assert p["well_kid"] == 1045079321 and p["producing_formation"] == "ARBUCKLE" and p["latlon_source"] == "from GPS"
    assert p["lat"] == pytest.approx(38.43746) and p["total_depth"] == 3845.0
    assert kgs.parse_well_page(EMPTY_WELL_HTML) == {}


def test_well_info_caches_and_rejects_empty(stub, tmp_path):
    base = stub["base_url"]
    Stub.hits.clear()
    w = kgs.well_info(1046105344, base_url=base, cache_dir=tmp_path)
    assert w["lat"] == pytest.approx(38.880121)
    assert (tmp_path / "wells" / "1046105344.json").exists()
    kgs.well_info(1046105344, base_url=base, cache_dir=tmp_path)  # served from cache
    assert len([h for h in Stub.hits if "DisplayWell" in h[1]]) == 1
    with pytest.raises(kgs.KGSError, match="empty"):
        kgs.well_info(1046139243, base_url=base)  # a LAS KID, not a well KID


def test_search_wells_with_coords_merges_well_pages(stub, tmp_path):
    rows = kgs.search_wells(13, 22, "W", 35, base_url=stub["base_url"], with_coords=True, cache_dir=tmp_path)
    pearson = next(r for r in rows if r["kid"] == 1046139243)
    assert pearson["lat"] == pytest.approx(38.880121) and pearson["elevation"] == 2395.0 and pearson["status"] == "Plugged and Abandoned"
    other = next(r for r in rows if r["kid"] == 1046427082)
    assert "lat" not in other  # no fixture for well 1046140288 -> row left without coords, no abort
    logs = []
    kgs.add_well_info([{"well_kid": 999}], base_url=stub["base_url"], log=logs.append)
    assert logs and "no header info" in logs[0]


@pytest.mark.skipif(not os.environ.get("LASANALYSIS_NETWORK"), reason="set LASANALYSIS_NETWORK=1 to hit KGS")
def test_live_kgs_head_only():
    rows = kgs.search_wells(13, 22, "W", 35)
    assert any(r["kid"] == 1046139243 for r in rows)
    url = kgs.resolve_las_url(1046139243)
    assert url == kgs.las_url(1046139243, 2016)
