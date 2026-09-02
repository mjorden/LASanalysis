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


class Stub(BaseHTTPRequestHandler):
    hits = []

    def log_message(self, *a):  # silence
        pass

    def _route(self):
        u = urlparse(self.path)
        q = {k: v[0] for k, v in parse_qs(u.query).items()}
        Stub.hits.append((self.command, u.path, q))
        if u.path.endswith("las.lasd5.SelectWells"):
            return 200, "text/html", (SEARCH_HTML if q.get("f_t") == "13" else EMPTY_HTML).encode()
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


@pytest.mark.skipif(not os.environ.get("LASANALYSIS_NETWORK"), reason="set LASANALYSIS_NETWORK=1 to hit KGS")
def test_live_kgs_head_only():
    rows = kgs.search_wells(13, 22, "W", 35)
    assert any(r["kid"] == 1046139243 for r in rows)
    url = kgs.resolve_las_url(1046139243)
    assert url == kgs.las_url(1046139243, 2016)
