import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from lasanalysis import viewer
from lasanalysis.viewer import build_viewer_html, viewer_data, write_viewer

DATA = Path(__file__).resolve().parents[1] / "data"
PEARSON = DATA / "1046139243.las"


def _embedded(html: str) -> dict:
    m = re.search(r"const D = (\{.*?\});\nconst P0", html, flags=re.S)
    assert m, "data blob not found"
    return json.loads(m.group(1).replace("<\\/", "</"))


@pytest.fixture(scope="module")
def pearson_df():
    from lasanalysis import curves, read_las, standardize

    return standardize(curves(read_las(PEARSON)))


def test_viewer_data_shape(pearson_df):
    d = viewer_data(pearson_df, {"well": "PEARSON"}, {"rw": 0.03}, depth_range=(3400, 4200))
    assert d["title"] == "PEARSON"
    assert len(d["depth"]) == 17444
    assert set(d["curves"]) == {"GR", "SP", "RT", "RM", "RXO", "RHOB", "NPHI", "DPHI", "DT", "CALI"}
    assert all(len(v) == 17444 for v in d["curves"].values())
    assert d["depth_range"] == [3400.0, 4200.0]
    assert d["full_range"] == [0.0, 4360.75]
    assert d["params"]["rw"] == 0.03 and d["params"]["rho_ma"] == 2.71 and d["params"]["matrix"] == "limestone"
    names = [c for t in d["tracks"] for c in t["curves"]]
    assert names[:2] == ["GR", "RT"] and "SW" in names and "VSH" in names and "PHID" in names
    # NaN -> null, never the JSON-invalid NaN literal
    assert None in d["curves"]["RHOB"]
    json.dumps(d, allow_nan=False)


def test_build_html_is_self_contained_and_parseable(pearson_df):
    html = build_viewer_html(pearson_df, {"well": "PEARSON FAMILY #1-35", "api": "15-195-23011"}, depth_range=(3400, 4200))
    assert html.startswith("<!DOCTYPE html>")
    assert "cdnjs.cloudflare.com/ajax/libs/plotly.js/" in html
    assert "NaN" not in html.split("const D = ")[1].split(";\nconst P0")[0]
    d = _embedded(html)
    assert d["meta"]["api"] == "15-195-23011"
    assert d["curves"]["GR"][0] is None or isinstance(d["curves"]["GR"][0], float)
    # ids the JS wires up must exist exactly once
    for el in ("logs", "pickett", "rw", "rw_s", "m_s", "matrix", "top", "base", "s_pay", "reset", "fullrange"):
        assert html.count(f'id="{el}"') == 1, el


def test_title_is_escaped_and_script_close_is_neutralised():
    df = pd.DataFrame({"GR": [10.0, 20.0]}, index=pd.Index([1.0, 2.0], name="DEPT"))
    html = build_viewer_html(df, {"well": "<b>x</b>", "note": "</script><script>alert(1)</script>"}, title="<T>")
    assert "<title>&lt;T&gt;</title>" in html
    assert "</script><script>alert" not in html
    assert _embedded(html)["meta"]["note"].startswith("</script>")


def test_minimal_frame_gets_minimal_tracks():
    df = pd.DataFrame({"GR": [10.0, 20.0, 30.0]}, index=pd.Index([1.0, 2.0, 3.0], name="DEPT"))
    d = viewer_data(df)
    assert list(d["curves"]) == ["GR"]
    assert [t["curves"] for t in d["tracks"]] == [["GR"], ["VSH"]]
    assert d["title"] == "Log viewer"
    empty = pd.DataFrame({"ABHV": [1.0]}, index=pd.Index([1.0], name="DEPT"))
    with pytest.raises(ValueError, match="no recognised curves"):  # #23: never a silent blank page
        viewer_data(empty)


def test_plotlyjs_is_attribute_escaped():
    # #24
    df = pd.DataFrame({"GR": [10.0, 20.0]}, index=pd.Index([1.0, 2.0], name="DEPT"))
    html = build_viewer_html(df, plotlyjs='x.js"></script><script>alert(1)</script>')
    assert '<script src="x.js&quot;&gt;&lt;/script&gt;&lt;script&gt;alert(1)&lt;/script&gt;"></script>' in html
    assert "<script>alert(1)</script>" not in html


def test_write_viewer_from_las_and_cli(tmp_path):
    out = write_viewer(PEARSON, tmp_path / "v.html", depth_range=(3400, 4200))
    d = _embedded(out.read_text(encoding="utf-8"))
    assert d["meta"]["well"] == "PEARSON FAMILY #1-35"
    assert d["meta"]["county"] == "TREGO" and d["meta"]["file"] == "1046139243.las"
    assert d["title"] == "PEARSON FAMILY #1-35"

    rc = viewer.main([str(PEARSON), "-o", str(tmp_path / "cli.html"), "--depth", "3400", "4200", "--param", "rw=0.04", "--param", "matrix=dolomite"])
    assert rc == 0
    d2 = _embedded((tmp_path / "cli.html").read_text(encoding="utf-8"))
    assert d2["params"]["rw"] == 0.04 and d2["params"]["matrix"] == "dolomite" and d2["params"]["rho_ma"] == 2.87
    with pytest.raises(SystemExit):
        viewer.main([str(PEARSON), "-o", str(tmp_path / "x.html"), "--param", "bogus=1"])
