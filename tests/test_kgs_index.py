"""Offline KGS LAS index (#30) against a 10-row cut of the real ks_las_files.zip (2026-09-03)."""

from pathlib import Path

import pandas as pd
import pytest

from lasanalysis import kgs

FIX = Path(__file__).parent / "fixtures" / "ks_las_files_sample.zip"


@pytest.fixture(scope="module")
def index():
    return kgs.load_index(FIX)


def test_load_index_shape_and_types(index):
    assert len(index) == 10
    assert list(index.columns[:4]) == ["kid", "well_kid", "well", "api"]
    pearson = index[index["well_kid"] == 1046105344].iloc[0]
    assert pearson["kid"] == 1046139243
    assert pearson["well"] == "PEARSON FAMILY 1-35" and pearson["api"] == "15-195-23011"
    assert pearson["lat_nad27"] == pytest.approx(38.8801099) and pearson["lon_nad27"] == pytest.approx(-99.7317284)
    assert (pearson["township"], pearson["range"], pearson["ew"], pearson["section"]) == (13, 22, "W", 35)
    assert (pearson["elevation"], pearson["elevation_datum"]) == (2395.0, "KB")
    assert (pearson["depth_start"], pearson["depth_stop"]) == (0.0, 4360.75)
    assert pearson["las_path"] == "kcc_logs_2016/1046139243.las"
    assert pearson["las_url"] == kgs.las_url(1046139243, 2016)  # index host is dead; path tail lives on the blob
    # a pre-KCC log filed under a township folder, and a blank API
    old = index.iloc[0]
    assert old["las_path"] == "01S02E/1020069094.las" and old["kid"] == 1020069094
    assert old["api"] == "" and pd.isna(old["township"]) is False


def test_index_las_path_handles_folderless_rows():
    assert kgs.index_las_path("https://www.kgs.ku.edu/b_1/WebDocs/WellLogs/01S02E/1020069094.las") == "01S02E/1020069094.las"
    assert kgs.index_las_path("https://www.kgs.ku.edu//1056890332.las") is None
    assert kgs.index_las_path("") is None


def test_search_index_filters(index):
    rows = kgs.search_index(index, township=13, range_=22, ew="W", section=35)
    assert [r["kid"] for r in rows] == [1046139243]
    assert rows[0]["las_url"].endswith("/kcc_logs_2016/1046139243.las") and rows[0]["header_url"].endswith("f_kid=1046139243")
    assert [r["well"] for r in kgs.search_index(index, lease="pbw")] == ["PBW 1-32"]
    assert len(kgs.search_index(index, operator="downing")) == 2
    assert [r["kid"] for r in kgs.search_index(index, api="15-165-22116")] == [1045399712]
    assert [r["kid"] for r in kgs.search_index(index, api="15165")] == [1045399712]
    near = kgs.search_index(index, within=(38.88, -99.73, 5))   # 5 km around Pearson
    assert [r["kid"] for r in near] == [1046139243] and near[0]["distance_km"] < 1
    far = kgs.search_index(index, within=(38.88, -99.73, 80))   # PBW is ~66 km away
    assert [r["kid"] for r in far][:2] == [1046139243, 1045399712]
    assert kgs.search_index(index, township=1, range_=1, ew="E", section=36) == []
    with pytest.raises(ValueError):
        kgs.search_index(index, township=99)


def test_load_index_rejects_unexpected_columns(tmp_path):
    bad = tmp_path / "x.txt"
    bad.write_text('"A","B"\n"1","2"\n', encoding="utf-8")
    with pytest.raises(kgs.KGSError, match="missing columns"):
        kgs.load_index(bad)


def test_fetch_index_caches_and_validates(tmp_path):
    class Resp:
        def __init__(self, content):
            self.content = content

        def raise_for_status(self):
            pass

    class S:
        calls = 0

        def get(self, *a, **k):
            S.calls += 1
            return Resp(FIX.read_bytes())

    p = kgs.fetch_index(tmp_path / "idx.zip", session=S())
    assert p.exists() and S.calls == 1
    kgs.fetch_index(tmp_path / "idx.zip", session=S())
    assert S.calls == 1  # fresh cache reused
    kgs.fetch_index(tmp_path / "idx.zip", session=S(), max_age_days=0)
    assert S.calls == 2

    class Bad(S):
        def get(self, *a, **k):
            return Resp(b"<html>nope</html>")

    with pytest.raises(kgs.KGSError, match="zip"):
        kgs.fetch_index(tmp_path / "bad.zip", session=Bad())
