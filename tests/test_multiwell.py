import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

from lasanalysis import multiwell
from lasanalysis.multiwell import DEFAULT_PARAMS, analyze, default_tracks, run_search, run_well, summarize

DATA = Path(__file__).resolve().parents[1] / "data"
PEARSON = DATA / "1046139243.las"


def _frame(**cols):
    depth = np.arange(1000.0, 1010.0, 0.5)
    return pd.DataFrame({k: np.full(depth.size, v, dtype=float) for k, v in cols.items()}, index=pd.Index(depth, name="DEPT"))


def test_analyze_adds_what_the_curves_allow():
    df = analyze(_frame(GR=65.0, RHOB=2.55, NPHI=10.0, RT=20.0))
    for c in ("VSH", "PHID", "PHIN", "PHIND", "SW", "PAY"):
        assert c in df
    assert df["PHIN"].iloc[0] == pytest.approx(0.10)
    assert df["PHID"].iloc[0] == pytest.approx((2.71 - 2.55) / 1.71)
    # limestone default, Rw=0.03, m=n=2: Sw = sqrt(0.03 / (phi^2 * 20))
    phi = df["PHIND"].iloc[0]
    assert df["SW"].iloc[0] == pytest.approx(min(1.0, np.sqrt(0.03 / (phi**2 * 20))))
    assert df["PAY"].dtype == bool

    partial = analyze(_frame(GR=65.0, RT=20.0))
    assert "VSH" in partial and "SW" not in partial and "PAY" not in partial
    only_density = analyze(_frame(RHOB=2.5, RT=20.0))
    assert "SW" in only_density  # falls back to PHID when there is no neutron


def test_default_tracks_only_uses_present_curves():
    t = default_tracks(["GR", "RT", "RHOB", "VSH", "SW"])
    names = [c for tr in t for c in tr["curves"]]
    assert names == ["GR", "RT", "RHOB", "VSH", "SW"]
    assert default_tracks(["ABHV"]) == []
    assert any(tr.get("twin") for tr in default_tracks(["RHOB", "NPHI"]))


def test_summarize_pay_feet_uses_the_sample_step():
    df = analyze(_frame(GR=30.0, RHOB=2.45, NPHI=15.0, RT=200.0))  # clean, porous, resistive -> pay
    assert df["PAY"].all()
    row = summarize(df)
    assert row["pay_ft"] == pytest.approx(len(df) * 0.5)
    assert (row["pay_top"], row["pay_base"]) == (1000.0, 1009.5)
    assert row["curves"] == "GR,RT,RHOB,NPHI"
    assert json.loads(row["params"])["rw"] == DEFAULT_PARAMS["rw"]
    wet = analyze(_frame(GR=30.0, RHOB=2.45, NPHI=15.0, RT=1.0))
    assert summarize(wet)["pay_ft"] == 0


def test_run_well_on_pearson(tmp_path):
    row = run_well(PEARSON, tmp_path, depth_range=(3400, 4200), meta={"kid": 1046139243})
    assert row["kid"] == 1046139243
    assert row["well_name"] == "PEARSON FAMILY #1-35"
    assert (row["depth_top"], row["depth_base"]) == (3400.0, 4200.0)
    assert row["step"] == 0.25
    assert row["n_samples"] == 3201
    assert 0 < row["phind_mean"] < 0.3
    assert 0 < row["sw_mean"] <= 1
    assert row["pay_ft"] >= 0
    assert json.loads(row["mask_report"])["RILD"]["off_scale"] == 268
    assert (tmp_path / "1046139243.png").exists() and row["png"] == "1046139243.png"
    derived = pd.read_csv(tmp_path / "1046139243_derived.csv", index_col=0)
    assert {"VSH", "PHIND", "SW", "PAY"} <= set(derived.columns)
    assert len(derived) == 3201


def test_run_well_without_plot(tmp_path):
    row = run_well(PEARSON, tmp_path, plot=False, params={"rw": 0.05, "matrix": "dolomite"})
    assert row["png"] == ""
    assert json.loads(row["params"]) == {**{k: DEFAULT_PARAMS[k] for k in ("a", "m", "n", "gr_clean", "gr_dirty")}, "rw": 0.05, "matrix": "dolomite"}
    assert not (tmp_path / "1046139243.png").exists()


def test_run_search_chains_search_fetch_run_and_survives_failures(tmp_path):
    hits = [
        {"kid": 1046139243, "well": "PEARSON FAMILY 1-35", "api": "15-195-23011", "las_url": "u1"},
        {"kid": 999, "well": "BROKEN", "api": "x", "las_url": "u2"},
    ]
    calls = []

    def search(**kw):
        calls.append(("search", kw))
        return hits

    def fetch(kid, cache_dir, url=None):
        calls.append(("fetch", kid, url))
        if kid == 999:
            raise RuntimeError("no such file")
        return PEARSON

    logs = []
    out = run_search({"township": 13, "range_": 22, "ew": "W", "section": 35}, tmp_path, "unused-cache",
                     depth_range=(3400, 4200), search=search, fetch=fetch, log=logs.append, plot=False)
    assert calls[0] == ("search", {"township": 13, "range_": 22, "ew": "W", "section": 35})
    assert ("fetch", 1046139243, "u1") in calls and ("fetch", 999, "u2") in calls
    assert list(out["kid"]) == [1046139243, 999]
    assert out.loc[0, "error"] == "" and out.loc[0, "well_name"] == "PEARSON FAMILY #1-35"
    assert out.loc[1, "error"].startswith("RuntimeError: no such file")
    assert (tmp_path / "summary.csv").exists()
    assert any("FAILED" in line for line in logs)


def test_cli_local_files(tmp_path, capsys):
    rc = multiwell.main(["--las", str(PEARSON), "--out", str(tmp_path), "--depth", "3400", "4200", "--no-plot", "--param", "rw=0.04"])
    assert rc == 0
    s = pd.read_csv(tmp_path / "summary.csv")
    assert len(s) == 1 and json.loads(s.loc[0, "params"])["rw"] == 0.04
    assert "summary.csv" in capsys.readouterr().out
    with pytest.raises(SystemExit):
        multiwell.main(["--out", str(tmp_path)])  # no source
    with pytest.raises(SystemExit):
        multiwell.main(["--las", str(PEARSON), "--out", str(tmp_path), "--param", "bogus=1"])
