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


def test_analyze_neutron_matrix_correction_and_sw_models():
    import warnings

    base = _frame(GR=65.0, RHOB=2.55, NPHI=10.0, RT=20.0)
    lime = analyze(base)  # default: neutron and density both limestone -> no shift
    assert lime["PHIN"].iloc[0] == pytest.approx(0.10) and lime.attrs["phin_corrected"] is True
    ss = analyze(base, {"matrix": "sandstone"})
    assert ss["PHIN"].iloc[0] == pytest.approx(0.14)  # limestone-scaled neutron -> sandstone
    assert ss["PHID"].iloc[0] == pytest.approx((2.65 - 2.55) / 1.65)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        salt = analyze(base, {"matrix": "salt"})
    assert salt["PHIN"].iloc[0] == pytest.approx(0.10) and salt.attrs["phin_corrected"] is False
    assert any("no lithology correction" in str(x.message) for x in w)

    # shaly-sand model with an explicit Rsh
    shaly = _frame(GR=90.0, RHOB=2.5, NPHI=15.0, RT=5.0)
    sim = analyze(shaly, {"sw_model": "simandoux", "rsh": 3.0})
    assert sim.attrs["sw_model"] == "simandoux" and sim.attrs["rsh"] == 3.0
    assert "SW_ARCHIE" in sim and (sim["SW"] < sim["SW_ARCHIE"]).all()
    # Rsh auto-pick needs shale samples; a clean frame falls back to Archie with a warning
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        clean = analyze(_frame(GR=30.0, RHOB=2.5, NPHI=15.0, RT=5.0), {"sw_model": "indonesia"})
    assert clean.attrs["sw_model"] == "archie" and "SW_ARCHIE" not in clean
    assert any("falling back to Archie" in str(x.message) for x in w)


def test_summarize_reports_both_rw_picks(tmp_path):
    row = run_well(PEARSON, tmp_path, depth_range=(3400, 4200), plot=False)
    assert 0.02 < row["rw_envelope"] < 0.05 and 1.8 < row["m_envelope"] < 2.2
    assert 0.01 < row["rw_rwa"] < 0.1
    assert "-" in row["rwa_interval"]
    assert row["sw_model"] == "archie" and row["phin_corrected"] is True
    assert json.loads(row["params"])["sw_model"] == "archie" and json.loads(row["params"])["rsh"] is None


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


def test_run_well_html(tmp_path):
    row = run_well(PEARSON, tmp_path, depth_range=(3400, 4200), plot=False, html=True, meta={"kid": 1046139243, "api": "15-195-23011"})
    assert row["html"] == "1046139243.html"
    html = (tmp_path / "1046139243.html").read_text(encoding="utf-8")
    assert '"well":"PEARSON FAMILY #1-35"' in html and '"api":"15-195-23011"' in html
    assert row["png"] == ""


def test_run_well_without_plot(tmp_path):
    row = run_well(PEARSON, tmp_path, plot=False, params={"rw": 0.05, "matrix": "dolomite"})
    assert row["png"] == ""
    got = json.loads(row["params"])
    assert got["rw"] == 0.05 and got["matrix"] == "dolomite" and got["m"] == DEFAULT_PARAMS["m"] and got["rsh"] is None
    assert not (tmp_path / "1046139243.png").exists()


def test_run_search_chains_search_fetch_run_and_survives_failures(tmp_path):
    hits = [
        {"kid": 1046139243, "well": "PEARSON FAMILY 1-35", "api": "15-195-23011", "las_url": "u1", "well_kid": 1046105344},
        {"kid": 999, "well": "BROKEN", "api": "x", "las_url": "u2", "well_kid": 1},
    ]
    calls = []

    def search(**kw):
        calls.append(("search", kw))
        return [dict(h) for h in hits]  # fresh rows: add_well_info mutates in place

    def well_info(rows, cache_dir=None, log=None):
        calls.append(("well_info", [r["well_kid"] for r in rows]))
        for r in rows:
            if r["well_kid"] == 1046105344:
                r.update({"lat": 38.880121, "lon": -99.7321228, "elevation": 2395.0, "status": "Plugged and Abandoned"})
        return rows

    def fetch(kid, cache_dir, url=None):
        calls.append(("fetch", kid, url))
        if kid == 999:
            raise RuntimeError("no such file")
        return PEARSON

    logs = []
    out = run_search({"township": 13, "range_": 22, "ew": "W", "section": 35}, tmp_path, "unused-cache",
                     depth_range=(3400, 4200), search=search, fetch=fetch, log=logs.append, plot=False, well_info=well_info)
    assert calls[0] == ("search", {"township": 13, "range_": 22, "ew": "W", "section": 35})
    assert calls[1] == ("well_info", [1046105344, 1])
    assert ("fetch", 1046139243, "u1") in calls and ("fetch", 999, "u2") in calls
    assert list(out["kid"]) == [1046139243, 999]
    assert out.loc[0, "error"] == "" and out.loc[0, "well_name"] == "PEARSON FAMILY #1-35"
    assert out.loc[0, "lat"] == pytest.approx(38.880121) and out.loc[0, "status"] == "Plugged and Abandoned"
    assert np.isnan(out.loc[1, "lat"])
    assert out.loc[1, "error"].startswith("RuntimeError: no such file")
    assert (tmp_path / "summary.csv").exists()
    assert (tmp_path / "wells.png").exists()  # #31: location map from the wells that have coordinates
    assert any("FAILED" in line for line in logs)

    # coords=False skips the lookup and the map entirely
    calls.clear()
    out2 = run_search({"lease": "PEARSON"}, tmp_path / "nocoords", "unused-cache", search=search, fetch=fetch, log=logs.append,
                      plot=False, coords=False, well_info=well_info)
    assert not any(c[0] == "well_info" for c in calls) and "lat" not in out2.columns
    assert not (tmp_path / "nocoords" / "wells.png").exists()


def test_cli_index_search_offline(tmp_path, monkeypatch):
    # --index: search the offline index; fetch/analysis stubbed so no network is touched
    from lasanalysis import kgs

    fixture = Path(__file__).parent / "fixtures" / "ks_las_files_sample.zip"
    monkeypatch.setattr(multiwell.kgs, "fetch_las", lambda kid, cache_dir, url=None: PEARSON)
    monkeypatch.setattr(multiwell.kgs, "add_well_info", lambda rows, cache_dir=None, log=None: rows)
    rc = multiwell.main(["--index", str(fixture), "--within", "38.88", "-99.73", "5", "--out", str(tmp_path), "--depth", "3400", "4200", "--no-plot"])
    assert rc == 0
    s = pd.read_csv(tmp_path / "summary.csv")
    assert list(s["kid"]) == [1046139243] and s.loc[0, "error"] != s.loc[0, "error"] or s.loc[0, "error"] == ""  # NaN or ""
    assert s.loc[0, "well"] == "PEARSON FAMILY 1-35"
    with pytest.raises(SystemExit):  # --within without --index
        multiwell.main(["--township", "13", "--within", "38.88", "-99.73", "5", "--out", str(tmp_path)])


def test_plot_wells(tmp_path):
    from lasanalysis.multiwell import plot_wells

    s = pd.DataFrame({"well": ["A", "B", "C"], "lat": [38.88, 38.44, np.nan], "lon": [-99.73, -99.23, np.nan], "pay_ft": [7.5, 120.5, 1.0]})
    p = plot_wells(s, tmp_path / "m.png")
    assert p.exists() and p.stat().st_size > 1000
    with pytest.raises(ValueError):
        plot_wells(s.iloc[2:], tmp_path / "n.png")


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
    with pytest.raises(SystemExit):  # #22: unknown matrix is rejected at parse time
        multiwell.main(["--las", str(PEARSON), "--out", str(tmp_path), "--param", "matrix=chalk"])


def test_parse_params():
    from lasanalysis.multiwell import parse_params

    assert parse_params(["rw=0.04", "matrix=Dolomite", "m=2.2"]) == {"rw": 0.04, "matrix": "dolomite", "m": 2.2}
    assert parse_params(["matrix=2.68"]) == {"matrix": 2.68}
    assert parse_params(["sw_model=Simandoux", "rsh=2.5", "neutron_matrix=sandstone"]) == {"sw_model": "simandoux", "rsh": 2.5, "neutron_matrix": "sandstone"}
    with pytest.raises(ValueError, match="unknown Sw model"):
        parse_params(["sw_model=waxman"])
    with pytest.raises(ValueError, match="no neutron lithology correction"):
        parse_params(["neutron_matrix=salt"])
    assert parse_params([]) == {}
    with pytest.raises(ValueError, match="unknown matrix 'chalk'"):
        parse_params(["matrix=chalk"])
    with pytest.raises(ValueError, match="must be a number"):
        parse_params(["rw=abc"])
    with pytest.raises(ValueError, match="unknown param"):
        parse_params(["rw"])
