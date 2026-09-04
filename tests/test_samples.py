from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from lasanalysis import samples as S

ROOT = Path(__file__).resolve().parents[1]
SRA = ROOT / "docs" / "templates" / "sra_template.csv"
XRD = ROOT / "docs" / "templates" / "xrd_template.csv"


def test_normalize_api():
    assert S.normalize_api("15-195-23011") == ("15-195-23011", "1519523011")
    assert S.normalize_api("15195230110000") == ("15195230110000", "1519523011")
    assert S.normalize_api(" 15-195-23011-00-00 ") == ("15-195-23011-00-00", "1519523011")
    assert S.normalize_api("1046139243") == ("1046139243", "1046139243")  # ten digits: accepted as given
    assert S.normalize_api("12-345") == ("12-345", None)
    assert S.normalize_api(None) == ("", None) and S.normalize_api(float("nan")) == ("", None)


def test_read_sra_template_to_long():
    long = S.read_samples(SRA)
    assert list(long.columns) == S.SAMPLE_COLUMNS
    assert len(long) == 11  # 6 analytes on row 1 + 5 on row 2 (Ro blank)
    assert set(long["analyte"]) == {"TOC", "S1", "S2", "S3", "Tmax", "Ro"}
    r = long[(long["analyte"] == "TOC") & (long["depth_top"] == 3612.0)].iloc[0]
    assert (r["api10"], r["depth_base"], r["sample_type"], r["lab"], r["method"], r["unit"]) == ("1519523011", 3612.5, "core", "Example Lab", "Rock-Eval 6", "wt%")
    assert r["value"] == 2.41 and r["source"] == "Report 2026-001"
    assert S.validate_samples(long) == []


def test_read_xrd_template_and_validate_totals():
    long = S.read_samples(XRD)
    assert set(long["analyte"]) >= {"quartz", "calcite", "illite", "total_clay"}
    assert S.validate_samples(long) == []
    bad = long.copy()
    bad.loc[bad["analyte"] == "quartz", "value"] += 20  # totals now off by 20
    probs = S.validate_samples(bad)
    assert any("XRD totals" in p for p in probs)
    with pytest.raises(ValueError, match="XRD totals"):
        S.validate_samples(bad, strict=True)


def test_read_samples_aliases_defaults_and_errors():
    wide = pd.DataFrame({"API": ["1519523011"], "Depth": [3600.0], "toc": [1.5], "kspar": [4.0], "GR_core": [55.0], "Notes": ["x"]})
    long = S.read_samples(wide, lab="L", method="M", source="R", sample_type="cuttings")
    assert set(long["analyte"]) == {"TOC", "k_feldspar", "core_GR"}  # aliases resolved, Notes ignored
    assert (long["depth_base"] == 3600.0).all() and (long["lab"] == "L").all() and (long["sample_type"] == "cuttings").all()
    only = S.read_samples(wide, analytes=["TOC"])
    assert list(only["analyte"]) == ["TOC"]
    with pytest.raises(ValueError, match="API column"):
        S.read_samples(pd.DataFrame({"depth": [1.0], "TOC": [1.0]}))
    with pytest.raises(ValueError, match="no analyte columns"):
        S.read_samples(pd.DataFrame({"api": ["1519523011"], "depth": [1.0], "foo": [1.0]}))


def test_validate_catches_ranges_apis_depths_duplicates():
    long = S.read_samples(SRA)
    rows = long.iloc[:3].copy()
    rows.loc[rows.index[0], "value"] = 750.0  # TOC out of range
    rows.loc[rows.index[1], "api10"] = None
    rows.loc[rows.index[2], "depth_base"] = 100.0  # inverted
    probs = S.validate_samples(pd.concat([long, rows, long.iloc[[5]]]))
    text = " | ".join(probs)
    assert "TOC: 1 values outside" in text and "without a usable API" in text and "inverted" in text and "duplicate" in text
    assert S.validate_samples(pd.DataFrame(columns=S.SAMPLE_COLUMNS)) == []
    assert "missing columns" in S.validate_samples(pd.DataFrame({"api": []}))[0]


def test_to_wide_round_trip():
    long = S.read_samples(SRA)
    wide = S.to_wide(long)
    assert len(wide) == 2 and {"TOC", "S2", "Tmax", "depth_mid"} <= set(wide.columns)
    assert wide.loc[0, "depth_mid"] == 3612.25 and wide.loc[1, "TOC"] == 1.05 and np.isnan(wide.loc[1, "Ro"])
    assert S.to_wide(long, analytes=["TOC"]).shape[0] == 2
    assert S.to_wide(long.iloc[0:0]).empty


@pytest.fixture
def log():
    depth = np.arange(3600.0, 3700.0, 0.5)
    gr = 40 + 30 * np.sin(depth / 4.0)  # a wiggly "log" so shifts are distinguishable
    return pd.DataFrame({"GR": gr, "RT": np.full(depth.size, 10.0)}, index=pd.Index(depth, name="DEPT"))


def test_join_to_log_points_and_intervals(log):
    wide = pd.DataFrame({"api10": ["1519523011"] * 3, "depth_top": [3610.0, 3620.0, 3699.9], "depth_base": [3610.0, 3630.0, 3699.9], "TOC": [1.0, 2.0, 3.0]})
    j = S.join_to_log(wide, log, ["GR", "RT"], tolerance=0.5)
    assert j.loc[0, "log_GR"] == pytest.approx(log.loc[3610.0, "GR"])          # point: nearest sample
    assert j.loc[1, "log_GR"] == pytest.approx(log.loc[3620.0:3630.0, "GR"].mean())  # interval: mean
    assert j.loc[2, "log_GR"] == pytest.approx(log.loc[3699.5, "GR"])          # 0.4 ft away, inside tolerance
    assert (j["log_RT"] == 10.0).all()
    far = S.join_to_log(pd.DataFrame({"depth_top": [3800.0], "depth_base": [3800.0]}), log, ["GR"])
    assert np.isnan(far.loc[0, "log_GR"])


def test_depth_shift_recovers_a_known_offset(log):
    # core GR read at depths 3 ft shallower than the log says they are
    true_shift = 3.0
    core_depths = np.arange(3615.0, 3685.0, 2.5)
    core_gr = np.interp(core_depths + true_shift, log.index.to_numpy(), log["GR"].to_numpy()) + np.random.default_rng(0).normal(0, 1, core_depths.size)
    wide = pd.DataFrame({"depth_top": core_depths, "depth_base": core_depths, "core_GR": core_gr})
    fit = S.depth_shift(wide, log, "core_GR", "GR", max_shift=10, step=0.5)
    assert fit["shift"] == pytest.approx(true_shift, abs=0.5)
    assert fit["r"] > 0.95 and fit["n"] == core_depths.size and fit["curve"] == "core_GR vs GR"
    shifted = S.apply_shift(wide, fit["shift"])
    assert shifted.loc[0, "depth_top"] == core_depths[0] + fit["shift"]
    with pytest.raises(ValueError, match="need 5"):
        S.depth_shift(wide.iloc[:3], log, "core_GR", "GR")
    with pytest.raises(ValueError):
        S.depth_shift(wide, log, "core_RHOB", "RHOB")


def test_sample_tracks_on_plot_and_viewer(log):
    import json
    import re

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from lasanalysis.plot import plot_tracks
    from lasanalysis.viewer import build_viewer_html

    long = S.read_samples(pd.DataFrame({"api": ["1519523011"] * 3, "depth": [3610.0, 3640.0, 3680.0], "TOC": [1.0, 2.5, 0.4], "Tmax": [430, 435, 440]}))
    tr = S.sample_tracks(long, ("TOC", "Tmax", "S2"), one_track=False, shift=1.0)
    assert [t["points"][0]["label"] for t in tr] == ["TOC [wt%]", "Tmax [degC]"]  # S2 absent -> skipped
    assert tr[0]["points"][0]["depth"] == [3611.0, 3641.0, 3681.0] and tr[0]["xlim"] == (0, pytest.approx(2.875))
    one = S.sample_tracks(long, ("TOC", "Tmax"))
    assert len(one) == 1 and len(one[0]["points"]) == 2

    fig = plot_tracks(log, [{"curves": ["GR"], "xlim": (0, 100)}] + tr, depth_range=(3600, 3700))
    try:
        assert len(fig.axes) == 3 and fig.axes[1].get_xlim() == (0, pytest.approx(2.875))
        assert len(fig.axes[1].collections) == 1  # the scatter
    finally:
        plt.close(fig)
    with pytest.raises(ValueError, match="curves or points"):
        plot_tracks(log, [{"xlim": (0, 1)}])

    html = build_viewer_html(log, samples=long, sample_analytes=("TOC",))
    d = json.loads(re.search(r"const D = (\{.*?\});\nconst P0", html, flags=re.S).group(1).replace("<\\/", "</"))
    st = [t for t in d["tracks"] if t.get("points")]
    assert len(st) == 1 and st[0]["points"][0]["label"] == "TOC [wt%]" and st[0]["points"][0]["value"] == [1.0, 2.5, 0.4]


def test_samples_on_grid_and_sample_points(log):
    wide = pd.DataFrame({"depth_top": [3610.0, 3620.0], "depth_base": [3610.0, 3630.0], "TOC": [1.0, 2.0]})
    s = S.samples_on_grid(wide, log.index, "TOC")
    assert s.loc[3610.0] == 1.0 and np.isnan(s.loc[3610.5])
    assert (s.loc[3620.0:3630.0] == 2.0).all() and np.isnan(s.loc[3631.0])
    assert s.notna().sum() == 1 + 21
    pts = S.sample_points(wide, "TOC", shift=1.0)
    assert pts == {"depth": [3611.0, 3626.0], "value": [1.0, 2.0]}
    assert S.sample_points(wide, "S2") == {"depth": [], "value": []}
