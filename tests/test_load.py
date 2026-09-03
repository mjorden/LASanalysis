import numpy as np
import pandas as pd
import pytest

from lasanalysis import curves, find_curve, read_las, standardize
from lasanalysis.load import RESISTIVITY_CEILING, clean_curves


def test_declared_null_is_masked_by_lasio(pearson):
    # -999.25 is declared in the header; lasio handles it before we do.
    for c in pearson.curves[1:]:
        assert not np.any(np.isclose(c.data, -999.25))


def test_undeclared_resistivity_sentinel_is_masked(pearson):
    # Raw file: 268 x 1e5 in RILD, 862 in RILM, none in RLL3 (1,130 total).
    rep = pearson.mask_report
    assert rep["RILD"]["off_scale"] == 268
    assert rep["RILM"]["off_scale"] == 862
    assert "RLL3" not in rep or "off_scale" not in rep["RLL3"]
    for m in ("RILD", "RILM", "RLL3"):
        assert np.nanmax(pearson[m]) < RESISTIVITY_CEILING


def test_undeclared_rhob_sentinel_is_masked(pearson):
    # RHOB carries -999.000 (not the declared -999.25) on some rows.
    assert pearson.mask_report["RHOB"]["extra_null"] > 0
    assert np.nanmin(pearson["RHOB"]) > 1.0
    # ...which drove DPOR to tens of thousands of pu on the same rows.
    assert np.nanmax(pearson["DPOR"]) <= 100.0
    assert np.nanmin(pearson["DPOR"]) >= -50.0


def test_density_correction_curve_survives(pearson_las_path, pearson):
    # #17: RHOC (density correction, G/CC, legitimately ~0) must not be range-masked.
    raw = read_las(pearson_las_path, clean=False)
    n_raw = int(np.isfinite(raw["RHOC"]).sum())
    assert n_raw == 3767
    assert int(np.isfinite(pearson["RHOC"]).sum()) == n_raw
    assert "RHOC" not in pearson.mask_report
    assert np.nanmin(pearson["RHOC"]) < 0.5  # real correction values, not a wiped curve


def test_density_rule_is_name_and_unit_gated():
    import lasio

    from lasanalysis.load import rule_for

    assert rule_for("RHOB", "G/CC") == ("density", 1.0)
    assert rule_for("DEN", "kg/m3") == ("density", 1000.0)
    assert rule_for("RHOC", "G/CC") == (None, 1.0)   # correction curve
    assert rule_for("DRHO", "G/CC") == (None, 1.0)
    assert rule_for("RILD", "OHM-M") == ("resistivity", 1.0)
    assert rule_for("RILD", "ohm - m") == ("resistivity", 1.0)  # whitespace-insensitive
    assert rule_for("CNPOR", "PU") == ("porosity_pu", 1.0)
    assert rule_for("SPOR", "P U") == ("porosity_pu", 1.0)
    assert rule_for("RHOB", "G/\tCC") == ("density", 1.0)     # tab in unit (gen-r1-005)
    assert rule_for("GR", "GAPI") == (None, 1.0)
    # no unit (CSV): the mnemonic decides
    assert rule_for("RHOB", None) == ("density", 1.0)
    assert rule_for("RHOC", None) == (None, 1.0)
    assert rule_for("SPOR", None) == ("porosity_pu", 1.0)
    assert rule_for("RILM", "") == ("resistivity", 1.0)

    las = lasio.LASFile()
    las.append_curve("DEPT", np.arange(4.0), unit="FT")
    las.append_curve("RHOC", np.array([0.02, -0.05, 0.31, 0.10]), unit="G/CC")
    las.append_curve("DEN", np.array([2450.0, 2650.0, 0.0, 9999.0]), unit="KG/M3")
    las.append_curve("SPOR", np.array([10.0, 500.0, -999.0, 20.0]), unit="PU")
    rep = clean_curves(las)
    assert "RHOC" not in rep
    assert rep["DEN"] == {"density": 2}          # 0 and 9999 kg/m3 masked, 2450/2650 kept
    assert rep["SPOR"] == {"extra_null": 1, "porosity": 1}
    np.testing.assert_array_equal(las["RHOC"], [0.02, -0.05, 0.31, 0.10])
    assert np.isnan(las["DEN"][2]) and np.isnan(las["DEN"][3]) and las["DEN"][0] == 2450.0


def test_depth_curve_is_never_touched(pearson_las_path):
    raw = read_las(pearson_las_path, clean=False)
    cleaned = read_las(pearson_las_path)
    np.testing.assert_array_equal(raw.curves[0].data, cleaned.curves[0].data)
    assert raw.mask_report == {}


def test_clean_is_idempotent(pearson_las_path):
    las = read_las(pearson_las_path)
    assert clean_curves(las) == {}


def test_clean_curves_rules_on_synthetic(tmp_path):
    import lasio

    las = lasio.LASFile()
    las.append_curve("DEPT", np.arange(5.0), unit="FT")
    las.append_curve("RILD", np.array([1.0, 1e5, 2e5, -999.0, np.nan]), unit="OHM-M")
    las.append_curve("NPHI", np.array([10.0, -999.0, 500.0, -60.0, 20.0]), unit="PU")
    las.append_curve("GR", np.array([1e5, -999.0, 30.0, 40.0, 50.0]), unit="GAPI")
    las.append_curve("RHOB", np.array([2.5, 0.0, -999.0, 4.0, 2.7]), unit="G/CC")
    rep = clean_curves(las)
    assert rep["RILD"] == {"extra_null": 1, "off_scale": 2}
    assert rep["NPHI"] == {"extra_null": 1, "porosity": 2}
    assert rep["GR"] == {"extra_null": 1}  # no ceiling for a non-resistivity curve
    assert rep["RHOB"] == {"extra_null": 1, "density": 2}
    assert np.isnan(las["RILD"]).sum() == 4
    assert np.isnan(las["GR"]).sum() == 1
    assert las["GR"][0] == 1e5


def test_clean_frame_uses_names_when_there_are_no_units():
    from lasanalysis import clean_frame

    df = pd.DataFrame(
        {
            "Depth": [1.0, 2.0, 3.0],
            "RILD": [1.0, 1e5, -999.25],
            "cnpor": [10.0, 500.0, 20.0],
            "RHOB": [2.5, 0.0, 2.6],
            "GR": [1e5, 30.0, -999.25],
            "NAME": ["a", "b", "c"],
        }
    )
    df["RHOC"] = [0.01, 0.2, -0.1]
    df["SPOR"] = [10.0, 300.0, 5.0]
    out, rep = clean_frame(df)
    assert rep == {
        "RILD": {"extra_null": 1, "off_scale": 1},
        "cnpor": {"porosity": 1},
        "RHOB": {"density": 1},
        "GR": {"extra_null": 1},
        "SPOR": {"porosity": 1},
    }
    assert out["RHOC"].tolist() == [0.01, 0.2, -0.1]
    assert out["GR"].iloc[0] == 1e5
    assert df["RILD"].iloc[1] == 1e5  # input untouched
    assert out["Depth"].tolist() == [1.0, 2.0, 3.0]


def test_read_log_csv_second_well():
    from pathlib import Path

    from lasanalysis import read_log_csv

    df = read_log_csv(Path(__file__).resolve().parents[1] / "data" / "1045399712.csv")
    assert df.index.name == "DEPT"
    assert len(df) == 7707
    assert not (df == -999.25).any().any()
    for col in ("RILD", "RILM"):
        assert df[col].max() < RESISTIVITY_CEILING
    assert df["DPOR"].max() <= 100
    assert df["RHOB"].min() >= 1.0
    assert df.attrs["mask_report"]["RILD"]["off_scale"] > 0


def test_curves_dataframe(pearson):
    df = curves(pearson)
    assert isinstance(df, pd.DataFrame)
    assert df.index.name == "DEPT"
    assert "RXORT" in df.columns  # lasio upper-cases RxoRt
    assert len(df) == 17444


def test_find_curve_and_standardize(pearson):
    df = curves(pearson)
    assert find_curve(df.columns, "RT") == "RILD"
    assert find_curve(df.columns, "NPHI") == "CNPOR"
    assert find_curve(df.columns, "CALI") == "DCAL"
    assert find_curve(["FOO"], "RT") is None
    with pytest.raises(KeyError):
        find_curve(df.columns, "NOPE")
    std = standardize(df)
    for name in ("GR", "RT", "RM", "RXO", "RHOB", "NPHI", "DPHI", "DT", "SP", "CALI"):
        assert name in std.columns
    assert "RILD" not in std.columns
    assert "ABHV" in std.columns  # unmatched columns are kept
    # lower-case aliases are matched too
    assert "RT" in standardize(pd.DataFrame({"rild": [1.0]})).columns
