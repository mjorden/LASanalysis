from pathlib import Path

import pytest

DATA = Path(__file__).resolve().parents[1] / "data"


@pytest.fixture(scope="session")
def pearson_las_path() -> Path:
    """KGS KID 1046139243 — Pearson Family #1-35, Trego County, KS."""
    return DATA / "1046139243.las"


@pytest.fixture(scope="session")
def pearson(pearson_las_path):
    from lasanalysis import read_las

    return read_las(pearson_las_path)
