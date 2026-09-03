"""Headless-browser smoke test of the viewer's JavaScript (#29).

Skipped unless Playwright and its Chromium build are installed
(``pip install playwright && python -m playwright install chromium``); CI does
both. Everything the Python tests cannot see is checked here: the page builds
its Plotly figures, a slider recomputes the derived curves, a zoom on one
track propagates to all of them, and the Sw-model switch changes the numbers.
"""

import json
from pathlib import Path

import pytest

pw = pytest.importorskip("playwright.sync_api")

from lasanalysis.viewer import write_viewer  # noqa: E402

PEARSON = Path(__file__).resolve().parents[1] / "data" / "1046139243.las"


@pytest.fixture(scope="module")
def page_url(tmp_path_factory):
    out = write_viewer(PEARSON, tmp_path_factory.mktemp("viewer") / "pearson.html", depth_range=(3400, 4200))
    # Plotly.js comes from cdnjs; the test needs network for that one script.
    return out.resolve().as_uri()


@pytest.fixture(scope="module")
def browser():
    with pw.sync_playwright() as p:
        try:
            b = p.chromium.launch()
        except Exception as e:  # noqa: BLE001 - browser build missing
            pytest.skip(f"chromium not available: {e}")
        yield b
        b.close()


def _read(page):
    return page.evaluate(
        """() => ({
            traces: document.getElementById('logs').data.length,
            pay: document.getElementById('s_pay').textContent,
            sw: document.getElementById('s_sw').textContent,
            n: document.getElementById('s_n').textContent,
            top: document.getElementById('top').value,
            base: document.getElementById('base').value,
            yranges: Object.entries(document.getElementById('logs')._fullLayout)
                .filter(([k, v]) => k.startsWith('yaxis') && v && v.range).map(([, v]) => v.range.map(Math.round)),
            pickett: document.querySelector('#pickett .gtitle') ? document.querySelector('#pickett .gtitle').textContent : ''
        })"""
    )


def test_viewer_renders_and_reacts(browser, page_url):
    page = browser.new_page()
    errors = []
    page.on("pageerror", lambda e: errors.append(str(e)))
    page.goto(page_url)
    page.wait_for_function("() => window.Plotly && document.getElementById('logs').data && document.getElementById('logs').data.length > 0", timeout=30000)
    s0 = _read(page)
    assert s0["traces"] >= 10
    assert s0["n"] == "3201" and s0["top"] == "3400.0" and s0["base"] == "4200.0"
    assert all(r == [4200, 3400] for r in s0["yranges"]), s0["yranges"]  # every track shares the inverted depth axis

    # m slider: Sw and pay must change
    page.evaluate("() => { const s = document.getElementById('m_s'); s.value = 2.6; s.dispatchEvent(new Event('input', {bubbles: true})); }")
    page.wait_for_timeout(500)
    s1 = _read(page)
    assert s1["pay"] != s0["pay"] and s1["sw"] != s0["sw"]
    assert "m=2.6" in s1["pickett"]

    # zoom one track: all y-axes follow and the depth inputs + window stats update
    page.evaluate("() => Plotly.relayout('logs', {'yaxis.range': [3700, 3560]})")
    page.wait_for_timeout(500)
    s2 = _read(page)
    assert all(r == [3700, 3560] for r in s2["yranges"])
    assert (s2["top"], s2["base"]) == ("3560.0", "3700.0") and s2["n"] == "561"

    # Sw model switch changes the numbers and the Pickett title
    page.select_option("#sw_model", "indonesia")
    page.wait_for_timeout(500)
    s3 = _read(page)
    assert s3["sw"] != s2["sw"] and "indonesia" in s3["pickett"]

    # reset restores the defaults
    page.click("#reset")
    page.wait_for_timeout(500)
    s4 = _read(page)
    assert "m=2," in s4["pickett"] and "archie" in s4["pickett"]

    assert errors == [], errors
    page.close()


def test_embedded_data_matches_page_state(browser, page_url):
    page = browser.new_page()
    page.goto(page_url)
    # D and P are top-level `const`/`let` in the page script: reachable by name, not as window properties
    page.wait_for_function("() => typeof D !== 'undefined' && typeof P !== 'undefined'", timeout=30000)
    d = page.evaluate("() => ({rw: D.params.rw, n: D.depth.length, curves: Object.keys(D.curves)})")
    assert d["rw"] == 0.03 and d["n"] == 17444 and "GR" in d["curves"]
    page.close()
