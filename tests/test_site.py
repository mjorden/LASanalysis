import json
import re

from lasanalysis import site


def _embedded(html: str) -> dict:
    m = re.search(r"const D = (\{.*?\});\nconst P0", html, flags=re.S)
    return json.loads(m.group(1).replace("<\\/", "</"))


def test_build_site_writes_both_wells_and_index(tmp_path):
    written = site.build(tmp_path)
    names = sorted(p.name for p in written)
    assert names == ["index.html", "pbw.html", "pearson.html"]
    assert (tmp_path / ".nojekyll").exists()
    index = (tmp_path / "index.html").read_text(encoding="utf-8")
    assert 'href="pearson.html"' in index and 'href="pbw.html"' in index
    assert "15-195-23011" in index and "15-165-22116" in index
    assert "Rw 0.03" in index and "Rw 0.06" in index  # per-well picks shown on the cards
    for page, well, rw in (("pearson.html", "Pearson Family #1-35", 0.03), ("pbw.html", "PBW #1-32", 0.06)):
        html = (tmp_path / page).read_text(encoding="utf-8")
        d = _embedded(html)
        assert d["title"] == well
        assert d["params"]["rw"] == rw and d["params"]["m"] == 2.0
        assert "RT" in d["curves"] and "RHOB" in d["curves"]
        assert any(t["curves"] == ["SW"] for t in d["tracks"])


def test_index_escapes_every_card_field(tmp_path):
    hostile = dict(site.WELLS[0])
    hostile.update({"page": "x.html", "well": "<b>x</b>", "api": '"><script>alert(1)</script>', "location": "a & b",
                    "logged": "<i>", "pick": "</a><a href=evil>"})
    written = site.build(tmp_path, wells=[hostile])
    index = (tmp_path / "index.html").read_text(encoding="utf-8")
    assert "<script>alert" not in index and "<b>x</b>" not in index and "</a><a href=evil>" not in index
    assert "&lt;b&gt;x&lt;/b&gt;" in index and "a &amp; b" in index
    assert len(written) == 2


def test_script_wrapper_delegates():
    import scripts.build_site as wrapper  # noqa: F401  (importable, delegates to lasanalysis.site)

    assert wrapper.main is site.main
