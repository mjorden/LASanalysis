import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import build_site  # noqa: E402


def test_build_site_writes_both_wells_and_index(tmp_path):
    written = build_site.build(tmp_path)
    names = sorted(p.name for p in written)
    assert names == ["index.html", "pbw.html", "pearson.html"]
    assert (tmp_path / ".nojekyll").exists()
    index = (tmp_path / "index.html").read_text(encoding="utf-8")
    assert 'href="pearson.html"' in index and 'href="pbw.html"' in index
    assert "15-195-23011" in index and "15-165-22116" in index
    for page, well in (("pearson.html", "Pearson Family #1-35"), ("pbw.html", "PBW #1-32")):
        html = (tmp_path / page).read_text(encoding="utf-8")
        m = re.search(r"const D = (\{.*?\});\nconst P0", html, flags=re.S)
        d = json.loads(m.group(1).replace("<\\/", "</"))
        assert d["title"] == well
        assert "RT" in d["curves"] and "RHOB" in d["curves"]
        assert any(t["curves"] == ["SW"] for t in d["tracks"])
