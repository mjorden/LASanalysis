"""Self-contained interactive log viewer: one HTML file, Plotly.js, no server.

The page carries the standard curves as JSON and does the petrophysics in
JavaScript, so Rw / a / m / n, the matrix density, the GR cutoffs and the pay
criteria are live sliders: Vsh, porosity, Sw, the pay shading and the Pickett
panel all update as you drag. The depth axis is shared across tracks (zoom one,
they all zoom) and a spike line gives the readout at the cursor depth.

    python -m lasanalysis.viewer data/1046139243.las -o output/pearson.html
    python -m lasanalysis.viewer data/1046139243.las -o output/pearson.html --depth 3400 4200

Plotly.js is loaded from cdnjs by default (``plotlyjs=``) so the file stays
small; pass a local path or another URL to change that.
"""

from __future__ import annotations

import argparse
import html
import json
import math
import sys
from pathlib import Path
from typing import Dict, Optional

import pandas as pd

from .load import curves, read_las, standardize
from .multiwell import DEFAULT_PARAMS, default_tracks
from .petro import MATRIX_DENSITY, NEUTRON_MATRIX_OFFSET, SW_MODELS, matrix_density

PLOTLY_CDN = "https://cdnjs.cloudflare.com/ajax/libs/plotly.js/3.1.0/plotly.min.js"

#: Raw curves shipped to the page, in this order, when present.
STANDARD_CURVES = ("GR", "SP", "RT", "RM", "RXO", "RHOB", "NPHI", "DPHI", "DT", "CALI")

#: Derived in JavaScript from the raw curves + parameters.
DERIVED = ("VSH", "PHID", "PHIN", "PHIND", "SW")


def _col(series: pd.Series, ndigits: int = 4) -> list:
    out = []
    for v in series.to_numpy(dtype=float):
        out.append(None if not math.isfinite(v) else round(float(v), ndigits))
    return out


def _well_meta(las) -> Dict[str, str]:
    meta = {}
    for key, label in (("WELL", "well"), ("UWI", "uwi"), ("API", "api"), ("CNTY", "county"), ("STAT", "state"),
                       ("COMP", "operator"), ("SRVC", "service"), ("DATE", "date"), ("FLD", "field"), ("LOC", "location")):
        try:
            v = str(las.well[key].value).strip()
        except (KeyError, AttributeError):
            v = ""
        if v:
            meta[label] = v
    return meta


def viewer_data(df: pd.DataFrame, meta: Optional[dict] = None, params: Optional[dict] = None,
                depth_range: Optional[tuple] = None, title: Optional[str] = None,
                samples=None, sample_analytes=("TOC",), sample_shift: float = 0.0) -> dict:
    """Everything the page needs, as plain JSON-able Python.

    ``samples`` (long or wide table from :mod:`lasanalysis.samples`) adds a
    track of lab-sample markers for ``sample_analytes``, shifted by
    ``sample_shift`` ft (see :func:`~lasanalysis.samples.depth_shift`).
    """
    p = {**DEFAULT_PARAMS, **(params or {})}
    present = [c for c in STANDARD_CURVES if c in df.columns and df[c].notna().any()]
    derived_present = []
    if "GR" in present:
        derived_present.append("VSH")
    if "RHOB" in present:
        derived_present.append("PHID")
    if "NPHI" in present:
        derived_present.append("PHIN")
    if "RHOB" in present and "NPHI" in present:
        derived_present.append("PHIND")
    if "RT" in present and ("RHOB" in present or "NPHI" in present):
        derived_present.append("SW")
    tracks = default_tracks(present + derived_present)
    if samples is not None:
        from .samples import sample_tracks

        tracks += sample_tracks(samples, sample_analytes, shift=sample_shift)
    if not tracks:
        raise ValueError(
            "no recognised curves to display; need at least one of " + ", ".join(STANDARD_CURVES) + f" (columns: {list(df.columns)})"
        )
    depth = df.index.to_numpy(dtype=float)
    top, base = (float(depth.min()), float(depth.max())) if depth_range is None else (float(depth_range[0]), float(depth_range[1]))
    matrix = p["matrix"]
    rho_ma = matrix_density(matrix)
    return {
        "title": title or (meta or {}).get("well") or "Log viewer",
        "meta": meta or {},
        "depth": _col(pd.Series(depth), 3),
        "curves": {c: _col(df[c]) for c in present},
        "tracks": tracks,
        "depth_range": [top, base],
        "full_range": [float(depth.min()), float(depth.max())],
        "params": {
            "gr_clean": float(p["gr_clean"]), "gr_dirty": float(p["gr_dirty"]),
            "matrix": matrix if isinstance(matrix, str) else "custom", "rho_ma": rho_ma, "rho_fluid": float(p["rho_fluid"]),
            "neutron_matrix": str(p["neutron_matrix"]).lower(),
            "sw_model": str(p["sw_model"]).lower(),
            "rsh": None if p["rsh"] is None or not math.isfinite(float(p["rsh"])) else float(p["rsh"]),
            "rw": float(p["rw"]), "a": float(p["a"]), "m": float(p["m"]), "n": float(p["n"]),
            "phi_cut": float(p["phi_cut"]), "vsh_cut": float(p["vsh_cut"]), "sw_cut": float(p["sw_cut"]),
        },
        "matrices": MATRIX_DENSITY,
        "neutron_offsets": NEUTRON_MATRIX_OFFSET,
        "sw_models": list(SW_MODELS),
    }


def build_viewer_html(df: pd.DataFrame, meta: Optional[dict] = None, params: Optional[dict] = None,
                      depth_range: Optional[tuple] = None, title: Optional[str] = None,
                      plotlyjs: str = PLOTLY_CDN, **sample_kwargs) -> str:
    """Render the viewer page for a standardised frame. Returns the HTML as a string.

    ``sample_kwargs`` (``samples``, ``sample_analytes``, ``sample_shift``) go to :func:`viewer_data`.
    """
    data = viewer_data(df, meta, params, depth_range, title, **sample_kwargs)
    blob = json.dumps(data, allow_nan=False, separators=(",", ":")).replace("</", "<\\/")
    return (
        _TEMPLATE.replace("__TITLE__", _esc(data["title"]))
        .replace("__PLOTLY__", html.escape(str(plotlyjs), quote=True))  # attribute context
        .replace("__DATA__", blob)
    )


def write_viewer(source, out_html, params: Optional[dict] = None, depth_range: Optional[tuple] = None,
                 title: Optional[str] = None, meta: Optional[dict] = None, plotlyjs: str = PLOTLY_CDN, **sample_kwargs) -> Path:
    """``source`` is a LAS path or an already-standardised DataFrame. Writes ``out_html`` and returns its path.

    Pass ``samples=`` (and optionally ``sample_analytes=``, ``sample_shift=``) to add lab-sample markers.
    """
    if isinstance(source, pd.DataFrame):
        df = source
        meta = meta or {}
    else:
        las = read_las(source)
        df = standardize(curves(las))
        meta = {**_well_meta(las), "file": Path(source).name, **(meta or {})}
    out = Path(out_html)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(build_viewer_html(df, meta, params, depth_range, title, plotlyjs, **sample_kwargs), encoding="utf-8")
    return out


def _esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Write a self-contained interactive log viewer (HTML).")
    ap.add_argument("las")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--depth", nargs=2, type=float, metavar=("TOP", "BASE"))
    ap.add_argument("--param", action="append", default=[], metavar="KEY=VALUE", help=f"override any of {sorted(DEFAULT_PARAMS)}")
    ap.add_argument("--plotlyjs", default=PLOTLY_CDN, help="URL or path of plotly.min.js")
    args = ap.parse_args(argv)
    from .multiwell import parse_params

    try:
        params = parse_params(args.param)
    except ValueError as e:
        ap.error(str(e))
    out = write_viewer(args.las, args.out, params, tuple(args.depth) if args.depth else None, plotlyjs=args.plotlyjs)
    print(f"wrote {out}")
    return 0


_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<script src="__PLOTLY__"></script>
<style>
  :root { --bg:#fafaf8; --fg:#222; --muted:#666; --line:#ddd; --accent:#2a6f97; }
  html, body { margin:0; background:var(--bg); color:var(--fg); font:14px/1.4 Georgia, "Times New Roman", serif; }
  header { padding:10px 16px; border-bottom:1px solid var(--line); display:flex; gap:24px; align-items:baseline; flex-wrap:wrap; }
  header h1 { font-size:18px; margin:0; font-weight:normal; }
  header .meta { color:var(--muted); font-size:12px; }
  main { display:grid; grid-template-columns:250px 1fr; gap:0; min-height:calc(100vh - 44px); }
  aside { border-right:1px solid var(--line); padding:12px; font-size:13px; overflow:auto; }
  aside h2 { font-size:12px; text-transform:uppercase; letter-spacing:.06em; color:var(--muted); margin:14px 0 6px; font-weight:normal; }
  aside h2:first-child { margin-top:0; }
  .row { display:grid; grid-template-columns:1fr 64px; align-items:center; gap:6px; margin:4px 0; }
  .row label { white-space:nowrap; }
  .row input[type=number] { width:100%; box-sizing:border-box; font:inherit; padding:2px 4px; }
  .row input[type=range] { grid-column:1 / span 2; width:100%; margin:0 0 2px; }
  select { width:100%; font:inherit; }
  .stat { display:flex; justify-content:space-between; padding:2px 0; }
  .stat b { font-weight:normal; color:var(--accent); }
  button { font:inherit; font-size:12px; padding:3px 8px; margin-top:6px; }
  #plots { display:grid; grid-template-rows:1fr auto; }
  #logs { width:100%; height:calc(100vh - 44px - 280px); min-height:520px; }
  #pickett { width:100%; height:280px; }
  .fine { color:var(--muted); font-size:11px; }
  @media (max-width:800px) { main { grid-template-columns:1fr; } aside { border-right:0; border-bottom:1px solid var(--line); } }
</style>
</head>
<body>
<header>
  <h1 id="title"></h1>
  <span class="meta" id="meta"></span>
</header>
<main>
<aside>
  <h2>Depth</h2>
  <div class="row"><label>Top</label><input type="number" id="top" step="10"></div>
  <div class="row"><label>Base</label><input type="number" id="base" step="10"></div>
  <button id="fullrange">Full log</button>

  <h2>Shale volume (Larionov, older rocks)</h2>
  <div class="row"><label>GR clean</label><input type="number" id="gr_clean" step="5"></div>
  <div class="row"><label>GR dirty</label><input type="number" id="gr_dirty" step="5"></div>

  <h2>Porosity</h2>
  <div class="row"><label>Matrix</label><select id="matrix"></select></div>
  <div class="row"><label>&rho;<sub>ma</sub> g/cc</label><input type="number" id="rho_ma" step="0.01"></div>
  <div class="row"><label>&rho;<sub>fluid</sub> g/cc</label><input type="number" id="rho_fluid" step="0.05"></div>
  <div class="row"><label>Neutron scale</label><select id="neutron_matrix"></select></div>
  <p class="fine" id="neutron_note"></p>

  <h2>Water saturation</h2>
  <div class="row"><label>Model</label><select id="sw_model"></select></div>
  <div class="row"><label>Rsh &Omega;m</label><input type="number" id="rsh" step="0.1" min="0.01" placeholder="auto"></div>
  <p class="fine" id="rsh_note"></p>
  <div class="row"><label>Rw &Omega;m</label><input type="number" id="rw" step="0.005" min="0.001"><input type="range" id="rw_s" min="-2.3" max="0" step="0.01"></div>
  <div class="row"><label>a</label><input type="number" id="a" step="0.05"><input type="range" id="a_s" min="0.5" max="1.5" step="0.01"></div>
  <div class="row"><label>m</label><input type="number" id="m" step="0.05"><input type="range" id="m_s" min="1.3" max="3" step="0.01"></div>
  <div class="row"><label>n</label><input type="number" id="n" step="0.05"><input type="range" id="n_s" min="1.3" max="3" step="0.01"></div>

  <h2>Pay criteria</h2>
  <div class="row"><label>&phi; &gt;</label><input type="number" id="phi_cut" step="0.01"></div>
  <div class="row"><label>Vsh &lt;</label><input type="number" id="vsh_cut" step="0.05"></div>
  <div class="row"><label>Sw &lt;</label><input type="number" id="sw_cut" step="0.05"></div>

  <h2>In window</h2>
  <div class="stat"><span>Samples</span><b id="s_n"></b></div>
  <div class="stat"><span>Mean &phi;<sub>ND</sub></span><b id="s_phi"></b></div>
  <div class="stat"><span>Mean Sw</span><b id="s_sw"></b></div>
  <div class="stat"><span>Pay</span><b id="s_pay"></b></div>
  <button id="reset">Reset parameters</button>
  <p class="fine">Drag a box or scroll on any track to zoom; double-click to reset. Pickett shows points with Vsh below the pay cutoff.</p>
</aside>
<div id="plots">
  <div id="logs"></div>
  <div id="pickett"></div>
</div>
</main>
<script>
const D = __DATA__;
const P0 = Object.assign({}, D.params);
let P = Object.assign({}, D.params);
const depth = D.depth, C = D.curves, N = depth.length;
const COLORS = {GR:"#2ca02c", SP:"#1f77b4", RT:"#1f77b4", RM:"#d62728", RXO:"#ff7f0e", RHOB:"#d62728", NPHI:"#1f77b4",
                DPHI:"#9467bd", DT:"#bcbd22", CALI:"#7f7f7f", VSH:"#8c564b", PHID:"#d62728", PHIN:"#1f77b4", PHIND:"#2ca02c", SW:"#1f77b4"};
const LOG_TRACES = [];   // {name, traceIndex}
let derived = {};

function derive(p) {
  const VSH=[], PHID=[], PHIN=[], PHIND=[], SW=[], PAY=[];
  for (let i=0;i<N;i++) {
    const gr = C.GR ? C.GR[i] : null; let vsh=null;
    if (gr!=null && p.gr_dirty>p.gr_clean) { let igr=(gr-p.gr_clean)/(p.gr_dirty-p.gr_clean); igr=Math.min(1,Math.max(0,igr)); vsh=Math.min(1,0.33*(Math.pow(2,2*igr)-1)); }
    const rhob = C.RHOB ? C.RHOB[i] : null; const phid = (rhob!=null && p.rho_ma>p.rho_fluid) ? (p.rho_ma-rhob)/(p.rho_ma-p.rho_fluid) : null;
    const nphi = C.NPHI ? C.NPHI[i] : null; let phin = nphi!=null ? nphi/100 : null;
    // #26: neutron is scaled on the logging matrix; shift it onto the density matrix when both are chart lithologies
    if (phin!=null && (p.neutron_matrix in D.neutron_offsets) && (p.matrix in D.neutron_offsets)) phin += D.neutron_offsets[p.matrix] - D.neutron_offsets[p.neutron_matrix];
    const phind = (phid!=null && phin!=null) ? (phid+phin)/2 : null;
    const phi = phind!=null ? phind : phid;
    const rt = C.RT ? C.RT[i] : null; let sw=null;
    if (phi!=null && rt!=null && phi>0 && rt>0) sw = Math.min(1, swModel(p, rt, phi, vsh));
    const pay = phi!=null && sw!=null && vsh!=null && phi>p.phi_cut && vsh<p.vsh_cut && sw<p.sw_cut;
    VSH.push(vsh); PHID.push(phid); PHIN.push(phin); PHIND.push(phind); SW.push(sw); PAY.push(pay);
  }
  return {VSH, PHID, PHIN, PHIND, SW, PAY};
}
function series(name) { return C[name] || derived[name]; }

// #28: Sw models mirroring petro.py. Archie when no Vsh is available.
function archie(p, rt, phi) { return Math.pow(p.a*p.rw/(Math.pow(phi,p.m)*rt), 1/p.n); }
function rshEff(p) { return (p.rsh!=null && p.rsh>0) ? p.rsh : P.rsh_auto; }
function swModel(p, rt, phi, vsh) {
  const model = p.sw_model, rsh = rshEff(p);
  if (model==="archie" || vsh==null || !(rsh>0)) return archie(p, rt, phi);
  const v = Math.min(1, Math.max(0, vsh));
  if (model==="indonesia") {
    const term = Math.pow(v, 1 - v/2)/Math.sqrt(rsh) + Math.pow(phi, p.m/2)/Math.sqrt(p.a*p.rw);
    return Math.pow((1/Math.sqrt(rt))/term, 2/p.n);
  }
  // modified Simandoux: c*Sw^n + d*Sw = 1/Rt
  const c = Math.pow(phi, p.m)/(p.a*p.rw), d = v/rsh;
  if (p.n===2) return (Math.sqrt(d*d + 4*c/rt) - d)/(2*c);
  let sw = Math.min(archie(p, rt, phi), 1) || 1;
  for (let k=0;k<30;k++) { const f = c*Math.pow(sw,p.n) + d*sw - 1/rt, fp = p.n*c*Math.pow(sw,p.n-1) + d; sw = Math.min(10, Math.max(1e-6, sw - f/fp)); }
  return sw;
}
function pickRshAuto() {
  // median Rt where Vsh >= 0.8 (petro.pick_rsh defaults); null when there is no shale
  const d0 = derive(Object.assign({}, P, {sw_model:"archie"}));
  const vals = [];
  for (let i=0;i<N;i++) { const v = d0.VSH[i], rt = C.RT ? C.RT[i] : null; if (v!=null && rt!=null && rt>0 && v>=0.8) vals.push(rt); }
  if (vals.length < 10) return null;
  vals.sort((a,b)=>a-b); return vals[Math.floor(vals.length/2)];
}

function buildLogs() {
  derived = derive(P);
  const tracks = D.tracks, n = tracks.length, gap = 0.012, w = (1 - gap*(n-1)) / n;
  const traces = [], layout = {
    margin: {l:60, r:10, t:60, b:30}, paper_bgcolor:"#fff", plot_bgcolor:"#fff", showlegend:false,
    hovermode:"y", font:{family:"Georgia, serif", size:11},
    yaxis: {range:[D.depth_range[1], D.depth_range[0]], title:{text:"Depth"}, showspikes:true, spikemode:"across", spikesnap:"cursor",
            spikethickness:1, spikecolor:"#999", spikedash:"solid", gridcolor:"#eee", domain:[0,1]},
  };
  let extraAxis = n;  // twin x-axes get ids beyond the track count
  tracks.forEach((t, k) => {
    const ax = k+1, xid = "x"+(ax>1?ax:""), yid = "y"+(ax>1?ax:"");
    const x0 = k*(w+gap);
    const xl = Array.isArray(t.xlim) && Array.isArray(t.xlim[0]) ? t.xlim : (t.curves.length ? t.curves : [null]).map(()=>t.xlim);
    const tcolor = t.twin ? (COLORS[t.curves[0]]||"#333") : "#333";
    layout["xaxis"+(ax>1?ax:"")] = {domain:[x0, x0+w], side:"top", title:{text:t.title || t.curves.join(" / "), font:{size:11}},
      type: t.log ? "log" : "linear", range: t.log ? [Math.log10(xl[0][0]), Math.log10(xl[0][1])] : xl[0], gridcolor:"#eee", zeroline:false,
      tickfont:{size:9, color:tcolor}, showspikes:false, fixedrange:false};
    if (ax>1) layout["yaxis"+ax] = {matches:"y", showticklabels:false, showspikes:true, spikemode:"across", spikesnap:"cursor", spikethickness:1, spikecolor:"#999", spikedash:"solid", gridcolor:"#eee", domain:[0,1]};
    t.curves.forEach((name, j) => {
      let xaxis = xid;
      if (t.twin && j>0) {
        extraAxis += 1;
        const tid = "x"+extraAxis;
        // Twin axis: ticks along the bottom in the curve's colour; the track title
        // (top) already names both curves, so no second title to collide with it.
        layout["xaxis"+extraAxis] = {overlaying:xid, side:"bottom", range: t.log?[Math.log10(xl[j][0]),Math.log10(xl[j][1])]:xl[j],
          type: t.log?"log":"linear", tickfont:{size:9, color:COLORS[name]||"#333"}, showgrid:false, zeroline:false, showspikes:false};
        xaxis = tid;
      }
      const tr = {x: series(name), y: depth, xaxis: xaxis, yaxis: yid, mode:"lines", name:name, line:{width:1, color:COLORS[name]||"#333"},
        hovertemplate: name+": %{x:.3g}<extra></extra>", connectgaps:false};
      if (t.fill && j===0) { tr.fill = "tozerox"; tr.fillcolor = hexA(COLORS[name]||"#333", 0.12); }
      LOG_TRACES.push({name:name, index:traces.length, derived: !!derived[name] && !C[name]});
      traces.push(tr);
    });
    (t.points || []).forEach(pt => {   // lab samples (core / cuttings) as markers
      traces.push({x: pt.value, y: pt.depth, xaxis: xid, yaxis: yid, mode:"markers", name: pt.label,
        marker:{size:7, color: pt.color || "#000", line:{color:"#fff", width:0.8}},
        hovertemplate: pt.label + ": %{x:.3g}<br>depth %{y:.1f}<extra></extra>"});
    });
    if (t.curves[0]==="SW") {
      // Plotly does not break a tozerox fill at null gaps, so pay is a 0/1 step
      // curve filled to x=0: zero width where there is no pay, full width where there is.
      traces.push({x: payX(), y: depth, xaxis: xid, yaxis: yid, mode:"lines", line:{width:0, shape:"vh"}, fill:"tozerox",
        fillcolor:"rgba(0,160,0,0.30)", hoverinfo:"skip", name:"pay"});
      LOG_TRACES.push({name:"PAY", index:traces.length-1, derived:true});
    }
  });
  Plotly.newPlot("logs", traces, layout, {responsive:true, displaylogo:false, scrollZoom:true});
  document.getElementById("logs").on("plotly_relayout", ev => {
    const r = ev["yaxis.range"] || (ev["yaxis.range[0]"]!==undefined ? [ev["yaxis.range[0]"], ev["yaxis.range[1]"]] : null);
    if (r) { setDepthInputs(Math.min(r[0],r[1]), Math.max(r[0],r[1])); stats(); }
    if (ev["yaxis.autorange"]) { setDepthInputs(D.full_range[0], D.full_range[1]); stats(); }
  });
}
function hexA(h, a) { const v=parseInt(h.slice(1),16); return `rgba(${v>>16},${(v>>8)&255},${v&255},${a})`; }
function payX() { return derived.PAY.map(v => v ? 1 : 0); }

function update() {
  derived = derive(P);
  const idx = [], xs = [];
  LOG_TRACES.forEach(t => { if (t.derived) { idx.push(t.index); xs.push(t.name==="PAY" ? payX() : derived[t.name]); } });
  if (idx.length) Plotly.restyle("logs", {x: xs}, idx);
  pickett(); stats();
}

function window_() { const t=+document.getElementById("top").value, b=+document.getElementById("base").value; return [Math.min(t,b), Math.max(t,b)]; }
function inWin(i, w) { return depth[i]>=w[0] && depth[i]<=w[1]; }
function stats() {
  const w = window_(); let n=0, sphi=0, nphi=0, ssw=0, nsw=0, npay=0;
  for (let i=0;i<N;i++) { if (!inWin(i,w)) continue; n++;
    const phi = derived.PHIND[i]!=null ? derived.PHIND[i] : derived.PHID[i]; if (phi!=null) { sphi+=phi; nphi++; }
    if (derived.SW[i]!=null) { ssw+=derived.SW[i]; nsw++; } if (derived.PAY[i]) npay++; }
  const step = N>1 ? Math.abs(depth[1]-depth[0]) : 0;
  document.getElementById("s_n").textContent = n;
  document.getElementById("s_phi").textContent = nphi ? (100*sphi/nphi).toFixed(1)+" pu" : "–";
  document.getElementById("s_sw").textContent = nsw ? (ssw/nsw).toFixed(2) : "–";
  document.getElementById("s_pay").textContent = (npay*step).toFixed(1)+" ft";
}

function pickett() {
  if (!C.RT || !(C.RHOB || C.NPHI)) { document.getElementById("pickett").style.display="none"; return; }
  const w = window_(), px=[], py=[], pc=[];
  for (let i=0;i<N;i++) { if (!inWin(i,w)) continue;
    const phi = derived.PHIND[i]!=null ? derived.PHIND[i] : derived.PHID[i], rt = C.RT[i], vsh = derived.VSH[i];
    if (phi!=null && rt!=null && phi>0 && rt>0 && (vsh==null || vsh<P.vsh_cut)) { px.push(phi); py.push(rt); pc.push(depth[i]); } }
  const traces = [{x:px, y:py, mode:"markers", marker:{size:3, color:pc, colorscale:"Viridis", opacity:0.6, colorbar:{title:{text:"Depth"}, thickness:8, len:0.9}}, name:"points",
                   hovertemplate:"phi %{x:.3f}<br>Rt %{y:.2f}<br>depth %{marker.color:.1f}<extra></extra>"}];
  const phis = []; for (let e=-2;e<=0.0001;e+=0.05) phis.push(Math.pow(10,e));
  [1,0.5,0.25,0.1].forEach(sw => traces.push({x:phis, y:phis.map(f=>P.a*P.rw/(Math.pow(f,P.m)*Math.pow(sw,P.n))), mode:"lines", line:{width:1}, name:"Sw="+sw, hoverinfo:"name"}));
  Plotly.react("pickett", traces, {margin:{l:50,r:10,t:26,b:36}, paper_bgcolor:"#fff", plot_bgcolor:"#fff", font:{family:"Georgia, serif", size:11},
    title:{text:`Pickett (Rw=${P.rw.toFixed(3)}, a=${P.a}, m=${P.m}, n=${P.n}; Sw track: ${P.sw_model})`, font:{size:12}}, showlegend:true, legend:{font:{size:9}, x:1.02},
    xaxis:{type:"log", range:[-2,0], title:{text:"phi (N-D avg or density)"}, gridcolor:"#eee"}, yaxis:{type:"log", range:[-0.5,3], title:{text:"Rt"}, gridcolor:"#eee"}},
    {responsive:true, displaylogo:false});
}

function setDepthInputs(t, b) { document.getElementById("top").value = t.toFixed(1); document.getElementById("base").value = b.toFixed(1); }
function applyDepth() { const w = window_(); Plotly.relayout("logs", {"yaxis.range":[w[1], w[0]]}); pickett(); stats(); }

const NUM = ["gr_clean","gr_dirty","rho_ma","rho_fluid","rw","a","m","n","phi_cut","vsh_cut","sw_cut"];
function fillControls() {
  NUM.forEach(k => document.getElementById(k).value = P[k]);
  ["a","m","n"].forEach(k => document.getElementById(k+"_s").value = P[k]);
  document.getElementById("rw_s").value = Math.log10(P.rw);
  const sel = document.getElementById("matrix"); sel.innerHTML = "";
  Object.entries(D.matrices).forEach(([k,v]) => { const o=document.createElement("option"); o.value=k; o.textContent=`${k} (${v})`; sel.appendChild(o); });
  const o=document.createElement("option"); o.value="custom"; o.textContent="custom"; sel.appendChild(o);
  sel.value = P.matrix in D.matrices ? P.matrix : "custom";
  const nsel = document.getElementById("neutron_matrix"); nsel.innerHTML = "";
  Object.keys(D.neutron_offsets).forEach(k => { const o=document.createElement("option"); o.value=k; o.textContent=k; nsel.appendChild(o); });
  nsel.value = P.neutron_matrix in D.neutron_offsets ? P.neutron_matrix : "limestone";
  const msel = document.getElementById("sw_model"); msel.innerHTML = "";
  D.sw_models.forEach(k => { const o=document.createElement("option"); o.value=k; o.textContent=k; msel.appendChild(o); });
  msel.value = P.sw_model;
  document.getElementById("rsh").value = (P.rsh!=null && P.rsh>0) ? P.rsh : "";
  notes();
}
function notes() {
  const nn = document.getElementById("neutron_note");
  nn.textContent = (P.matrix in D.neutron_offsets) ? `Neutron shifted ${((D.neutron_offsets[P.matrix]-D.neutron_offsets[P.neutron_matrix])*100).toFixed(0)} pu from the ${P.neutron_matrix} scale to ${P.matrix}.`
                                                   : `No neutron lithology correction for matrix "${P.matrix}" — neutron left on the ${P.neutron_matrix} scale.`;
  const rn = document.getElementById("rsh_note");
  if (P.sw_model==="archie") rn.textContent = "Archie ignores shale conductivity; pay flag relies on the Vsh cutoff.";
  else rn.textContent = (P.rsh!=null && P.rsh>0) ? `${P.sw_model} with Rsh = ${P.rsh}` : (P.rsh_auto ? `${P.sw_model}; Rsh auto-picked as ${P.rsh_auto.toFixed(2)} Ωm (median Rt where Vsh ≥ 0.8)` : `${P.sw_model} needs Rsh — no shale in this log to pick it from; Archie used`);
}
function wire() {
  NUM.forEach(k => document.getElementById(k).addEventListener("input", e => { const v=+e.target.value; if (!isFinite(v)) return; P[k]=v;
    if (k==="rho_ma") { document.getElementById("matrix").value="custom"; P.matrix="custom"; notes(); }
    if (k==="rw") document.getElementById("rw_s").value=Math.log10(v); if (["a","m","n"].includes(k)) document.getElementById(k+"_s").value=v; update(); }));
  ["a","m","n"].forEach(k => document.getElementById(k+"_s").addEventListener("input", e => { P[k]=+e.target.value; document.getElementById(k).value=P[k].toFixed(2); update(); }));
  document.getElementById("rw_s").addEventListener("input", e => { P.rw=Math.pow(10,+e.target.value); document.getElementById("rw").value=P.rw.toFixed(3); update(); });
  document.getElementById("matrix").addEventListener("change", e => { const k=e.target.value; P.matrix=k; if (k in D.matrices) { P.rho_ma=D.matrices[k]; document.getElementById("rho_ma").value=P.rho_ma; } notes(); update(); });
  document.getElementById("neutron_matrix").addEventListener("change", e => { P.neutron_matrix=e.target.value; notes(); update(); });
  document.getElementById("sw_model").addEventListener("change", e => { P.sw_model=e.target.value; notes(); update(); });
  document.getElementById("rsh").addEventListener("input", e => { const v=+e.target.value; P.rsh = (isFinite(v) && v>0) ? v : null; notes(); update(); });
  document.getElementById("top").addEventListener("change", applyDepth);
  document.getElementById("base").addEventListener("change", applyDepth);
  document.getElementById("fullrange").addEventListener("click", () => { setDepthInputs(D.full_range[0], D.full_range[1]); applyDepth(); });
  document.getElementById("reset").addEventListener("click", () => { P = Object.assign({}, P0); fillControls(); update(); });
}

document.getElementById("title").textContent = D.title;
document.getElementById("meta").textContent = Object.entries(D.meta).filter(([k])=>k!=="well").map(([k,v])=>`${k}: ${v}`).join("  ·  ");
setDepthInputs(D.depth_range[0], D.depth_range[1]);
P.rsh_auto = pickRshAuto(); P0.rsh_auto = P.rsh_auto;
fillControls(); wire(); buildLogs(); pickett(); stats();
</script>
</body>
</html>
"""


if __name__ == "__main__":
    sys.exit(main())
