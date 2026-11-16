#!/usr/bin/env python3
"""fig3 render - the three-dimensional reimaging UMAP as a rotatable, self-contained HTML page.

The interactive counterpart of `reimagingUmap3d.py`: same bundled table, same two colorings, but you
drag it around instead of picking an azimuth. A 3-D scatter really wants to be rotated - occlusion and
depth are what carry the structure - so this is the honest way to read the manifold, and the static
PNGs are the version that goes on paper.

Self-contained by the package's rules: one HTML file, no CDN, no sibling assets, no network. It opens
offline by double-clicking. Rendering is plain canvas 2-D with the points depth-sorted every frame
(painter's algorithm), which is ample for 3,669 dots and needs no WebGL.

  drag            rotate (turntable: horizontal = azimuth, vertical = elevation)
  shift-drag      pan
  scroll          zoom
  space           start/stop the spin
  double-click    reset the view
  hover           identify a replicate

Colorings are BOTH in the one file and switch in the toolbar - functional groups (the Fig 3B look) or
WT-normalized peak biomass on the plasma ramp (the retired-3A look). Not Panel 3D; see ../README.md.

Reads:  data/reimagingLandscape3d_nn10_md0.10_coords.csv
Writes: figures/reimagingUmap3d_interactive.html
"""
import sys
import json
import argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig3/ for figlib
import matplotlib as mpl
import numpy as np
import pandas as pd
from figlib import config, plotting

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--cmap', default='plasma', help='colormap for the biomass coloring (default plasma)')
ap.add_argument('--out', default='reimagingUmap3d_interactive.html', help='output filename')
args = ap.parse_args()

HS, FC, BG = plotting.HIGHLIGHT_SETS, plotting.FUNCTION_COLORS, plotting.BACKGROUND_COLOR

emb = pd.read_csv(config.TABLES / 'reimagingLandscape3d_nn10_md0.10_coords.csv')
emb['geneLocus'] = emb['geneLocus'].fillna('')
emb['function'] = emb['function'].fillna('')
emb['mutant'] = emb['mutant'].astype(str)

# ── group assignment: 0 = unclassified, 1..6 = functional groups, 7 = WT ──────
groupNames = ['Unclassified'] + [plotting.functionLabel(k) for k in HS] + ['WT']
groupColors = [BG] + [FC[k] for k in HS] + ['#000000']
locusGroup = {locus: i for i, genes in enumerate(HS.values(), start=1) for locus in genes}
group = emb['geneLocus'].map(locusGroup).fillna(0).astype(int)
group[emb['mutant'] == 'WT'] = len(groupNames) - 1

# ── dedupe the long strings: 158 mutants and 48 plates, not 3669 copies each ──
mutantKeys = emb['mutant'].drop_duplicates().tolist()
mutantIdx = {m: i for i, m in enumerate(mutantKeys)}
byMutant = emb.drop_duplicates('mutant').set_index('mutant')
mutantMeta = [[m, str(byMutant.at[m, 'geneLocus']), str(byMutant.at[m, 'function'])[:240]]
              for m in mutantKeys]
plateKeys = emb['plateId'].astype(str).drop_duplicates().tolist()
plateIdx = {p: i for i, p in enumerate(plateKeys)}

# ── the biomass ramp, taken from matplotlib so it matches the static panel ────
cmap = mpl.colormaps[args.cmap]
ramp = [mpl.colors.to_hex(cmap(i / 255)) for i in range(256)]
vmin = float(emb['peakBiomassNorm'].quantile(0.02))
vmax = float(emb['peakBiomassNorm'].quantile(0.98))
bIdx = np.clip((emb['peakBiomassNorm'] - vmin) / (vmax - vmin), 0, 1).mul(255).round().astype(int)

payload = {
    'x': [round(v, 3) for v in emb['umap1']],
    'y': [round(v, 3) for v in emb['umap2']],
    'z': [round(v, 3) for v in emb['umap3']],
    'g': group.tolist(),
    'b': bIdx.tolist(),
    'bv': [round(v, 3) for v in emb['peakBiomassNorm']],
    'mi': [mutantIdx[m] for m in emb['mutant']],
    'pi': [plateIdx[str(p)] for p in emb['plateId']],
    'w': emb['wellId'].astype(str).tolist(),
    'mutants': mutantMeta,
    'plates': plateKeys,
    'groupNames': groupNames,
    'groupColors': groupColors,
    'ramp': ramp,
    'vmin': round(vmin, 3), 'vmax': round(vmax, 3),
    'n': int(len(emb)),
}

TEMPLATE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>Reimaging UMAP - 3D</title>
<style>
  * { box-sizing:border-box; }
  html,body { margin:0; height:100%; }
  body { font-family:'Gillius ADF','Helvetica Neue',Arial,sans-serif; background:#1a1a2e; color:#eee;
         display:flex; flex-direction:column; overflow:hidden; }
  #controls { display:flex; align-items:center; gap:18px; padding:10px 18px; background:#16213e;
              border-bottom:1px solid #334; flex-shrink:0; flex-wrap:wrap; }
  #controls .grp { display:flex; align-items:center; gap:6px; }
  .btn { background:#2a2a4a; border:1px solid #445; border-radius:5px; color:#ccc; padding:5px 13px;
         font-size:14px; cursor:pointer; font-family:inherit; }
  .btn:hover { background:#3a3a5a; }
  .btn.active { background:#4a6fa5; border-color:#6a9fd5; color:#fff; }
  label.lab { font-size:13px; color:#aab; }
  input[type=range] { width:110px; accent-color:#6a9fd5; }
  #legend { display:flex; align-items:center; gap:8px; padding:7px 18px; background:#1a2040;
            border-bottom:1px solid #334; flex-shrink:0; flex-wrap:wrap; font-size:13px; }
  .sw { display:flex; align-items:center; gap:6px; cursor:pointer; padding:2px 8px; border-radius:11px;
        border:1px solid transparent; user-select:none; }
  .sw:hover { border-color:#445; }
  .sw.off { opacity:.34; }
  .dot { width:12px; height:12px; border-radius:50%; border:1px solid #000; display:inline-block; }
  #cbar { display:none; align-items:center; gap:10px; font-size:13px; color:#ccd; }
  #cbar .bar { width:230px; height:13px; border:1px solid #445; border-radius:2px; }
  #stage { flex:1; position:relative; min-height:0; }
  canvas { display:block; width:100%; height:100%; cursor:grab; }
  canvas.drag { cursor:grabbing; }
  #tip { position:absolute; pointer-events:none; display:none; background:#0d1024ee; border:1px solid #4a6fa5;
         border-radius:6px; padding:8px 11px; font-size:13px; max-width:330px; line-height:1.45;
         box-shadow:0 4px 14px #0008; }
  #tip b { color:#8fc0ff; }
  #tip .fn { color:#aab; font-size:12px; display:block; margin-top:3px; }
  #foot { padding:7px 18px; background:#16213e; border-top:1px solid #334; font-size:12px; color:#8890a8;
          flex-shrink:0; }
</style></head><body>

<div id="controls">
  <div class="grp"><span class="lab">Color</span>
    <button class="btn active" id="cFunc">Gene function</button>
    <button class="btn" id="cBio">Peak biomass</button></div>
  <div class="grp"><button class="btn" id="spin">Spin</button>
    <button class="btn" id="reset">Reset view</button></div>
  <div class="grp"><label class="lab" for="size">Point size</label>
    <input type="range" id="size" min="1" max="9" step="0.5" value="3.4"></div>
  <div class="grp"><label class="lab" for="fade">Depth fade</label>
    <input type="range" id="fade" min="0" max="100" step="5" value="55"></div>
  <div id="cbar"><span id="cbarLo"></span><div class="bar" id="cbarBar"></div><span id="cbarHi"></span>
    <span class="lab">normalized peak biofilm biomass</span></div>
</div>

<div id="legend"></div>
<div id="stage"><canvas id="cv"></canvas><div id="tip"></div></div>
<div id="foot">
  <b>Drag</b> rotate &middot; <b>shift-drag</b> pan &middot; <b>scroll</b> zoom &middot;
  <b>space</b> spin &middot; <b>double-click</b> reset &middot; hover a dot to identify it.
  &nbsp;|&nbsp; <span id="count"></span> replicates from a separate <i>n_components</i>&thinsp;=&thinsp;3
  UMAP fit (nn&thinsp;=&thinsp;10, min_dist&thinsp;=&thinsp;0.1) &mdash; not the published 2-D coordinates.
</div>

<script>
const D = __PAYLOAD__;
const cv = document.getElementById('cv'), ctx = cv.getContext('2d'), tip = document.getElementById('tip');
const N = D.n;

// centre + scale the cloud once, so the view maths works in unit-ish space
const mid = a => { let lo = Math.min.apply(null, a), hi = Math.max.apply(null, a); return [(lo+hi)/2, hi-lo]; };
const [cx, sx] = mid(D.x), [cy, sy] = mid(D.y), [cz, sz] = mid(D.z);
const span = Math.max(sx, sy, sz) || 1;
const X = new Float32Array(N), Y = new Float32Array(N), Z = new Float32Array(N);
for (let i = 0; i < N; i++) { X[i] = (D.x[i]-cx)/span; Y[i] = (D.y[i]-cy)/span; Z[i] = (D.z[i]-cz)/span; }

const HOME = { az: -1.047, el: 0.384, zoom: 1, panX: 0, panY: 0 };   // = the static default, elev 22 / azim -60
let view = Object.assign({}, HOME);
let mode = 'func', spinning = false, hoverIdx = -1;
let ptSize = 3.4, fade = 0.55;
const shown = D.groupNames.map(() => true);

// ── legend ───────────────────────────────────────────────────────────────────
const legend = document.getElementById('legend');
D.groupNames.forEach((name, gi) => {
  const el = document.createElement('div');
  el.className = 'sw'; el.dataset.g = gi;
  el.innerHTML = '<span class="dot" style="background:' + D.groupColors[gi] + '"></span>' + name;
  el.onclick = () => { shown[gi] = !shown[gi]; el.classList.toggle('off', !shown[gi]); draw(); };
  legend.appendChild(el);
});
const cbar = document.getElementById('cbar');
document.getElementById('cbarLo').textContent = D.vmin.toFixed(1);
document.getElementById('cbarHi').textContent = D.vmax.toFixed(1);
document.getElementById('cbarBar').style.background =
  'linear-gradient(to right,' + D.ramp.filter((_, i) => i % 8 === 0).join(',') + ')';
document.getElementById('count').textContent = N.toLocaleString();

// ── projection ───────────────────────────────────────────────────────────────
let W = 0, H = 0, dpr = 1, unit = 1;
const px = new Float32Array(N), py = new Float32Array(N), pd = new Float32Array(N), ps = new Float32Array(N);
const order = new Int32Array(N);

function resize() {
  dpr = window.devicePixelRatio || 1;
  const r = cv.getBoundingClientRect();
  W = r.width; H = r.height;
  cv.width = Math.round(W*dpr); cv.height = Math.round(H*dpr);
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  draw();
}

function project() {
  const ca = Math.cos(view.az), sa = Math.sin(view.az);
  const ce = Math.cos(view.el), se = Math.sin(view.el);
  unit = Math.min(W, H) * 0.78 * view.zoom;
  const ox = W/2 + view.panX, oy = H/2 + view.panY;
  for (let i = 0; i < N; i++) {
    const x1 =  X[i]*ca - Y[i]*sa, y1 = X[i]*sa + Y[i]*ca, z1 = Z[i];
    const depth = y1*ce - z1*se;          // toward the viewer
    const up    = y1*se + z1*ce;          // screen vertical
    const persp = 1 / (1 + depth*0.28);   // weak perspective: a depth cue, not a wide lens
    px[i] = ox + x1*unit*persp;
    py[i] = oy - up*unit*persp;
    pd[i] = depth;
    ps[i] = persp;
  }
}

function colorOf(i, alpha) {
  const hex = mode === 'func' ? D.groupColors[D.g[i]] : D.ramp[D.b[i]];
  const n = parseInt(hex.slice(1), 16);
  return 'rgba(' + ((n>>16)&255) + ',' + ((n>>8)&255) + ',' + (n&255) + ',' + alpha.toFixed(3) + ')';
}

function draw() {
  project();
  let m = 0;
  for (let i = 0; i < N; i++) if (shown[D.g[i]]) order[m++] = i;
  const sub = order.subarray(0, m);
  Array.prototype.sort.call(sub, (a, b) => pd[b] - pd[a]);   // far first

  ctx.clearRect(0, 0, W, H);
  let dLo = Infinity, dHi = -Infinity;
  for (let k = 0; k < m; k++) { const d = pd[sub[k]]; if (d < dLo) dLo = d; if (d > dHi) dHi = d; }
  const dRange = (dHi - dLo) || 1;

  for (let k = 0; k < m; k++) {
    const i = sub[k];
    const near = (pd[i] - dLo) / dRange;                 // 0 far, 1 near
    const dim = 1 - fade * (1 - near);
    const isWT = D.g[i] === D.groupNames.length - 1;
    const plain = mode === 'func' && D.g[i] === 0;
    let a = (plain ? 0.20 : isWT ? 0.9 : 0.78) * dim;
    let r = ptSize * ps[i] * (plain ? 0.82 : isWT ? 1.25 : 1.05);
    if (i === hoverIdx) { a = 1; r *= 2.1; }
    ctx.beginPath();
    ctx.arc(px[i], py[i], r, 0, 6.283185);
    ctx.fillStyle = colorOf(i, a);
    ctx.fill();
    if (r > 2.2 && !plain) {
      ctx.lineWidth = Math.min(1, r*0.28);
      ctx.strokeStyle = i === hoverIdx ? '#fff' : 'rgba(0,0,0,' + (0.55*dim).toFixed(2) + ')';
      ctx.stroke();
    }
  }
}

// ── interaction ──────────────────────────────────────────────────────────────
let dragging = false, panning = false, lastX = 0, lastY = 0, moved = 0;
cv.addEventListener('mousedown', e => {
  dragging = true; panning = e.shiftKey; moved = 0;
  lastX = e.clientX; lastY = e.clientY; cv.classList.add('drag'); tip.style.display = 'none';
});
window.addEventListener('mouseup', () => { dragging = false; cv.classList.remove('drag'); });
window.addEventListener('mousemove', e => {
  if (dragging) {
    const dx = e.clientX - lastX, dy = e.clientY - lastY;
    lastX = e.clientX; lastY = e.clientY; moved += Math.abs(dx) + Math.abs(dy);
    if (panning) { view.panX += dx; view.panY += dy; }
    else {
      view.az += dx * 0.0085;
      view.el = Math.max(-1.53, Math.min(1.53, view.el + dy * 0.0085));
    }
    draw();
    return;
  }
  hover(e);
});
cv.addEventListener('wheel', e => {
  e.preventDefault();
  view.zoom = Math.max(0.25, Math.min(14, view.zoom * (e.deltaY < 0 ? 1.11 : 0.9)));
  draw();
}, { passive: false });
cv.addEventListener('dblclick', () => { view = Object.assign({}, HOME); draw(); });
cv.addEventListener('mouseleave', () => { hoverIdx = -1; tip.style.display = 'none'; draw(); });

function hover(e) {
  const r = cv.getBoundingClientRect(), mx = e.clientX - r.left, my = e.clientY - r.top;
  let best = -1, bestDepth = -Infinity; const PICK2 = 13*13;
  for (let i = 0; i < N; i++) {
    if (!shown[D.g[i]]) continue;
    const dx = px[i] - mx, dy = py[i] - my, d2 = dx*dx + dy*dy;
    if (d2 < PICK2 && pd[i] > bestDepth) { best = i; bestDepth = pd[i]; }   // nearest to the viewer wins
  }
  if (best !== hoverIdx) { hoverIdx = best; draw(); }
  if (best < 0) { tip.style.display = 'none'; return; }
  const mt = D.mutants[D.mi[best]];
  tip.innerHTML = '<b>' + mt[0] + '</b>' + (mt[1] ? ' &middot; ' + mt[1] : '') +
    '<br>' + D.plates[D.pi[best]].split(' ')[0] + ' &middot; well ' + D.w[best] +
    '<br>peak biomass ' + D.bv[best].toFixed(2) + '&times; WT' +
    (mt[2] ? '<span class="fn">' + mt[2] + '</span>' : '');
  tip.style.display = 'block';
  const tw = tip.offsetWidth, th = tip.offsetHeight;
  tip.style.left = Math.min(Math.max(4, mx + 16), r.width - tw - 4) + 'px';
  tip.style.top = Math.min(Math.max(4, my + 16), r.height - th - 4) + 'px';
}

// ── toolbar ──────────────────────────────────────────────────────────────────
const bFunc = document.getElementById('cFunc'), bBio = document.getElementById('cBio');
function setMode(m) {
  mode = m;
  bFunc.classList.toggle('active', m === 'func');
  bBio.classList.toggle('active', m === 'bio');
  legend.style.display = m === 'func' ? 'flex' : 'none';
  cbar.style.display = m === 'bio' ? 'flex' : 'none';
  draw();
}
bFunc.onclick = () => setMode('func');
bBio.onclick = () => setMode('bio');
document.getElementById('size').oninput = e => { ptSize = +e.target.value; draw(); };
document.getElementById('fade').oninput = e => { fade = +e.target.value / 100; draw(); };
document.getElementById('reset').onclick = () => { view = Object.assign({}, HOME); draw(); };

const spinBtn = document.getElementById('spin');
let raf = null, last = 0;
function tick(t) {
  if (!spinning) return;
  if (last) { view.az += (t - last) * 0.00035; draw(); }
  last = t; raf = requestAnimationFrame(tick);
}
function setSpin(on) {
  spinning = on; spinBtn.classList.toggle('active', on); last = 0;
  if (on) raf = requestAnimationFrame(tick); else if (raf) cancelAnimationFrame(raf);
}
spinBtn.onclick = () => setSpin(!spinning);
window.addEventListener('keydown', e => {
  if (e.code === 'Space') { e.preventDefault(); setSpin(!spinning); }
});

window.addEventListener('resize', resize);
setMode('func');
resize();
</script></body></html>
"""

html = TEMPLATE.replace('__PAYLOAD__', json.dumps(payload, separators=(',', ':')))
out = config.ensure(config.FIGURES) / args.out
out.write_text(html, encoding='utf-8')
print(f'Saved: {out}  ({out.stat().st_size/1e6:.2f} MB, {payload["n"]} replicates)')
