#!/usr/bin/env python3
"""interactive build — Interactive Plot 3: dendrogram + heatmap timelapse of the reimaging atlas.

The interactive counterpart of Fig. S6A: the same Ward dendrogram over PCA-50 per-mutant centroids
and the same z-scored heatmap, but the heatmap animates across the imaging timecourse and every
mutant is clickable, opening its replicates' images — which FOLLOW the time slider, so a genotype can
be watched developing next to its features.

Self-contained single-file HTML (thumbnails inlined as base64 data URIs; no server, no asset
folder). Cell colors are precomputed with the same RdBu_r +-3 scale as the static figure, so the
interactive and printed heatmaps match exactly.

Carrying every frame of every well is not shippable (~330 MB of JPEG, ~440 MB base64), so the
timelapse pack holds a few representative replicates per mutant at every frame — the wells whose peak
biomass is closest to the mutant's median peak. Regenerate it at a different trade-off with
`v2/reimaging/umap/exportTimelapsePack.py --reps N --px P --quality Q`.

Reads:  data/fullAtlas_pcaLinkage_{linkage.npy, cluster_order.csv},
        data/fullAtlas_peakBiomass_featureMatrix.csv, timelapse thumbnail pack (config.input)
Writes: figures/interactivePlot3_dendrogramTimelapse.html

Usage:  python build/IP3_interactiveDendrogramTimelapse.py
"""
import sys
import json
import base64
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> interactive/ for figlib
import numpy as np
import pandas as pd
import matplotlib as mpl
from scipy.cluster.hierarchy import dendrogram
from figlib import config, plotting

VMIN, VMAX = -3, 3
outHtml = config.ensure(config.FIGURES / 'interactivePlot3_dendrogramTimelapse.html')

# ── clustering + heatmap, exactly as the static panel ─────────────────────────
order = pd.read_csv(config.TABLES / 'fullAtlas_pcaLinkage_cluster_order.csv')
linkage = np.load(config.TABLES / 'fullAtlas_pcaLinkage_linkage.npy')
fm = pd.read_csv(config.TABLES / 'fullAtlas_peakBiomass_featureMatrix.csv').rename(
    columns={'Unnamed: 0': 'feature'}).set_index('feature')

featureNames = list(fm.index)
# heatmap columns are keyed by `mutant` (the gene name where one exists, else the locus), which is
# also the leaf label -- NOT by `gene`, which is the locus and only matches for named genes
mat = fm[list(order['mutant'])].to_numpy(dtype=float).T    # -> mutants (leaf order) x features

cmap = mpl.colormaps['RdBu_r']
# 255 color steps over [-3, 3] plus a reserved index 255 for "no data" (colony features before t5).
# Cells ship as base64 palette indices -- one byte each -- so all 32 timepoints cost ~180 kB rather
# than ~1.4 MB of hex strings, while still matching the static figure's colors exactly.
PALETTE = [mpl.colors.to_hex(cmap(i / 254)) for i in range(255)] + ['#2a2a3a']


def encodeFrame(m):
    idx = np.full(m.shape, 255, dtype=np.uint8)
    ok = np.isfinite(m)
    idx[ok] = np.clip(np.round((m[ok] - VMIN) / (VMAX - VMIN) * 254), 0, 254).astype(np.uint8)
    return base64.b64encode(idx.tobytes()).decode()


# frames 0-30 from step 0, then the Fig. S6A peak-biomass matrix as the final slider position
npz = np.load(config.TABLES / 'fullAtlas_frameMatrices.npz', allow_pickle=False)
frameMats = npz['matrices'].astype(np.float32)
framesData = {str(int(t)): encodeFrame(frameMats[i]) for i, t in enumerate(npz['frames'])}
framesData['peak'] = encodeFrame(mat)
frameKeys = [str(int(t)) for t in npz['frames']] + ['peak']
frameLabels = [f'{int(t)} h' for t in npz['frames']] + ['Peak biomass']

# dendrogram geometry, computed here so the page needs no clustering code
dd = dendrogram(linkage, no_plot=True, no_labels=True)
dmax = max((max(d) for d in dd['dcoord']), default=1.0) or 1.0
segments = [{'i': [x / 10.0 - 0.5 for x in ic], 'd': [y / dmax for y in dc]}
            for ic, dc in zip(dd['icoord'], dd['dcoord'])]

# ── replicate thumbnails (every frame of a few representative wells), grouped under their mutant ──
thumbList, byMutant, missing = [], {}, 0
packDir = Path(config.input('reimaging/thumbnailsTimelapse/'))
man = pd.read_csv(packDir / 'manifest.csv')
# the pack is built from the wide table, which carries more mutants than survive the manifold
# filters; inlining images for leaves the dendrogram does not show is pure dead payload
leaves = set(order['mutant'].astype(str))
man = man[man['mutant'].astype(str).isin(leaves)]
for mutant, grp in man.groupby('mutant'):
    idxs = []
    for (plate, well), wg in grp.groupby(['plateId', 'wellId'], sort=True):
        wg = wg.sort_values('frame')
        byFrame = {}
        for _, r in wg.iterrows():
            jpg = packDir / str(r['file'])
            if not jpg.exists():
                missing += 1
                continue
            byFrame[int(r['frame'])] = len(thumbList)
            thumbList.append('data:image/jpeg;base64,' + base64.b64encode(jpg.read_bytes()).decode())
        if not byFrame:
            continue
        pf = int(wg['peakFrame'].iloc[0])
        peakIdx = byFrame.get(pf, next(iter(byFrame.values())))
        idxs.append({'t': peakIdx, 'w': f'{plate} · {well}', 'pf': pf, 'bio': None,
                     'fr': {str(k): v for k, v in byFrame.items()}})
    byMutant[str(mutant)] = idxs

rows = []
for i, r in order.iterrows():
    m = str(r['mutant'])
    reps = byMutant.get(m, [])
    rows.append({'mutant': m, 'display': str(r['display']), 'locus': str(r['gene']),
                 'anno': str(r['annotation']), 'color': str(r['color']),
                 'pf': int(r['peakFrame']), 'n': len(reps), 'reps': reps})
print(f'{len(rows)} mutants, {len(thumbList)} replicate thumbnails inlined ({missing} missing); '
      f'heatmap {mat.shape[0]}x{mat.shape[1]}')

# feature-group brackets, matching the static figure's grouping
def groupOf(f):
    if f.startswith('biomass'):
        return 'Biofilm Biomass'
    return 'Colony Features' if f.startswith(('colony_', 'nColonies')) else 'Whole-Image Haralick Features'


groups, start = [], 0
for i, f in enumerate(featureNames):
    if i and groupOf(f) != groupOf(featureNames[i - 1]):
        groups.append({'name': groupOf(featureNames[i - 1]), 'a': start, 'b': i})
        start = i
groups.append({'name': groupOf(featureNames[-1]), 'a': start, 'b': len(featureNames)})

payload = {'rows': rows, 'features': featureNames, 'segments': segments, 'groups': groups,
           'palette': PALETTE, 'framesData': framesData, 'frameKeys': frameKeys,
           'frameLabels': frameLabels, 'vmin': VMIN, 'vmax': VMAX}
dataJson = json.dumps(payload)
thumbsJson = json.dumps(thumbList)
legendJson = json.dumps({**plotting.FUNCTION_COLORS, 'WT': '#000000'})

html = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8">
<title>Reimaging atlas - interactive dendrogram and heatmap</title>
<style>
  * {{ margin:0; padding:0; box-sizing:border-box; }}
  body {{ font-family:'Gillius ADF','Helvetica Neue',Arial,sans-serif; background:#1a1a2e; color:#eee; display:flex; flex-direction:column; height:100vh; }}
  #controls {{ display:flex; align-items:center; gap:16px; padding:12px 20px; background:#16213e; border-bottom:1px solid #334; flex-shrink:0; flex-wrap:wrap; }}
  #controls label {{ font-size:15px; color:#99a; }}
  #search {{ background:#2a2a4a; border:1px solid #445; border-radius:5px; color:#eee; font-size:15px; padding:6px 12px; width:260px; outline:none; }}
  #search::placeholder {{ color:#667; }}
  .btn {{ background:#2a2a4a; border:1px solid #445; border-radius:5px; color:#ccc; padding:6px 14px; font-size:14px; cursor:pointer; }}
  .btn:hover {{ background:#3a3a5a; }}
  input[type=range] {{ width:170px; accent-color:#4a6fa5; cursor:pointer; }}
  #fslider {{ width:90px; }}
  .ro {{ font-size:14px; color:#ccd; min-width:74px; display:inline-block; }}
  #cbar {{ display:flex; align-items:center; gap:8px; margin-left:auto; font-size:13px; color:#99a; }}
  #cbar .ramp {{ width:160px; height:12px; border:1px solid #556; border-radius:2px;
                 background:linear-gradient(to right,#053061,#f7f7f7,#67001f); }}
  #main {{ flex:1; display:flex; overflow:hidden; }}
  #plot {{ flex:1; position:relative; overflow:hidden; }}
  canvas {{ display:block; width:100%; height:100%; }}
  #tooltip {{ display:none; position:absolute; pointer-events:none; background:rgba(20,20,40,0.96); border:1px solid #556; border-radius:6px; padding:8px 12px; font-size:14px; line-height:1.45; z-index:10; }}
  #tooltip .g {{ font-weight:bold; color:#7ecfff; font-size:16px; }}
  #hint {{ position:absolute; top:10px; right:14px; font-size:13px; color:#667; z-index:5; }}
  #sidebar {{ width:520px; background:#16213e; border-left:1px solid #334; display:flex; flex-direction:column; overflow:hidden; }}
  #sidebar.hidden {{ display:none; }}
  #sb-head {{ padding:12px 16px; border-bottom:1px solid #334; position:relative; }}
  #sb-head .g {{ font-weight:bold; font-size:19px; color:#7ecfff; font-style:italic; }}
  #sb-head .l {{ font-size:13px; color:#9bd; }}
  #sb-head .a {{ font-size:13px; color:#dda; font-style:italic; }}
  #sb-head .m {{ font-size:13px; color:#aab; }}
  #sb-close {{ position:absolute; top:8px; right:12px; background:none; border:none; color:#889; font-size:22px; cursor:pointer; }}
  #sb-close:hover {{ color:#eee; }}
  #gallery {{ flex:1; overflow-y:auto; padding:10px; display:grid; grid-template-columns:repeat(auto-fill,minmax(112px,1fr)); gap:8px; align-content:start; }}
  .cell {{ background:#000; border:1px solid #334; border-radius:4px; overflow:hidden; }}
  .cell img {{ display:block; width:100%; }}
  .cell .cap {{ font-size:10px; color:#889; padding:3px 4px; line-height:1.3; }}
  #empty {{ color:#667; font-size:14px; padding:20px; }}
</style></head><body>
<div id="controls">
  <button class="btn" id="play">&#9654; Play</button>
  <label>Time</label><input id="tslider" type="range" min="0" step="1" />
  <span id="tlabel" class="ro"></span>
  <label>fps</label><input id="fslider" type="range" min="1" max="24" step="1" value="8" />
  <span id="flabel" class="ro">8</span>
  <input id="search" type="text" placeholder="Search mutant, locus, or pathway..." />
  <span id="hits" style="font-size:13px;color:#889"></span>
  <button class="btn" id="fit">Fit</button>
  <div id="cbar"><span>-3</span><span class="ramp"></span><span>+3</span><span style="margin-left:6px">Z-score</span></div>
</div>
<div id="main">
  <div id="plot">
    <canvas id="cv"></canvas>
    <div id="tooltip"></div>
    <div id="hint">scroll=pan · ctrl+scroll=zoom · drag=pan · dblclick=fit</div>
  </div>
  <div id="sidebar" class="hidden">
    <div id="sb-head"></div><button id="sb-close">&times;</button>
    <div id="gallery"></div>
  </div>
</div>
<script>
const D = {dataJson};
const THUMBS = {thumbsJson};
const FUNC = {legendJson};
const ROWS = D.rows.length, NF = D.features.length;
const ROW_H = 15, CELL_W = 26, DEND_W = 210, STRIP_W = 14, LAB_W = 150;
const HEAT_X = DEND_W + STRIP_W + 6;

const cv = document.getElementById('cv'), ctx = cv.getContext('2d');
// Feature labels run at 45 degrees above AND below the heatmap, so the header and footer have to be
// as tall as the longest label's diagonal rise (width / sqrt2). A fixed height clipped the long
// names ("Nearest-Neighbor Distance Variability") and let them collide with the cells.
const FEAT_FONT = '10px sans-serif';
ctx.font = FEAT_FONT;
const FEAT_RISE = Math.ceil(Math.max(...D.features.map(f => ctx.measureText(f).width)) * Math.SQRT1_2);
const TOP_H = FEAT_RISE + 46;   // + room for the group brackets, which sit between labels and cells
const BOT_H = FEAT_RISE + 20;
const WORLD_W = HEAT_X + NF * CELL_W + 8 + Math.max(LAB_W, FEAT_RISE + 20);
const WORLD_H = TOP_H + ROWS * ROW_H + BOT_H;
const tip = document.getElementById('tooltip'), sb = document.getElementById('sidebar');
const gallery = document.getElementById('gallery'), sbHead = document.getElementById('sb-head');
let dpr = window.devicePixelRatio || 1, W, H, scale = 1, ox = 0, oy = 0, hover = -1, hits = null;

// ── timepoints: one byte per cell, decoded once ───────────────────────────────
const PAL = D.palette, KEYS = D.frameKeys, LABELS = D.frameLabels;
function decode(b64) {{
  const s = atob(b64), a = new Uint8Array(s.length);
  for (let i = 0; i < s.length; i++) a[i] = s.charCodeAt(i);
  return a;
}}
const CELLS = {{}};
for (const k of KEYS) CELLS[k] = decode(D.framesData[k]);
let curIdx = KEYS.length - 1;          // opens on the peak-biomass frame, matching the static figure
let playing = false, timer = null, fps = 8;
const tsl = document.getElementById('tslider'), fsl = document.getElementById('fslider');
const tlab = document.getElementById('tlabel'), flab = document.getElementById('flabel');
const playBtn = document.getElementById('play');
tsl.max = KEYS.length - 1; tsl.value = curIdx;
function setFrame(i) {{ curIdx = i; tsl.value = i; tlab.textContent = LABELS[i]; draw(); refreshGallery(); }}
function stop() {{ playing = false; if (timer) clearInterval(timer); timer = null; playBtn.innerHTML = '&#9654; Play'; }}
function start() {{
  playing = true; playBtn.innerHTML = '&#10073;&#10073; Pause';
  if (timer) clearInterval(timer);
  timer = setInterval(() => setFrame((curIdx + 1) % KEYS.length), 1000 / fps);
}}
playBtn.addEventListener('click', () => playing ? stop() : start());
tsl.addEventListener('input', e => {{ stop(); setFrame(+e.target.value); }});
fsl.addEventListener('input', e => {{ fps = +e.target.value; flab.textContent = fps; if (playing) start(); }});
document.addEventListener('keydown', e => {{
  if (e.target.matches('input')) return;
  if (e.key === 'ArrowRight') {{ stop(); setFrame(Math.min(curIdx + 1, KEYS.length - 1)); }}
  if (e.key === 'ArrowLeft') {{ stop(); setFrame(Math.max(curIdx - 1, 0)); }}
  if (e.key === ' ') {{ e.preventDefault(); playing ? stop() : start(); }}
}});

function fit() {{
  const s = Math.min(W / WORLD_W, H / WORLD_H);
  scale = s; ox = (W - WORLD_W * s) / 2; oy = 10;
}}
const wx = x => x * scale + ox, wy = y => y * scale + oy;
const sxToWorld = x => (x - ox) / scale, syToWorld = y => (y - oy) / scale;

function draw() {{
  if (!W) return;
  ctx.clearRect(0, 0, W, H);
  ctx.save(); ctx.translate(ox, oy); ctx.scale(scale, scale);

  // dendrogram: root at left, leaves at right
  ctx.strokeStyle = '#8fa'; ctx.lineWidth = 1 / scale * 1.2; ctx.globalAlpha = 0.75;
  for (const s of D.segments) {{
    ctx.beginPath();
    for (let k = 0; k < 4; k++) {{
      const X = DEND_W - s.d[k] * (DEND_W - 6);
      const Y = TOP_H + (s.i[k] + 0.5) * ROW_H;
      k === 0 ? ctx.moveTo(X, Y) : ctx.lineTo(X, Y);
    }}
    ctx.stroke();
  }}
  ctx.globalAlpha = 1;

  // feature-group brackets + rotated feature labels
  ctx.fillStyle = '#aab'; ctx.font = `${{Math.max(9, 12)}}px sans-serif`;
  for (const g of D.groups) {{
    const x0 = HEAT_X + g.a * CELL_W, x1 = HEAT_X + g.b * CELL_W;
    ctx.strokeStyle = '#889'; ctx.lineWidth = 1.2 / scale;
    ctx.beginPath(); ctx.moveTo(x0 + 2, TOP_H - 12); ctx.lineTo(x1 - 2, TOP_H - 12); ctx.stroke();
    ctx.textAlign = 'center'; ctx.fillText(g.name, (x0 + x1) / 2, TOP_H - 18);
  }}
  // Feature names on BOTH x axes. Above the heatmap they rise to the right (textAlign left); below it
  // they fall to the left (textAlign right) -- the static figure's rotation=45 / ha='right'. Getting
  // the alignment wrong is what previously sent the top labels down into the cells.
  ctx.fillStyle = '#99a'; ctx.font = FEAT_FONT;
  const botY = TOP_H + ROWS * ROW_H + 14;
  for (let j = 0; j < NF; j++) {{
    const cx = HEAT_X + j * CELL_W + CELL_W / 2;
    ctx.save(); ctx.translate(cx, TOP_H - 26); ctx.rotate(-Math.PI / 4);
    ctx.textAlign = 'left'; ctx.fillText(D.features[j], 0, 0); ctx.restore();
    ctx.save(); ctx.translate(cx, botY); ctx.rotate(-Math.PI / 4);
    ctx.textAlign = 'right'; ctx.fillText(D.features[j], 0, 0); ctx.restore();
  }}

  // annotation strip + heatmap + labels
  for (let i = 0; i < ROWS; i++) {{
    const y = TOP_H + i * ROW_H, r = D.rows[i];
    const dim = hits && !hits.has(i) ? 0.25 : 1;
    ctx.globalAlpha = dim;
    ctx.fillStyle = r.color; ctx.fillRect(DEND_W + 2, y, STRIP_W, ROW_H - 1);
    const buf = CELLS[KEYS[curIdx]], base = i * NF;
    for (let j = 0; j < NF; j++) {{ ctx.fillStyle = PAL[buf[base + j]]; ctx.fillRect(HEAT_X + j * CELL_W, y, CELL_W, ROW_H - 1); }}
    ctx.fillStyle = i === hover ? '#fff' : '#ccd'; ctx.textAlign = 'left';
    ctx.font = `${{i === hover ? 'bold ' : ''}}italic 11px sans-serif`;
    ctx.fillText(r.display, HEAT_X + NF * CELL_W + 8, y + ROW_H - 4);
    ctx.globalAlpha = 1;
    if (i === hover) {{
      ctx.strokeStyle = '#7ecfff'; ctx.lineWidth = 1.6 / scale;
      ctx.strokeRect(DEND_W + 2, y - 0.5, STRIP_W + 4 + NF * CELL_W, ROW_H);
    }}
  }}
  ctx.restore();
}}

function rowAt(sx, sy) {{
  const x = sxToWorld(sx), y = syToWorld(sy);
  if (x < DEND_W || x > HEAT_X + NF * CELL_W + 8 + LAB_W) return -1;
  const i = Math.floor((y - TOP_H) / ROW_H);
  return (i >= 0 && i < ROWS) ? i : -1;
}}
// In timelapse builds each replicate carries an index per frame ("fr"), so the gallery can follow
// the slider; otherwise "fr" is null and the peak-biomass image stays put.
let curMutant = -1;
function frameNumber() {{
  const k = KEYS[curIdx];
  return k === 'peak' ? null : parseInt(k, 10);      // null => each well's own peak frame
}}
function refreshGallery() {{
  if (curMutant < 0) return;
  const reps = D.rows[curMutant].reps, fn = frameNumber();
  gallery.querySelectorAll('.cell').forEach((cell, k) => {{
    const v = reps[k]; if (!v || !v.fr) return;
    const t = fn === null ? v.pf : fn;
    const ti = v.fr[String(t)];
    const img = cell.querySelector('img');
    if (ti !== undefined) img.src = THUMBS[ti];
    cell.querySelector('.cap').innerHTML = `${{v.w}}<br>t=${{t}}h${{t === v.pf ? ' (peak)' : ''}}`;
  }});
}}
function openMutant(i) {{
  curMutant = i;
  const r = D.rows[i], fn = frameNumber();
  sb.classList.remove('hidden');
  sbHead.innerHTML = `<div class="g">${{r.display}}</div><div class="l">${{r.locus}}</div>` +
    `<div class="a">${{r.anno}}</div>` +
    `<div class="m">${{r.n}} replicates · peak frame t=${{r.pf}}h</div>`;
  gallery.innerHTML = r.reps.length
    ? r.reps.map(v => {{
        const t = v.fr ? (fn === null ? v.pf : fn) : v.pf;
        const ti = v.fr ? (v.fr[String(t)] !== undefined ? v.fr[String(t)] : v.t) : v.t;
        return `<div class="cell"><img src="${{THUMBS[ti]}}" alt="biofilm image">` +
               `<div class="cap">${{v.w}}<br>t=${{t}}h${{t === v.pf ? ' (peak)' : ''}}</div></div>`;
      }}).join('')
    : '<div id="empty">No replicate images available</div>';
  resize();
}}
document.getElementById('sb-close').addEventListener('click', () => {{
  sb.classList.add('hidden'); gallery.innerHTML = ''; curMutant = -1; resize();
}});
cv.addEventListener('mousemove', e => {{
  const b = cv.getBoundingClientRect(), sx = e.clientX - b.left, sy = e.clientY - b.top;
  if (panning) {{ ox += e.clientX - px; oy += e.clientY - py; px = e.clientX; py = e.clientY; draw(); return; }}
  const i = rowAt(sx, sy);
  if (i !== hover) {{ hover = i; draw(); }}
  if (i >= 0) {{
    const r = D.rows[i];
    tip.innerHTML = `<div class="g">${{r.display}}</div><div>${{r.locus}} · ${{r.anno}}</div>` +
                    `<div style="color:#889">${{r.n}} replicates · click to view images</div>`;
    tip.style.display = 'block';
    let tx = sx + 14, ty = sy - 10;
    if (tx + tip.offsetWidth > W) tx = sx - tip.offsetWidth - 10;
    if (ty + tip.offsetHeight > H) ty = H - tip.offsetHeight - 4;
    tip.style.left = tx + 'px'; tip.style.top = Math.max(0, ty) + 'px'; cv.style.cursor = 'pointer';
  }} else {{ tip.style.display = 'none'; cv.style.cursor = 'default'; }}
}});
let panning = false, px, py, downX, downY;
cv.addEventListener('mousedown', e => {{ panning = true; px = downX = e.clientX; py = downY = e.clientY; }});
window.addEventListener('mouseup', e => {{
  if (panning && Math.abs(e.clientX - downX) + Math.abs(e.clientY - downY) < 4) {{
    const b = cv.getBoundingClientRect(); const i = rowAt(e.clientX - b.left, e.clientY - b.top);
    if (i >= 0) openMutant(i);
  }}
  panning = false;
}});
cv.addEventListener('wheel', e => {{
  e.preventDefault();
  const b = cv.getBoundingClientRect(), sx = e.clientX - b.left, sy = e.clientY - b.top;
  if (e.ctrlKey || e.metaKey) {{
    // ctrl/cmd + wheel: zoom about the cursor (also what a trackpad pinch sends)
    const wxp = sxToWorld(sx), wyp = syToWorld(sy), f = e.deltaY > 0 ? 1 / 1.12 : 1.12;
    scale *= f; ox = sx - wxp * scale; oy = sy - wyp * scale;
  }} else {{
    // plain wheel: scroll the plot. deltaX covers horizontal/shift-wheel and trackpad two-finger pans.
    ox -= e.deltaX; oy -= e.deltaY;
  }}
  draw();
}}, {{ passive: false }});
cv.addEventListener('dblclick', () => {{ fit(); draw(); }});
document.getElementById('fit').addEventListener('click', () => {{ fit(); draw(); }});
document.getElementById('search').addEventListener('input', e => {{
  const q = e.target.value.trim().toLowerCase();
  if (!q) {{ hits = null; document.getElementById('hits').textContent = ''; draw(); return; }}
  hits = new Set();
  D.rows.forEach((r, i) => {{
    if (r.mutant.toLowerCase().includes(q) || r.locus.toLowerCase().includes(q) ||
        r.anno.toLowerCase().includes(q)) hits.add(i);
  }});
  document.getElementById('hits').textContent = hits.size + ' hits';
  draw();
}});
function resize() {{
  const b = cv.parentElement.getBoundingClientRect(); W = b.width; H = b.height;
  cv.width = W * dpr; cv.height = H * dpr; cv.style.width = W + 'px'; cv.style.height = H + 'px';
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  if (!scale || scale === 1) fit();
  draw();
}}
window.addEventListener('resize', resize);
setFrame(curIdx); resize(); fit(); draw();
</script></body></html>
"""

outHtml.write_text(html)
print(f'Wrote {outHtml}  ({len(html.encode("utf-8")) / 1e6:.1f} MB)')
