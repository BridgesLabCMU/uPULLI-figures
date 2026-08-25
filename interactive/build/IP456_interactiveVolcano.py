#!/usr/bin/env python3
"""interactive build — Interactive Plots 4-6: RNA-seq volcano explorers for the clean deletions.

One self-contained single-file HTML per mutant (ΔbioD, ΔpdhE2, ΔmanA vs WT). Every quantified gene is
a point; hovering gives its name, locus, description, log2FC and adjusted p; the sidebar toggles
pathway highlights; the search box matches gene name, either locus form, or description text.

Port of the lab's `Interactive_Volcano_Explorer_light_V2.py` (RNAseq/Scripts on the lab share), with
the repo's requirements applied: reads only the BUNDLED tables in data/ (no absolute paths, no
network, nothing to download), and the six shared functional groups keep the loci and colors they
carry in every other figure, so a pathway reads the same here as in Fig 3 or Fig 4.

Significant genes (past BOTH |log2FC| > 2 and q < 0.05) are drawn in that mutant's own color -- orange
for bioD, green for pdhE2, blue for manA -- the same colors they carry in Figs 3-5 and Fig S8B-D.

Navigation matches Interactive Plots 1-3: scroll pans, ctrl/cmd + scroll zooms about the cursor,
drag pans, shift-drag box-zooms, double-click resets.

Reads:  data/rnaseq_volcano_{BioD,PdhE2,ManA}.csv   (bundled; produced by figS8/build/S8BCD_volcano.py)
Writes: figures/interactivePlot{4,5,6}_volcano<mutant>.html

Usage:
  python build/IP456_interactiveVolcano.py                 # all three
  python build/IP456_interactiveVolcano.py --mutant BioD
"""
import sys
import json
import argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> interactive/ for figlib
import numpy as np
import pandas as pd
from figlib import (config, VOLCANO_PANELS, VOLCANO_HIGHLIGHTS, VOLCANO_COLORS,
                    VOLCANO_SIG_COLORS, VOLCANO_SEARCH_COLOR)

# Dashed guides. |log2FC| > 2 AND q < 0.05 matches Fig S8B-D exactly, so the interactive and printed
# volcanoes call the same genes significant (the lab's original explorer used 1, which disagreed).
FC_LINE, Q_LINE = 2.0, 0.05

ap = argparse.ArgumentParser()
ap.add_argument('--mutant', default='all', choices=['all'] + list(VOLCANO_PANELS))
args = ap.parse_args()


def normalizeId(s):
    """Match ids across the three id styles in the tables: VC_RS00005 / VC0002; VC_0002 / gene name."""
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return ''
    return str(s).replace('_', '').replace(' ', '').upper()


normHighlights = {k: [normalizeId(v) for v in vals] for k, vals in VOLCANO_HIGHLIGHTS.items()}


def loadTable(mutant):
    df = pd.read_csv(config.TABLES / f'rnaseq_volcano_{mutant}.csv').dropna(subset=['qvalue'])
    df['y'] = -np.log10(df['qvalue'])
    return df


def rowIds(r):
    """Every id form a highlight set might be keyed by: new locus, gene name, each old locus."""
    oldLocus = '' if pd.isna(r['oldLocus']) else str(r['oldLocus'])
    gene = '' if pd.isna(r['gene']) else str(r['gene'])
    return {i for i in [normalizeId(r['locus']), normalizeId(gene)] + [normalizeId(o) for o in oldLocus.split(';')] if i}


# Shared axes across Plots 4-6: one view range over all three contrasts, so the panels are directly
# comparable and "reset" lands in the same place in each. Computed once, injected into every page.
_tables = {m: loadTable(m) for m in VOLCANO_PANELS}
_allX = pd.concat([t['logFC'] for t in _tables.values()])
_allY = pd.concat([t['y'] for t in _tables.values()])
SHARED_VIEW = {'x0': float(np.floor(_allX.min())), 'x1': float(np.ceil(_allX.max())),
               'y0': 0.0, 'y1': float(np.ceil(_allY.max()))}
print(f"shared axes: log2FC [{SHARED_VIEW['x0']:g}, {SHARED_VIEW['x1']:g}]  "
      f"-log10 q [0, {SHARED_VIEW['y1']:g}]")

# ── page template. Plain string + placeholders (not an f-string) so the JS braces stay untouched. ──
HTML = """<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8">
<title>TITLE_TEXT_PLACEHOLDER</title>
<style>
  body { font-family: 'Gillius ADF', -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
         background:#f6f8fa; color:#24292f; margin:0; display:flex; flex-direction:column; height:100vh; overflow:hidden; }
  #header { background:#fff; padding:10px 20px; display:flex; align-items:center; border-bottom:1px solid #d0d7de; flex-shrink:0; }
  #title { font-size:18px; font-weight:600; margin:0; }
  #search { background:#fff; border:1px solid #d0d7de; color:#24292f; padding:6px 12px; border-radius:6px; width:250px; outline:none; }
  #search:focus { border-color:#0969da; box-shadow:0 0 0 3px rgba(9,105,218,0.3); }
  #sidebar { background:#fff; width:250px; padding:15px; border-right:1px solid #d0d7de; overflow-y:auto; flex-shrink:0; }
  #main { flex:1; position:relative; background:#fff; }
  canvas { width:100%; height:100%; cursor:crosshair; }
  #tooltip { position:absolute; background:rgba(255,255,255,0.98); border:1px solid #d0d7de; padding:12px; display:none;
             pointer-events:none; border-radius:8px; z-index:1000; max-width:320px; box-shadow:0 8px 24px rgba(140,149,159,0.2); }
  .t-gene { color:#0969da; font-size:16px; font-weight:bold; margin-bottom:2px; }
  .t-meta { color:#57606a; font-size:12px; margin-bottom:8px; border-bottom:1px solid #d0d7de; padding-bottom:4px; }
  .t-desc { font-size:13px; line-height:1.4; }
  .pathway-item { display:flex; align-items:center; gap:8px; margin-bottom:8px; font-size:13px; cursor:pointer; font-weight:500; }
  .pathway-item .n { margin-left:auto; color:#57606a; font-weight:400; font-size:11px; }
  .mini { background:#f6f8fa; border:1px solid #d0d7de; border-radius:5px; color:#57606a; font-size:11px; padding:2px 8px; cursor:pointer; }
  .mini:hover { background:#eaeef2; }
  #tooltip.pinned { pointer-events:auto; border-color:#0969da; box-shadow:0 8px 24px rgba(9,105,218,0.25); }
  #tooltip .pin-note { font-size:10px; color:#57606a; margin-top:6px; }
</style></head><body>

<div id="header">
  <div id="title">PLOT_TITLE_PLACEHOLDER</div>
  <div style="margin-left:auto; display:flex; align-items:center; gap:10px;">
    <input type="text" id="search" placeholder="Search genes, locus, description...">
    <span id="hit-count" style="font-size:12px; color:#57606a; width:60px;"></span>
  </div>
</div>

<div style="display:flex; flex:1; min-height:0;">
  <div id="sidebar">
    <div style="display:flex; align-items:center; margin-bottom:10px;">
      <span style="font-weight:bold; color:#57606a; font-size:11px; text-transform:uppercase; letter-spacing:0.5px;">Highlight pathways</span>
      <span style="margin-left:auto; display:flex; gap:6px;">
        <button id="sel-all" class="mini">all</button><button id="sel-none" class="mini">none</button>
      </span>
    </div>
    <div id="pathway-list"></div>
    <div style="margin-top:24px; font-size:11px; color:#57606a; line-height:1.7; border-top:1px solid #d0d7de; padding-top:14px;">
      <b>Key:</b><br>
      <span style="color:SIG_COLOR_PLACEHOLDER; font-size:15px; -webkit-text-stroke:0.5px rgba(0,0,0,0.55);">&#9679;</span> |log<sub>2</sub>FC| &gt; FC_LINE_PLACEHOLDER and adjusted <i>p</i> &lt; Q_LINE_PLACEHOLDER<br>
      <span style="color:#afb8c1; font-size:15px;">&#9679;</span> all other genes<br>
      <span style="color:SEARCH_COLOR_PLACEHOLDER; font-size:15px;">&#9679;</span> search hit<br>
      <span style="color:#57606a;">Dashed guides mark both cutoffs.</span>
    </div>
    <div style="margin-top:18px; font-size:11px; color:#57606a; line-height:1.6; border-top:1px solid #d0d7de; padding-top:14px;">
      <b>Controls:</b><br>
      &bull; Click a gene: pin its card<br>
      &bull; Scroll: pan<br>
      &bull; Ctrl/&#8984; + scroll: zoom at cursor<br>
      &bull; Click &amp; drag: pan<br>
      &bull; Shift + drag: box zoom<br>
      &bull; Double-click: reset view
    </div>
  </div>
  <div id="main"><canvas id="canvas"></canvas><div id="tooltip"></div></div>
</div>

<script>
const points = DATA_JSON;
const highlightSets = HIGHLIGHT_JSON;
const colors = COLOR_JSON;
const FC_LINE = FC_LINE_PLACEHOLDER, Q_LINE = Q_LINE_PLACEHOLDER;
const SIG_COLOR = 'SIG_COLOR_PLACEHOLDER';        // this mutant's color, as in Figs 3-5 and S8
const SEARCH_COLOR = 'SEARCH_COLOR_PLACEHOLDER';
const SHARED_VIEW = SHARED_VIEW_JSON;   // same axes in Plots 4-6, so the three are comparable
const PATH_COUNTS = COUNTS_JSON;        // pathway -> [genes matched, of which significant]
let pinned = -1;

let activePathways = new Set();
const canvas = document.getElementById('canvas');
const ctx = canvas.getContext('2d');
const tooltip = document.getElementById('tooltip');
let W, H, vx, vy, vw, vh;
let dpr = window.devicePixelRatio || 1;
let mode = 'none', dragX, dragY, boxX, boxY, downX = 0, downY = 0;
let searchTerm = '', searchHits = null;
const margin = { top: 25, right: 30, bottom: 55, left: 65 };

const idToPathway = {};
for (const [name, ids] of Object.entries(highlightSets)) { ids.forEach(id => idToPathway[id] = name); }

function resetView() {
  const padX = (SHARED_VIEW.x1 - SHARED_VIEW.x0) * 0.04, padY = (SHARED_VIEW.y1 - SHARED_VIEW.y0) * 0.06;
  vx = SHARED_VIEW.x0 - padX; vw = (SHARED_VIEW.x1 - SHARED_VIEW.x0) + padX * 2;
  vy = SHARED_VIEW.y0 - padY; vh = (SHARED_VIEW.y1 - SHARED_VIEW.y0) + padY * 2;
}
function plotW() { return W - margin.left - margin.right; }
function plotH() { return H - margin.top - margin.bottom; }
function worldToScreen(wx, wy) {
  return [margin.left + (wx - vx) / vw * plotW(), margin.top + (1 - (wy - vy) / vh) * plotH()];
}
function screenToWorld(sx, sy) {
  return [(sx - margin.left) / plotW() * vw + vx, (1 - (sy - margin.top) / plotH()) * vh + vy];
}
function getTicks(min, max, targetCount) {
  const span = max - min; if (span <= 0) return [];
  const rawStep = span / targetCount, mag = Math.pow(10, Math.floor(Math.log10(rawStep)));
  const residual = rawStep / mag;
  const step = residual < 1.5 ? mag : residual < 3 ? 2 * mag : residual < 7 ? 5 * mag : 10 * mag;
  const ticks = [];
  for (let v = Math.ceil(min / step) * step; v <= max; v += step) ticks.push(Number(v.toFixed(4)));
  return ticks;
}

function draw() {
  if (!W) return;
  ctx.clearRect(0, 0, W, H);
  const pw = plotW(), ph = plotH();
  ctx.save(); ctx.beginPath(); ctx.rect(margin.left, margin.top, pw, ph); ctx.clip();

  ctx.strokeStyle = '#d0d7de'; ctx.setLineDash([5, 5]);
  const [lx1] = worldToScreen(-FC_LINE, 0), [lx2] = worldToScreen(FC_LINE, 0);
  const [, ly] = worldToScreen(0, -Math.log10(Q_LINE));
  ctx.beginPath(); ctx.moveTo(margin.left, ly); ctx.lineTo(margin.left + pw, ly); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(lx1, margin.top); ctx.lineTo(lx1, margin.top + ph); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(lx2, margin.top); ctx.lineTo(lx2, margin.top + ph); ctx.stroke();
  ctx.setLineDash([]);

  // Three passes so a highlighted gene is never hidden under the 3,500-point cloud: plain points
  // first, then pathway highlights, then search hits on top. (Drawing in index order would let any
  // later gene paint over a highlighted one.)
  const r = vw > 8 ? 2 : 3;
  const plain = [], highlighted = [], hits = [];
  for (let i = 0; i < points.length; i++) {
    if (searchHits && searchHits.has(i)) { hits.push(i); continue; }
    let activePath = null;
    for (const id of points[i].ids) { if (idToPathway[id] && activePathways.has(idToPathway[id])) { activePath = idToPathway[id]; break; } }
    if (activePath) highlighted.push([i, activePath]); else plain.push(i);
  }
  for (const i of plain) {
    const p = points[i], [sx, sy] = worldToScreen(p.x, p.y);
    const isSig = Math.abs(p.x) > FC_LINE && p.y > -Math.log10(Q_LINE);
    ctx.beginPath(); ctx.arc(sx, sy, r, 0, 6.3);
    ctx.fillStyle = isSig ? SIG_COLOR : '#afb8c1';
    ctx.fill();
    if (isSig) {  // thin dark edge, as the printed panels use -- the bright colors need it on white
      ctx.strokeStyle = 'rgba(0,0,0,0.55)'; ctx.lineWidth = 0.8; ctx.stroke();
    }
  }
  for (const [i, path] of highlighted) {
    const p = points[i], [sx, sy] = worldToScreen(p.x, p.y);
    ctx.strokeStyle = colors[path]; ctx.lineWidth = 2.5;
    ctx.beginPath(); ctx.arc(sx, sy, r + 2, 0, 6.3); ctx.stroke();
    ctx.fillStyle = colors[path]; ctx.beginPath(); ctx.arc(sx, sy, r, 0, 6.3); ctx.fill();
  }
  for (const i of hits) {
    const p = points[i], [sx, sy] = worldToScreen(p.x, p.y);
    ctx.fillStyle = SEARCH_COLOR; ctx.beginPath(); ctx.arc(sx, sy, r + 3, 0, 6.3); ctx.fill();
  }
  if (mode === 'box') {
    ctx.strokeStyle = '#0969da'; ctx.lineWidth = 1; ctx.strokeRect(dragX, dragY, boxX - dragX, boxY - dragY);
    ctx.fillStyle = 'rgba(9,105,218,0.15)'; ctx.fillRect(dragX, dragY, boxX - dragX, boxY - dragY);
  }
  ctx.restore();

  ctx.fillStyle = '#57606a'; ctx.strokeStyle = '#d0d7de'; ctx.lineWidth = 1;
  ctx.font = '11px -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif';
  ctx.textAlign = 'center'; ctx.textBaseline = 'top';
  getTicks(vx, vx + vw, 8).forEach(val => {
    const [sx] = worldToScreen(val, 0);
    if (sx >= margin.left && sx <= margin.left + pw) {
      ctx.beginPath(); ctx.moveTo(sx, margin.top + ph); ctx.lineTo(sx, margin.top + ph + 5); ctx.stroke();
      ctx.fillText(val.toString(), sx, margin.top + ph + 8);
    }
  });
  ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
  getTicks(vy, vy + vh, 6).forEach(val => {
    if (val < 0) return;
    const [, sy] = worldToScreen(0, val);
    if (sy >= margin.top && sy <= margin.top + ph) {
      ctx.beginPath(); ctx.moveTo(margin.left - 5, sy); ctx.lineTo(margin.left, sy); ctx.stroke();
      ctx.fillText(val.toString(), margin.left - 8, sy);
    }
  });
  ctx.fillStyle = '#24292f';
  ctx.font = '500 13px -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif';
  ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
  ctx.fillText('log\u2082 fold change vs WT', margin.left + pw / 2, H - 15);
  ctx.save(); ctx.translate(18, margin.top + ph / 2); ctx.rotate(-Math.PI / 2);
  ctx.fillText('\u2212log\u2081\u2080 adjusted p', 0, 0); ctx.restore();
  if (pinned >= 0) placePinnedCard();
}

function showCard(i, sx, sy, isPinned) {
  const p = points[i];
  tooltip.style.display = 'block';
  tooltip.classList.toggle('pinned', !!isPinned);
  tooltip.innerHTML = '<div class="t-gene">' + p.gene + '</div>' +
    '<div class="t-meta">' + p.locus + ' | ' + p.old_locus + '</div>' +
    '<div class="t-desc">' + p.desc + '</div>' +
    '<div style="font-size:11px; margin-top:8px; color:#57606a; font-weight:600;">log\u2082FC: ' +
    p.x.toFixed(2) + ' | adj. p: ' + p.padj.toExponential(2) + '</div>' +
    (isPinned ? '<div class="pin-note">pinned \u2014 click empty space to dismiss</div>' : '');
  // clamp inside the plot so a gene near an edge is not cut off
  const mw = document.getElementById('main').clientWidth, mh = document.getElementById('main').clientHeight;
  let tx = sx + 15, ty = sy - 15;
  if (tx + tooltip.offsetWidth > mw) tx = sx - tooltip.offsetWidth - 15;
  if (ty + tooltip.offsetHeight > mh) ty = mh - tooltip.offsetHeight - 6;
  tooltip.style.left = Math.max(4, tx) + 'px';
  tooltip.style.top = Math.max(4, ty) + 'px';
}
function placePinnedCard() {
  if (pinned < 0) { tooltip.style.display = 'none'; return; }
  const [px, py] = worldToScreen(points[pinned].x, points[pinned].y);
  showCard(pinned, px, py, true);
}
function hitTest(sx, sy) {
  let found = -1, best = 12, bestTier = -1;
  for (let i = 0; i < points.length; i++) {
    const [px, py] = worldToScreen(points[i].x, points[i].y);
    const d = Math.sqrt((px - sx) ** 2 + (py - sy) ** 2);
    if (d >= 12) continue;
    let tier = 0;
    if (searchHits && searchHits.has(i)) tier = 2;
    else { for (const id of points[i].ids) { if (idToPathway[id] && activePathways.has(idToPathway[id])) { tier = 1; break; } } }
    if (tier > bestTier || (tier === bestTier && d < best)) { bestTier = tier; best = d; found = i; }
  }
  return found;
}

canvas.addEventListener('mousedown', e => {
  const rect = canvas.getBoundingClientRect();
  dragX = e.clientX - rect.left; dragY = e.clientY - rect.top;
  downX = dragX; downY = dragY;
  mode = e.shiftKey ? 'box' : 'pan';
});
window.addEventListener('mousemove', e => {
  const rect = canvas.getBoundingClientRect();
  const sx = e.clientX - rect.left, sy = e.clientY - rect.top;
  if (mode === 'pan') {
    vx -= (sx - dragX) / plotW() * vw; vy += (sy - dragY) / plotH() * vh;
    dragX = sx; dragY = sy; draw();
  } else if (mode === 'box') {
    boxX = sx; boxY = sy; draw();
  } else {
    // Match the draw order: a search hit wins over a pathway highlight, which wins over a plain
    // point, so the gene you can see on top is the one the tooltip reports.
    let found = -1, best = 12, bestTier = -1;
    for (let i = 0; i < points.length; i++) {
      const [px, py] = worldToScreen(points[i].x, points[i].y);
      const d = Math.sqrt((px - sx) ** 2 + (py - sy) ** 2);
      if (d >= 12) continue;
      let tier = 0;
      if (searchHits && searchHits.has(i)) tier = 2;
      else { for (const id of points[i].ids) { if (idToPathway[id] && activePathways.has(idToPathway[id])) { tier = 1; break; } } }
      if (tier > bestTier || (tier === bestTier && d < best)) { bestTier = tier; best = d; found = i; }
    }
    if (found >= 0) showCard(found, sx, sy, false);
    else if (pinned < 0) tooltip.style.display = 'none';
    else placePinnedCard();
  }
});
window.addEventListener('mouseup', e => {
  if (mode === 'pan') {
    const rect = canvas.getBoundingClientRect();
    const sx = e.clientX - rect.left, sy = e.clientY - rect.top;
    if (Math.abs(sx - downX) + Math.abs(sy - downY) < 4) {   // a click, not a drag
      const i = hitTest(sx, sy);
      pinned = (i >= 0 && i !== pinned) ? i : -1;
      mode = 'none';
      draw(); placePinnedCard();
      return;
    }
  }
  if (mode === 'box') {
    const [wx0, wy0] = screenToWorld(dragX, dragY), [wx1, wy1] = screenToWorld(boxX, boxY);
    vx = Math.min(wx0, wx1); vw = Math.abs(wx1 - wx0);
    vy = Math.min(wy0, wy1); vh = Math.abs(wy1 - wy0);
  }
  mode = 'none'; draw();
});
canvas.addEventListener('wheel', e => {
  e.preventDefault();
  const rect = canvas.getBoundingClientRect();
  const sx = e.clientX - rect.left, sy = e.clientY - rect.top;
  if (e.ctrlKey || e.metaKey) {
    // ctrl/cmd + wheel: zoom about the cursor (also what a trackpad pinch sends)
    const [wx, wy] = screenToWorld(sx, sy);
    const factor = e.deltaY > 0 ? 1.2 : 0.8;
    vw *= factor; vh *= factor;
    vx = wx - ((sx - margin.left) / plotW()) * vw;
    vy = wy - (1 - (sy - margin.top) / plotH()) * vh;
  } else {
    // plain wheel: pan. World y is up, so scrolling down lowers vy and the points move up on screen.
    vy -= e.deltaY / plotH() * vh; vx += e.deltaX / plotW() * vw;
  }
  draw();
}, { passive: false });
canvas.addEventListener('dblclick', () => { resetView(); draw(); placePinnedCard(); });

document.getElementById('search').addEventListener('input', e => {
  searchTerm = e.target.value.toLowerCase();
  if (!searchTerm) { searchHits = null; }
  else {
    searchHits = new Set();
    points.forEach((p, i) => {
      if (p.gene.toLowerCase().includes(searchTerm) || p.locus.toLowerCase().includes(searchTerm) ||
          p.old_locus.toLowerCase().includes(searchTerm) || p.desc.toLowerCase().includes(searchTerm)) searchHits.add(i);
    });
  }
  document.getElementById('hit-count').innerText = searchTerm ? searchHits.size + ' hits' : '';
  draw();
});

const list = document.getElementById('pathway-list');
const boxes = [];
for (const [name, color] of Object.entries(colors)) {
  const div = document.createElement('div');
  div.className = 'pathway-item';
  const c = PATH_COUNTS[name] || [0, 0];
  div.innerHTML = '<input type="checkbox"> <span style="color:' + color + '; font-size:16px;">&#9632;</span> ' +
                  name + '<span class="n" title="genes in this set / of those, past both cutoffs">' +
                  c[0] + ' \u00b7 ' + c[1] + ' sig</span>';
  const box = div.querySelector('input');
  box.onchange = ev => {
    if (ev.target.checked) activePathways.add(name); else activePathways.delete(name);
    draw();
  };
  boxes.push([box, name]);
  list.appendChild(div);
}
document.getElementById('sel-all').onclick = () => {
  boxes.forEach(([b, n]) => { b.checked = true; activePathways.add(n); }); draw();
};
document.getElementById('sel-none').onclick = () => {
  boxes.forEach(([b, n]) => { b.checked = false; activePathways.delete(n); }); draw();
};

function resize() {
  W = canvas.clientWidth; H = canvas.clientHeight;
  canvas.width = W * dpr; canvas.height = H * dpr;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  draw();
}
window.addEventListener('resize', resize);
resetView(); resize();
</script></body></html>
"""


def build(mutant):
    plotNo, titleHtml = VOLCANO_PANELS[mutant]
    df = _tables[mutant]
    sig = (df['logFC'].abs() > FC_LINE) & (df['qvalue'] < Q_LINE)

    pts, idSets = [], []
    for _, r in df.iterrows():
        oldLocus = '' if pd.isna(r['oldLocus']) else str(r['oldLocus'])
        gene = '' if pd.isna(r['gene']) else str(r['gene'])
        ids = rowIds(r)
        idSets.append(ids)
        pts.append({'x': round(float(r['logFC']), 4), 'y': round(float(r['y']), 4),
                    'gene': gene or str(r['locus']), 'locus': str(r['locus']), 'old_locus': oldLocus,
                    'desc': '' if pd.isna(r['description']) else str(r['description']),
                    'padj': float(r['qvalue']), 'ids': sorted(ids)})

    # per-pathway counts for the sidebar: genes present, and how many pass both cutoffs
    sigFlags = sig.tolist()
    counts = {}
    for name, members in normHighlights.items():
        want = set(members)
        idx = [i for i, s_ in enumerate(idSets) if s_ & want]
        counts[name] = [len(idx), int(sum(1 for i in idx if sigFlags[i]))]

    page = (HTML
            .replace('DATA_JSON', json.dumps(pts))
            .replace('HIGHLIGHT_JSON', json.dumps(normHighlights))
            .replace('COLOR_JSON', json.dumps(VOLCANO_COLORS))
            .replace('SHARED_VIEW_JSON', json.dumps(SHARED_VIEW))
            .replace('COUNTS_JSON', json.dumps(counts))
            .replace('PLOT_TITLE_PLACEHOLDER', titleHtml)
            .replace('TITLE_TEXT_PLACEHOLDER', f'Volcano explorer - {mutant} vs WT')
            .replace('SIG_COLOR_PLACEHOLDER', VOLCANO_SIG_COLORS[mutant])
            .replace('SEARCH_COLOR_PLACEHOLDER', VOLCANO_SEARCH_COLOR)
            .replace('FC_LINE_PLACEHOLDER', f'{FC_LINE:g}')
            .replace('Q_LINE_PLACEHOLDER', f'{Q_LINE:g}'))

    out = config.ensure(config.FIGURES / f'interactivePlot{plotNo}_volcano{mutant}.html')
    out.write_text(page, encoding='utf-8')
    print(f'Plot {plotNo} {mutant:6s}: {len(pts)} genes ({int(sig.sum())} past both guides) -> {out.name}'
          f'  ({out.stat().st_size / 1e6:.1f} MB)')


for m in (list(VOLCANO_PANELS) if args.mutant == 'all' else [args.mutant]):
    build(m)
