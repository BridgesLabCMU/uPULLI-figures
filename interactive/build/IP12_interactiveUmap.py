#!/usr/bin/env python3
"""interactive build — Interactive Plots 1 and 2: UMAP explorers of the reimaging atlas.

Self-contained single-file HTML: dark-theme canvas, functional-annotation highlights, gene coloring,
search, zoom/pan, an n_neighbors/min_dist selector, and click-a-replicate -> its peak biofilm biomass
image. Every thumbnail is inlined as a base64 data URI, so the file works offline by double-clicking
with no server and no sibling asset folder -- what a journal needs for an interactive supplement.

  --track numerical   Interactive Plot 1: the quantitative uPULLI-I feature manifold
  --track embedding   Interactive Plot 2: the PCA-50 DINOv2 CLS embedding manifold (uPULLI-DL)

Thumbnails are taken pre-encoded from the deposited pack (112px grayscale JPEG, keyed by
<plate-well>_t<peakFrame>.jpg), so this needs no image-processing dependency and no access to the
image tree -- only the pinned core environment. Identical wells and peak frames underlie both
tracks, so the two plots differ solely in the coordinates.

Reads:  data/wells.csv, data/umapCoords_<track>.parquet, and the thumbnail pack (config.input)
Writes: figures/interactivePlot{1,2}_<track>Umap.html

Usage:  python build/IP12_interactiveUmap.py --track numerical
        python build/IP12_interactiveUmap.py --track embedding
"""
import sys
import json
import base64
import argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> interactive/ for figlib
import pandas as pd
from figlib import config, plotting

TRACKS = {
    'numerical': (1, 'umapCoords_numerical.parquet',
                  'Reimaging atlas - quantitative features (uPULLI-I)'),
    'embedding': (2, 'umapCoords_embedding.parquet',
                  'Reimaging atlas - DINOv2 CLS embeddings, PCA-50 (uPULLI-DL)'),
}

ap = argparse.ArgumentParser()
ap.add_argument('--track', choices=sorted(TRACKS), required=True)
ap.add_argument('--nn', type=int, default=10, help='n_neighbors selected on load')
ap.add_argument('--md', type=float, default=0.1, help='min_dist selected on load')
args = ap.parse_args()

plotNo, coordsFile, title = TRACKS[args.track]
outHtml = config.ensure(config.FIGURES / f'interactivePlot{plotNo}_{args.track}Umap.html')

WT_LABEL = 'WT'
highlightSets = {'WT': [], **plotting.HIGHLIGHT_SETS}
functionColors = {'WT': '#ffffff', **plotting.FUNCTION_COLORS}

# ── wells: labels, peak frame/biomass, and the pre-encoded thumbnail for each ──
packDir = Path(config.input('reimaging/thumbnails112/'))
wells = pd.read_csv(config.TABLES / 'wells.csv')
for c in ('mutant', 'geneLocus', 'function'):
    wells[c] = wells[c].fillna('')

thumbList, thumbIdx, missing = [], {}, 0
for _, r in wells.iterrows():
    jpg = packDir / str(r['file'])
    key = (str(r['plateId']), str(r['wellId']))
    if not jpg.exists():
        missing += 1
        continue
    thumbIdx[key] = len(thumbList)
    thumbList.append('data:image/jpeg;base64,' + base64.b64encode(jpg.read_bytes()).decode())
print(f'{len(thumbList)} thumbnails inlined from the pack ({missing} missing)')

meta = {(str(r['plateId']), str(r['wellId'])): r for _, r in wells.iterrows()}

# ── coordinates: the full (n_neighbors, min_dist) grid, so the selector works ──
emb = pd.read_parquet(config.TABLES / coordsFile)
for c in ('geneLocus', 'function', 'mutant'):
    emb[c] = emb[c].fillna('') if c in emb.columns else ''
paramCombos = (emb[['n_neighbors', 'min_dist']].drop_duplicates()
               .sort_values(['n_neighbors', 'min_dist']).values.tolist())

pointsByParam = {}
for nn, md in paramCombos:
    sub = emb[(emb['n_neighbors'] == nn) & (emb['min_dist'] == md)]
    pts = []
    for _, r in sub.iterrows():
        plate, well = str(r['plateId']), str(r['wellId'])
        m = meta.get((plate, well))
        pts.append({'x': round(float(r['umap1']), 4), 'y': round(float(r['umap2']), 4),
                    'gene': str(r['mutant']), 'locus': str(r['geneLocus']), 'func': str(r['function']),
                    'plate': plate, 'well': well, 'isWt': str(r['mutant']) == WT_LABEL,
                    'proj': False, 'pcolor': '',
                    'pf': int(m['peakFrame']) if m is not None else -1,
                    'bio': round(float(m['peakBiomass']), 5) if m is not None else None,
                    'ti': thumbIdx.get((plate, well), -1)})
    pointsByParam[f'{int(nn)}_{md}'] = pts

projPoints = []                      # plots 1 and 2 draw the atlas alone; kept so the viewer JS is unchanged
dataJson = json.dumps(pointsByParam)
paramJson = json.dumps(paramCombos)
hsJson = json.dumps(highlightSets)
fcJson = json.dumps(functionColors)
thumbsJson = json.dumps(thumbList)

html = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8">
<title>{title}</title>
<style>
  * {{ margin:0; padding:0; box-sizing:border-box; }}
  body {{ font-family:'Gillius ADF','Helvetica Neue',Arial,sans-serif; background:#1a1a2e; color:#eee; display:flex; flex-direction:column; height:100vh; }}
  #controls {{ display:flex; align-items:center; gap:20px; padding:12px 20px; background:#16213e; border-bottom:1px solid #334; flex-shrink:0; flex-wrap:wrap; }}
  #controls label {{ font-size:16px; color:#99a; }}
  #search-box {{ display:flex; align-items:center; gap:8px; margin-left:auto; }}
  #search-box input {{ background:#2a2a4a; border:1px solid #445; border-radius:5px; color:#eee; font-size:15px; padding:6px 12px; width:240px; outline:none; }}
  #search-box input::placeholder {{ color:#667; }}
  #search-box .count {{ font-size:14px; color:#889; }}
  #annotation-bar {{ display:flex; align-items:center; gap:10px; padding:8px 20px; background:#1a2040; border-bottom:1px solid #334; flex-shrink:0; flex-wrap:nowrap; overflow-x:auto; }}
  #annotation-bar label.section-label {{ font-size:15px; color:#778; margin-right:6px; white-space:nowrap; }}
  #anno-checks {{ display:flex; flex-direction:row; gap:8px; flex-wrap:nowrap; }}
  .anno-check {{ display:flex; align-items:center; gap:6px; cursor:pointer; font-size:14px; padding:4px 10px; border-radius:5px; border:1px solid transparent; user-select:none; white-space:nowrap; }}
  .anno-check:hover {{ background:rgba(255,255,255,0.05); }}
  .anno-check.active {{ border-color:currentColor; background:rgba(255,255,255,0.08); }}
  .anno-check input {{ display:none; }}
  .anno-swatch {{ width:12px; height:12px; border-radius:50%; flex-shrink:0; border:2px solid transparent; }}
  #toggle-all {{ background:#2a2a4a; border:1px solid #445; border-radius:5px; color:#aab; padding:4px 12px; font-size:13px; cursor:pointer; white-space:nowrap; }}
  #toggle-all:hover {{ background:#3a3a5a; }}
  .param-btn {{ background:#2a2a4a; border:1px solid #445; border-radius:5px; color:#ccc; padding:6px 14px; font-size:15px; cursor:pointer; margin:0 2px; }}
  .param-btn.active {{ background:#4a6fa5; border-color:#6a9fd5; color:#fff; }}
  #main {{ flex:1; display:flex; overflow:hidden; }}
  #plot-container {{ flex:1; position:relative; overflow:hidden; }}
  canvas {{ display:block; width:100%; height:100%; }}
  #tooltip {{ display:none; position:absolute; pointer-events:none; background:rgba(20,20,40,0.95); border:1px solid #556; border-radius:6px; padding:10px 14px; font-size:15px; line-height:1.5; max-width:400px; z-index:10; }}
  #tooltip .gene {{ font-weight:bold; font-size:18px; color:#7ecfff; }}
  #tooltip .locus {{ font-size:15px; color:#9bd; }}
  #tooltip .meta {{ color:#aab; font-size:14px; }}
  #tooltip .pathway {{ color:#dda; font-style:italic; font-size:14px; margin-top:3px; }}
  #legend {{ display:none; position:absolute; bottom:14px; left:14px; z-index:5; background:rgba(20,20,40,0.9); border:1px solid #556; border-radius:6px; padding:10px 14px; max-height:50vh; overflow-y:auto; font-size:14px; }}
  .legend-item {{ display:flex; align-items:center; gap:8px; padding:2px 0; }}
  .swatch {{ width:12px; height:12px; border-radius:50%; flex-shrink:0; }}
  .swatch-outline {{ width:12px; height:12px; border-radius:50%; flex-shrink:0; background:transparent; border:2px solid; }}
  #info {{ position:absolute; top:10px; right:14px; font-size:13px; color:#667; z-index:5; }}
  #sidebar {{ width:460px; background:#16213e; border-left:1px solid #334; display:flex; flex-direction:column; overflow:hidden; }}
  #sidebar.hidden {{ display:none; }}
  #thumb-header {{ padding:12px 16px; border-bottom:1px solid #334; font-size:14px; position:relative; }}
  #thumb-header .gene {{ font-weight:bold; font-size:16px; color:#7ecfff; }}
  #thumb-header .locus {{ font-size:13px; color:#9bd; }}
  #thumb-header .meta {{ color:#aab; font-size:13px; }}
  #thumb-header .pathway {{ color:#dda; font-style:italic; font-size:12px; }}
  #thumb-close {{ position:absolute; top:8px; right:12px; background:none; border:none; color:#889; font-size:20px; cursor:pointer; }}
  #thumb-close:hover {{ color:#eee; }}
  #thumb-wrap {{ flex:1; display:flex; align-items:center; justify-content:center; background:#000; }}
  #thumb-wrap img {{ max-width:100%; max-height:100%; image-rendering:auto; }}
  #no-thumb {{ color:#667; font-size:14px; text-align:center; padding:20px; }}
</style></head><body>
<div id="controls">
  <div><label>n_neighbors:</label><span id="nn-btns"></span></div>
  <div><label>min_dist:</label><span id="md-btns"></span></div>
  <div id="search-box"><input id="search" type="text" placeholder="Search gene, locus, function..." /><span class="count" id="search-count"></span></div>
</div>
<div id="annotation-bar">
  <label class="section-label">Highlights:</label>
  <button id="toggle-all">All on/off</button>
  <div id="anno-checks"></div>
</div>
<div id="main">
  <div id="plot-container">
    <canvas id="canvas"></canvas>
    <div id="tooltip"></div>
    <div id="legend"></div>
    <div id="info">scroll=pan · ctrl+scroll=zoom · drag=pan · shift-drag=box · dblclick=reset · click=peak-biomass thumbnail</div>
  </div>
  <div id="sidebar" class="hidden">
    <div id="thumb-header"></div>
    <button id="thumb-close">&times;</button>
    <div id="thumb-wrap"><div id="no-thumb">Click a point to view its peak-biomass frame</div></div>
  </div>
</div>
<script>
const allData = {dataJson};
const paramCombos = {paramJson};
const highlightSets = {hsJson};
const functionColors = {fcJson};
const THUMBS = {thumbsJson};
const nnValues = [...new Set(paramCombos.map(p => p[0]))].sort((a,b)=>a-b);
const mdValues = [...new Set(paramCombos.map(p => p[1]))].sort((a,b)=>a-b);
let currentNn = {args.nn}, currentMd = {args.md};
if (!nnValues.includes(currentNn)) currentNn = nnValues[0];
if (!mdValues.includes(currentMd)) currentMd = mdValues[0];
const locusToFunc = {{}};
for (const [f, loci] of Object.entries(highlightSets)) {{ if (f==='WT') continue; for (const l of loci) locusToFunc[l]=f; }}

const activeAnnotations = new Set();
const canvas = document.getElementById('canvas'), ctx = canvas.getContext('2d');
const tooltip = document.getElementById('tooltip');
const sidebar = document.getElementById('sidebar');
const thumbHeader = document.getElementById('thumb-header'), thumbWrap = document.getElementById('thumb-wrap');
const searchInput = document.getElementById('search'), searchCount = document.getElementById('search-count');
let dpr = window.devicePixelRatio || 1, W, H, vx, vy, vw, vh, searchTerm='', highlightSet=null, points=[];
const GENE_COLOR_LIMIT = 80;
function loadPoints() {{ points = allData[currentNn + '_' + currentMd] || []; resetView(); applySearch(); draw(); }}

function geneColor(i,total,a) {{ const hue=(i*360/total)%360; return `hsla(${{hue}},${{70+(i%3)*10}}%,${{55+(i%2)*10}}%,${{a}})`; }}
function hexToRgba(hex,a) {{ return `rgba(${{parseInt(hex.slice(1,3),16)}},${{parseInt(hex.slice(3,5),16)}},${{parseInt(hex.slice(5,7),16)}},${{a}})`; }}
function pointHl(p) {{ if (activeAnnotations.has('WT') && p.isWt) return 'WT'; const f=locusToFunc[p.locus]; return (f&&activeAnnotations.has(f))?f:null; }}
function resetView() {{
  if (!points.length) return;
  const xs=points.map(p=>p.x), ys=points.map(p=>p.y);
  const x0=Math.min(...xs), x1=Math.max(...xs), y0=Math.min(...ys), y1=Math.max(...ys);
  const dw=(x1-x0)||1, dh=(y1-y0)||1, pad=0.05;
  vx=x0-pad*dw; vy=y0-pad*dh; vw=dw*(1+2*pad); vh=dh*(1+2*pad);
}}
function w2s(wx,wy) {{ return [(wx-vx)/vw*W,(1-(wy-vy)/vh)*H]; }}
function s2w(sx,sy) {{ return [sx/W*vw+vx,(1-sy/H)*vh+vy]; }}

function draw() {{
  if (!W) return;
  ctx.clearRect(0,0,W,H);
  const bg = points.filter(p=>!p.proj), proj = points.filter(p=>p.proj);
  const r = Math.max(1.5, Math.min(4.5, 2200/Math.sqrt(bg.length||1)));
  const annoActive = activeAnnotations.size>0, searchActive = highlightSet!==null;
  const hasOverlay = annoActive || searchActive || proj.length>0;
  const visibleGenes = new Set(), visibleIdx = [];
  for (let i=0;i<points.length;i++) {{ const p=points[i]; if (p.proj) continue;
    const [sx,sy]=w2s(p.x,p.y); if (sx>=-5&&sx<=W+5&&sy>=-5&&sy<=H+5) {{ visibleIdx.push(i); visibleGenes.add(p.gene); }} }}
  const colorByGene = !hasOverlay && visibleGenes.size<=GENE_COLOR_LIMIT && visibleGenes.size>1;
  let geneList=[], geneToIdx={{}};
  if (colorByGene) {{ geneList=[...visibleGenes].sort(); geneList.forEach((g,i)=>geneToIdx[g]=i); }}
  for (const i of visibleIdx) {{
    const p=points[i]; if (pointHl(p) || (searchActive&&highlightSet.has(i))) continue;
    const [sx,sy]=w2s(p.x,p.y);
    ctx.beginPath(); ctx.arc(sx,sy,r,0,Math.PI*2);
    if (colorByGene) ctx.fillStyle=geneColor(geneToIdx[p.gene],geneList.length,0.75);
    else if (hasOverlay) ctx.fillStyle='rgba(150,160,200,0.14)';
    else if (p.isWt) ctx.fillStyle='rgba(255,140,0,0.8)';
    else ctx.fillStyle='rgba(100,180,255,0.5)';
    ctx.fill();
  }}
  if (annoActive) for (const i of visibleIdx) {{ const p=points[i], hl=pointHl(p); if (!hl) continue;
    const [sx,sy]=w2s(p.x,p.y); ctx.beginPath(); ctx.arc(sx,sy,r*1.3,0,Math.PI*2);
    ctx.strokeStyle=hexToRgba(functionColors[hl],0.95); ctx.lineWidth=2.5; ctx.stroke(); }}
  if (searchActive) for (const i of highlightSet) {{ const p=points[i]; if (p.proj) continue; const [sx,sy]=w2s(p.x,p.y);
    ctx.beginPath(); ctx.arc(sx,sy,r*1.5,0,Math.PI*2); ctx.fillStyle=p.isWt?'rgba(255,140,0,0.95)':'rgba(50,220,120,0.9)'; ctx.fill(); }}
  for (let i=0;i<points.length;i++) {{ const p=points[i]; if (!p.proj) continue;
    const [sx,sy]=w2s(p.x,p.y); if (sx<-8||sx>W+8||sy<-8||sy>H+8) continue;
    ctx.beginPath(); ctx.arc(sx,sy,r*1.7,0,Math.PI*2);
    ctx.fillStyle=(p.pcolor||'#ffffff'); ctx.fill(); ctx.lineWidth=1.6; ctx.strokeStyle='rgba(0,0,0,0.9)'; ctx.stroke(); }}
  const legend=document.getElementById('legend');
  if (annoActive) {{ let h=''; for (const f of activeAnnotations) h+=`<div class="legend-item"><span class="swatch-outline" style="border-color:${{functionColors[f]}}"></span>${{f}}</div>`;
    legend.innerHTML=h; legend.style.display='block'; }}
  else if (colorByGene) {{ let h=''; for (let i=0;i<geneList.length;i++) h+=`<div class="legend-item"><span class="swatch" style="background:${{geneColor(i,geneList.length,1)}}"></span>${{geneList[i]}}</div>`;
    legend.innerHTML=h; legend.style.display='block'; }}
  else legend.style.display='none';
}}
function nearest(sx,sy,maxD) {{ let best=-1, bd2=maxD*maxD;
  for (let i=0;i<points.length;i++) {{ const [px,py]=w2s(points[i].x,points[i].y); const d2=(px-sx)**2+(py-sy)**2;
    if (d2<bd2) {{ bd2=d2; best=i; }} }} return best; }}
function applySearch() {{
  if (!searchTerm) {{ highlightSet=null; searchCount.textContent=''; return; }}
  highlightSet=new Set();
  for (let i=0;i<points.length;i++) {{ const p=points[i];
    if (p.gene.toLowerCase().includes(searchTerm)||p.locus.toLowerCase().includes(searchTerm)||p.func.toLowerCase().includes(searchTerm)) highlightSet.add(i); }}
  searchCount.textContent=highlightSet.size+' hits';
}}
function headerHtml(p) {{
  let h=`<div class="gene">${{p.gene||'(unlabeled)'}}</div>`;
  if (p.locus) h+=`<div class="locus">${{p.locus}}</div>`;
  h+=`<div class="meta">Plate: ${{p.plate}} · Well: ${{p.well}}</div>`;
  h+=`<div class="meta">peak frame t=${{p.pf}}h${{p.bio!=null?` · biomass ${{p.bio}}`:''}}</div>`;
  const pw=locusToFunc[p.locus]||''; if (pw) h+=`<div class="pathway">${{pw}}</div>`;
  else if (p.func) h+=`<div class="pathway">${{p.func}}</div>`;
  return h;
}}
function showThumb(p) {{
  sidebar.classList.remove('hidden');
  thumbHeader.innerHTML=headerHtml(p);
  const uri = (p.ti>=0 && THUMBS[p.ti]) ? THUMBS[p.ti] : '';
  if (uri) thumbWrap.innerHTML=`<img src="${{uri}}" alt="peak biomass frame">`;
  else thumbWrap.innerHTML='<div id="no-thumb">No thumbnail available for this well</div>';
  resize();
}}
document.getElementById('thumb-close').addEventListener('click',()=>{{ sidebar.classList.add('hidden'); thumbWrap.innerHTML=''; resize(); }});
function buildAnnotations() {{
  const c=document.getElementById('anno-checks');
  for (const f of Object.keys(highlightSets)) {{
    const label=document.createElement('label'); label.className='anno-check'; label.style.color=functionColors[f];
    const cb=document.createElement('input'); cb.type='checkbox'; cb.dataset.func=f;
    cb.addEventListener('change',()=>{{ if (cb.checked) {{ activeAnnotations.add(f); label.classList.add('active'); }}
      else {{ activeAnnotations.delete(f); label.classList.remove('active'); }} draw(); }});
    const sw=document.createElement('span'); sw.className='anno-swatch'; sw.style.borderColor=functionColors[f];
    label.appendChild(cb); label.appendChild(sw); label.appendChild(document.createTextNode(f)); c.appendChild(label);
  }}
  document.getElementById('toggle-all').addEventListener('click',()=>{{
    const allOn=activeAnnotations.size===Object.keys(highlightSets).length;
    c.querySelectorAll('input').forEach(cb=>{{ cb.checked=!allOn; const f=cb.dataset.func;
      if (!allOn) {{ activeAnnotations.add(f); cb.parentElement.classList.add('active'); }}
      else {{ activeAnnotations.delete(f); cb.parentElement.classList.remove('active'); }} }}); draw();
  }});
}}
function buildButtons() {{
  const nnC=document.getElementById('nn-btns'), mdC=document.getElementById('md-btns');
  nnValues.forEach(v=>{{ const b=document.createElement('button'); b.className='param-btn'+(v===currentNn?' active':''); b.textContent=v;
    b.onclick=()=>{{ currentNn=v; updButtons(); loadPoints(); }}; nnC.appendChild(b); }});
  mdValues.forEach(v=>{{ const b=document.createElement('button'); b.className='param-btn'+(v===currentMd?' active':''); b.textContent=v;
    b.onclick=()=>{{ currentMd=v; updButtons(); loadPoints(); }}; mdC.appendChild(b); }});
}}
function updButtons() {{
  document.querySelectorAll('#nn-btns .param-btn').forEach((b,i)=>b.classList.toggle('active',nnValues[i]===currentNn));
  document.querySelectorAll('#md-btns .param-btn').forEach((b,i)=>b.classList.toggle('active',mdValues[i]===currentMd));
}}
let mode='none', dx, dy, bx0, by0, bx1, by1, downX, downY;
canvas.addEventListener('mousemove', e=>{{
  const rect=canvas.getBoundingClientRect(), sx=e.clientX-rect.left, sy=e.clientY-rect.top;
  if (mode==='box') {{ bx1=sx; by1=sy; draw(); ctx.strokeStyle='rgba(120,210,255,0.8)'; ctx.lineWidth=1.5; ctx.setLineDash([6,3]);
    ctx.strokeRect(Math.min(bx0,bx1),Math.min(by0,by1),Math.abs(bx1-bx0),Math.abs(by1-by0)); ctx.setLineDash([]); return; }}
  if (mode==='pan') {{ vx-=(e.clientX-dx)/W*vw; vy+=(e.clientY-dy)/H*vh; dx=e.clientX; dy=e.clientY; draw(); return; }}
  const idx=nearest(sx,sy,13);
  if (idx>=0) {{ const p=points[idx]; tooltip.innerHTML=headerHtml(p).concat('<div class="meta">click to view frame</div>');
    tooltip.style.display='block'; let tx=sx+14, ty=sy-10;
    if (tx+tooltip.offsetWidth>W) tx=sx-tooltip.offsetWidth-10; if (ty<0) ty=4; if (ty+tooltip.offsetHeight>H) ty=H-tooltip.offsetHeight-4;
    tooltip.style.left=tx+'px'; tooltip.style.top=ty+'px'; canvas.style.cursor='pointer';
  }} else {{ tooltip.style.display='none'; canvas.style.cursor=e.shiftKey?'crosshair':'default'; }}
}});
canvas.addEventListener('mousedown', e=>{{ if (e.button!==0) return; const rect=canvas.getBoundingClientRect();
  const sx=e.clientX-rect.left, sy=e.clientY-rect.top; downX=e.clientX; downY=e.clientY;
  if (e.shiftKey) {{ mode='box'; bx0=sx; by0=sy; bx1=sx; by1=sy; }} else {{ mode='pan'; dx=e.clientX; dy=e.clientY; }} }});
window.addEventListener('mouseup', e=>{{
  if (mode==='box') {{ const bw=Math.abs(bx1-bx0), bh=Math.abs(by1-by0);
    if (bw>5&&bh>5) {{ const [wx0,wy0]=s2w(Math.min(bx0,bx1),Math.max(by0,by1)); const [wx1,wy1]=s2w(Math.max(bx0,bx1),Math.min(by0,by1));
      vx=wx0; vy=wy0; vw=wx1-wx0; vh=wy1-wy0; }} mode='none'; draw(); return; }}
  if (mode==='pan') {{ const moved=Math.abs(e.clientX-downX)+Math.abs(e.clientY-downY);
    if (moved<4) {{ const rect=canvas.getBoundingClientRect(); const idx=nearest(e.clientX-rect.left,e.clientY-rect.top,13); if (idx>=0) showThumb(points[idx]); }}
    mode='none'; }}
}});
canvas.addEventListener('wheel', e=>{{ e.preventDefault(); const rect=canvas.getBoundingClientRect();
  const sx=e.clientX-rect.left, sy=e.clientY-rect.top;
  if (e.ctrlKey || e.metaKey) {{
    // ctrl/cmd + wheel: zoom about the cursor (also what a trackpad pinch sends)
    const [wx,wy]=s2w(sx,sy); const f=e.deltaY>0?1.15:1/1.15;
    vw*=f; vh*=f; vx=wx-sx/W*vw; vy=wy-(1-sy/H)*vh;
  }} else {{
    // plain wheel: scroll the view. World y is up, so scrolling down (deltaY>0) lowers vy, which
    // moves the content up on screen -- the same direction as scrolling a page.
    vy-=e.deltaY/H*vh; vx+=e.deltaX/W*vw;
  }}
  draw(); }}, {{passive:false}});
canvas.addEventListener('dblclick', ()=>{{ resetView(); draw(); }});
searchInput.addEventListener('input', e=>{{ searchTerm=e.target.value.trim().toLowerCase(); applySearch(); draw(); }});
function resize() {{ const rect=canvas.parentElement.getBoundingClientRect(); W=rect.width; H=rect.height;
  canvas.width=W*dpr; canvas.height=H*dpr; canvas.style.width=W+'px'; canvas.style.height=H+'px'; ctx.setTransform(dpr,0,0,dpr,0,0); draw(); }}
window.addEventListener('resize', resize);
buildButtons(); buildAnnotations(); loadPoints(); resize();
</script></body></html>
"""

outHtml.write_text(html)
mb = len(html.encode('utf-8')) / 1e6
selKey = f'{args.nn}_{args.md}'
sel = pointsByParam.get(selKey, next(iter(pointsByParam.values()), []))
print(f'Wrote {outHtml}  ({mb:.1f} MB)')
print(f'  combos: {list(pointsByParam)} | {len(thumbList)} thumbnails | '
      f'selected {selKey}: {len(sel)} replicates')
