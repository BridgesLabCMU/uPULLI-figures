#!/usr/bin/env python3
"""interactive build — Interactive Plots 1 and 2: UMAP explorers of the reimaging atlas.

Self-contained single-file HTML: dark-theme canvas, functional-annotation highlights, gene coloring,
search, zoom/pan, an n_neighbors/min_dist selector, and click-a-replicate -> its peak biofilm biomass
image. Every thumbnail is inlined as a base64 data URI, so the file works offline by double-clicking
with no server and no sibling asset folder -- what a journal needs for an interactive supplement.

  --track numerical   Interactive Plot 1: the quantitative uPULLI-I feature manifold
  --track embedding   Interactive Plot 2: the PCA-50 DINOv2 CLS embedding manifold (uPULLI-DL)

  --timelapse         Plots 1T/2T: the same explorer, but a clicked replicate PLAYS its whole
                      timecourse instead of showing one still (slider, play/pause, adjustable fps).

Thumbnails are taken pre-encoded from the deposited packs, so this needs no image-processing
dependency and no access to the image tree -- only the pinned core environment. Identical wells and
peak frames underlie both tracks, so the two plots differ solely in the coordinates.

The timelapse mode's cost is coverage, not detail: every frame of every well is ~330 MB of JPEG, so
the deposited timelapse pack carries all 31 frames for a few representative replicates per mutant
(471 of the 3,669 manifold wells, spanning 157 of 158 mutants) instead of one frame for all of them.
Those wells are drawn with a ring and play; the rest keep their peak-biomass still from the
112 px pack, so no point loses the image it had in Plots 1/2. `--stills none` drops those extra
stills, halving the file at the cost of leaving non-representative points imageless.

Reads:  data/wells.csv, data/umapCoords_<track>.parquet, and the thumbnail pack(s) (config.input)
Writes: figures/interactivePlot{1,2}[T]_<track>Umap[Timelapse].html

Usage:  python build/IP12_interactiveUmap.py --track numerical
        python build/IP12_interactiveUmap.py --track embedding
        python build/IP12_interactiveUmap.py --track numerical --timelapse
        python build/IP12_interactiveUmap.py --track embedding --timelapse
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
ap.add_argument('--timelapse', action='store_true',
                help='play a clicked replicate through the timecourse (plots 1T/2T)')
ap.add_argument('--stills', choices=['all', 'none'], default='all',
                help='timelapse mode: also inline peak stills for wells outside the timelapse pack')
# ── lab-only: reference the image trees instead of inlining them (NOT for deposit) ──
ap.add_argument('--nas', choices=['mp4', 'frames'], default=None,
                help='LAB BUILD: link images from the mounted trees by relative path instead of '
                     'inlining them. mp4 = each well\'s processed.mp4; frames = the all-wells frame '
                     'pack. ~1 MB HTML covering every replicate, but it only opens where the '
                     'mount is present, so it is not depositable.')
ap.add_argument('--video-base', default='../../../../Sehna/multiphenotype-data-06-25-26/reimagingData',
                help='--nas mp4: first candidate path to the reimaging image tree, RELATIVE to the '
                     'built HTML as it sits in the mirror. Only a starting guess -- the page probes '
                     'this and then searches upward for the tree (see --asset-tail), so the file '
                     'still works if it is moved.')
ap.add_argument('--frame-base', default='../../../data/v2/reimaging/framesPack',
                help='--nas frames: first candidate path to the frame pack, relative to the built '
                     'HTML. Same story: a starting guess, not a requirement.')
ap.add_argument('--asset-tail', default='',
                help='extra directory tails to search for when hunting upward from the page, comma '
                     'separated (e.g. "phenotyper/Sehna/..."). Added ahead of the built-in tails.')
ap.add_argument('--search-depth', type=int, default=8,
                help='how many parent directories up from the page to search for the assets')
ap.add_argument('--mag', default='02', help='--nas mp4: magnification suffix on the image filenames')
ap.add_argument('--nas-frames-count', type=int, default=31,
                help='--nas frames: frames per well assumed present in the pack')
args = ap.parse_args()

plotNo, coordsFile, title = TRACKS[args.track]
if args.nas:
    # the mode is in the filename: the two lab builds differ in what they play and where it comes
    # from, so they must not overwrite each other
    title += f' - linked images ({args.nas})'
    outHtml = config.ensure(
        config.FIGURES / f'interactivePlot{plotNo}N_{args.track}UmapLinked_{args.nas}.html')
elif args.timelapse:
    title += ' - timelapse'
    outHtml = config.ensure(config.FIGURES / f'interactivePlot{plotNo}T_{args.track}UmapTimelapse.html')
else:
    outHtml = config.ensure(config.FIGURES / f'interactivePlot{plotNo}_{args.track}Umap.html')

WT_LABEL = 'WT'
highlightSets = {'WT': [], **plotting.HIGHLIGHT_SETS}
functionColors = {'WT': '#ffffff', **plotting.FUNCTION_COLORS}


def inline(jpg):
    thumbList.append('data:image/jpeg;base64,' + base64.b64encode(jpg.read_bytes()).decode())
    return len(thumbList) - 1


# ── wells: labels, peak frame/biomass, and the pre-encoded thumbnail(s) for each ──
wells = pd.read_csv(config.TABLES / 'wells.csv')
for c in ('mutant', 'geneLocus', 'function'):
    wells[c] = wells[c].fillna('')

thumbList, thumbIdx, frameIdx, missing = [], {}, {}, 0

# --nas: nothing is inlined at all. The page builds each well's URL from its plate/well at click
# time, so neither the images nor even a per-point URL list touch the payload.
if args.nas:
    print(f'--nas {args.nas}: no images inlined; the page links '
          f'{args.video_base if args.nas == "mp4" else args.frame_base}')

# timelapse mode first: wells in that pack get every frame, and their peak frame comes from it too,
# so nothing is inlined twice
if args.timelapse and not args.nas:
    tlDir = Path(config.input('reimaging/thumbnailsTimelapse/'))
    onManifold = set(zip(wells['plateId'].astype(str), wells['wellId'].astype(str)))
    man = pd.read_csv(tlDir / 'manifest.csv')
    man = man[[(p, w) in onManifold for p, w in
               zip(man['plateId'].astype(str), man['wellId'].astype(str))]]
    for (plate, well), wg in man.groupby(['plateId', 'wellId'], sort=True):
        key = (str(plate), str(well))
        byFrame = {}
        for _, r in wg.sort_values('frame').iterrows():
            jpg = tlDir / str(r['file'])
            if not jpg.exists():
                missing += 1
                continue
            byFrame[int(r['frame'])] = inline(jpg)
        if not byFrame:
            continue
        frameIdx[key] = byFrame
        pf = int(wg['peakFrame'].iloc[0])
        thumbIdx[key] = byFrame.get(pf, next(iter(byFrame.values())))
    print(f'{len(frameIdx)} wells with a full timelapse ({len(thumbList)} frames inlined)')

# every remaining well keeps its single peak-biomass still from the 112px pack
if not args.nas and (not args.timelapse or args.stills == 'all'):
    packDir = Path(config.input('reimaging/thumbnails112/'))
    nStill = 0
    for _, r in wells.iterrows():
        key = (str(r['plateId']), str(r['wellId']))
        if key in thumbIdx:
            continue
        jpg = packDir / str(r['file'])
        if not jpg.exists():
            missing += 1
            continue
        thumbIdx[key] = inline(jpg)
        nStill += 1
    print(f'{nStill} peak-biomass stills inlined from the pack')
print(f'{len(thumbList)} thumbnails inlined in total ({missing} missing)')

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
        fr = frameIdx.get((plate, well))
        pt = {'x': round(float(r['umap1']), 4), 'y': round(float(r['umap2']), 4),
              'gene': str(r['mutant']), 'locus': str(r['geneLocus']), 'func': str(r['function']),
              'plate': plate, 'well': well, 'isWt': str(r['mutant']) == WT_LABEL,
              'proj': False, 'pcolor': '',
              'pf': int(m['peakFrame']) if m is not None else -1,
              'bio': round(float(m['peakBiomass']), 5) if m is not None else None,
              'ti': thumbIdx.get((plate, well), -1)}
        if fr:      # frame -> thumbnail index; carried only for the timelapse pack's wells
            pt['fr'] = {str(k): v for k, v in fr.items()}
        pts.append(pt)
    pointsByParam[f'{int(nn)}_{md}'] = pts

projPoints = []                      # plots 1 and 2 draw the atlas alone; kept so the viewer JS is unchanged
# Candidate asset locations, tried in order at load time. All RELATIVE: an absolute path would both
# break the package's no-absolute-paths rule and pin the file to one machine's mount spelling. The
# page walks up from wherever it happens to sit looking for these tails, so moving the HTML (to
# another share, a USB stick, a web server) keeps it working as long as the assets are somewhere
# above it. Anything stranger than that is handled by the in-page folder box.
ASSET_TAILS = {
    'mp4': ['Sehna/multiphenotype-data-06-25-26/reimagingData',
            'multiphenotype-data-06-25-26/reimagingData', 'reimagingData'],
    'frames': ['data/v2/reimaging/framesPack', 'v2/reimaging/framesPack', 'framesPack'],
}
if args.nas:
    tails = [t.strip() for t in args.asset_tail.split(',') if t.strip()] + ASSET_TAILS[args.nas]
    firstBase = args.video_base if args.nas == 'mp4' else args.frame_base
    nasJson = json.dumps({'mode': args.nas, 'mag': args.mag, 'nFrames': args.nas_frames_count,
                          'firstBase': firstBase, 'tails': tails, 'depth': args.search_depth})
else:
    nasJson = json.dumps(None)
# --nas frames gives every well a timecourse, so the player is live there too
tlJson = json.dumps(bool(args.timelapse) or args.nas == 'frames')
# The timelapse pack is encoded at 96px, so a frame at native size is a stamp in a 460px sidebar --
# and here the image IS the panel. It is scaled to fit with the browser's default smoothing: at 4.8x,
# nearest-neighbour turns the microcolony texture into noise, while the interpolated version still
# reads as a biofilm. The footer says what the source resolution is, so nobody reads it as full
# detail. Plots 1/2 are untouched and keep native size.
thumbFillCss = ('  #thumb-wrap img, #thumb-wrap video { width:100%; height:auto; }'
                if (args.timelapse or args.nas) else '')
clickHint = 'peak-biomass thumbnail'
emptyHint = 'Click a point to view its peak-biomass frame'
if args.timelapse or args.nas:
    clickHint = 'play the timecourse'
    emptyHint = 'Click a point to play its timecourse'
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
{thumbFillCss}
  #no-thumb {{ color:#667; font-size:14px; text-align:center; padding:20px; }}
  #player {{ display:flex; align-items:center; gap:8px; padding:10px 14px; border-top:1px solid #334; flex-wrap:nowrap; }}
  #player.hidden {{ display:none; }}
  #player .btn {{ background:#2a2a4a; border:1px solid #445; border-radius:5px; color:#ccc; padding:5px 12px; font-size:14px; cursor:pointer; white-space:nowrap; }}
  #player .btn:hover {{ background:#3a3a5a; }}
  #player input[type=range] {{ accent-color:#4a6fa5; cursor:pointer; }}
  #tl-slider {{ flex:1; min-width:60px; }}
  #tl-fps {{ width:64px; flex-shrink:0; }}
  #tl-label {{ font-size:13px; color:#ccd; min-width:88px; white-space:nowrap; }}
  #tl-note {{ font-size:12px; color:#778; padding:8px 14px; border-top:1px solid #334; }}
  #tl-note.hidden {{ display:none; }}
  .tl-toggle {{ display:flex; align-items:center; gap:6px; font-size:14px; color:#aab; cursor:pointer; user-select:none; }}
  .tl-toggle.hidden {{ display:none; }}
  #asset-wrap {{ display:flex; align-items:center; gap:8px; font-size:13px; }}
  #asset-wrap.hidden, #asset-box.hidden {{ display:none; }}
  #asset-box {{ display:flex; align-items:center; gap:6px; }}
  #asset-input {{ background:#2a2a4a; border:1px solid #445; border-radius:5px; color:#eee; font-size:13px; padding:4px 8px; width:280px; outline:none; }}
  #asset-input::placeholder {{ color:#667; }}
</style></head><body>
<div id="controls">
  <div><label>n_neighbors:</label><span id="nn-btns"></span></div>
  <div><label>min_dist:</label><span id="md-btns"></span></div>
  <label class="tl-toggle" id="tl-only-wrap"><input id="tl-only" type="checkbox" />Timelapse wells only</label>
  <span id="asset-wrap" class="hidden"><span id="asset-status"></span>
    <span id="asset-box" class="hidden"><input id="asset-input" type="text" placeholder="path to the image folder" /><button class="param-btn" id="asset-set">Use</button></span></span>
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
    <div id="info">scroll=pan · ctrl+scroll=zoom · drag=pan · shift-drag=box · dblclick=reset · click={clickHint}</div>
  </div>
  <div id="sidebar" class="hidden">
    <div id="thumb-header"></div>
    <button id="thumb-close">&times;</button>
    <div id="thumb-wrap"><div id="no-thumb">{emptyHint}</div></div>
    <div id="player" class="hidden">
      <button class="btn" id="tl-play">&#9654; Play</button>
      <input id="tl-slider" type="range" min="0" step="1" />
      <span id="tl-label"></span>
      <label style="font-size:13px;color:#99a">fps</label>
      <input id="tl-fps" type="range" min="1" max="24" step="1" value="8" />
      <span id="tl-fpslabel" style="font-size:13px;color:#ccd">8</span>
    </div>
    <div id="tl-note" class="hidden"></div>
  </div>
</div>
<script>
const HAS_TL = {tlJson};
// Lab builds only: images are referenced from the mounted trees rather than inlined. Every URL is
// derived from the point's own plate/well, so no per-point URL list is carried. The slug must match
// the exporters' thumbKey() exactly: every run of non-alphanumerics collapses to a single '-'.
const NAS = {nasJson};
const allData = {dataJson};
function nasSlug(s) {{ return String(s).replace(/[^A-Za-z0-9]+/g,'-'); }}
// ── where are the assets? resolved at load, not baked in ──────────────────────
// The build can only guess where this file will end up. So it ships a list of directory TAILS and
// the page hunts for them: the build's own guess first, then '', '../', '../../' ... up to depth,
// then whatever the user last typed into the folder box (remembered). That makes the file portable
// to any location that has the assets somewhere at or above it.
let nasBase = null;
function nasCandidates() {{
  const out = [];
  const saved = (()=>{{ try {{ return localStorage.getItem('upulliAssetBase'); }} catch(e) {{ return null; }} }})();
  if (saved) out.push(saved);
  if (NAS.firstBase) out.push(NAS.firstBase);
  for (let k = 0; k <= NAS.depth; k++) {{
    const up = '../'.repeat(k);
    for (const tail of NAS.tails) out.push(up + tail);
  }}
  return [...new Set(out)];
}}
function videoUrlAt(base,p) {{
  return base + '/' + encodeURI(p.plate) + '/processedImages/' +
         encodeURI(p.well + '_' + NAS.mag + '_processed.mp4');
}}
function frameUrlAt(base,p,t) {{
  return base + '/' + nasSlug(p.plate) + '/' + nasSlug(p.plate + '__' + p.well) +
         '_t' + String(t).padStart(2,'0') + '.jpg';
}}
const nasVideoUrl = p => videoUrlAt(nasBase,p);
const nasFrameUrl = (p,t) => frameUrlAt(nasBase,p,t);

function probeBase(base,p) {{
  return new Promise(res => {{
    let done=false; const fin=ok=>{{ if(!done){{done=true; res(ok);}} }};
    if (NAS.mode === 'frames') {{
      const im=new Image(); im.onload=()=>fin(true); im.onerror=()=>fin(false);
      im.src=frameUrlAt(base,p,0); setTimeout(()=>fin(false),4000);
    }} else {{
      const v=document.createElement('video'); v.preload='metadata';
      v.onloadedmetadata=()=>fin(true); v.onerror=()=>fin(false);
      v.src=videoUrlAt(base,p); setTimeout(()=>fin(false),4000);
    }}
  }});
}}
async function resolveAssets() {{
  if (!NAS) return;
  const probe = points.find(q=>q.pf>=0) || points[0];
  if (!probe) return;
  setAssetStatus('searching', 'Locating images...');
  for (const base of nasCandidates()) {{
    if (await probeBase(base, probe)) {{
      nasBase = base;
      setAssetStatus('ok', 'Images found');
      if (curPoint) showThumb(curPoint);          // redraw whatever is open with real URLs
      return;
    }}
  }}
  setAssetStatus('fail', 'Images not found - set folder');
}}
function setAssetStatus(state,text) {{
  const el=document.getElementById('asset-status'); if(!el) return;
  el.textContent=text;
  el.style.color = state==='ok' ? '#7ecfff' : state==='fail' ? '#ffa0a0' : '#889';
  document.getElementById('asset-box').classList.toggle('hidden', state!=='fail');
}}
// frames available for a point: the inlined map, or the whole pack range in --nas frames
function frameList(p) {{
  if (NAS && NAS.mode === 'frames') return Array.from({{length: NAS.nFrames}}, (_,i) => i);
  return p.fr ? Object.keys(p.fr).map(Number).sort((a,b)=>a-b) : [];
}}
function frameSrc(p,t) {{
  if (NAS && NAS.mode === 'frames') return nasFrameUrl(p,t);
  const ti = p.fr ? p.fr[String(t)] : undefined;
  return ti === undefined ? '' : THUMBS[ti];
}}
function peakSrc(p) {{
  if (NAS && NAS.mode === 'frames') return nasFrameUrl(p, p.pf);
  return (p.ti>=0 && THUMBS[p.ti]) ? THUMBS[p.ti] : '';
}}
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
// Timelapse builds carry every frame for the representative replicates only (p.fr), so those points
// are ringed and can be isolated; every other point keeps its single peak-biomass still.
// In --nas builds every well has a timecourse, so there is nothing to ring or to filter down to.
let tlOnly = false;
const playable = p => NAS ? NAS.mode !== null : !!p.fr;
const shown = p => !tlOnly || playable(p);
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
  for (let i=0;i<points.length;i++) {{ const p=points[i]; if (p.proj || !shown(p)) continue;
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
  // wells that carry a full timecourse: a faint ring, drawn under the annotation/search overlays so
  // it marks what is playable without competing with them
  if (HAS_TL && !tlOnly && !NAS) for (const i of visibleIdx) {{ const p=points[i]; if (!p.fr) continue;
    const [sx,sy]=w2s(p.x,p.y); ctx.beginPath(); ctx.arc(sx,sy,r*1.7,0,Math.PI*2);
    ctx.strokeStyle='rgba(255,255,255,0.5)'; ctx.lineWidth=1.2; ctx.stroke(); }}
  if (annoActive) for (const i of visibleIdx) {{ const p=points[i], hl=pointHl(p); if (!hl) continue;
    const [sx,sy]=w2s(p.x,p.y); ctx.beginPath(); ctx.arc(sx,sy,r*1.3,0,Math.PI*2);
    ctx.strokeStyle=hexToRgba(functionColors[hl],0.95); ctx.lineWidth=2.5; ctx.stroke(); }}
  if (searchActive) for (const i of highlightSet) {{ const p=points[i]; if (p.proj || !shown(p)) continue; const [sx,sy]=w2s(p.x,p.y);
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
  for (let i=0;i<points.length;i++) {{ if (!shown(points[i])) continue;
    const [px,py]=w2s(points[i].x,points[i].y); const d2=(px-sx)**2+(py-sy)**2;
    if (d2<bd2) {{ bd2=d2; best=i; }} }} return best; }}
function applySearch() {{
  if (!searchTerm) {{ highlightSet=null; searchCount.textContent=''; return; }}
  highlightSet=new Set();
  for (let i=0;i<points.length;i++) {{ const p=points[i]; if (!shown(p)) continue;
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
// ── the timecourse player (inert unless the open well carries frames) ─────────
const player=document.getElementById('player'), tlNote=document.getElementById('tl-note');
const tlSlider=document.getElementById('tl-slider'), tlLabel=document.getElementById('tl-label');
const tlPlayBtn=document.getElementById('tl-play'), tlFpsSl=document.getElementById('tl-fps');
const tlFpsLab=document.getElementById('tl-fpslabel');
let curPoint=null, tlFrames=[], tlIdx=0, tlPlaying=false, tlTimer=null, tlFps=8;

// ── frame preloading (linked builds only) ────────────────────────────────────
// Inlined frames are data URIs and are decoded the moment they are assigned. LINKED frames are one
// network fetch each, so naively swapping img.src on a fixed timer just outruns the source: the
// slider advances (it is only a number) while the visible image stays on the last frame that
// actually decoded -- the peak frame the sidebar opened on. So on open we pull the whole well into
// the browser cache, and playback below waits for a frame rather than assuming it is there.
let preImgs=[], preErr=[], preKey=null, preDone=0;
function frameReady(i) {{
  if (!NAS) return true;                       // data URIs: always ready
  if (preErr[i]) return true;                  // a frame that failed counts as ready, or waiting for
  const im=preImgs[i];                         // it would stall playback on that index forever
  return !!im && im.complete && im.naturalWidth>0;
}}
function preloadNote() {{
  if (!NAS || !curPoint) return;
  const n=tlFrames.length;
  tlNote.textContent = preDone<n
    ? `Loading frames from the image tree... ${{preDone}}/${{n}}`
    : `${{n}} hourly frames of this replicate, linked from the frame pack.`;
}}
function preloadFrames(p, frames) {{
  preImgs=[]; preErr=[]; preDone=0;
  if (!NAS) return;
  const key=p.plate+'|'+p.well; preKey=key;
  frames.forEach((t,i)=>{{
    const im=new Image();
    const settle=bad=>{{ if (preKey!==key) return; if (bad) preErr[i]=true; preDone++; preloadNote(); }};
    im.onload=()=>settle(false); im.onerror=()=>settle(true);
    im.src=frameSrc(p,t);                      // same URL tlSet uses -> a cache hit later
    preImgs[i]=im;
  }});
  preloadNote();
}}
function tlStop() {{ tlPlaying=false; if (tlTimer) clearInterval(tlTimer); tlTimer=null; tlPlayBtn.innerHTML='&#9654; Play'; }}
function tlStart() {{ if (tlFrames.length<2) return; tlPlaying=true; tlPlayBtn.innerHTML='&#10073;&#10073; Pause';
  if (tlTimer) clearInterval(tlTimer);
  tlTimer=setInterval(()=>{{
    const next=(tlIdx+1)%tlFrames.length;
    if (!frameReady(next)) return;             // hold, don't advance past what has arrived
    tlSet(next);
  }},1000/tlFps); }}
function tlSet(i) {{
  if (!curPoint || !tlFrames.length) return;
  tlIdx=i; tlSlider.value=i;
  const t=tlFrames[i], src=frameSrc(curPoint,t);
  const img=thumbWrap.querySelector('img');
  if (img && src) img.src=src;
  tlLabel.textContent=`t=${{t}}h`+(t===curPoint.pf?' (peak)':'');
}}
tlPlayBtn.addEventListener('click',()=>tlPlaying?tlStop():tlStart());
tlSlider.addEventListener('input',e=>{{ tlStop(); tlSet(+e.target.value); }});
tlFpsSl.addEventListener('input',e=>{{ tlFps=+e.target.value; tlFpsLab.textContent=tlFps; if (tlPlaying) tlStart(); }});
document.addEventListener('keydown',e=>{{
  if (e.target.matches('input')||player.classList.contains('hidden')) return;
  if (e.key==='ArrowRight') {{ tlStop(); tlSet(Math.min(tlIdx+1,tlFrames.length-1)); }}
  if (e.key==='ArrowLeft') {{ tlStop(); tlSet(Math.max(tlIdx-1,0)); }}
  if (e.key===' ') {{ e.preventDefault(); tlPlaying?tlStop():tlStart(); }}
}});

function showThumb(p) {{
  tlStop(); curPoint=p;
  sidebar.classList.remove('hidden');
  thumbHeader.innerHTML=headerHtml(p);
  // clicked before the asset hunt finished (or it failed): say so instead of showing a broken frame.
  // resolveAssets() re-opens whatever is on screen once it finds them.
  if (NAS && !nasBase) {{
    thumbWrap.innerHTML='<div id="no-thumb">Locating the image folder...</div>';
    tlFrames=[]; player.classList.add('hidden'); tlNote.classList.add('hidden');
    resize(); return;
  }}
  // --nas mp4: the well's own full-resolution movie, streamed from the mount by the <video> element,
  // which brings its own transport controls -- so the frame player stays hidden here.
  if (NAS && NAS.mode === 'mp4') {{
    thumbWrap.innerHTML=`<video src="${{nasVideoUrl(p)}}" controls autoplay loop muted playsinline></video>`;
    tlFrames=[]; player.classList.add('hidden'); tlNote.classList.remove('hidden');
    tlNote.textContent='Full-resolution processed.mp4, linked from the image tree. If nothing plays, '
      + 'the page is open somewhere the mount is not visible.';
    resize(); return;
  }}
  const uri = peakSrc(p);
  if (uri) thumbWrap.innerHTML=`<img src="${{uri}}" alt="biofilm image">`;
  else thumbWrap.innerHTML='<div id="no-thumb">No image available for this well</div>';
  const frames = HAS_TL ? frameList(p) : [];
  tlFrames=frames;
  if (frames.length) {{
    tlSlider.max=frames.length-1; player.classList.remove('hidden');
    tlNote.classList.remove('hidden');
    preloadFrames(p, frames);                     // linked builds: warm the cache before Play is hit
    if (!NAS) tlNote.textContent=`${{frames.length}} hourly frames of this replicate, 96 px previews scaled to the panel.`;
    tlSet(Math.max(0,frames.indexOf(p.pf)));      // opens on the well's own peak, as plots 1/2 show
  }} else {{
    player.classList.add('hidden');
    // in a timelapse build, say why this well is a still rather than leaving it looking broken
    if (HAS_TL) {{ tlNote.classList.remove('hidden');
      tlNote.textContent = uri
        ? 'Peak-biomass frame. The whole timecourse is carried only for the representative replicates of each mutant, ringed on the plot.'
        : 'No image is carried for this replicate.'; }}
    else tlNote.classList.add('hidden');
  }}
  resize();
}}
document.getElementById('thumb-close').addEventListener('click',()=>{{ tlStop(); curPoint=null; sidebar.classList.add('hidden'); thumbWrap.innerHTML=''; player.classList.add('hidden'); tlNote.classList.add('hidden'); resize(); }});
if (!HAS_TL || NAS) document.getElementById('tl-only-wrap').classList.add('hidden');
document.getElementById('tl-only').addEventListener('change',e=>{{ tlOnly=e.target.checked; applySearch(); draw(); }});
if (NAS) {{
  document.getElementById('asset-wrap').classList.remove('hidden');
  document.getElementById('asset-set').addEventListener('click', async ()=>{{
    const v=document.getElementById('asset-input').value.trim().replace(/\\/$/,'');
    if (!v) return;
    setAssetStatus('searching','Checking...');
    const probe=points.find(q=>q.pf>=0)||points[0];
    if (await probeBase(v,probe)) {{
      nasBase=v;
      try {{ localStorage.setItem('upulliAssetBase',v); }} catch(e) {{}}   // file:// may refuse; harmless
      setAssetStatus('ok','Images found');
      if (curPoint) showThumb(curPoint);
    }} else setAssetStatus('fail','Not found there - check the path');
  }});
}}
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
  if (idx>=0) {{ const p=points[idx];
    tooltip.innerHTML=headerHtml(p).concat(`<div class="meta">click to ${{(NAS||(HAS_TL&&p.fr))?'play its timecourse':'view frame'}}</div>`);
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
if (NAS) resolveAssets();      // async: the plot is usable while the image folder is being found
</script></body></html>
"""

outHtml.write_text(html)
mb = len(html.encode('utf-8')) / 1e6
selKey = f'{args.nn}_{args.md}'
sel = pointsByParam.get(selKey, next(iter(pointsByParam.values()), []))
print(f'Wrote {outHtml}  ({mb:.1f} MB)')
print(f'  combos: {list(pointsByParam)} | {len(thumbList)} thumbnails | '
      f'selected {selKey}: {len(sel)} replicates')
if args.nas:
    print(f'  all {len(sel)} replicates play a timecourse, linked (not inlined). '
          f'LAB BUILD -- opens only where the image tree is mounted; do not deposit.')
elif args.timelapse:
    nTl = sum(1 for p in sel if p.get('fr'))
    print(f'  {nTl} of {len(sel)} replicates play a timecourse '
          f'({len({p["gene"] for p in sel if p.get("fr")})} mutants covered); '
          f'{sum(1 for p in sel if not p.get("fr") and p["ti"] >= 0)} keep a peak-biomass still')
