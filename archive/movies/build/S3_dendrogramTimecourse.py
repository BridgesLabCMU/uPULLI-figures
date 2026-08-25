#!/usr/bin/env python3
"""movies build — Movie S3: the Figure 3D dendrogram + heatmap animated across the timecourse.

Figure 3D shows the functional-annotation subset (50 mutants) at each mutant's peak-biomass frame.
This renders the identical figure once per hourly timepoint instead, so the reader can watch each
feature class evolve as the biofilms develop. Nothing moves except the heatmap: the dendrogram, the
leaf order, the labels and the annotation strip are fixed, so change on screen is change in the
features, not in the clustering.

Layout, colors, brackets and fonts are a direct port of fig3/render/3D_dendrogramHeatmap.py. The only
difference is the title, which per spec is just "Time = X h" -- there is no other text.

Encoded twice from the same PNG frames: MP4 (H.264, yuv420p, +faststart) for general use, and AVI
(MJPEG) for journals that require it. Frames are z-scored across the full 158-mutant atlas -- the
same normalization behind Fig. 3D -- and then subset to the 50 functional leaves, so a color means
the same thing in the movie and in the printed panel.

Reads:  data/functional_pcaLinkage_{linkage.npy, cluster_order.csv}, data/fullAtlas_frameMatrices.npz
Writes: videos/movieS3_dendrogramTimecourse.{mp4,avi}

Usage:  python build/S3_dendrogramTimecourse.py [--fps 8] [--width 1920] [--keep-frames]
"""
import sys
import re
import shutil
import argparse
import subprocess
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> movies/ for figlib
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from matplotlib.transforms import Bbox
from figlib import config, plotting

ap = argparse.ArgumentParser()
ap.add_argument('--fps', type=int, default=8)
ap.add_argument('--width', type=int, default=1920, help='output width in px (height keeps aspect)')
ap.add_argument('--dpi', type=int, default=60, help='render dpi before the ffmpeg rescale')
ap.add_argument('--keep-frames', action='store_true', help='leave the per-frame PNGs on disk')
args = ap.parse_args()

plotting.setStyle(extra={'font.size': 44, 'axes.titlesize': 44, 'axes.labelsize': 44,
                         'xtick.labelsize': 44, 'ytick.labelsize': 44, 'axes.linewidth': 2})
WT = 'WT'
cellSize, titleFontsize = 0.7, 52
NODATA = '#dcdcdc'          # colony features do not exist before t5 -- shown, not silently zeroed

canonicalHaralickOrder = ['energy', 'contrast', 'correlation', 'variance', 'inverse_difference_moment',
                          'sum_average', 'sum_variance', 'sum_entropy', 'entropy',
                          'difference_variance', 'difference_entropy', 'imc1', 'imc2']
haralickPretty = {'energy': 'Energy', 'contrast': 'Contrast', 'correlation': 'Correlation', 'variance': 'Variance',
                  'inverse_difference_moment': 'Inv Diff Moment', 'sum_average': 'Sum Average', 'sum_variance': 'Sum Variance',
                  'sum_entropy': 'Sum Entropy', 'entropy': 'Entropy', 'difference_variance': 'Diff Variance',
                  'difference_entropy': 'Diff Entropy', 'imc1': 'IMC 1', 'imc2': 'IMC 2'}
prettyMap = {'biomass': 'Biofilm Biomass (a.u.)', 'nColonies': 'Colonies (count)', 'colony_area_um2_mean': r'Area ($\mu$m$^2$)',
             'colony_area_um2_std': r'Area Variability ($\mu$m$^2$)', 'colony_bgCV': 'Intensity Variability (CV)',
             'colony_centroidOffset_um_mean': r'Radial Offset ($\mu$m)', 'colony_eccentricity_mean': 'Eccentricity',
             'colony_majorAxisLength_um_mean': r'Major Axis Length ($\mu$m)', 'colony_meanIntensity_mean': 'Intensity (a.u.)',
             'colony_meanIntensity_kurtosis': 'Intensity Kurtosis', 'colony_mstEdgeMax_um_mean': r'Max Distance ($\mu$m)',
             'colony_nnDistance1_um_mean': r'Nearest Neighbor ($\mu$m)', 'colony_nnDistance1_um_std': r'NN Variability ($\mu$m)'}


def featureGroup(base):
    if base == 'biomass':
        return 'biomass'
    if 'haralick' in base:
        return 'haralick'
    return 'entropy' if 'entropy' in base else 'colony'


def prettyFeatureName(base):
    if base in prettyMap:
        return prettyMap[base]
    if base.startswith('whole_entropy'):
        return 'Global Entropy'
    if base.startswith('whole_haralick_'):
        key = re.sub(r'_(mean|std|var)$', '', base.replace('whole_haralick_', ''))
        if key.isdigit() and int(key) < len(canonicalHaralickOrder):
            key = canonicalHaralickOrder[int(key)]
        return haralickPretty.get(key, key)
    return base


def cleanName(s):
    s = str(s)
    if ';' in s:
        parts = [p.strip().replace('_', '') for p in s.split(';')]
        if len(parts) >= 2 and parts[0].startswith('VC') and parts[1].startswith('VC'):
            return f'{parts[0]}-{parts[1][-2:]}'
    return s.replace('_', '')


def formatMutantLabel(display, annotation):
    if annotation == 'WT':
        return 'WT'
    return cleanName(display) if annotation == 'Compound' else f'${cleanName(display)}$'


def drawBracket(ax, y0, y1, label):
    ax.plot([0.80, 0.80], [y0, y1], lw=3, color='black', clip_on=False)
    ax.plot([0.80, 0.97], [y0, y0], lw=3, color='black', clip_on=False)
    ax.plot([0.80, 0.97], [y1, y1], lw=3, color='black', clip_on=False)
    ax.text(0.42, (y0 + y1) / 2, label, rotation=90, va='center', ha='center', fontsize=44, clip_on=False)


# ── clustering + per-frame matrices ───────────────────────────────────────────
Z = np.load(config.TABLES / 'functional_pcaLinkage_linkage.npy')
order = pd.read_csv(config.TABLES / 'functional_pcaLinkage_cluster_order.csv')
npz = np.load(config.TABLES / 'fullAtlas_frameMatrices.npz', allow_pickle=False)
allMutants = [str(m) for m in npz['mutants']]
featureBases = [str(f) for f in npz['features']]
frames = [int(t) for t in npz['frames']]
mats = npz['matrices'].astype(np.float32)          # frames x 158 mutants x 27 features

ordered = order['mutant'].astype(str).tolist()
leafColors = order['color'].tolist()
labels = [formatMutantLabel(d, a) for d, a in zip(order['display'].astype(str), order['annotation'].astype(str))]
sel = [allMutants.index(m) for m in ordered]       # 50 functional leaves, in leaf order
dend = dendrogram(Z, no_plot=True)
nF, nM = len(featureBases), len(ordered)
print(f'{nM} leaves x {nF} features x {len(frames)} frames (t{frames[0]}-t{frames[-1]}h)')

cmap = mpl.colormaps['RdBu_r'].copy()
cmap.set_bad(NODATA)

# ── layout, identical to fig3/render/3D_dendrogramHeatmap.py ──────────────────
scale = min(2.0, 60000.0 / (300.0 * max(nM * cellSize + 7, 1)))
fw, fh = (nM * cellSize + 7) * scale, (nF * cellSize + 6) * scale
_vis = lambda l: re.sub(r'\$|\\mathit|\{|\}|\\Delta|\^|_', '', str(l))
maxFeatChars = max((len(prettyFeatureName(b)) for b in featureBases), default=10)
maxMutChars = min(max((len(_vis(l)) for l in labels), default=6), 14)
yLabelIn = maxFeatChars * (44 / 72.0) * 0.40
bracketIn, stripIn, padIn = 1.5, 0.5, 0.3
mutLabelIn = maxMutChars * (44 / 72.0) * 0.55 * 0.71
leftIn = bracketIn + yLabelIn + padIn
bottomIn = mutLabelIn + stripIn + 1.5 * padIn
left, bottom = leftIn / fw, bottomIn / fh
heatW, heatH = (nM * cellSize) / fw, (nF * cellSize) / fh

icoord, dcoord = dend['icoord'], dend['dcoord']
xMin = min(min(i) for i in icoord)
scaleX = (nM - 1) / (max(max(i) for i in icoord) - xMin)
maxD = max(max(d) for d in dcoord)
groups = [featureGroup(b) for b in featureBases]


def rangeOf(name):
    idx = [i for i, g in enumerate(groups) if g == name]
    return (min(idx), max(idx)) if idx else None


def buildFig(t, heat):
    fig = plt.figure(figsize=(fw, fh), dpi=args.dpi)
    axH = fig.add_axes([left, bottom, heatW, heatH])
    axD = fig.add_axes([left, bottom + heatH + 0.003, heatW, 0.07])
    axD.set_title(f'Time = {t} h', pad=50, fontsize=titleFontsize, fontweight='bold')
    axC = fig.add_axes([left + heatW + 0.015, bottom, 0.012, heatH])
    stripBottom = max((bottomIn - mutLabelIn - padIn - stripIn) / fh, 0.004)
    axFunc = fig.add_axes([left, stripBottom, heatW, stripIn / fh])
    axBracket = fig.add_axes([0.4 / fw, bottom, bracketIn / fw, heatH])

    for xs, ys in zip(icoord, dcoord):
        axD.plot([(x - xMin) * scaleX for x in xs], ys, color='black', linewidth=3)
    for i, m in enumerate(ordered):
        axD.plot(i, -0.05 * maxD, 'o', color=leafColors[i], markersize=(24 if m == WT else 20), clip_on=False)
    axD.set_xlim(-0.5, nM - 0.5); axD.set_ylim(-0.1 * maxD, maxD); axD.axis('off')

    im = axH.imshow(np.ma.masked_invalid(heat), cmap=cmap, vmin=-3, vmax=3,
                    interpolation='nearest', aspect='auto')
    axH.set_yticks(np.arange(nF)); axH.set_yticklabels([prettyFeatureName(b) for b in featureBases])
    axH.set_xticks(np.arange(nM)); axH.set_xticklabels(labels, rotation=45, ha='right', fontsize=44)
    axH.tick_params(axis='x', pad=4)
    fig.colorbar(im, cax=axC).set_label(r'Z-score')

    axFunc.imshow([range(nM)], aspect='auto', cmap=mpl.colors.ListedColormap(leafColors)); axFunc.axis('off')
    axBracket.set_xlim(0, 1); axBracket.set_ylim(nF - 0.5, -0.5); axBracket.axis('off')
    if (r := rangeOf('colony')):
        drawBracket(axBracket, r[0], r[1], 'Colony Segmentation-\nDerived Features')
    if (r := rangeOf('haralick')):
        drawBracket(axBracket, r[0], r[1], 'Whole Image\nHaralick Features')

    return fig


# The static panel crops with bbox_inches='tight', but a per-frame crop would let the frame size
# drift with the title width ("Time = 0 h" vs "Time = 30 h") and every frame in a video must be the
# same size. So: measure the tight box on the shortest- and longest-titled frames, union them, and
# save every frame against that one fixed box -- cropped like the figure, yet pixel-identical.
boxes = []
for t, i in ((frames[0], 0), (frames[-1], len(frames) - 1)):
    f = buildFig(t, mats[i][sel].T)
    f.canvas.draw()
    boxes.append(f.get_tightbbox(f.canvas.get_renderer()))
    plt.close(f)
bbox = Bbox.union(boxes).expanded(1.01, 1.01)
print(f'frame box {bbox.width:.1f} x {bbox.height:.1f} in -> '
      f'{round(bbox.width * args.dpi)} x {round(bbox.height * args.dpi)} px before rescale')

tmp = Path(tempfile.mkdtemp(prefix='movieS3_'))
for i, t in enumerate(frames):
    fig = buildFig(t, mats[i][sel].T)                          # -> features x leaves
    fig.savefig(tmp / f'f{i:04d}.png', dpi=args.dpi, facecolor='white', bbox_inches=bbox)
    plt.close(fig)
    if (i + 1) % 10 == 0 or i == len(frames) - 1:
        print(f'  rendered {i + 1}/{len(frames)} frames')

from PIL import Image
sizes = {Image.open(im).size
         for im in sorted(tmp.glob("f*.png"))}
if len(sizes) != 1:
    raise SystemExit(f'frames differ in size, ffmpeg would reject them: {sizes}')
print(f'  all {len(frames)} frames {sizes.pop()} px')

outDir = Path(config.ensure(config.VIDEOS))
common = ['-y', '-loglevel', 'error', '-framerate', str(args.fps), '-i', str(tmp / 'f%04d.png'),
          '-vf', f'scale={args.width}:-2:flags=lanczos']
mp4 = outDir / 'movieS3_dendrogramTimecourse.mp4'
avi = outDir / 'movieS3_dendrogramTimecourse.avi'
subprocess.run(['ffmpeg', *common, '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-crf', '18',
                '-movflags', '+faststart', '-r', str(args.fps), str(mp4)], check=True)
subprocess.run(['ffmpeg', *common, '-c:v', 'mjpeg', '-q:v', '3', '-r', str(args.fps), str(avi)], check=True)

if args.keep_frames:
    dest = outDir / 'movieS3_frames'
    shutil.rmtree(dest, ignore_errors=True); shutil.copytree(tmp, dest)
    print(f'  frames kept in {dest}')
shutil.rmtree(tmp, ignore_errors=True)
for p in (mp4, avi):
    print(f'Wrote {p}  ({p.stat().st_size / 1e6:.1f} MB)')
print(f'  {len(frames)} frames at {args.fps} fps = {len(frames) / args.fps:.1f} s')
