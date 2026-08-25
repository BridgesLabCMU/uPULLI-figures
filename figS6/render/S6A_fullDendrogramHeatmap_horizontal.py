#!/usr/bin/env python3
"""figS6 render: FULL reimaging atlas dendrogram + peak-biomass feature heatmap.

The full-atlas (all 158 mutants) counterpart of Fig 3D (which shows only the functional subset). Renders
FROM the saved PCA-Ward linkage + leaf order + z-scored peak-biomass feature matrix, no re-clustering.

Reads:  data/fullAtlas_pcaLinkage_linkage.npy, data/fullAtlas_pcaLinkage_cluster_order.csv,
        data/fullAtlas_peakBiomass_featureMatrix.csv
Writes: figures/S6A_fullDendrogramHeatmap_horizontal.{png,svg}
"""
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS6/ for figlib
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from figlib import config, plotting

plotting.setStyle(extra={'font.size': 44, 'axes.titlesize': 44, 'axes.labelsize': 44,
                         'xtick.labelsize': 44, 'ytick.labelsize': 44, 'axes.linewidth': 2})
WT = 'WT'
cellSize, titleFontsize = 0.7, 52

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
    if 'entropy' in base:
        return 'entropy'
    return 'colony'


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
    if annotation == 'Compound':
        return cleanName(display)
    return f'${cleanName(display)}$'


def drawBracket(ax, y0, y1, label):
    ax.plot([0.80, 0.80], [y0, y1], lw=3, color='black', clip_on=False)
    ax.plot([0.80, 0.97], [y0, y0], lw=3, color='black', clip_on=False)
    ax.plot([0.80, 0.97], [y1, y1], lw=3, color='black', clip_on=False)
    ax.text(0.42, (y0 + y1) / 2, label, rotation=90, va='center', ha='center', fontsize=44, clip_on=False)


# ── load the saved clustering + tables (full atlas) ────────────────────────────
Z = np.load(config.TABLES / 'fullAtlas_pcaLinkage_linkage.npy')
order = pd.read_csv(config.TABLES / 'fullAtlas_pcaLinkage_cluster_order.csv')
featureMatrix = pd.read_csv(config.TABLES / 'fullAtlas_peakBiomass_featureMatrix.csv', index_col=0)

ordered = order['mutant'].astype(str).tolist()
leafColors = order['color'].tolist()
labels = [formatMutantLabel(d, a) for d, a in zip(order['display'].astype(str), order['annotation'].astype(str))]
dend = dendrogram(Z, no_plot=True)
heat = featureMatrix[ordered]
nF, nM = heat.shape

# ── layout (inch-based; matches the v2 dendrogram figure) ──────────────────────
scale = min(2.0, 60000.0 / (300.0 * max(nM * cellSize + 7, 1)))
fw = (nM * cellSize + 7) * scale
fh = (nF * cellSize + 6) * scale
fig = plt.figure(figsize=(fw, fh), dpi=150)

_vis = lambda l: re.sub(r'\$|\\mathit|\{|\}|\\Delta|\^|_', '', str(l))
maxFeatChars = max((len(prettyFeatureName(b)) for b in heat.index), default=10)
maxMutChars = min(max((len(_vis(l)) for l in labels), default=6), 14)
yLabelIn = maxFeatChars * (44 / 72.0) * 0.40
bracketIn, stripIn, padIn = 1.5, 0.5, 0.3
mutLabelIn = maxMutChars * (44 / 72.0) * 0.55 * 0.71
leftIn = bracketIn + yLabelIn + padIn
bottomIn = mutLabelIn + stripIn + 1.5 * padIn
left, bottom = leftIn / fw, bottomIn / fh
heatW, heatH = (nM * cellSize) / fw, (nF * cellSize) / fh

axH = fig.add_axes([left, bottom, heatW, heatH])
axD = fig.add_axes([left, bottom + heatH + 0.003, heatW, 0.07])
axD.set_title('Peak Biofilm Biomass Frame', pad=50, fontsize=titleFontsize, fontweight='bold')
axC = fig.add_axes([left + heatW + 0.015, bottom, 0.012, heatH])
stripBottom = max((bottomIn - mutLabelIn - padIn - stripIn) / fh, 0.004)
axFunc = fig.add_axes([left, stripBottom, heatW, stripIn / fh])
axBracket = fig.add_axes([0.4 / fw, bottom, bracketIn / fw, heatH])

icoord, dcoord = dend['icoord'], dend['dcoord']
xMin = min(min(i) for i in icoord)
scaleX = (nM - 1) / (max(max(i) for i in icoord) - xMin)
maxD = max(max(d) for d in dcoord)
for xs, ys in zip(icoord, dcoord):
    axD.plot([(x - xMin) * scaleX for x in xs], ys, color='black', linewidth=3)
for i, m in enumerate(ordered):
    axD.plot(i, -0.05 * maxD, 'o', color=leafColors[i], markersize=(24 if m == WT else 20), clip_on=False)
axD.set_xlim(-0.5, nM - 0.5); axD.set_ylim(-0.1 * maxD, maxD); axD.axis('off')

im = axH.imshow(heat.values, cmap='RdBu_r', vmin=-3, vmax=3, interpolation='nearest', aspect='auto')
im.set_rasterized(True)   # keep the (very wide) heatmap raster under the SVG 32767px image limit
axH.set_yticks(np.arange(nF)); axH.set_yticklabels([prettyFeatureName(b) for b in heat.index])
axH.set_xticks(np.arange(nM)); axH.set_xticklabels(labels, rotation=45, ha='right', fontsize=38)
axH.tick_params(axis='x', pad=4)
fig.colorbar(im, cax=axC).set_label(r'Z-score')

axFunc.imshow([range(nM)], aspect='auto', cmap=mpl.colors.ListedColormap(leafColors)).set_rasterized(True); axFunc.axis('off')
axBracket.set_xlim(0, 1); axBracket.set_ylim(nF - 0.5, -0.5); axBracket.axis('off')
groups = [featureGroup(b) for b in heat.index]


def getRange(name):
    idx = [i for i, g in enumerate(groups) if g == name]
    return (min(idx), max(idx)) if idx else None


if (r := getRange('colony')):
    drawBracket(axBracket, r[0], r[1], 'Colony Segmentation-\nDerived Features')
if (r := getRange('haralick')):
    drawBracket(axBracket, r[0], r[1], 'Whole Image\nHaralick Features')

out = config.ensure(config.FIGURES) / 'S6A_fullDendrogramHeatmap_horizontal.png'
fig.savefig(out, dpi=300, bbox_inches='tight')
fig.savefig(str(out).replace('.png', '.svg'), dpi=200, bbox_inches='tight')  # dpi=200 -> raster < 32767px wide
plt.close(fig)
print(f'Saved: {out}  ({nF} features x {nM} leaves)')
