#!/usr/bin/env python3
"""figS6 render (VERTICAL): full-atlas dendrogram + peak-biomass heatmap, transposed.

Portrait version of S7A: mutants on the y-axis (leaf order, dendrogram on the left), features on the
x-axis (labels rotated 30 deg), feature-group brackets along the top, functional-annotation color strip
between the dendrogram and the heatmap. Large tick fonts.

Reads:  data/fullAtlas_pcaLinkage_linkage.npy, data/fullAtlas_pcaLinkage_cluster_order.csv,
        data/fullAtlas_peakBiomass_featureMatrix.csv
Writes: figures/S6A_fullDendrogramHeatmap_vertical.{png,svg}
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

plotting.setStyle(extra={'font.size': 40, 'axes.linewidth': 2})
WT = 'WT'
cellSize = 0.7
FONT_MUT, FONT_FEAT = 38, 38

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
    return 'biomass' if base == 'biomass' else 'haralick' if 'haralick' in base else 'entropy' if 'entropy' in base else 'colony'


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
    return 'WT' if annotation == 'WT' else f'${cleanName(display)}$'


Z = np.load(config.TABLES / 'fullAtlas_pcaLinkage_linkage.npy')
order = pd.read_csv(config.TABLES / 'fullAtlas_pcaLinkage_cluster_order.csv')
featureMatrix = pd.read_csv(config.TABLES / 'fullAtlas_peakBiomass_featureMatrix.csv', index_col=0)

ordered = order['mutant'].astype(str).tolist()
leafColors = order['color'].tolist()
mutLabels = [formatMutantLabel(d, a) for d, a in zip(order['display'].astype(str), order['annotation'].astype(str))]
heat = featureMatrix[ordered].T                              # mutants (rows, leaf order) x features (cols)
featNames = [prettyFeatureName(b) for b in heat.columns]
groups = [featureGroup(b) for b in heat.columns]
nM, nF = heat.shape

# ── layout (inches) ──
padIn, dendIn, stripIn, bracketIn, cbarIn = 0.4, 6.0, 0.6, 2.2, 1.8
mutLabelIn = 3.0
maxFeatChars = max(len(f) for f in featNames)
featLabelIn = maxFeatChars * (FONT_FEAT / 72.0) * 0.55 + 0.5
heatW, heatH = nF * cellSize, nM * cellSize
fwPre = padIn + dendIn + stripIn + heatW + mutLabelIn + cbarIn + padIn
fhPre = padIn + featLabelIn + heatH + bracketIn + padIn
scale = min(1.0, 60000.0 / (300.0 * fhPre))   # never upscale (158 rows is already very tall)
fig = plt.figure(figsize=(fwPre * scale, fhPre * scale), dpi=150)

xDend = padIn / fwPre
xStrip = (padIn + dendIn) / fwPre
xHeat = (padIn + dendIn + stripIn) / fwPre
yHeat = (padIn + featLabelIn) / fhPre
wDend, wStrip, wHeat = dendIn / fwPre, stripIn / fwPre, heatW / fwPre
hHeat = heatH / fhPre

axDend = fig.add_axes([xDend, yHeat, wDend, hHeat])
axStrip = fig.add_axes([xStrip, yHeat, wStrip, hHeat])
axH = fig.add_axes([xHeat, yHeat, wHeat, hHeat])
axBr = fig.add_axes([xHeat, yHeat + hHeat + 0.002, wHeat, bracketIn / fhPre]); axBr.set_xlim(-0.5, nF - 0.5); axBr.set_ylim(0, 1); axBr.axis('off')
axC = fig.add_axes([(padIn + dendIn + stripIn + heatW + mutLabelIn + 0.4) / fwPre, yHeat + hHeat * 0.25, 0.22 / fwPre, hHeat * 0.5])

# heatmap (rows = mutants, cols = features)
im = axH.imshow(heat.values, cmap='RdBu_r', vmin=-3, vmax=3, interpolation='nearest', aspect='auto')
im.set_rasterized(True)
axH.set_xticks(np.arange(nF)); axH.set_xticklabels(featNames, rotation=30, ha='right', rotation_mode='anchor', fontsize=FONT_FEAT)
axH.yaxis.set_ticks_position('right'); axH.yaxis.set_label_position('right')
axH.set_yticks(np.arange(nM)); axH.set_yticklabels(mutLabels, fontsize=FONT_MUT)
cb = fig.colorbar(im, cax=axC); cb.set_label('Z-score', fontsize=34); cb.ax.tick_params(labelsize=28)

# functional strip (one color cell per mutant row)
axStrip.imshow([[i] for i in range(nM)], aspect='auto', cmap=mpl.colors.ListedColormap(leafColors)).set_rasterized(True)
axStrip.axis('off')

# left dendrogram (leaves = mutant rows; distance increases leftward)
dend = dendrogram(Z, no_plot=True)
maxD = max(max(d) for d in dend['dcoord'])
for xs, ys in zip(dend['icoord'], dend['dcoord']):
    rows = [(v - 5.0) / 10.0 for v in xs]
    axDend.plot(ys, rows, color='black', linewidth=2)
for i in range(nM):
    axDend.plot(0, i, 'o', color=leafColors[i], markersize=(11 if ordered[i] == WT else 9), clip_on=False)
axDend.set_xlim(maxD * 1.02, -0.02 * maxD); axDend.set_ylim(nM - 0.5, -0.5); axDend.axis('off')


def bracket(x0, x1, label):
    axBr.plot([x0, x1], [0.18, 0.18], lw=3, color='black', clip_on=False)
    axBr.plot([x0, x0], [0.05, 0.18], lw=3, color='black', clip_on=False)
    axBr.plot([x1, x1], [0.05, 0.18], lw=3, color='black', clip_on=False)
    axBr.text((x0 + x1) / 2, 0.30, label, ha='center', va='bottom', fontsize=34)


def rng(name):
    idx = [i for i, g in enumerate(groups) if g == name]
    return (min(idx), max(idx)) if idx else None


if (r := rng('colony')):
    bracket(r[0], r[1], 'Colony Features')
if (r := rng('haralick')):
    bracket(r[0], r[1], 'Whole-Image Haralick Features')

out = config.ensure(config.FIGURES) / 'S6A_fullDendrogramHeatmap_vertical'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
fig.savefig(str(out) + '.svg', dpi=200, bbox_inches='tight')
plt.close(fig)
print(f'Saved: {out}.png  ({nM} mutants x {nF} features)')
