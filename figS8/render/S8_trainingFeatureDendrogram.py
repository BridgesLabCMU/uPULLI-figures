#!/usr/bin/env python3
"""figS8 render — Fig. 2C with the FEATURES clustered instead of held in family order.

Figure 2C fixes the 27 features in their family order (biomass, colony-segmentation, whole-image)
and asks how the eight training mutants differ. This asks the transposed question: which
*measurements* behave alike across those mutants. The matrix, colors and scale are identical to
2C -- only the row order changes, and a dendrogram is added.

Rows are clustered by Ward linkage on Euclidean distances between feature profiles (each row is a
feature's z-scores across the 8 mutants, already standardized in the bundled matrix, so no further
scaling is applied). Columns stay in the fixed strain order -- mutants are deliberately NOT
re-clustered, so this panel is directly comparable to 2C column by column.

Clustering reorders the rows, so 2C's contiguous family brackets no longer apply; a feature-family
color strip on the right carries that information instead, using the same three classes and colors as
Fig. S2D. Where families interleave, that is the result.

Reads:  data/trainingHeatmap_peakBiomass_featureMatrix.csv
Writes: figures/S8_trainingFeatureDendrogram.{png,svg}
"""
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS8/ for figlib
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
from figlib import (config, plotting, STRAIN_ORDER, DISPLAY_NAMES,
                    FAMILY_COLORS, FAMILY_LABELS, FAMILY_ORDER, featureFamily)

plotting.setStyle(extra={'font.size': 30, 'axes.titlesize': 34, 'axes.labelsize': 30,
                         'xtick.labelsize': 36, 'ytick.labelsize': 34, 'axes.linewidth': 2})
cellSize = 0.7

canonHaralick = ['energy', 'contrast', 'correlation', 'variance', 'inverse_difference_moment', 'sum_average',
                 'sum_variance', 'sum_entropy', 'entropy', 'difference_variance', 'difference_entropy', 'imc1', 'imc2']
haralickPretty = {'energy': 'Energy', 'contrast': 'Contrast', 'correlation': 'Correlation', 'variance': 'Variance',
                  'inverse_difference_moment': 'Inv Diff Moment', 'sum_average': 'Sum Average', 'sum_variance': 'Sum Variance',
                  'sum_entropy': 'Sum Entropy', 'entropy': 'Entropy', 'difference_variance': 'Diff Variance',
                  'difference_entropy': 'Diff Entropy', 'imc1': 'IMC 1', 'imc2': 'IMC 2'}
shortColony = {'biomass': r'Biofilm Biomass (a.u.)', 'nColonies': r'Colonies (count)', 'colony_area_um2_mean': r'Area ($\mu$m$^2$)',
               'colony_area_um2_std': r'Area Variability ($\mu$m$^2$)', 'colony_bgCV': r'Intensity Variability (CV)',
               'colony_centroidOffset_um_mean': r'Radial Offset ($\mu$m)', 'colony_eccentricity_mean': r'Eccentricity',
               'colony_majorAxisLength_um_mean': r'Major Axis Length ($\mu$m)', 'colony_meanIntensity_mean': r'Intensity (a.u.)',
               'colony_meanIntensity_kurtosis': r'Intensity Kurtosis', 'colony_mstEdgeMax_um_mean': r'Max Distance ($\mu$m)',
               'colony_nnDistance1_um_mean': r'Nearest Neighbor ($\mu$m)', 'colony_nnDistance1_um_std': r'NN Variability ($\mu$m)'}


def prettyFeatureName(b):
    if b in shortColony:
        return shortColony[b]
    if b.startswith('whole_entropy'):
        return 'Global Entropy'
    if b.startswith('whole_haralick_'):
        k = re.sub(r'_(mean|std|var)$', '', b.replace('whole_haralick_', ''))
        if k.isdigit() and int(k) < len(canonHaralick):
            k = canonHaralick[int(k)]
        return haralickPretty.get(k, k)
    return b


# ── cluster the FEATURES (rows); columns keep the fixed strain order ──────────
heat = pd.read_csv(config.TABLES / 'trainingHeatmap_peakBiomass_featureMatrix.csv', index_col=0)
heat = heat[[m for m in STRAIN_ORDER if m in heat.columns]]
# Cluster on the available information: undefined cells (vpsL's colony features) are treated as
# the row mean, which is 0 in z-space. They are still DRAWN as missing (grey), not as 0.
Z = linkage(pdist(np.nan_to_num(heat.to_numpy(dtype=float), nan=0.0), metric='euclidean'),
            method='ward')
dend = dendrogram(Z, no_plot=True)
rowOrder = dend['leaves']
heat = heat.iloc[rowOrder]
nF, nM = heat.shape
families = [featureFamily(b) for b in heat.index]
print(f'{nF} features clustered (Ward/euclidean) x {nM} mutants')

# ── layout: dendrogram | labels | heatmap | family strip | colorbar | legend ──
dendIn, stripIn, padIn = 2.2, 0.35, 0.3
fw = nM * cellSize + 9
fh = nF * cellSize + 4
def widestTextIn(labels, fontsize):
    """Widest rendered label, in inches. Measured rather than estimated from character counts, so
    mathtext (Delta, italics, superscripts) is sized correctly and the dendrogram can sit right up
    against the label column without a slack gap."""
    probe = plt.figure(figsize=(1, 1), dpi=200)
    renderer = probe.canvas.get_renderer()
    widest = max((probe.text(0, 0, s, fontsize=fontsize).get_window_extent(renderer=renderer).width
                  for s in labels), default=0.0) / probe.dpi
    plt.close(probe)
    return widest

yLabelIn = widestTextIn([prettyFeatureName(b) for b in heat.index], mpl.rcParams['ytick.labelsize'])
mutLabelIn = widestTextIn([DISPLAY_NAMES[m] for m in heat.columns], mpl.rcParams['xtick.labelsize']) * 0.71
# dendGapIn = clearance between the dendrogram leaves and the longest feature label
dendGapIn = 0.25
leftIn = dendIn + dendGapIn + yLabelIn + padIn
left, bottom = leftIn / fw, (mutLabelIn + padIn) / fh
heatW, heatH = (nM * cellSize) / fw, (nF * cellSize) / fh

fig = plt.figure(figsize=(fw, fh), dpi=200)
axH = fig.add_axes([left, bottom, heatW, heatH])
axD = fig.add_axes([(leftIn - yLabelIn - dendGapIn - dendIn) / fw, bottom, dendIn / fw, heatH])
axS = fig.add_axes([left + heatW + padIn / fw, bottom, stripIn / fw, heatH])
cbarX = left + heatW + (2 * padIn + stripIn) / fw
cbarW = 0.035
axC = fig.add_axes([cbarX, bottom, cbarW, heatH])

# dendrogram with the LEAVES on the inside edge, so each leaf points at its heatmap row
# (merge distance increases leftwards; the root is at the far left).
icoord, dcoord = np.array(dend['icoord']), np.array(dend['dcoord'])
yMin, yMax = icoord.min(), icoord.max()
scaleY = (nF - 1) / (yMax - yMin)
maxD = dcoord.max()
for xs, ys in zip(dcoord, icoord):
    axD.plot(xs, [(y - yMin) * scaleY for y in ys], color='black', linewidth=2.5)
axD.set_xlim(maxD * 1.02, -maxD * 0.02)
axD.set_ylim(nF - 0.5, -0.5)
axD.axis('off')

axS.imshow(np.array([[FAMILY_ORDER.index(f)] for f in families]), aspect='auto',
           cmap=mpl.colors.ListedColormap([FAMILY_COLORS[k] for k in FAMILY_ORDER]),
           vmin=0, vmax=len(FAMILY_ORDER) - 1)
axS.axis('off')

cmap = mpl.colormaps['RdBu_r'].copy(); cmap.set_bad('#b0b0b0')
im = axH.imshow(np.ma.masked_invalid(heat.values), cmap=cmap, vmin=-3, vmax=3,
                interpolation='nearest', aspect='auto')
axH.set_yticks(np.arange(nF)); axH.set_yticklabels([prettyFeatureName(b) for b in heat.index])
axH.yaxis.tick_right()
axH.set_xticks(np.arange(nM)); axH.set_xticklabels([DISPLAY_NAMES[m] for m in heat.columns], rotation=45, ha='right')
axH.tick_params(axis='x', pad=4); axH.tick_params(axis='y', pad=6)
axH.set_title('Peak Biofilm Biomass Frame', pad=24, fontweight='bold')
fig.colorbar(im, cax=axC).set_label('Z-score')

present = [k for k in FAMILY_ORDER if k in set(families)]
legendX = 1.0 + ((cbarX + cbarW + 0.012) - (left + heatW)) / heatW   # just right of the colorbar
axH.legend(handles=[Line2D([0], [0], marker='s', linestyle='none', markersize=16,
                           markerfacecolor=FAMILY_COLORS[k], markeredgecolor='black', label=FAMILY_LABELS[k])
                    for k in present],
           frameon=False, fontsize=24, loc='upper left', bbox_to_anchor=(legendX, 1.0))

out = config.ensure(config.FIGURES) / 'S8_trainingFeatureDendrogram'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
fig.savefig(str(out) + '.svg', bbox_inches='tight')
plt.close(fig)
print(f'Saved: {out}.png')
print('  cluster order: ' + ', '.join(prettyFeatureName(b) for b in heat.index[:6]) + ' ...')
