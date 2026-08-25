#!/usr/bin/env python3
"""fig2 render — Panel 2C, family-ordered variant (superseded by the clustered 2C).

Renders FROM the bundled z-scored feature matrix (features x mutants). RdBu_r (-3..3), fixed mutant
order, feature-group brackets (colony / whole-image haralick) on the left. No dendrogram.

Reads:  data/trainingHeatmap_peakBiomass_featureMatrix.csv
Writes: figures/2C_trainingHeatmap_familyOrder.{png,svg}
"""
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig2/ for figlib
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from figlib import config, plotting, STRAIN_ORDER, DISPLAY_NAMES

plotting.setStyle(extra={'font.size': 30, 'axes.titlesize': 34, 'axes.labelsize': 30,
                         'xtick.labelsize': 30, 'ytick.labelsize': 28, 'axes.linewidth': 2})
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


def featureGroup(b):
    return 'biomass' if b == 'biomass' else 'haralick' if 'haralick' in b else 'entropy' if 'entropy' in b else 'colony'


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


def drawBracket(ax, y0, y1, label):
    ax.plot([0.80, 0.80], [y0, y1], lw=3, color='black', clip_on=False)
    ax.plot([0.80, 0.97], [y0, y0], lw=3, color='black', clip_on=False)
    ax.plot([0.80, 0.97], [y1, y1], lw=3, color='black', clip_on=False)
    ax.text(0.38, (y0 + y1) / 2, label, rotation=90, va='center', ha='center', fontsize=30, transform=ax.transData, clip_on=False)


heat = pd.read_csv(config.TABLES / 'trainingHeatmap_peakBiomass_featureMatrix.csv', index_col=0)
heat = heat[[m for m in STRAIN_ORDER if m in heat.columns]]
nF, nM = heat.shape
fw = nM * cellSize + 9
fh = nF * cellSize + 4
fig = plt.figure(figsize=(fw, fh), dpi=200)

maxFeatChars = max((len(re.sub(r'\$|\\mathit|\{|\}|\^|_', '', prettyFeatureName(b))) for b in heat.index), default=10)
bracketIn, padIn = 1.4, 0.3
yLabelIn = maxFeatChars * (30 / 72.0) * 0.46
maxMutChars = max((len(re.sub(r'\$|\\mathit|\{|\}|\\Delta|\^|_', '', DISPLAY_NAMES[m])) for m in heat.columns), default=6)
mutLabelIn = maxMutChars * (30 / 72.0) * 0.55 * 0.71
left, bottom = (bracketIn + yLabelIn + padIn) / fw, (mutLabelIn + padIn) / fh
heatW, heatH = (nM * cellSize) / fw, (nF * cellSize) / fh

axH = fig.add_axes([left, bottom, heatW, heatH])
axC = fig.add_axes([left + heatW + 0.02, bottom, 0.02, heatH])
axBracket = fig.add_axes([0.35 / fw, bottom, bracketIn / fw, heatH])

# undefined measurements (vpsL forms no microcolonies) are drawn grey, not as z = 0
cmap = mpl.colormaps['RdBu_r'].copy(); cmap.set_bad('#b0b0b0')
im = axH.imshow(np.ma.masked_invalid(heat.values), cmap=cmap, vmin=-3, vmax=3,
                interpolation='nearest', aspect='auto')
axH.set_yticks(np.arange(nF)); axH.set_yticklabels([prettyFeatureName(b) for b in heat.index])
axH.set_xticks(np.arange(nM)); axH.set_xticklabels([DISPLAY_NAMES[m] for m in heat.columns], rotation=45, ha='right')
axH.tick_params(axis='x', pad=4)
axH.set_title('Peak Biofilm Biomass Frame', pad=24, fontweight='bold')
fig.colorbar(im, cax=axC).set_label('Z-score')

axBracket.set_xlim(0, 1); axBracket.set_ylim(nF - 0.5, -0.5); axBracket.axis('off')
groups = [featureGroup(b) for b in heat.index]


def rng(name):
    idx = [i for i, g in enumerate(groups) if g == name]
    return (min(idx), max(idx)) if idx else None


if (r := rng('colony')):
    drawBracket(axBracket, r[0], r[1], 'Colony Segmentation-\nDerived Features')
if (r := rng('haralick')):
    drawBracket(axBracket, r[0], r[1], 'Whole Image\nHaralick Features')

out = config.ensure(config.FIGURES) / '2C_trainingHeatmap_familyOrder'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
fig.savefig(str(out) + '.svg', bbox_inches='tight')
print(f'Saved: {out}.png  ({nF} features x {nM} mutants)')
