"""fig5 render: Panel 5A right - biotin barcode heatmap (vsWT, vertical).

Renders FROM the bundled matrix (features x conditions). RdBu_r (-3..3), feature-group brackets on the
left; WT+DMSO is the zero-baseline column. Column headers in the CSV are already the display labels.

Reads:  data/compounds_biotinBarcode_vsWT_matrix.csv
Writes: figures/5A_right_biotinBarcode.{png,svg}
"""
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig5/ for figlib
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from figlib import config, plotting

plotting.setStyle(extra={'font.size': 30, 'axes.labelsize': 30, 'xtick.labelsize': 30,
                         'ytick.labelsize': 28, 'axes.linewidth': 2, 'mathtext.fontset': 'stixsans'})
cellSize = 0.7
canonHaralick = ['energy', 'contrast', 'correlation', 'variance', 'inverse_difference_moment', 'sum_average',
                 'sum_variance', 'sum_entropy', 'entropy', 'difference_variance', 'difference_entropy', 'imc1', 'imc2']
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


def featGroup(b):
    return 'biomass' if b == 'biomass' else 'haralick' if 'haralick' in b else 'entropy' if 'entropy' in b else 'colony'


def pretty(b):
    if b in prettyMap:
        return prettyMap[b]
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
    ax.text(0.38, (y0 + y1) / 2, label, rotation=90, va='center', ha='center', fontsize=28, clip_on=False)


heat = pd.read_csv(config.TABLES / 'compounds_biotinBarcode_vsWT_matrix.csv', index_col=0)
colLabels = list(heat.columns)   # already display labels
nF, nM = heat.shape
fw = nM * cellSize + 9
fh = nF * cellSize + 4
fig = plt.figure(figsize=(fw, fh), dpi=200)
maxFeatChars = max((len(re.sub(r'\$|\\mathit|\{|\}|\^|_', '', pretty(b))) for b in heat.index), default=10)
bracketIn, padIn = 1.4, 0.3
yLabelIn = maxFeatChars * (30 / 72.0) * 0.46
maxMutChars = max((len(re.sub(r'\$|\\mathit|\{|\}|\\Delta|\^|_|\\mu', '', c)) for c in colLabels), default=6)
mutLabelIn = maxMutChars * (30 / 72.0) * 0.55 * 0.71
left, bottom = (bracketIn + yLabelIn + padIn) / fw, (mutLabelIn + padIn) / fh
heatW, heatH = (nM * cellSize) / fw, (nF * cellSize) / fh
axH = fig.add_axes([left, bottom, heatW, heatH])
axC = fig.add_axes([left + heatW + 0.02, bottom, 0.02, heatH])
axBr = fig.add_axes([0.35 / fw, bottom, bracketIn / fw, heatH])

im = axH.imshow(heat.values, cmap='RdBu_r', vmin=-3, vmax=3, interpolation='nearest', aspect='auto')
axH.set_yticks(np.arange(nF)); axH.set_yticklabels([pretty(b) for b in heat.index])
axH.set_xticks(np.arange(nM)); axH.set_xticklabels(colLabels, rotation=45, ha='right')
axH.tick_params(axis='x', pad=4)
cbar = fig.colorbar(im, cax=axC); cbar.set_label(r'$\sigma$ from WT', fontsize=26); cbar.ax.tick_params(labelsize=22)

axBr.set_xlim(0, 1); axBr.set_ylim(nF - 0.5, -0.5); axBr.axis('off')
groups = [featGroup(b) for b in heat.index]


def rng(name):
    idx = [i for i, g in enumerate(groups) if g == name]
    return (min(idx), max(idx)) if idx else None


if (r := rng('colony')):
    drawBracket(axBr, r[0], r[1], 'Colony Segmentation-\nDerived Features')
if (r := rng('haralick')):
    drawBracket(axBr, r[0], r[1], 'Whole Image\nHaralick Features')

out = config.ensure(config.FIGURES) / '5A_right_biotinBarcode'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
fig.savefig(str(out) + '.svg', bbox_inches='tight')
print(f'Saved: {out}.png  ({nF} features x {nM} conditions)')
