"""fig4 render: Panel 4C right - WT-normalized clean-deletion feature heatmap (horizontal).

Renders FROM the bundled matrix (conditions x features). Landscape orientation: rows = clean-deletion
conditions, columns = features, feature-group brackets along the top. Each cell = (condition - WT) /
reimaging-atlas sigma at the condition's peak-biomass frame (WT already dropped).

Reads:  data/cleanDel_vsWT_horizontal_matrix.csv
Writes: figures/4C_right_vsWTheatmap.{png,svg}
"""
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig4/ for figlib
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from figlib import config, plotting, CLEANDEL_DISPLAY

plotting.setStyle()
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


H = pd.read_csv(config.TABLES / 'cleanDel_vsWT_horizontal_matrix.csv', index_col=0)
condLabels = [CLEANDEL_DISPLAY.get(m, f'${m}$') for m in H.index]
nM, nF = H.shape

yLabelIn, featLabelIn, bracketTopIn, cbarIn, padIn = 1.8, 4.2, 1.6, 1.4, 0.3
fw = nF * cellSize + yLabelIn + cbarIn + padIn
fh = nM * cellSize + featLabelIn + bracketTopIn + padIn
fig = plt.figure(figsize=(fw, fh), dpi=200)
heatW, heatH = (nF * cellSize) / fw, (nM * cellSize) / fh
left, bottom = yLabelIn / fw, featLabelIn / fh
axH = fig.add_axes([left, bottom, heatW, heatH])
axBr = fig.add_axes([left, bottom + heatH + 0.01, heatW, bracketTopIn / fh]); axBr.set_xlim(-0.5, nF - 0.5); axBr.set_ylim(0, 1); axBr.axis('off')
axC = fig.add_axes([left + heatW + 0.02, bottom, 0.018, heatH])

im = axH.imshow(H.values, cmap='RdBu_r', vmin=-3, vmax=3, aspect='auto', interpolation='nearest')
axH.set_yticks(range(nM)); axH.set_yticklabels(condLabels, fontsize=30)
axH.set_xticks(range(nF)); axH.set_xticklabels([pretty(b) for b in H.columns], rotation=45, ha='right', fontsize=24)
cb = fig.colorbar(im, cax=axC); cb.set_label(r'$\sigma$ from WT', fontsize=24); cb.ax.tick_params(labelsize=20)

groups = [featGroup(b) for b in H.columns]


def rng(name):
    idx = [i for i, g in enumerate(groups) if g == name]
    return (min(idx), max(idx)) if idx else None


def bracket(x0, x1, label):
    axBr.plot([x0, x1], [0.25, 0.25], lw=3, color='black', clip_on=False)
    axBr.plot([x0, x0], [0.10, 0.25], lw=3, color='black', clip_on=False)
    axBr.plot([x1, x1], [0.10, 0.25], lw=3, color='black', clip_on=False)
    axBr.text((x0 + x1) / 2, 0.45, label, ha='center', va='bottom', fontsize=24)


if (r := rng('colony')):
    bracket(r[0], r[1], 'Colony Features')
if (r := rng('haralick')):
    bracket(r[0], r[1], 'Haralick Features')

out = config.ensure(config.FIGURES) / '4C_right_vsWTheatmap'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
fig.savefig(str(out) + '.svg', bbox_inches='tight')
print(f'Saved: {out}.png  ({nM} conditions x {nF} features)')
