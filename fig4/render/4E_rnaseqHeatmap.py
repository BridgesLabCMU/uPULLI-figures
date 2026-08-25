"""fig4 render: Panel 4E - RNA-seq log2FC heatmap for the three clean-deletion mutants.

Matches the clean-deletion feature-heatmap conventions (Panel 4C right): RdBu_r, vmin/vmax = ±3,
inch-based square cells, operon-group brackets along the top, house style. Genes on the x-axis
(italic), mutants on the y-axis ordered ΔbioD (top), ΔpdhE2, ΔmanA (bottom).

Reads:  data/rnaseq_logFC_matrix.csv
Writes: figures/4E_rnaseqHeatmap.{png,svg}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig4/ for figlib
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from figlib import config, plotting

plotting.setStyle()
cellSize = 0.7
MUT_DISPLAY = {'bioD': r'$\Delta\mathit{bioD}$', 'pdhE2': r'$\Delta\mathit{pdhE2}$', 'manA': r'$\Delta\mathit{manA}$'}


def operon(g):
    for p in ('vps', 'rbm', 'msh', 'lap'):
        if g.startswith(p):
            return p.upper()
    return 'Other'


H = pd.read_csv(config.TABLES / 'rnaseq_logFC_matrix.csv', index_col=0)
mutLabels = [MUT_DISPLAY.get(m, f'${m}$') for m in H.index]
nM, nF = H.shape

yLabelIn, geneLabelIn, bracketTopIn, cbarIn, padIn = 2.0, 2.6, 1.4, 1.4, 0.3
TICK_FONT = 34   # same size for x (gene) and y (mutant) tick labels
fw = nF * cellSize + yLabelIn + cbarIn + padIn
fh = nM * cellSize + geneLabelIn + bracketTopIn + padIn
fig = plt.figure(figsize=(fw, fh), dpi=200)
heatW, heatH = (nF * cellSize) / fw, (nM * cellSize) / fh
left, bottom = yLabelIn / fw, geneLabelIn / fh
axH = fig.add_axes([left, bottom, heatW, heatH])
axBr = fig.add_axes([left, bottom + heatH + 0.01, heatW, bracketTopIn / fh]); axBr.set_xlim(-0.5, nF - 0.5); axBr.set_ylim(0, 1); axBr.axis('off')
axC = fig.add_axes([left + heatW + 0.02, bottom, 0.018, heatH])

im = axH.imshow(H.values, cmap='RdBu_r', vmin=-3, vmax=3, aspect='auto', interpolation='nearest')
axH.set_yticks(range(nM)); axH.set_yticklabels(mutLabels, fontsize=TICK_FONT)
axH.set_xticks(range(nF)); axH.set_xticklabels([f'${g}$' for g in H.columns], rotation=45, ha='right',
                                               rotation_mode='anchor', fontsize=TICK_FONT, fontstyle='italic')
axH.tick_params(axis='x', length=5, width=1, color='black'); axH.tick_params(axis='y', length=5, width=1, color='black')
cb = fig.colorbar(im, cax=axC); cb.set_label(r'log$_2$(FC)', fontsize=24); cb.ax.tick_params(labelsize=20)

groups = [operon(g) for g in H.columns]


def bracket(x0, x1, label):
    axBr.plot([x0, x1], [0.25, 0.25], lw=3, color='black', clip_on=False)
    axBr.plot([x0, x0], [0.10, 0.25], lw=3, color='black', clip_on=False)
    axBr.plot([x1, x1], [0.10, 0.25], lw=3, color='black', clip_on=False)
    axBr.text((x0 + x1) / 2, 0.45, label, ha='center', va='bottom', fontsize=24, fontstyle='italic')


start = 0
for i in range(1, nF + 1):
    if i == nF or groups[i] != groups[start]:
        if groups[start] != 'Other':
            bracket(start, i - 1, groups[start].lower())
        start = i

out = config.ensure(config.FIGURES) / '4E_rnaseqHeatmap'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
fig.savefig(str(out) + '.svg', bbox_inches='tight')
print(f'Saved: {out}.png  ({nM} mutants x {nF} genes)')
