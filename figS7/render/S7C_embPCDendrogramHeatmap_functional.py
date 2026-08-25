#!/usr/bin/env python3
"""figS7 render — Panel C: top-15 embedding-PC dendrogram + heatmap, FUNCTIONAL-ANNOTATION subset.

The functional-subset counterpart of Panel A (the embedding-PC dendrogram): only the functionally-
annotated reimaging mutants (the six highlight pathways + WT), clustered in the top-15 DINOv2-embedding
PC space, with a z-scored heatmap of those 15 PCs. Renders FROM the saved linkage + leaf order + matrix
(build with build/S7C_embPC15_functional.py).

Reads:  data/embPC15func_linkage.npy, data/embPC15func_cluster_order.csv, data/embPC15func_heatmap_matrix.csv
Writes: figures/S7C_embPCDendrogramHeatmap_functional.{png,svg}
"""
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS7/ for figlib
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram
from figlib import config, plotting

plotting.setStyle(extra={'font.size': 44, 'axes.titlesize': 44, 'axes.labelsize': 44,
                         'xtick.labelsize': 30, 'ytick.labelsize': 44, 'axes.linewidth': 2})
WT = 'WT'
cellSize, titleFontsize = 0.7, 52
FC = plotting.FUNCTION_COLORS


def formatLabel(mutant, annotation):
    m = str(mutant).replace('_', '')
    return 'WT' if annotation == 'WT' else f'${m}$'


Z = np.load(config.TABLES / 'embPC15func_linkage.npy')
order = pd.read_csv(config.TABLES / 'embPC15func_cluster_order.csv')
heat = pd.read_csv(config.TABLES / 'embPC15func_heatmap_matrix.csv', index_col=0)

ordered = order['mutant'].astype(str).tolist()
heat = heat[ordered]
leafColors = order['color'].tolist()
labels = [formatLabel(m, a) for m, a in zip(order['mutant'], order['annotation'].astype(str))]
dend = dendrogram(Z, no_plot=True)
nF, nM = heat.shape

scale = min(2.0, 60000.0 / (300.0 * max(nM * cellSize + 7, 1)))
fw = (nM * cellSize + 7) * scale
fh = (nF * cellSize + 7.5) * scale
fig = plt.figure(figsize=(fw, fh), dpi=150)

_vis = lambda l: re.sub(r'\$|_', '', str(l))
maxMutChars = min(max((len(_vis(l)) for l in labels), default=6), 14)
yLabelIn, padIn, stripIn = 1.4, 0.3, 0.5
mutLabelIn = maxMutChars * (44 / 72.0) * 0.55 * 0.71
leftIn = yLabelIn + padIn
bottomIn = mutLabelIn + stripIn + 1.5 * padIn
left, bottom = leftIn / fw, bottomIn / fh
heatW, heatH = (nM * cellSize) / fw, (nF * cellSize) / fh

axH = fig.add_axes([left, bottom, heatW, heatH])
axD = fig.add_axes([left, bottom + heatH + 0.003, heatW, 0.07])
axD.set_title('Top 15 Embedding PCs — Functional Subset', pad=50, fontsize=titleFontsize, fontweight='bold')
axC = fig.add_axes([left + heatW + 0.015, bottom, 0.012, heatH])
stripBottom = max((bottomIn - mutLabelIn - padIn - stripIn) / fh, 0.004)
axFunc = fig.add_axes([left, stripBottom, heatW, stripIn / fh])

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
axH.set_yticks(np.arange(nF)); axH.set_yticklabels(list(heat.index))
axH.set_xticks(np.arange(nM)); axH.set_xticklabels(labels, rotation=45, ha='right', fontsize=38)
axH.tick_params(axis='x', pad=4)
fig.colorbar(im, cax=axC).set_label(r'Z-score')

axFunc.imshow([range(nM)], aspect='auto', cmap=mpl.colors.ListedColormap(leafColors)); axFunc.axis('off')

# functional-annotation legend (groups present, in canonical order)
present = [g for g in FC if g in set(order['annotation'])]
handles = [Patch(facecolor=FC[g], edgecolor='black', label=g) for g in present]
if 'WT' in set(order['annotation']):
    handles.append(Patch(facecolor='#000000', edgecolor='black', label='WT'))
fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(left, 0.99), frameon=False,
           ncol=3, fontsize=40, title='Gene Function', title_fontsize=44)

out = config.ensure(config.FIGURES) / 'S7C_embPCDendrogramHeatmap_functional.png'
fig.savefig(out, dpi=300, bbox_inches='tight')
fig.savefig(str(out).replace('.png', '.svg'), bbox_inches='tight')
plt.close(fig)
print(f'Saved: {out}  ({nF} PCs x {nM} leaves)')
