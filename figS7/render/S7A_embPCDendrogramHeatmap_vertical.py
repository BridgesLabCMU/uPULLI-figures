#!/usr/bin/env python3
"""figS7 render — Panel A (VERTICAL): full-atlas embedding-PC dendrogram + heatmap, transposed.

Portrait version of S7A: mutants on the y-axis (leaf order, dendrogram on the left), the top-15 DINOv2
embedding PCs on the x-axis (labels rotated 30 deg), functional-annotation color strip between the
dendrogram and the heatmap. Large tick fonts.

Reads:  data/embPC15_linkage.npy, data/embPC15_cluster_order.csv, data/embPC15_heatmap_matrix.csv
Writes: figures/S7A_embPCDendrogramHeatmap_vertical.{png,svg}
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
from scipy.cluster.hierarchy import dendrogram
from figlib import config, plotting

plotting.setStyle(extra={'font.size': 40, 'axes.linewidth': 2})
WT = 'WT'
cellSize = 0.7
FONT_MUT, FONT_PC = 38, 38


def formatLabel(mutant, annotation):
    return 'WT' if annotation == 'WT' else f'${str(mutant).replace("_", "")}$'


Z = np.load(config.TABLES / 'embPC15_linkage.npy')
order = pd.read_csv(config.TABLES / 'embPC15_cluster_order.csv')
heatMat = pd.read_csv(config.TABLES / 'embPC15_heatmap_matrix.csv', index_col=0)   # PC (rows) x mutants

ordered = order['mutant'].astype(str).tolist()
leafColors = order['color'].tolist()
mutLabels = [formatLabel(m, a) for m, a in zip(order['mutant'], order['annotation'].astype(str))]
heat = heatMat[ordered].T                                    # mutants (rows, leaf order) x PCs (cols)
pcNames = list(heat.columns)
nM, nF = heat.shape

# ── layout (inches) ──
padIn, dendIn, stripIn, titleIn, cbarIn = 0.4, 6.0, 0.6, 1.4, 1.8
mutLabelIn = 3.0
pcLabelIn = max(len(p) for p in pcNames) * (FONT_PC / 72.0) * 0.55 + 0.5
heatW, heatH = nF * cellSize, nM * cellSize
fwPre = padIn + dendIn + stripIn + heatW + mutLabelIn + cbarIn + padIn
fhPre = padIn + pcLabelIn + heatH + titleIn + padIn
scale = min(1.0, 60000.0 / (300.0 * fhPre))
fig = plt.figure(figsize=(fwPre * scale, fhPre * scale), dpi=150)

xDend = padIn / fwPre
xStrip = (padIn + dendIn) / fwPre
xHeat = (padIn + dendIn + stripIn) / fwPre
yHeat = (padIn + pcLabelIn) / fhPre
wDend, wStrip, wHeat = dendIn / fwPre, stripIn / fwPre, heatW / fwPre
hHeat = heatH / fhPre

axDend = fig.add_axes([xDend, yHeat, wDend, hHeat])
axStrip = fig.add_axes([xStrip, yHeat, wStrip, hHeat])
axH = fig.add_axes([xHeat, yHeat, wHeat, hHeat])
axC = fig.add_axes([(padIn + dendIn + stripIn + heatW + mutLabelIn + 0.4) / fwPre, yHeat + hHeat * 0.25, 0.22 / fwPre, hHeat * 0.5])

im = axH.imshow(heat.values, cmap='RdBu_r', vmin=-3, vmax=3, interpolation='nearest', aspect='auto')
im.set_rasterized(True)
axH.set_xticks(np.arange(nF)); axH.set_xticklabels(pcNames, rotation=30, ha='right', rotation_mode='anchor', fontsize=FONT_PC)
axH.yaxis.set_ticks_position('right'); axH.yaxis.set_label_position('right')
axH.set_yticks(np.arange(nM)); axH.set_yticklabels(mutLabels, fontsize=FONT_MUT)
cb = fig.colorbar(im, cax=axC); cb.set_label('Z-score', fontsize=34); cb.ax.tick_params(labelsize=28)

axStrip.imshow([[i] for i in range(nM)], aspect='auto', cmap=mpl.colors.ListedColormap(leafColors)).set_rasterized(True)
axStrip.axis('off')

dend = dendrogram(Z, no_plot=True)
maxD = max(max(d) for d in dend['dcoord'])
for xs, ys in zip(dend['icoord'], dend['dcoord']):
    rows = [(v - 5.0) / 10.0 for v in xs]
    axDend.plot(ys, rows, color='black', linewidth=2)
for i in range(nM):
    axDend.plot(0, i, 'o', color=leafColors[i], markersize=(11 if ordered[i] == WT else 9), clip_on=False)
axDend.set_xlim(maxD * 1.02, -0.02 * maxD); axDend.set_ylim(nM - 0.5, -0.5); axDend.axis('off')

fig.text((xHeat + wHeat / 2), yHeat + hHeat + titleIn * 0.5 / fhPre, 'Top 15 Embedding Principal Components',
         ha='center', va='center', fontsize=46, fontweight='bold')

out = config.ensure(config.FIGURES) / 'S7A_embPCDendrogramHeatmap_vertical'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
fig.savefig(str(out) + '.svg', dpi=200, bbox_inches='tight')
plt.close(fig)
print(f'Saved: {out}.png  ({nM} mutants x {nF} PCs)')
