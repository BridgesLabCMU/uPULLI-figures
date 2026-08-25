#!/usr/bin/env python3
"""figS7 render — Panel A: full reimaging dendrogram + heatmap on the top-15 embedding PCs.

Embedding-PC counterpart of Fig S6: all 158 mutants clustered (Ward) in the top-15 DINOv2-embedding
principal-component space, with a heatmap of those 15 PCs (z-scored across mutants). Renders FROM the
saved linkage + leaf order + PC heatmap matrix (build with build/S7C_embPC15.py).

Reads:  data/embPC15_linkage.npy, data/embPC15_cluster_order.csv, data/embPC15_heatmap_matrix.csv
Writes: figures/S7A_embPCDendrogramHeatmap.{png,svg}
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

plotting.setStyle(extra={'font.size': 44, 'axes.titlesize': 44, 'axes.labelsize': 44,
                         'xtick.labelsize': 30, 'ytick.labelsize': 44, 'axes.linewidth': 2})
WT = 'WT'
cellSize, titleFontsize = 0.7, 52


def formatLabel(mutant, annotation):
    m = str(mutant).replace('_', '')
    return 'WT' if annotation == 'WT' else f'${m}$'


Z = np.load(config.TABLES / 'embPC15_linkage.npy')
order = pd.read_csv(config.TABLES / 'embPC15_cluster_order.csv')
heat = pd.read_csv(config.TABLES / 'embPC15_heatmap_matrix.csv', index_col=0)

ordered = order['mutant'].astype(str).tolist()
heat = heat[ordered]                                          # PC (rows) x mutants (leaf order)
leafColors = order['color'].tolist()
labels = [formatLabel(m, a) for m, a in zip(order['mutant'], order['annotation'].astype(str))]
dend = dendrogram(Z, no_plot=True)
nF, nM = heat.shape

# ── layout (inch-based; mirrors Panel A) ──
scale = min(2.0, 60000.0 / (300.0 * max(nM * cellSize + 7, 1)))
fw = (nM * cellSize + 7) * scale
fh = (nF * cellSize + 6) * scale
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
axD.set_title('Top 15 Embedding Principal Components', pad=50, fontsize=titleFontsize, fontweight='bold')
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
im.set_rasterized(True)   # keep the (very wide) heatmap raster under the SVG 32767px image limit
axH.set_yticks(np.arange(nF)); axH.set_yticklabels(list(heat.index))
axH.set_xticks(np.arange(nM)); axH.set_xticklabels(labels, rotation=45, ha='right', fontsize=38)
axH.tick_params(axis='x', pad=4)
fig.colorbar(im, cax=axC).set_label(r'Z-score')

axFunc.imshow([range(nM)], aspect='auto', cmap=mpl.colors.ListedColormap(leafColors)).set_rasterized(True); axFunc.axis('off')

out = config.ensure(config.FIGURES) / 'S7A_embPCDendrogramHeatmap.png'
fig.savefig(out, dpi=300, bbox_inches='tight')
fig.savefig(str(out).replace('.png', '.svg'), dpi=200, bbox_inches='tight')  # dpi=200 -> raster < 32767px wide
plt.close(fig)
print(f'Saved: {out}  ({nF} PCs x {nM} leaves)')
