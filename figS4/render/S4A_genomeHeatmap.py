#!/usr/bin/env python3
"""figS4 render: Panel S4A - genome-wide transposon biomass screen (horizontal / transposed).

Wide-banner overview of ALL ~2,900 transposon mutants. Genome on the x-axis (ordered Chr I VC_#### ->
Chr II VC_A####, sparsely labeled with the gene name where one exists, otherwise the locus number, in a
large font), imaging time on the y-axis, color = biofilm biomass normalized to the WT peak mean
(RdYlBu_r 0..3). A tall strip along the bottom marks the reimaging-selected mutants (High Biofilm or
Dispersal Defect). The biomass colorbar is vertical, at the left of the plot.

Reads:  data/tn_biomass_matrix.csv, data/tn_locus_meta.csv, data/reimagingGeneNames.csv
Writes: figures/S4A_genomeHeatmap.{png,svg}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS4/ for figlib
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from figlib import config, plotting, PHENO_COLORS

plotting.setStyle(extra={'axes.linewidth': 1.2})

mat = pd.read_csv(config.TABLES / 'tn_biomass_matrix.csv', index_col=0)
meta = pd.read_csv(config.TABLES / 'tn_locus_meta.csv')
mat = mat.reindex(meta['geneLocus'])                 # align matrix rows to genome-ordered meta
frames = [int(c) for c in mat.columns]
nFrames, nLoci = len(frames), len(meta)
H = np.ma.masked_invalid(mat.values.T.astype(float))  # (nFrames, nLoci)
cmap = plt.cm.RdYlBu_r.copy(); cmap.set_bad('white')

# gene names (locus -> name; fall back to locus)
ri = pd.read_csv(config.TABLES / 'reimagingGeneNames.csv').dropna(subset=['geneLocus'])
nameMap = dict(zip(ri['geneLocus'], ri['geneName'].astype(str)))


def dispName(locus):
    n = nameMap.get(locus, '')
    return n.strip() if (isinstance(n, str) and n.strip() and n.strip() != str(locus)) else str(locus)


NOHIT = '#e8e8e8'
stripColors = [PHENO_COLORS.get(p, NOHIT) if hit else NOHIT for p, hit in zip(meta['phenotype'], meta['isHit'])]
strip = np.array([[mpl.colors.to_rgb(c) for c in stripColors]])   # (1, nLoci, 3)

# ── layout (inches) ──
cellW, cellH = 0.015, 0.34       # gene column width, time row height
heatW, heatH = nLoci * cellW, nFrames * cellH
STRIP_H = 0.75                   # reimaging indicator strip (tall / "much wider")
GAP = 0.08
LEFT, RIGHT, TOP = 2.6, 0.6, 2.0
LABELS = 2.8                     # rotated gene-name labels below the strip
FONT_GENE = 18
figW = LEFT + heatW + RIGHT
figH = LABELS + STRIP_H + GAP + heatH + TOP

xHM, wHM = LEFT / figW, heatW / figW
yHeat, hHeat = (LABELS + STRIP_H + GAP) / figH, heatH / figH
yStrip, hStrip = LABELS / figH, STRIP_H / figH

fig = plt.figure(figsize=(figW, figH))
ax = fig.add_axes([xHM, yHeat, wHM, hHeat])
axS = fig.add_axes([xHM, yStrip, wHM, hStrip])

im = ax.imshow(H, aspect='auto', cmap=cmap, vmin=0, vmax=3, interpolation='nearest', origin='lower',
               extent=[-0.5, nLoci - 0.5, -0.5, nFrames - 0.5])
ax.set_xticks([])
ax.set_yticks(np.arange(0, nFrames, 5)); ax.set_yticklabels(np.arange(0, nFrames, 5), fontsize=18)
ax.set_ylabel('Time (h)', fontsize=24)
ax.set_xlim(-0.5, nLoci - 0.5)

axS.imshow(strip, aspect='auto', interpolation='nearest', extent=[-0.5, nLoci - 0.5, 0, 1])
axS.set_yticks([]); axS.set_ylabel('Reimaged', fontsize=16, rotation=0, ha='right', va='center', labelpad=14)
axS.set_xlim(-0.5, nLoci - 0.5)
labelStep = max(1, nLoci // 70)
xpos = list(range(0, nLoci, labelStep))
axS.set_xticks(xpos)
axS.set_xticklabels([dispName(meta['geneLocus'].iloc[i]) for i in xpos], rotation=90, fontsize=FONT_GENE, va='top')
axS.tick_params(axis='x', length=3)

# vertical biomass colorbar at the far left
cbW, cbH = 0.34, heatH * 0.55
cbax = fig.add_axes([0.6 / figW, yHeat + (heatH * 0.22) / figH, cbW / figW, cbH / figH])
cbar = fig.colorbar(im, cax=cbax, orientation='vertical')
cbar.set_label('Biofilm Biomass (a.u.)', fontsize=18, labelpad=8)
cbar.ax.yaxis.set_ticks_position('left'); cbar.ax.yaxis.set_label_position('left'); cbar.ax.tick_params(labelsize=15)

# reimaging-strip legend + title (top)
fig.legend(handles=[Patch(facecolor=PHENO_COLORS['High Biofilm'], edgecolor='black', label='High Biofilm'),
                    Patch(facecolor=PHENO_COLORS['Dispersal Defect'], edgecolor='black', label='Dispersal Defect'),
                    Patch(facecolor=NOHIT, edgecolor='black', label='Not reimaged')],
           loc='center', bbox_to_anchor=(0.5, 1 - 0.9 / figH), ncol=3, frameon=False, fontsize=16,
           title='Reimaging-selected (bottom strip)', title_fontsize=17)
fig.text(0.5, 1 - 0.35 / figH, 'Genome-Wide Transposon Biomass Screen', ha='center', fontsize=28, fontweight='bold')

out = config.ensure(config.FIGURES) / 'S4A_genomeHeatmap'
fig.savefig(str(out) + '.png', dpi=200, bbox_inches='tight')
fig.savefig(str(out) + '.svg', bbox_inches='tight')
print(f'Saved: {out}.png  ({nLoci} loci x {nFrames} frames, {figW:.0f}x{figH:.0f} in)')
