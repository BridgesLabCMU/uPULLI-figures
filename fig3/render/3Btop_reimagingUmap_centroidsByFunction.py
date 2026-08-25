#!/usr/bin/env python3
"""fig3 render — Figure 3 Panel B top inset: per-mutant UMAP centroids colored by functional group.

Renders FROM the bundled tables. Faint grey replicate dots (from the landscape coords) + one larger
centroid per mutant colored by functional group; unclassified dark grey; WT black on top.

Reads:  data/centroidsByFunction_nn10_md0.10_centroids.csv  (centroids + color)
        data/reimagingLandscape_nn10_md0.10_coords.csv       (replicate background)
Writes: figures/3Btop_centroidsByFunction.{png,svg}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig3/ for figlib
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from figlib import config, plotting

plotting.setStyle(extra={'font.size': 56, 'axes.linewidth': 2})
FC = plotting.FUNCTION_COLORS
REPLICATE_COLOR, UNCLASSIFIED_COLOR, WT_COLOR = '#e3e0e0', '#939090', '#000000'

reps = pd.read_csv(config.TABLES / 'reimagingLandscape_nn10_md0.10_coords.csv')
cen = pd.read_csv(config.TABLES / 'centroidsByFunction_nn10_md0.10_centroids.csv')

wtC = cen[cen['mutant'] == 'WT']
hiC = cen[(cen['mutant'] != 'WT') & (cen['color'] != UNCLASSIFIED_COLOR)]
unC = cen[(cen['mutant'] != 'WT') & (cen['color'] == UNCLASSIFIED_COLOR)]
pad = 0.5
xlim = (reps.umap1.min() - pad, reps.umap1.max() + pad)
ylim = (reps.umap2.min() - pad, reps.umap2.max() + pad)

fig, ax = plt.subplots(figsize=(25, 25)); ax.set_box_aspect(1)
ax.scatter(reps.umap1, reps.umap2, c=REPLICATE_COLOR, s=250, alpha=0.1, edgecolors='none', zorder=1)
if not unC.empty:
    ax.scatter(unC.umap1, unC.umap2, c=UNCLASSIFIED_COLOR, s=900, alpha=0.7,
               edgecolors='black', linewidth=1.0, zorder=3, label='Unclassified')
for func, color in FC.items():
    sub = hiC[hiC['color'] == color]
    if not sub.empty:
        ax.scatter(sub.umap1, sub.umap2, c=color, s=1100, alpha=0.85,
                   edgecolors='black', linewidth=1.0, zorder=4, label=func)
if not wtC.empty:
    ax.scatter(wtC.umap1, wtC.umap2, c=WT_COLOR, s=1200, alpha=0.85,
               edgecolors='black', linewidth=1.0, zorder=5, label='WT')
ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
ax.legend(title='Gene Function', frameon=False, bbox_to_anchor=(1, 0.8), loc='upper left')
plt.tight_layout()
out = config.ensure(config.FIGURES) / '3Btop_centroidsByFunction'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
fig.savefig(str(out) + '.svg', bbox_inches='tight')
print(f'Saved: {out}.png')
