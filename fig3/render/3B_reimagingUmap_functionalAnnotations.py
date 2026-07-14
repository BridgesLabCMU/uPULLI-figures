#!/usr/bin/env python3
"""fig3 render — Figure 3 Panel B: reimaging UMAP functional-annotation landscape.

Renders the panel FROM the bundled source-data table (no recompute). Grey = unclassified, six
functional-group colors, black = WT; solid legend swatches.

Reads:  data/reimagingLandscape_nn10_md0.10_coords.csv
Writes: figures/3B_landscape.{png,svg}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig3/ for figlib
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
from figlib import config, plotting

plotting.setStyle(extra={'font.size': 28, 'axes.labelsize': 32, 'xtick.labelsize': 28, 'ytick.labelsize': 28,
                         'legend.fontsize': 32, 'legend.title_fontsize': 32, 'axes.linewidth': 2})
HS, FC, BG = plotting.HIGHLIGHT_SETS, plotting.FUNCTION_COLORS, plotting.BACKGROUND_COLOR

emb = pd.read_csv(config.TABLES / 'reimagingLandscape_nn10_md0.10_coords.csv')
emb['geneLocus'] = emb['geneLocus'].fillna('')
wt = emb['mutant'].astype(str) == 'WT'
pad = 0.5
xlim = (emb.umap1.min() - pad, emb.umap1.max() + pad)
ylim = (emb.umap2.min() - pad, emb.umap2.max() + pad)
allHi = set(g for gs in HS.values() for g in gs)
unknown = ~emb['geneLocus'].isin(allHi) & ~wt


def solid(label, face):
    return Line2D([0], [0], marker='o', linestyle='none', markersize=18, markerfacecolor=face,
                  markeredgecolor='black', markeredgewidth=0.8, alpha=1.0, label=label)


fig, ax = plt.subplots(figsize=(15, 15))
ax.set_box_aspect(1)
ax.scatter(emb.loc[unknown, 'umap1'], emb.loc[unknown, 'umap2'], c=BG, s=400, alpha=0.15,
           edgecolors='black', linewidth=0.5, zorder=1)
handles = []
for func, genes in HS.items():
    sub = emb[emb['geneLocus'].isin(genes)]
    if not sub.empty:
        ax.scatter(sub.umap1, sub.umap2, c=FC[func], s=400, alpha=0.3, edgecolors='black', linewidth=0.5, zorder=3)
        handles.append(solid(func, FC[func]))
if wt.any():
    ax.scatter(emb.loc[wt, 'umap1'], emb.loc[wt, 'umap2'], c='black', s=500, alpha=0.6,
               edgecolors='black', linewidth=0.5, zorder=5)
    handles.append(solid('WT', 'black'))
ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
ax.legend(handles=handles, title='Gene Function', frameon=False, bbox_to_anchor=(1, 0.8), loc='upper left')
plt.tight_layout()
out = config.ensure(config.FIGURES) / '3B_landscape'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
fig.savefig(str(out) + '.svg', bbox_inches='tight')
print(f'Saved: {out}.png')
