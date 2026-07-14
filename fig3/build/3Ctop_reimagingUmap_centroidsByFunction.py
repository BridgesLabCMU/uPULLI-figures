#!/usr/bin/env python3
"""fig3 copy — v2 reimaging per-mutant centroids colored by functional annotation.

Copy of v2/reimaging/umap/umap_centroidsByFunction.py for the figure package, with the import path
fixed for fig3/ and a CENTROIDS CSV added (for publication: the plotted per-mutant centroid position,
replicate count, functional group, and color).

Reads:  data/v2/reimaging/reimaging_umapEmbeddings.parquet
Writes: results/v2/reimaging/umaps/centroidsByFunction_nn{nn}_md{md}.{png,svg}
        results/v2/reimaging/umaps/centroidsByFunction_nn{nn}_md{md}_centroids.csv
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig3/ for figlib
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from figlib import config, plotting

plotting.setStyle(extra={'font.size': 56, 'axes.linewidth': 2})

nn, md = 10, 0.1
embeddingsParquet = config.EMB
outRoot = config.ensure(config.FIGURES)
outBase = outRoot / f'centroidsByFunction_nn{nn}_md{md:.2f}'

highlightSets = plotting.HIGHLIGHT_SETS
functionColors = plotting.FUNCTION_COLORS
REPLICATE_COLOR, UNCLASSIFIED_COLOR, WT_COLOR = '#e3e0e0', '#939090', '#000000'

emb = pd.read_parquet(embeddingsParquet)
emb = emb[(emb['n_neighbors'] == nn) & (emb['min_dist'] == md)].reset_index(drop=True)
emb['geneLocus'] = emb['geneLocus'].fillna('')
print(f'Filtered nn={nn}, md={md}: {len(emb)} rows')

centroids = (emb.groupby('mutant', sort=False)
             .agg(umap1=('umap1', 'mean'), umap2=('umap2', 'mean'),
                  geneLocus=('geneLocus', 'first'), nReps=('umap1', 'count')).reset_index())


def centroidFuncName(row):
    if row['mutant'] == 'WT':
        return 'WT'
    for func, loci in highlightSets.items():
        if row['geneLocus'] in loci:
            return func
    return 'Unclassified'


def centroidColor(row):
    if row['mutant'] == 'WT':
        return WT_COLOR
    for func, loci in highlightSets.items():
        if row['geneLocus'] in loci:
            return functionColors[func]
    return UNCLASSIFIED_COLOR


centroids['functionalGroup'] = centroids.apply(centroidFuncName, axis=1)
centroids['color'] = centroids.apply(centroidColor, axis=1)

# ── publication CSV: the plotted per-mutant centroids ──
centroids[['mutant', 'geneLocus', 'nReps', 'umap1', 'umap2', 'functionalGroup', 'color']].to_csv(
    str(config.TABLES / outBase.name) + '_centroids.csv', index=False)
print(f'Saved centroids CSV: {outBase}_centroids.csv ({len(centroids)} mutants)')

wtC = centroids[centroids['mutant'] == 'WT']
hiC = centroids[(centroids['mutant'] != 'WT') & (centroids['color'] != UNCLASSIFIED_COLOR)]
unC = centroids[(centroids['mutant'] != 'WT') & (centroids['color'] == UNCLASSIFIED_COLOR)]
print(f'Centroids: {len(centroids)} ({len(hiC)} highlighted, {len(unC)} unclassified, {len(wtC)} WT)')

pad = 0.5
xlim = (emb['umap1'].min() - pad, emb['umap1'].max() + pad)
ylim = (emb['umap2'].min() - pad, emb['umap2'].max() + pad)
fig, ax = plt.subplots(figsize=(25, 25)); ax.set_box_aspect(1)
ax.scatter(emb['umap1'], emb['umap2'], c=REPLICATE_COLOR, s=250, alpha=0.1, edgecolors='none', zorder=1)
if not unC.empty:
    ax.scatter(unC['umap1'], unC['umap2'], c=UNCLASSIFIED_COLOR, s=900, alpha=0.7,
               edgecolors='black', linewidth=1.0, zorder=3, label='Unclassified')
for func, color in functionColors.items():
    sub = hiC[hiC['color'] == color]
    if not sub.empty:
        ax.scatter(sub['umap1'], sub['umap2'], c=color, s=1100, alpha=0.85,
                   edgecolors='black', linewidth=1.0, zorder=4, label=func)
if not wtC.empty:
    ax.scatter(wtC['umap1'], wtC['umap2'], c=WT_COLOR, s=1200, alpha=0.85,
               edgecolors='black', linewidth=1.0, zorder=5, label='WT')
ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
ax.legend(title='Gene Function', frameon=False, bbox_to_anchor=(1, 0.8), loc='upper left')
plt.tight_layout()
plt.savefig(str(outBase) + '.png', dpi=300, bbox_inches='tight')
plt.savefig(str(outBase) + '.svg', bbox_inches='tight')
plt.close(fig)
print(f'Saved: {outBase}.png')
