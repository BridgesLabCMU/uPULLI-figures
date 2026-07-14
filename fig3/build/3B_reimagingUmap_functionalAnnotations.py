#!/usr/bin/env python3
"""Figure 3 — reimaging UMAP landscape for the paper (nn=10, min_dist=0.1), self-contained.

The narrowed-down single script behind Figure 3: it builds the reimaging UMAP manifold at ONLY
nn=10 / min_dist=0.1 (the version used in the paper), writes the per-replicate coordinates CSV, and
renders the square functional-landscape plot. Nothing else — no nn/md grid, no per-function panels,
no per-gene PDF, no saved scaler/models.

It reproduces the two things that make the paper figure:
  1. COORDINATES — built exactly as v2/reimaging/umap/buildReimagingUmaps.py builds the canonical
     manifold (drop-NA-mutant -> EXCLUDE_LOCI -> select_umap_feature_columns -> fillna(0) ->
     drop zero-variance -> growth filter -> min-replicate filter -> StandardScaler ->
     UMAP(nn=10, md=0.1, metric='euclidean', random_state=42, low_memory=True)), all via
     v2/common/features, so it matches data/v2/reimaging/reimaging_umapEmbeddings.parquet at (10, 0.1).
  2. PLOT — the "overview" panel of v2/reimaging/umap/reimagingLandscape.py: square axes, grey
     unclassified dots + functional-group highlights (low alpha) + WT in black, solid legend swatches.

Reads:  data/v2/reimaging/reimaging_collapsedWide.parquet
Writes: results/v2/reimaging/umaps/reimagingLandscape_nn10_md0.10.{png,svg}
        results/v2/reimaging/umaps/reimagingLandscape_nn10_md0.10_coords.csv
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig3/ for figlib
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.preprocessing import StandardScaler
import umap.umap_ as umap
from figlib import config, features

NN, MD, RANDOM_STATE = 10, 0.1, 42

# functional groups + colors for the paper figure (same as reimagingLandscape.py / common.plotting)
HIGHLIGHT_SETS = {
    'Motility': ['VC_2059','VC_2066','VC_2067','VC_2069','VC_2120','VC_2121','VC_2122','VC_2123',
                 'VC_2129','VC_2130','VC_2134','VC_2136','VC_2137','VC_2138','VC_2140','VC_2188','VC_2191',
                 'VC_2196','VC_2197','VC_2198','VC_2200','VC_2203','VC_2204','VC_2206','VC_2207','VC_2208'],
    'O-Antigen Biosynthesis': ['VC_0212','VC_0223','VC_0239','VC_0241','VC_0242','VC_0245','VC_0247','VC_0249',
                               'VC_0250','VC_0251','VC_0259','VC_0269'],
    'Polyamine Import': ['VC_1424','VC_1426','VC_1427','VC_1428'],
    'Biotin Biosynthesis': ['VC_1111','VC_1113','VC_1114','VC_1115'],
    'Pyruvate Flux': ['VC_2413','VC_0943'],
    'Vibriobactin Biosynthesis': ['VC_0771','VC_0772'],
}
FUNCTION_COLORS = {
    'Motility': '#ff0004', 'O-Antigen Biosynthesis': '#0096ff', 'Polyamine Import': '#14f7f0',
    'Biotin Biosynthesis': '#ff9f1c', 'Pyruvate Flux': '#39ff14', 'Vibriobactin Biosynthesis': '#ba17f6',
}
BACKGROUND_COLOR = '#d0d0d0'

outDir = config.ensure(config.FIGURES)
stem = f'reimagingLandscape_nn{NN}_md{MD:.2f}'

# build coordinates (nn=10, md=0.1 only)
wide = pd.read_parquet(config.WIDE)
wide = wide[wide['mutant'].notna()].reset_index(drop=True)
if 'geneLocus' in wide.columns:
    wide = wide[~wide['geneLocus'].isin(features.EXCLUDE_LOCI)].reset_index(drop=True)

featureCols = features.select_umap_feature_columns(wide)
X = wide[featureCols].copy().fillna(0)
zv = X.columns[X.nunique(dropna=False) <= 1].tolist()
if zv:
    X = X.drop(columns=zv)

gmask = features.growth_mask(X)
wide, X = wide[gmask].reset_index(drop=True), X[gmask].reset_index(drop=True)

low = features.low_replicate_genes(wide['mutant'])
if low:
    keep = ~wide['mutant'].isin(low)
    wide, X = wide[keep].reset_index(drop=True), X[keep].reset_index(drop=True)

print(f'Feature matrix: {X.shape}  (WT replicates: {int((wide["mutant"].astype(str) == "WT").sum())})')

XScaled = StandardScaler().fit_transform(X).astype(np.float32)
reducer = umap.UMAP(n_neighbors=NN, min_dist=MD, n_components=2, metric='euclidean',
                    random_state=RANDOM_STATE, low_memory=True)
embedding = reducer.fit_transform(XScaled)

metaCols = [c for c in ['plateId', 'wellId', 'mutant', 'geneLocus', 'function'] if c in wide.columns]
coords = wide[metaCols].copy()
coords['n_neighbors'], coords['min_dist'] = NN, MD
coords['umap1'], coords['umap2'] = embedding[:, 0], embedding[:, 1]
coordsCsv = config.TABLES / f'{stem}_coords.csv'
coords.to_csv(coordsCsv, index=False)
print(f'Saved coordinates: {coordsCsv}  ({len(coords)} replicates)')

# plot (square functional-landscape overview)
mpl.rcParams.update({'font.family': 'Gillius ADF', 'mathtext.fontset': 'stixsans', 'font.size': 28,
                     'axes.labelsize': 32, 'xtick.labelsize': 28, 'ytick.labelsize': 28,
                     'legend.fontsize': 32, 'legend.title_fontsize': 32, 'axes.linewidth': 2,
                     'savefig.dpi': 300})


def solidLegendHandle(label, faceColor, edgeColor='black'):
    """Fully-opaque legend marker so swatches read solid even though the plotted dots are low-alpha."""
    return Line2D([0], [0], marker='o', linestyle='none', markersize=18, markerfacecolor=faceColor,
                  markeredgecolor=edgeColor, markeredgewidth=0.8, alpha=1.0, label=label)


emb = coords
emb['geneLocus'] = emb['geneLocus'].fillna('') if 'geneLocus' in emb.columns else ''
wtMask = emb['mutant'].astype(str) == 'WT'
pad = 0.5
xlim = (emb['umap1'].min() - pad, emb['umap1'].max() + pad)
ylim = (emb['umap2'].min() - pad, emb['umap2'].max() + pad)
allHighlightLoci = set(g for genes in HIGHLIGHT_SETS.values() for g in genes)
unknownMask = ~emb['geneLocus'].isin(allHighlightLoci) & ~wtMask

fig, ax = plt.subplots(figsize=(15, 15))
ax.set_box_aspect(1)
ax.scatter(emb.loc[unknownMask, 'umap1'], emb.loc[unknownMask, 'umap2'],
           c=BACKGROUND_COLOR, s=400, alpha=0.15, edgecolors='black', linewidth=0.5, zorder=1)
legendHandles = []
for func, genes in HIGHLIGHT_SETS.items():
    subset = emb[emb['geneLocus'].isin(genes)]
    if not subset.empty:
        ax.scatter(subset['umap1'], subset['umap2'], c=FUNCTION_COLORS[func], s=400, alpha=0.3,
                   edgecolors='black', linewidth=0.5, zorder=3)
        legendHandles.append(solidLegendHandle(func, FUNCTION_COLORS[func]))
if wtMask.any():
    ax.scatter(emb.loc[wtMask, 'umap1'], emb.loc[wtMask, 'umap2'], c='black', s=500, alpha=0.6,
               edgecolors='black', linewidth=0.5, zorder=5)
    legendHandles.append(solidLegendHandle('WT', 'black'))
ax.set_xlim(xlim); ax.set_ylim(ylim)
ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
ax.legend(handles=legendHandles, title='Gene Function', frameon=False,
          bbox_to_anchor=(1, 0.8), loc='upper left')
plt.tight_layout()
fig.savefig(outDir / f'{stem}.png', dpi=300, bbox_inches='tight')
fig.savefig(outDir / f'{stem}.svg', bbox_inches='tight')
plt.close(fig)
print(f'Saved: {outDir / f"{stem}.png"}')
