#!/usr/bin/env python3
"""fig3 copy — v2 reimaging UMAP with every replicate colored by peak biofilm biomass.

Copy of v2/reimaging/umap/umap_coloredByBiomass.py for the figure package, with the import path fixed
for fig3/ and a per-replicate COORDINATES CSV added (for publication: every plotted dot's umap position,
raw peak biomass, and peak biomass normalized to the WT mean, i.e. the value driving the color).

Peak biomass = max over frames; normalized to WT mean; plasma cmap, 2nd-98th pct range; WT drawn on top
with a white outline. nn=10, md=0.1.

Reads:  data/v2/reimaging/reimaging_umapEmbeddings.parquet + reimaging_collapsedWide.parquet
Writes: results/v2/reimaging/umaps/coloredByBiomassNormWT_plasma_nn{nn}_md{md}.{png,svg}
        results/v2/reimaging/umaps/coloredByBiomassNormWT_plasma_nn{nn}_md{md}_coords.csv
"""
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig3/ for figlib
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from figlib import config, plotting

plotting.setStyle(extra={'font.size': 44, 'axes.linewidth': 2})

nn, md, cmap = 10, 0.1, 'plasma'
embeddingsParquet = config.input('reimaging/umapEmbeddings.parquet')
wideParquet = config.input('reimaging/collapsedWide.parquet')
outRoot = config.ensure(config.FIGURES)
outBase = outRoot / f'coloredByBiomassNormWT_plasma_nn{nn}_md{md:.2f}'

emb = pd.read_parquet(embeddingsParquet)
emb = emb[(emb['n_neighbors'] == nn) & (emb['min_dist'] == md)].reset_index(drop=True)

wide = pd.read_parquet(wideParquet)
biomassCols = [c for c in wide.columns if re.match(r'^biomass_t\d+$', c)]
wide['peakBiomass'] = wide[biomassCols].max(axis=1)
emb = emb.merge(wide[['plateId', 'wellId', 'peakBiomass']], on=['plateId', 'wellId'], how='left')

wtMask = emb['mutant'].astype(str) == 'WT'
wtMeanPeak = emb.loc[wtMask, 'peakBiomass'].mean()
emb['peakBiomassNorm'] = emb['peakBiomass'] / wtMeanPeak
print(f'WT reps: {int(wtMask.sum())} | WT mean peak biomass: {wtMeanPeak:.4f}')

# CSV: every plotted dot's coordinates + peak biomass (raw + WT-normalized color value)
coordCols = [c for c in ['plateId', 'wellId', 'mutant', 'geneLocus', 'function'] if c in emb.columns]
emb[coordCols + ['umap1', 'umap2', 'peakBiomass', 'peakBiomassNorm']].to_csv(
    str(config.TABLES / outBase.name) + '_coords.csv', index=False)
print(f'Saved coordinates CSV: {outBase}_coords.csv ({len(emb)} replicates)')

vmin = emb['peakBiomassNorm'].quantile(0.02)
vmax = emb['peakBiomassNorm'].quantile(0.98)
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
pad = 0.5
xlim = (emb['umap1'].min() - pad, emb['umap1'].max() + pad)
ylim = (emb['umap2'].min() - pad, emb['umap2'].max() + pad)

fig = plt.figure(figsize=(18, 20))
gs = GridSpec(2, 1, height_ratios=[0.06, 1], hspace=0.08, figure=fig)
axCbar = fig.add_subplot(gs[0]); ax = fig.add_subplot(gs[1]); ax.set_box_aspect(1)

nonWt = emb[~wtMask]
ax.scatter(nonWt['umap1'], nonWt['umap2'], c=nonWt['peakBiomassNorm'], cmap=cmap, norm=norm,
           s=900, alpha=0.3, edgecolors='#8f8f8f', linewidth=1.0, zorder=2)
wtDf = emb[wtMask]
ax.scatter(wtDf['umap1'], wtDf['umap2'], c=wtDf['peakBiomassNorm'], cmap=cmap, norm=norm,
           s=1000, alpha=0.8, edgecolors='#FFFFFF', linewidth=3.5, zorder=5, label='WT')
ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
ax.legend(frameon=False, loc='upper right')

sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
cbar = fig.colorbar(sm, cax=axCbar, orientation='horizontal')
cbar.set_label('Normalized Peak \n Biofilm Biomass (a.u.)', labelpad=40, fontsize=56)
axCbar.xaxis.set_ticks_position('top'); axCbar.xaxis.set_label_position('top')
plt.savefig(str(outBase) + '.png', dpi=300, bbox_inches='tight')
plt.savefig(str(outBase) + '.svg', bbox_inches='tight')
plt.close(fig)
print(f'Saved: {outBase}.png')
