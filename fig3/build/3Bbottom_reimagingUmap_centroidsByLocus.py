#!/usr/bin/env python3
"""fig3 copy v2 reimaging per-mutant centroids colored by gene-locus number.

Copy of v2/reimaging/umap/umap_centroidsByLocus.py for the figure package, with the import path fixed
for fig3/ and a CENTROIDS CSV added (for publication: the plotted per-mutant centroid position, the raw
and Chr-II-shifted locus number, and the chromosome flag driving the colormaps).

Plot 1: trimmed-viridis by locus number (Chr II shifted up by Chr I max). Plot 2: YlOrRd (Chr I) +
PuBuGn (Chr II) split. WT black on top. nn=10, md=0.1.

Reads:  data/v2/reimaging/reimaging_umapEmbeddings.parquet
Writes: results/v2/reimaging/umaps/centroidsByLocus_nn{nn}_md{md}.{png,svg} (+ _byChr, cbars)
        results/v2/reimaging/umaps/centroidsByLocus_nn{nn}_md{md}_centroids.csv
"""
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig3/ for figlib
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from figlib import config, plotting

plotting.setStyle(extra={'font.size': 56, 'axes.linewidth': 2})

nn, md = 10, 0.1
embeddingsParquet = config.input('reimaging/umapEmbeddings.parquet')
outRoot = config.ensure(config.FIGURES)
outBase = outRoot / f'centroidsByLocus_nn{nn}_md{md:.2f}'
outBase2 = outRoot / f'centroidsByLocus_byChr_nn{nn}_md{md:.2f}'

REPLICATE_COLOR, WT_COLOR, OUTLINE_COLOR = '#e3e0e0', '#000000', '#000000'
_colors = mpl.colormaps['viridis'](np.linspace(0.0, 1.0, 1500))
locusColormap = LinearSegmentedColormap.from_list('viridis_trimmed', _colors, N=1500)


def parseLocus(s):
    if not isinstance(s, str):
        return False, np.nan
    m = re.search(r'VC_([A-Z]?)(\d+)', s)
    if not m:
        return False, np.nan
    return (m.group(1) != ''), float(m.group(2))


emb = pd.read_parquet(embeddingsParquet)
emb = emb[(emb['n_neighbors'] == nn) & (emb['min_dist'] == md)].reset_index(drop=True)
emb['geneLocus'] = emb['geneLocus'].fillna('')

centroids = (emb.groupby('mutant', sort=False)
             .agg(umap1=('umap1', 'mean'), umap2=('umap2', 'mean'),
                  geneLocus=('geneLocus', 'first'), nReps=('umap1', 'count')).reset_index())
parsed = centroids['geneLocus'].apply(parseLocus)
centroids['isChrII'] = parsed.apply(lambda x: x[0])
centroids['locusNum'] = parsed.apply(lambda x: x[1])
centroids['locusNumRaw'] = centroids['locusNum'].copy()
chrIMax = centroids.loc[~centroids['isChrII'], 'locusNum'].max()
centroids.loc[centroids['isChrII'], 'locusNum'] += chrIMax

# CSV: plotted per-mutant centroids + locus/chromosome info
centroids[['mutant', 'geneLocus', 'nReps', 'umap1', 'umap2', 'isChrII', 'locusNumRaw', 'locusNum']].to_csv(
    str(config.TABLES / outBase.name) + '_centroids.csv', index=False)
print(f'Saved centroids CSV: {outBase}_centroids.csv ({len(centroids)} mutants)')

wtC = centroids[centroids['mutant'] == 'WT']
mutC = centroids[centroids['mutant'] != 'WT']
validMut = mutC[mutC['locusNum'].notna()]
invalidMut = mutC[mutC['locusNum'].isna()]
norm = mpl.colors.Normalize(vmin=validMut['locusNum'].min(), vmax=validMut['locusNum'].max())

pad = 0.5
xlim = (emb['umap1'].min() - pad, emb['umap1'].max() + pad)
ylim = (emb['umap2'].min() - pad, emb['umap2'].max() + pad)

fig, ax = plt.subplots(figsize=(15, 15)); ax.set_box_aspect(1)
ax.scatter(emb['umap1'], emb['umap2'], c=REPLICATE_COLOR, s=350, alpha=0.1, edgecolors='none', zorder=1)
if not validMut.empty:
    ax.scatter(validMut['umap1'], validMut['umap2'], c=validMut['locusNum'], cmap=locusColormap, norm=norm,
               s=1400, alpha=0.85, edgecolors=OUTLINE_COLOR, linewidth=1.0, zorder=3)
if not invalidMut.empty:
    ax.scatter(invalidMut['umap1'], invalidMut['umap2'], c='#939090', s=1400, alpha=0.7,
               edgecolors=OUTLINE_COLOR, linewidth=1.0, zorder=3)
if not wtC.empty:
    ax.scatter(wtC['umap1'], wtC['umap2'], c=WT_COLOR, s=1500, alpha=0.85,
               edgecolors=OUTLINE_COLOR, linewidth=1.0, zorder=5)
ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
plt.tight_layout()
plt.savefig(str(outBase) + '.png', dpi=300, bbox_inches='tight')
plt.savefig(str(outBase) + '.svg', bbox_inches='tight')
plt.close(fig)
print(f'Saved: {outBase}.png')


def saveCbar(cmap, norm, label, path):
    fig_cb, ax_cb = plt.subplots(figsize=(2, 10))
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    fig_cb.colorbar(sm, cax=ax_cb).set_label(label)
    fig_cb.savefig(str(path) + '.png', dpi=300, bbox_inches='tight')
    fig_cb.savefig(str(path) + '.svg', bbox_inches='tight')
    plt.close(fig_cb)


saveCbar(locusColormap, norm, 'Gene Locus Number', Path(str(outBase) + '_cbar'))

# plot 2: Chr I (YlOrRd) + Chr II (PuBuGn) 
chrIMut = validMut[~validMut['isChrII']]
chrIIMut = validMut[validMut['isChrII']]
chrI_cmap = LinearSegmentedColormap.from_list('YlOrRd_trimmed', mpl.colormaps['YlOrRd'](np.linspace(0.2, 1.0, 1500)), N=1500)
chrII_cmap = LinearSegmentedColormap.from_list('PuBuGn_trimmed', mpl.colormaps['PuBuGn'](np.linspace(0.4, 1.0, 1500)), N=1500)
chrI_norm = mpl.colors.Normalize(vmin=chrIMut['locusNumRaw'].min(), vmax=chrIMut['locusNumRaw'].max())
chrII_norm = mpl.colors.Normalize(vmin=chrIIMut['locusNumRaw'].min(), vmax=chrIIMut['locusNumRaw'].max())

fig2, ax2 = plt.subplots(figsize=(15, 15)); ax2.set_box_aspect(1)
ax2.scatter(emb['umap1'], emb['umap2'], c=REPLICATE_COLOR, s=350, alpha=0.1, edgecolors='none', zorder=1)
if not invalidMut.empty:
    ax2.scatter(invalidMut['umap1'], invalidMut['umap2'], c='#939090', s=1400, alpha=0.7, edgecolors=OUTLINE_COLOR, linewidth=1.0, zorder=3)
if not chrIMut.empty:
    ax2.scatter(chrIMut['umap1'], chrIMut['umap2'], c=chrIMut['locusNumRaw'], cmap=chrI_cmap, norm=chrI_norm, s=1400, alpha=0.85, edgecolors=OUTLINE_COLOR, linewidth=1.0, zorder=3)
if not chrIIMut.empty:
    ax2.scatter(chrIIMut['umap1'], chrIIMut['umap2'], c=chrIIMut['locusNumRaw'], cmap=chrII_cmap, norm=chrII_norm, s=1400, alpha=0.85, edgecolors=OUTLINE_COLOR, linewidth=1.0, zorder=4)
if not wtC.empty:
    ax2.scatter(wtC['umap1'], wtC['umap2'], c=WT_COLOR, s=1500, alpha=0.85, edgecolors=OUTLINE_COLOR, linewidth=1.5, zorder=5)
ax2.set_xlim(xlim); ax2.set_ylim(ylim); ax2.set_xlabel('UMAP 1'); ax2.set_ylabel('UMAP 2')
plt.tight_layout()
plt.savefig(str(outBase2) + '.png', dpi=300, bbox_inches='tight')
plt.savefig(str(outBase2) + '.svg', bbox_inches='tight')
plt.close(fig2)
print(f'Saved: {outBase2}.png')
saveCbar(chrI_cmap, chrI_norm, 'Chr I locus number', Path(str(outBase2) + '_cbar_chrI'))
saveCbar(chrII_cmap, chrII_norm, 'Chr II locus number', Path(str(outBase2) + '_cbar_chrII'))
