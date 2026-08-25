#!/usr/bin/env python3
"""fig3 render — Figure 3 Panel B bottom inset: per-mutant UMAP centroids colored by genome position.

Renders FROM the bundled tables. Plot 1: continuous viridis by locus number (Chr II shifted above Chr I).
Plot 2: split Chr I (YlOrRd) + Chr II (PuBuGn). WT black on top. Faint grey replicate background.

Reads:  data/centroidsByLocus_nn10_md0.10_centroids.csv  (centroids + locus numbers + chromosome flag)
        data/reimagingLandscape_nn10_md0.10_coords.csv     (replicate background)
Writes: figures/3Bbottom_centroidsByLocus.{png,svg} (+ _byChr + colorbars)
"""
import sys
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
REPLICATE_COLOR, WT_COLOR, OUTLINE_COLOR = '#e3e0e0', '#000000', '#000000'
locusColormap = LinearSegmentedColormap.from_list('viridis_trimmed', mpl.colormaps['viridis'](np.linspace(0, 1, 1500)), N=1500)

reps = pd.read_csv(config.TABLES / 'reimagingLandscape_nn10_md0.10_coords.csv')
cen = pd.read_csv(config.TABLES / 'centroidsByLocus_nn10_md0.10_centroids.csv')

wtC = cen[cen['mutant'] == 'WT']
mutC = cen[cen['mutant'] != 'WT']
validMut = mutC[mutC['locusNum'].notna()]
invalidMut = mutC[mutC['locusNum'].isna()]
norm = mpl.colors.Normalize(vmin=validMut['locusNum'].min(), vmax=validMut['locusNum'].max())
pad = 0.5
xlim = (reps.umap1.min() - pad, reps.umap1.max() + pad)
ylim = (reps.umap2.min() - pad, reps.umap2.max() + pad)
outDir = config.ensure(config.FIGURES)
outBase = outDir / '3Bbottom_centroidsByLocus'
outBase2 = outDir / '3Bbottom_centroidsByLocus_byChr'

fig, ax = plt.subplots(figsize=(15, 15)); ax.set_box_aspect(1)
ax.scatter(reps.umap1, reps.umap2, c=REPLICATE_COLOR, s=350, alpha=0.1, edgecolors='none', zorder=1)
if not validMut.empty:
    ax.scatter(validMut.umap1, validMut.umap2, c=validMut['locusNum'], cmap=locusColormap, norm=norm,
               s=1400, alpha=0.85, edgecolors=OUTLINE_COLOR, linewidth=1.0, zorder=3)
if not invalidMut.empty:
    ax.scatter(invalidMut.umap1, invalidMut.umap2, c='#939090', s=1400, alpha=0.7,
               edgecolors=OUTLINE_COLOR, linewidth=1.0, zorder=3)
if not wtC.empty:
    ax.scatter(wtC.umap1, wtC.umap2, c=WT_COLOR, s=1500, alpha=0.85, edgecolors=OUTLINE_COLOR, linewidth=1.0, zorder=5)
ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
plt.tight_layout()
fig.savefig(str(outBase) + '.png', dpi=300, bbox_inches='tight')
fig.savefig(str(outBase) + '.svg', bbox_inches='tight')
plt.close(fig)
print(f'Saved: {outBase}.png')


def saveCbar(cmap, nrm, label, path):
    f, a = plt.subplots(figsize=(2, 10))
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=nrm); sm.set_array([])
    f.colorbar(sm, cax=a).set_label(label)
    f.savefig(str(path) + '.png', dpi=300, bbox_inches='tight')
    f.savefig(str(path) + '.svg', bbox_inches='tight')
    plt.close(f)


saveCbar(locusColormap, norm, 'Gene Locus Number', Path(str(outBase) + '_cbar'))

# plot 2: split Chr I / Chr II
chrIMut = validMut[~validMut['isChrII']]
chrIIMut = validMut[validMut['isChrII']]
chrI_cmap = LinearSegmentedColormap.from_list('YlOrRd_trimmed', mpl.colormaps['YlOrRd'](np.linspace(0.2, 1.0, 1500)), N=1500)
chrII_cmap = LinearSegmentedColormap.from_list('PuBuGn_trimmed', mpl.colormaps['PuBuGn'](np.linspace(0.4, 1.0, 1500)), N=1500)
chrI_norm = mpl.colors.Normalize(vmin=chrIMut['locusNumRaw'].min(), vmax=chrIMut['locusNumRaw'].max())
chrII_norm = mpl.colors.Normalize(vmin=chrIIMut['locusNumRaw'].min(), vmax=chrIIMut['locusNumRaw'].max())

fig2, ax2 = plt.subplots(figsize=(15, 15)); ax2.set_box_aspect(1)
ax2.scatter(reps.umap1, reps.umap2, c=REPLICATE_COLOR, s=350, alpha=0.1, edgecolors='none', zorder=1)
if not invalidMut.empty:
    ax2.scatter(invalidMut.umap1, invalidMut.umap2, c='#939090', s=1400, alpha=0.7, edgecolors=OUTLINE_COLOR, linewidth=1.0, zorder=3)
if not chrIMut.empty:
    ax2.scatter(chrIMut.umap1, chrIMut.umap2, c=chrIMut['locusNumRaw'], cmap=chrI_cmap, norm=chrI_norm, s=1400, alpha=0.85, edgecolors=OUTLINE_COLOR, linewidth=1.0, zorder=3)
if not chrIIMut.empty:
    ax2.scatter(chrIIMut.umap1, chrIIMut.umap2, c=chrIIMut['locusNumRaw'], cmap=chrII_cmap, norm=chrII_norm, s=1400, alpha=0.85, edgecolors=OUTLINE_COLOR, linewidth=1.0, zorder=4)
if not wtC.empty:
    ax2.scatter(wtC.umap1, wtC.umap2, c=WT_COLOR, s=1500, alpha=0.85, edgecolors=OUTLINE_COLOR, linewidth=1.5, zorder=5)
ax2.set_xlim(xlim); ax2.set_ylim(ylim); ax2.set_xlabel('UMAP 1'); ax2.set_ylabel('UMAP 2')
plt.tight_layout()
fig2.savefig(str(outBase2) + '.png', dpi=300, bbox_inches='tight')
fig2.savefig(str(outBase2) + '.svg', bbox_inches='tight')
plt.close(fig2)
print(f'Saved: {outBase2}.png')
saveCbar(chrI_cmap, chrI_norm, 'Chr I locus number', Path(str(outBase2) + '_cbar_chrI'))
saveCbar(chrII_cmap, chrII_norm, 'Chr II locus number', Path(str(outBase2) + '_cbar_chrII'))
