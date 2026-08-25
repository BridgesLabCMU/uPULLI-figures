#!/usr/bin/env python3
"""fig3 render — Figure 3 retired panel (was 3A): reimaging UMAP colored by peak biofilm biomass.

Renders FROM the bundled source-data table. Each replicate is colored by its WT-normalized peak
biomass (plasma, 2nd-98th pct); WT drawn on top with a white outline.

Reads:  data/coloredByBiomassNormWT_plasma_nn10_md0.10_coords.csv
Writes: figures/3A_coloredByBiomass.{png,svg}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig3/ for figlib
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd
from figlib import config, plotting

plotting.setStyle(extra={'font.size': 44, 'axes.linewidth': 2})
CMAP = 'plasma'

emb = pd.read_csv(config.TABLES / 'coloredByBiomassNormWT_plasma_nn10_md0.10_coords.csv')
wt = emb['mutant'].astype(str) == 'WT'
vmin, vmax = emb['peakBiomassNorm'].quantile(0.02), emb['peakBiomassNorm'].quantile(0.98)
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
pad = 0.5
xlim = (emb.umap1.min() - pad, emb.umap1.max() + pad)
ylim = (emb.umap2.min() - pad, emb.umap2.max() + pad)

fig = plt.figure(figsize=(18, 20))
gs = GridSpec(2, 1, height_ratios=[0.06, 1], hspace=0.08, figure=fig)
axCbar = fig.add_subplot(gs[0]); ax = fig.add_subplot(gs[1]); ax.set_box_aspect(1)

nonWt = emb[~wt]
ax.scatter(nonWt.umap1, nonWt.umap2, c=nonWt['peakBiomassNorm'], cmap=CMAP, norm=norm,
           s=900, alpha=0.3, edgecolors='#8f8f8f', linewidth=1.0, zorder=2)
wtDf = emb[wt]
ax.scatter(wtDf.umap1, wtDf.umap2, c=wtDf['peakBiomassNorm'], cmap=CMAP, norm=norm,
           s=1000, alpha=0.8, edgecolors='#FFFFFF', linewidth=3.5, zorder=5, label='WT')
ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
ax.legend(frameon=False, loc='upper right')

sm = mpl.cm.ScalarMappable(norm=norm, cmap=CMAP); sm.set_array([])
cbar = fig.colorbar(sm, cax=axCbar, orientation='horizontal')
cbar.set_label('Normalized Peak \n Biofilm Biomass (a.u.)', labelpad=40, fontsize=56)
axCbar.xaxis.set_ticks_position('top'); axCbar.xaxis.set_label_position('top')
out = config.ensure(config.FIGURES) / 'reimagingUmap_coloredByBiomass'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
fig.savefig(str(out) + '.svg', bbox_inches='tight')
print(f'Saved: {out}.png')
