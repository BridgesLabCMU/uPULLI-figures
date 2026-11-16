#!/usr/bin/env python3
"""fig3 render - the reimaging UMAP landscape in three dimensions.

Renders FROM the bundled source-data table (no recompute), in the two colorings the atlas is drawn
with elsewhere:

  --color function   grey unclassified + the six functional groups + WT in black   (the Fig 3B look)
  --color biomass    every replicate by WT-normalized peak biofilm biomass, plasma (the retired 3A look)

This is an alternative VIEW, not a re-labelled Figure 3B: `reimagingUmap3d.py` in build/ is a separate
UMAP fit at n_components=3, and UMAP does not nest, so these umap1/umap2 are not the published 2-D
coordinates. Not Panel 3D either - that is the dendrogram + heatmap.

A 3-D scatter is only as good as its viewpoint, so `--elev/--azim` set it and `--grid` renders four
azimuths in one figure, which is the readable option on paper. Depth cueing is by draw order between
groups (matplotlib does not z-sort points across collections), so WT and the highlighted groups always
sit in front of the grey cloud.

Reads:  data/reimagingLandscape3d_nn10_md0.10_coords.csv
Writes: figures/reimagingUmap3d_<color>_e<elev>a<azim>.{png,svg}   (or ..._grid.{png,svg})
"""
import sys
import argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig3/ for figlib
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers the 3d projection)
import pandas as pd
from figlib import config, plotting

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--color', choices=['function', 'biomass'], default='function',
                help='functional-group highlights (default) or WT-normalized peak biomass')
ap.add_argument('--elev', type=float, default=22, help='elevation angle in degrees (default 22)')
ap.add_argument('--azim', type=float, default=-60, help='azimuth angle in degrees (default -60)')
ap.add_argument('--grid', action='store_true', help='four azimuths in one figure instead of one view')
ap.add_argument('--cmap', default='plasma', help='colormap for --color biomass (default plasma)')
args = ap.parse_args()

plotting.setStyle(extra={'font.size': 20, 'axes.labelsize': 24, 'xtick.labelsize': 16,
                         'ytick.labelsize': 16, 'legend.fontsize': 22, 'legend.title_fontsize': 22,
                         'axes.linewidth': 1.5})
HS, FC, BG = plotting.HIGHLIGHT_SETS, plotting.FUNCTION_COLORS, plotting.BACKGROUND_COLOR

emb = pd.read_csv(config.TABLES / 'reimagingLandscape3d_nn10_md0.10_coords.csv')
emb['geneLocus'] = emb['geneLocus'].fillna('') if 'geneLocus' in emb.columns else ''
wt = emb['mutant'].astype(str) == 'WT'
allHi = set(g for gs in HS.values() for g in gs)
unknown = ~emb['geneLocus'].isin(allHi) & ~wt

pad = 0.5
lims = {a: (emb[a].min() - pad, emb[a].max() + pad) for a in ('umap1', 'umap2', 'umap3')}
norm = None
if args.color == 'biomass':
    norm = mpl.colors.Normalize(vmin=emb['peakBiomassNorm'].quantile(0.02),
                                vmax=emb['peakBiomassNorm'].quantile(0.98))


def solid(label, face, edge='black'):
    return Line2D([0], [0], marker='o', linestyle='none', markersize=13, markerfacecolor=face,
                  markeredgecolor=edge, markeredgewidth=0.8, alpha=1.0, label=label)


def drawAxes(ax, elev, azim):
    """One 3-D view. Returns the legend handles (functional coloring) or the mappable (biomass)."""
    handles, sm = [], None
    if args.color == 'function':
        ax.scatter(emb.loc[unknown, 'umap1'], emb.loc[unknown, 'umap2'], emb.loc[unknown, 'umap3'],
                   c=BG, s=55, alpha=0.13, edgecolors='black', linewidth=0.3, depthshade=False)
        for func, genes in HS.items():
            sub = emb[emb['geneLocus'].isin(genes)]
            if sub.empty:
                continue
            ax.scatter(sub['umap1'], sub['umap2'], sub['umap3'], c=FC[func], s=75, alpha=0.75,
                       edgecolors='black', linewidth=0.4, depthshade=False)
            handles.append(solid(plotting.functionLabel(func), FC[func]))
        if wt.any():
            ax.scatter(emb.loc[wt, 'umap1'], emb.loc[wt, 'umap2'], emb.loc[wt, 'umap3'], c='black',
                       s=95, alpha=0.9, edgecolors='black', linewidth=0.4, depthshade=False)
            handles.append(solid('WT', 'black'))
    else:
        nonWt = emb[~wt]
        sm = ax.scatter(nonWt['umap1'], nonWt['umap2'], nonWt['umap3'], c=nonWt['peakBiomassNorm'],
                        cmap=args.cmap, norm=norm, s=65, alpha=0.5, edgecolors='#8f8f8f',
                        linewidth=0.3, depthshade=False)
        wtDf = emb[wt]
        ax.scatter(wtDf['umap1'], wtDf['umap2'], wtDf['umap3'], c=wtDf['peakBiomassNorm'],
                   cmap=args.cmap, norm=norm, s=95, alpha=0.95, edgecolors='#FFFFFF', linewidth=1.6,
                   depthshade=False, label='WT')
        handles.append(solid('WT', 'none', edge='#FFFFFF'))

    ax.set_xlim(*lims['umap1']); ax.set_ylim(*lims['umap2']); ax.set_zlim(*lims['umap3'])
    ax.set_xlabel('UMAP 1', labelpad=14)
    ax.set_ylabel('UMAP 2', labelpad=14)
    ax.set_zlabel('UMAP 3', labelpad=14)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=elev, azim=azim)
    # keep the box light so the dots carry the figure
    for pane in (ax.xaxis, ax.yaxis, ax.zaxis):
        pane.pane.set_facecolor('white')
        pane.pane.set_edgecolor('#cccccc')
        pane._axinfo['grid'].update(color='#e8e8e8', linewidth=0.8)
    return handles, sm


if args.grid:
    azims = [args.azim + d for d in (0, 90, 180, 270)]
    fig = plt.figure(figsize=(22, 22))
    for i, az in enumerate(azims, start=1):
        ax = fig.add_subplot(2, 2, i, projection='3d')
        handles, sm = drawAxes(ax, args.elev, az)
        ax.set_title(f'azim = {az:g}°', pad=4)
    stem = f'reimagingUmap3d_{args.color}_e{args.elev:g}_grid'
else:
    fig = plt.figure(figsize=(14, 13))
    ax = fig.add_subplot(111, projection='3d')
    handles, sm = drawAxes(ax, args.elev, args.azim)
    stem = f'reimagingUmap3d_{args.color}_e{args.elev:g}a{args.azim:g}'

if args.color == 'function':
    fig.legend(handles=handles, title='Gene Function', frameon=False,
               loc='center left', bbox_to_anchor=(0.86, 0.5))
    fig.subplots_adjust(right=0.85)
else:
    # 3-D axes reserve wide internal margins, so the bar is shrunk and pulled in tight rather than
    # taking the default fraction of a box that is mostly empty.
    mappable = mpl.cm.ScalarMappable(norm=norm, cmap=args.cmap); mappable.set_array([])
    cbar = fig.colorbar(mappable, ax=fig.axes, orientation='horizontal', location='top',
                        fraction=0.022, shrink=0.55, aspect=32, pad=-0.02)
    cbar.set_label('Normalized Peak Biofilm Biomass (a.u.)', labelpad=12)

out = config.ensure(config.FIGURES) / stem
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
fig.savefig(str(out) + '.svg', bbox_inches='tight')
plt.close(fig)
print(f'Saved: {out}.png / .svg  ({len(emb)} replicates)')
