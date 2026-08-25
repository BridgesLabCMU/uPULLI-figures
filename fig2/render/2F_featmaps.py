#!/usr/bin/env python3
"""fig2 render — Panel 2F: single-feature mutant x frame heatmaps (four features).

Renders FROM the bundled per-feature matrices + meta. plasma colormap, masked (low-biomass colony)
cells black, 2-line title: feature name (+ unit) and the single-feature RF balanced accuracy.

Reads:  data/featmaps_meta.csv + data/featmap_<feat>.csv (x4)
Writes: figures/2F_<feat>.{png,svg}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig2/ for figlib
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from figlib import config, plotting, STRAIN_ORDER, DISPLAY_NAMES

plotting.setStyle(extra={'font.size': 36, 'axes.titlesize': 40, 'axes.labelsize': 36,
                         'xtick.labelsize': 30, 'ytick.labelsize': 30, 'axes.linewidth': 2})

meta = pd.read_csv(config.TABLES / 'featmaps_meta.csv')
outDir = config.ensure(config.FIGURES)

# per-family colormap (consistent with Fig S3): biomass = red-blue diverging, whole-image = viridis, colony = plasma
CMAP_BY_GROUP = {'biomass': 'coolwarm', 'whole': 'viridis', 'colony': 'plasma'}

for _, r in meta.iterrows():
    feat, label, unit = r['feature'], r['label'], ('' if pd.isna(r['unit']) else str(r['unit']))
    mat = pd.read_csv(config.TABLES / f'featmap_{feat}.csv', index_col=0)
    mat = mat.loc[[m for m in STRAIN_ORDER if m in mat.index]]
    frames = [int(c) for c in mat.columns]
    heatmap = np.ma.masked_invalid(mat.values.astype(float))
    vmin, vmax = float(np.nanmin(heatmap)), float(np.nanmax(heatmap))
    cmap = plt.get_cmap(CMAP_BY_GROUP.get(r['group'], 'plasma')).copy(); cmap.set_bad('black')

    fig = plt.figure(figsize=(9, 8))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 14], left=0.25, right=0.95, bottom=0.10, top=0.86, hspace=0.15)
    cax = fig.add_subplot(gs[0]); ax = fig.add_subplot(gs[1])
    im = ax.imshow(heatmap, aspect='auto', cmap=cmap, interpolation='nearest', vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
    cbar.set_label(f'({unit})' if unit else '(unitless)', fontsize=24, labelpad=8)
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.tick_params(labelsize=20)
    unitTxt = f' ({unit})' if unit else ''
    rfTxt = '' if pd.isna(r['rfAccuracy']) else f'\nClassification Accuracy = {float(r["rfAccuracy"]):.3f}'
    cax.set_title(f'{label}{unitTxt}{rfTxt}', fontsize=26, pad=14)

    tickFrames = [t for t in frames if t >= 10 and t % 5 == 0]
    ax.set_xticks([frames.index(t) for t in tickFrames]); ax.set_xticklabels(tickFrames)
    ax.set_xlabel('Time (h)', fontsize=26)
    ax.set_yticks(range(len(mat.index))); ax.set_yticklabels([DISPLAY_NAMES[m] for m in mat.index])

    out = outDir / f'2F_{feat}'
    fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
    fig.savefig(str(out) + '.svg', bbox_inches='tight')
    plt.close(fig)
    print(f'Saved: {out}.png')
