#!/usr/bin/env python3
"""fig2 render — Panel 2B: per-replicate biofilm-biomass traces heatmap.

Renders FROM the bundled source table. One row per replicate, grouped by mutant (STRAIN_ORDER); within a
mutant, rows are in (plateId, wellId) order. Values are WT-peak-normalized; RdYlBu_r, vmin 0 / vmax 3,
white gap rows between mutant blocks, horizontal colorbar on top.

Reads:  data/biomassTraces_normWTpeak.csv
Writes: figures/2B_biomassTraces.{png,svg}
"""
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig2/ for figlib
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from figlib import config, plotting, STRAIN_ORDER, DISPLAY_NAMES

plotting.setStyle()
GAP_ROWS = 4

df = pd.read_csv(config.TABLES / 'biomassTraces_normWTpeak.csv')
bcols = sorted([c for c in df.columns if re.match(r'^biomass_t\d+$', c)], key=lambda c: int(c.split('_t')[1]))
frames = [int(c.split('_t')[1]) for c in bcols]
nFrames = len(frames)

stacked, yticks, yticklabels, y = [], [], [], 0
for m in STRAIN_ORDER:
    sub = df[df['mutant'] == m].sort_values(['plateId', 'wellId'])
    if sub.empty:
        continue
    mat = sub[bcols].astype(float).ffill(axis=1).bfill(axis=1).values
    stacked.append(mat)
    yticks.append(y + mat.shape[0] / 2)
    yticklabels.append(DISPLAY_NAMES[m])
    y += mat.shape[0]
    stacked.append(np.full((GAP_ROWS, nFrames), np.nan))
    y += GAP_ROWS

heatmap = np.vstack(stacked)
cmap = plt.cm.RdYlBu_r.copy(); cmap.set_bad('white')

fig = plt.figure(figsize=(8, 12))
gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[1.2, 20], left=0.22, right=0.98, bottom=0.08, top=0.84, hspace=0.15)
cax = fig.add_subplot(gs[0]); ax = fig.add_subplot(gs[1])
im = ax.imshow(heatmap, aspect='auto', cmap=cmap, vmin=0, vmax=3, interpolation='nearest')

cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
cbar.set_label('Biofilm Biomass (a.u.)', labelpad=8, fontsize=44)
cbar.ax.xaxis.set_label_position('top'); cbar.ax.xaxis.set_ticks_position('top')
cbar.ax.tick_params(labelsize=36)

ax.set_xticks([i for i, t in enumerate(frames) if t % 5 == 0])
ax.set_xticklabels([t for t in frames if t % 5 == 0])
ax.set_xlabel('Time (h)', fontsize=36)
ax.tick_params(axis='x', labelsize=36)
ax.set_yticks(yticks); ax.set_yticklabels(yticklabels, fontsize=36)
ax.tick_params(axis='y', length=0)

out = config.ensure(config.FIGURES) / '2B_biomassTraces'
fig.savefig(str(out) + '.png', dpi=400, bbox_inches='tight')
fig.savefig(str(out) + '.svg', bbox_inches='tight')
print(f'Saved: {out}.png')
