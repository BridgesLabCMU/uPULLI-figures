#!/usr/bin/env python3
"""fig2 render — Panel 2D (frames 9-30 variant): training-set UMAP over a single 9-30 window for all
feature modalities, colored by mutant. Does NOT overwrite the standard 9-27 version.

Reads:  data/trainingUmap_all_three_frames9-30_coords.csv
Writes: figures/2D_trainingUmap_frames9-30.{png,svg}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig2/ for figlib
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from figlib import config, plotting, STRAIN_ORDER, DISPLAY_NAMES, STRAIN_COLORS

plotting.setStyle(extra={'font.size': 28, 'axes.titlesize': 30, 'axes.labelsize': 30,
                         'xtick.labelsize': 24, 'ytick.labelsize': 24, 'legend.fontsize': 26, 'axes.linewidth': 2})

df = pd.read_csv(config.TABLES / 'trainingUmap_all_three_frames9-30_coords.csv')
fig, ax = plt.subplots(figsize=(11, 8))
for m in STRAIN_ORDER:
    sub = df[df['mutant'] == m]
    if sub.empty:
        continue
    ax.scatter(sub['umap1'], sub['umap2'], s=90, color=STRAIN_COLORS[m], label=DISPLAY_NAMES[m],
               edgecolor='black', linewidth=1.0, alpha=0.7)
ax.set_title('Biofilm Biomass + Whole-image + Colony-level Features (frames 9-30)')
ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)
plt.tight_layout()
out = config.ensure(config.FIGURES) / '2D_trainingUmap_frames9-30'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
fig.savefig(str(out) + '.svg', bbox_inches='tight')
print(f'Saved: {out}.png')
