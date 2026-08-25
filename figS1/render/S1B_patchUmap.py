#!/usr/bin/env python3
"""figS1 render: Panel S1A - training-set UMAP from DINOv2 patch tokens, colored by mutant.

Renders FROM the bundled coordinates table (euclid representation; the CSV also holds cosine).
Compare to fig1 Panel 1D (same pipeline on the CLS token).

Reads:  data/trainingPatchUmap_coords.csv
Writes: figures/S1B_patchUmap.{png,svg}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS1/ for figlib
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from figlib import config, plotting, STRAIN_ORDER, DISPLAY_NAMES, STRAIN_COLORS

METRIC = 'euclid'
plotting.setStyle(extra={'font.size': 28, 'axes.titlesize': 30, 'axes.labelsize': 30,
                         'xtick.labelsize': 24, 'ytick.labelsize': 24, 'legend.fontsize': 26, 'axes.linewidth': 2})

df = pd.read_csv(config.TABLES / 'trainingPatchUmap_coords.csv')
df = df[df['metric'] == METRIC]
fig, ax = plt.subplots(figsize=(11, 8))
for m in STRAIN_ORDER:
    sub = df[df['mutant'] == m]
    if sub.empty:
        continue
    ax.scatter(sub['umap1'], sub['umap2'], s=90, color=STRAIN_COLORS[m], label=DISPLAY_NAMES[m],
               edgecolor='black', linewidth=1.0, alpha=0.7)
ax.set_title('DINOv2 Patch Tokens (euclidean)')
ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)
plt.tight_layout()
out = config.ensure(config.FIGURES) / 'S1B_patchUmap'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
fig.savefig(str(out) + '.svg', bbox_inches='tight')
print(f'Saved: {out}.png')
