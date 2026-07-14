"""fig5 render: Panel 5C left - multispecies embedding UMAP (100% LB, 10X), colored by species.

Reads:  data/multispecies_100pctLB_10X_umap_coords.csv   Writes: figures/5C_left_multispeciesUmap.{png,svg}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig5/ for figlib
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from figlib import config, plotting, SPECIES, SPECIES_COLORS, SPECIES_DISPLAY

plotting.setStyle(extra={'font.size': 28, 'axes.titlesize': 30, 'axes.labelsize': 30,
                         'xtick.labelsize': 24, 'ytick.labelsize': 24, 'legend.fontsize': 26, 'axes.linewidth': 2})
df = pd.read_csv(config.TABLES / 'multispecies_100pctLB_10X_umap_coords.csv')
fig, ax = plt.subplots(figsize=(11, 8))
for s in SPECIES:
    sub = df[df['species'] == s]
    if sub.empty:
        continue
    ax.scatter(sub['umap1'], sub['umap2'], s=90, color=SPECIES_COLORS[s], label=SPECIES_DISPLAY[s],
               edgecolor='black', linewidth=1.0, alpha=0.7)
ax.set_title('Multispecies (Embeddings), 100% LB, 10X')
ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)
plt.tight_layout()
out = config.ensure(config.FIGURES) / '5C_left_multispeciesUmap'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
fig.savefig(str(out) + '.svg', bbox_inches='tight')
print(f'Saved: {out}.png')
