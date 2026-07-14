"""fig5 render: Panel 5B left - K. pneumoniae embedding UMAP (euclid), colored by mutant.

Reads:  data/kleb_embeddingUmap_coords.csv   Writes: figures/5B_left_klebUmap.{png,svg}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig5/ for figlib
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from figlib import config, plotting, KLEB_ORDER, KLEB_DISPLAY, KLEB_COLORS

plotting.setStyle(extra={'font.size': 28, 'axes.titlesize': 30, 'axes.labelsize': 30,
                         'xtick.labelsize': 24, 'ytick.labelsize': 24, 'legend.fontsize': 26, 'axes.linewidth': 2})
df = pd.read_csv(config.TABLES / 'kleb_embeddingUmap_coords.csv')
df = df[df['metric'] == 'euclid']
fig, ax = plt.subplots(figsize=(11, 8))
for m in KLEB_ORDER:
    sub = df[df['mutant'] == m]
    if sub.empty:
        continue
    ax.scatter(sub['umap1'], sub['umap2'], s=90, color=KLEB_COLORS[m], label=KLEB_DISPLAY[m],
               edgecolor='black', linewidth=1.0, alpha=0.7)
ax.set_title(r'$\mathit{K.\ pneumoniae}$ (Euclidean)')
ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)
plt.tight_layout()
out = config.ensure(config.FIGURES) / '5B_left_klebUmap'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
fig.savefig(str(out) + '.svg', bbox_inches='tight')
print(f'Saved: {out}.png')
