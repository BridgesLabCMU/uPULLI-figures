#!/usr/bin/env python3
"""figS6 render — Panel B: within- vs between-mutant pairwise-distance distributions.

Renders FROM the bundled source-data table (no recompute). Deliberately spare: the pair counts,
AUC, Cohen's d and the permutation test belong in the figure legend, not on the axes (see
../README.md).

Reads:  data/replicate_distanceHistogram.csv
Writes: figures/S6B_left_distanceDistributions.{png,svg}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS6/ for figlib
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
from figlib import config, plotting

WITHIN_COLOR, BETWEEN_COLOR = '#2b6cb0', '#9aa0a6'

plotting.setStyle()
h = pd.read_csv(config.TABLES / 'replicate_distanceHistogram.csv')

fig, ax = plt.subplots(figsize=(11, 11))
ax.set_box_aspect(1)
for col, color in (('betweenDensity', BETWEEN_COLOR), ('withinDensity', WITHIN_COLOR)):
    ax.fill_between(h.binCenter, h[col], step='mid', color=color, alpha=0.55, linewidth=0, zorder=2)
    ax.step(h.binCenter, h[col], where='mid', color=color, linewidth=3, zorder=3)

ax.set_xlabel('Pairwise Distance', fontsize=32)
ax.set_ylabel('Density', fontsize=32)
ax.tick_params(labelsize=28)
ax.set_xlim(h.binLeft.min(), h.binRight.max())
ax.set_ylim(bottom=0)

# plural "Different Mutants": each grey observation is a PAIR of wells from two different mutants
handles = [Line2D([0], [0], color=WITHIN_COLOR, linewidth=10, alpha=0.8, label='Same Mutant'),
           Line2D([0], [0], color=BETWEEN_COLOR, linewidth=10, alpha=0.8, label='Different Mutants')]
ax.legend(handles=handles, frameon=False, fontsize=28, loc='upper right')
fig.tight_layout()

out = config.ensure(config.FIGURES / 'S6B_left_distanceDistributions.png')
fig.savefig(out, dpi=300, bbox_inches='tight')
fig.savefig(str(out).replace('.png', '.svg'), bbox_inches='tight')
plt.close(fig)
print(f'Saved: {out}')
