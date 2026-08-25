#!/usr/bin/env python3
"""figS6 render — Panel C: observed vs null distribution of per-mutant within-replicate distance.

Both curves are distributions of the SAME quantity — one mutant's mean within-replicate distance —
so they sit on one unbroken axis: blue = the 158 observed mutants, grey = the same quantity under
the null pooled over all permutations (free shuffle and within-plate shuffle).

Reads:  data/replicate_perMutantNullHistogram.csv
Writes: figures/S6C_left_perMutantNullOverlay.{png,svg}
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

OBS_COLOR, NULL_COLOR, NULL_STRAT_COLOR = '#2b6cb0', '#9aa0a6', '#4d4d4d'

plotting.setStyle()
h = pd.read_csv(config.TABLES / 'replicate_perMutantNullHistogram.csv')

fig, ax = plt.subplots(figsize=(11, 11))
ax.set_box_aspect(1)
for col, color in (('nullGlobalDensity', NULL_COLOR), ('nullWithinPlateDensity', NULL_STRAT_COLOR),
                   ('observedDensity', OBS_COLOR)):
    ax.fill_between(h.binCenter, h[col], step='mid', color=color, alpha=0.55, linewidth=0, zorder=2)
    ax.step(h.binCenter, h[col], where='mid', color=color, linewidth=3, zorder=3)

ax.set_xlabel('Mean Within-Replicate Distance', fontsize=32)
ax.set_ylabel('Density', fontsize=32)
ax.tick_params(labelsize=28)
ax.set_xlim(h.binLeft.min(), h.binRight.max())
ax.set_ylim(bottom=0)

handles = [Line2D([0], [0], color=OBS_COLOR, linewidth=10, alpha=0.8, label='Observed'),
           Line2D([0], [0], color=NULL_COLOR, linewidth=10, alpha=0.8, label='Shuffled Labels'),
           Line2D([0], [0], color=NULL_STRAT_COLOR, linewidth=10, alpha=0.8,
                  label='Within-Plate Shuffled')]
ax.legend(handles=handles, frameon=False, fontsize=26, loc='upper right')
fig.tight_layout()

out = config.ensure(config.FIGURES / 'S6C_left_perMutantNullOverlay.png')
fig.savefig(out, dpi=300, bbox_inches='tight')
fig.savefig(str(out).replace('.png', '.svg'), bbox_inches='tight')
plt.close(fig)
print(f'Saved: {out}')
