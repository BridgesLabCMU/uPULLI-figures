#!/usr/bin/env python3
"""figS6 render — Panel D: every mutant's within-replicate distance against the permutation null.

One point per mutant, sorted tightest-to-loosest, colored by functional group. Dashed line = the
null expectation; grey band = +-2 SD of a single mutant's null (drawn at the across-mutant mean SD
as a visual guide — each mutant's actual test uses its own null, reported as BH-FDR q in the source
table). Counts and thresholds live in the figure legend, not on the axes.

Reads:  data/replicate_perMutant.csv
Writes: figures/S6C_right_perMutantVsNull.{png,svg}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS6/ for figlib
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from figlib import config, plotting

WT_COLOR, OTHER_COLOR = '#000000', '#bdbdbd'

plotting.setStyle()
pm = pd.read_csv(config.TABLES / 'replicate_perMutant.csv')

locusToFunc = {l: f for f, loci in plotting.HIGHLIGHT_SETS.items() for l in loci}
pm['funcGroup'] = np.where(pm.mutant.astype(str) == 'WT', 'WT',
                           pm.geneLocus.fillna('').map(locusToFunc).fillna('Unclassified'))
colorOf = dict(plotting.FUNCTION_COLORS, WT=WT_COLOR, Unclassified=OTHER_COLOR)
pm['color'] = pm.funcGroup.map(colorOf)
pm['z'] = pm.funcGroup.map(lambda f: {'Unclassified': 2, 'WT': 6}.get(f, 4))

order = [f for f in plotting.HIGHLIGHT_SETS if f in set(pm.funcGroup)] + ['WT', 'Unclassified']
handles = [Line2D([0], [0], color='black', linestyle='--', linewidth=2.5, label='Chance')] + [
    Line2D([0], [0], marker='o', linestyle='none', markersize=16, markerfacecolor=colorOf[f],
           markeredgecolor='black', markeredgewidth=1.0, label=plotting.functionLabel(f))
    for f in order]

s = pm.sort_values('meanWithin').reset_index(drop=True)
nullMean, nullSd = float(s.nullMeanWithin.mean()), float(s.nullSdWithin.mean())

fig, ax = plt.subplots(figsize=(13, 10))
ax.axhspan(nullMean - 2 * nullSd, nullMean + 2 * nullSd, color='#c8c8c8', alpha=0.6, zorder=1)
ax.axhline(nullMean, color='black', linestyle='--', linewidth=2.5, zorder=3)
for z in sorted(set(s.z)):
    sub = s[s.z == z]
    ax.scatter(sub.index, sub.meanWithin, s=200, c=sub.color, edgecolors='black',
               linewidth=1.0, alpha=0.9, zorder=int(z))

ax.set_xlabel('Ranked Mutants', fontsize=32)
ax.set_ylabel('Mean Within-Replicate Distance', fontsize=32)
ax.tick_params(labelsize=28)
ax.set_xlim(-3, len(s) + 2)
ax.legend(handles=handles, frameon=False, fontsize=24, loc='center left', bbox_to_anchor=(1.02, 0.5))
fig.tight_layout()

out = config.ensure(config.FIGURES / 'S6C_right_perMutantVsNull.png')
fig.savefig(out, dpi=300, bbox_inches='tight')
fig.savefig(str(out).replace('.png', '.svg'), bbox_inches='tight')
plt.close(fig)
print(f'Saved: {out}')
