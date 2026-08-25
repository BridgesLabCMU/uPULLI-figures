#!/usr/bin/env python3
"""figS6 render — Panel E: paired within- vs between-mutant distance, one line per mutant.

The statistically clean form of the claim: the unit of analysis is the MUTANT (n = 158), so unlike
pairwise distances these observations are approximately independent and a paired test is valid.

The bracket prints an inequality rather than the raw Wilcoxon p. That p (5.6e-28) is correct
arithmetic — the normal approximation with W = 0, n = 158 — and is even conservative against the
exact signed-rank floor of 2^-158, but it assumes the 158 paired differences are independent, and
each mutant's distance to "other mutants" is computed against nearly the same pool of wells. Two
significant figures would imply precision the design does not support. Use --sig stars for the
conventional star code instead.

Reads:  data/replicate_perMutant.csv, data/replicate_summary.csv
Writes: figures/S6B_right_pairedWithinVsBetween.{png,svg}
"""
import sys
import argparse
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
P_FLOOR = 1e-20      # do not print a p below this -- see the docstring

ap = argparse.ArgumentParser()
ap.add_argument('--sig', choices=['value', 'stars', 'none'], default='value')
args = ap.parse_args()

plotting.setStyle()
pm = pd.read_csv(config.TABLES / 'replicate_perMutant.csv')
sm = pd.read_csv(config.TABLES / 'replicate_summary.csv').set_index('statistic')['value']

locusToFunc = {l: f for f, loci in plotting.HIGHLIGHT_SETS.items() for l in loci}
pm['funcGroup'] = np.where(pm.mutant.astype(str) == 'WT', 'WT',
                           pm.geneLocus.fillna('').map(locusToFunc).fillna('Unclassified'))
colorOf = dict(plotting.FUNCTION_COLORS, WT=WT_COLOR, Unclassified=OTHER_COLOR)
pm['color'] = pm.funcGroup.map(colorOf)
pm['z'] = pm.funcGroup.map(lambda f: {'Unclassified': 2, 'WT': 6}.get(f, 4))

order = [f for f in plotting.HIGHLIGHT_SETS if f in set(pm.funcGroup)] + ['WT', 'Unclassified']
handles = [Line2D([0], [0], marker='o', linestyle='none', markersize=16, markerfacecolor=colorOf[f],
                  markeredgecolor='black', markeredgewidth=1.0, label=plotting.functionLabel(f))
           for f in order]

fig, ax = plt.subplots(figsize=(10, 11))
for _, r in pm.iterrows():
    ax.plot([0, 1], [r.meanWithin, r.meanToOtherMutants], color=r.color, linewidth=1.6,
            alpha=0.55, zorder=int(r.z), solid_capstyle='round')
for z in sorted(set(pm.z)):
    sub = pm[pm.z == z]
    ax.scatter(np.zeros(len(sub)), sub.meanWithin, s=150, c=sub.color, edgecolors='black',
               linewidth=1.0, alpha=0.95, zorder=int(z) + 10)
    ax.scatter(np.ones(len(sub)), sub.meanToOtherMutants, s=150, c=sub.color, edgecolors='black',
               linewidth=1.0, alpha=0.95, zorder=int(z) + 10)

ax.set_xticks([0, 1])
# same vocabulary as panel B's legend, so one wording covers both panels
ax.set_xticklabels(['Same\nMutant', 'Different\nMutants'])
ax.set_xlim(-0.35, 1.35)
ax.set_ylabel('Mean Distance', fontsize=32)
ax.tick_params(labelsize=28)

if args.sig != 'none':
    p = float(sm['wilcoxon_p'])
    if args.sig == 'stars':
        txt = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
    elif p < P_FLOOR:
        txt = f'$p$ < 10$^{{{int(np.log10(P_FLOOR))}}}$'
    else:
        txt = f'$p$ = {p:.3g}'
    top = max(pm.meanWithin.max(), pm.meanToOtherMutants.max())
    span = top - min(pm.meanWithin.min(), pm.meanToOtherMutants.min())
    ax.set_ylim(top=top + 0.20 * span)
    y, h = top + 0.05 * span, 0.02 * span
    ax.plot([0, 0, 1, 1], [y, y + h, y + h, y], color='black', linewidth=2.0, zorder=20)
    ax.annotate(txt, xy=(0.5, y + h), xytext=(0, 6), textcoords='offset points',
                ha='center', va='bottom', fontsize=26, zorder=20)

ax.legend(handles=handles, frameon=False, fontsize=22, loc='center left', bbox_to_anchor=(1.02, 0.5))
fig.tight_layout()

out = config.ensure(config.FIGURES / 'S6B_right_pairedWithinVsBetween.png')
fig.savefig(out, dpi=300, bbox_inches='tight')
fig.savefig(str(out).replace('.png', '.svg'), bbox_inches='tight')
plt.close(fig)
print(f'Saved: {out}')
