"""fig5 render: Panel 5A left - compound projection onto the reimaging landscape (biotin view).

Renders FROM bundled tables: the reimaging landscape (background, Biotin Biosynthesis + WT highlighted)
and the compound projected coords. Overlay = open markers for WT+DMSO, WT+anti-biotin, ΔbioD+biotin.

Reads:  data/reimaging_landscape_coords.csv, data/compounds_projectedCoords.csv
Writes: figures/5A_left_compoundProjection.{png,svg}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig5/ for figlib
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from figlib import config, plotting, CMPD_MARKERS, CMPD_PRETTY

plotting.setStyle()
KEEP_FUNCS = ['Biotin Biosynthesis']
CONDS = ['WT_DMSO', 'WT_antiBio', 'bioD_biotin']

reim = pd.read_csv(config.TABLES / 'reimaging_landscape_coords.csv')
reim['geneLocus'] = reim['geneLocus'].fillna('')
isWT = reim['mutant'].astype(str) == 'WT'
keptLoci = {l: f for f in KEEP_FUNCS for l in plotting.HIGHLIGHT_SETS[f]}
funcOfRow = reim['geneLocus'].map(keptLoci)
isGrey = funcOfRow.isna() & ~isWT
proj = pd.read_csv(config.TABLES / 'compounds_projectedCoords.csv')

fig, ax = plt.subplots(figsize=(15, 14)); ax.set_box_aspect(1)
ax.scatter(reim.loc[isGrey, 'umap1'], reim.loc[isGrey, 'umap2'], c=plotting.BACKGROUND_COLOR,
           s=300, alpha=0.15, edgecolors='black', linewidth=0.5, zorder=1)
funcH = []
for f in KEEP_FUNCS:
    mm = funcOfRow == f
    if mm.any():
        ax.scatter(reim.loc[mm, 'umap1'], reim.loc[mm, 'umap2'], c=plotting.FUNCTION_COLORS[f], s=300, alpha=0.3,
                   edgecolors='black', linewidth=0.5, zorder=3)
        funcH.append(Line2D([0], [0], marker='o', linestyle='none', markersize=16,
                            markerfacecolor=plotting.FUNCTION_COLORS[f], markeredgecolor='black', markeredgewidth=0.6, label=f))
if isWT.any():
    ax.scatter(reim.loc[isWT, 'umap1'], reim.loc[isWT, 'umap2'], c='black', s=380, alpha=0.6,
               edgecolors='black', linewidth=0.5, zorder=5)
    funcH.append(Line2D([0], [0], marker='o', linestyle='none', markersize=16, markerfacecolor='black',
                        markeredgecolor='black', label='WT (reimaging)'))

projH = []
for c in [c for c in CONDS if c in set(proj['mutant'].dropna())]:
    sub = proj[proj['mutant'] == c]; mk = CMPD_MARKERS.get(c, 'X')
    ax.scatter(sub['umap1'], sub['umap2'], marker=mk, s=260, facecolors='none', edgecolors='black',
               linewidths=2.2, alpha=1.0, zorder=8)
    projH.append(Line2D([0], [0], marker=mk, linestyle='none', markersize=16, markerfacecolor='none',
                        markeredgecolor='black', markeredgewidth=2.2, label=CMPD_PRETTY.get(c, c)))

pad = 0.5
ax.set_xlim(reim['umap1'].min() - pad, reim['umap1'].max() + pad)
ax.set_ylim(reim['umap2'].min() - pad, reim['umap2'].max() + pad)
ax.set_xlabel('UMAP 1', fontsize=32); ax.set_ylabel('UMAP 2', fontsize=32); ax.tick_params(labelsize=28)
leg1 = ax.legend(handles=projH, title='Compound treatment', frameon=False, fontsize=28, title_fontsize=30,
                 loc='upper left', bbox_to_anchor=(1.01, 1.0)); ax.add_artist(leg1)
leg2 = ax.legend(handles=funcH, title='Reimaging function', frameon=False, fontsize=26, title_fontsize=28,
                 loc='upper left', bbox_to_anchor=(1.01, 0.5))
fig.tight_layout()
out = config.ensure(config.FIGURES) / '5A_left_compoundProjection'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight', bbox_extra_artists=(leg1, leg2))
fig.savefig(str(out) + '.svg', bbox_inches='tight', bbox_extra_artists=(leg1, leg2))
print(f'Saved: {out}.png')
