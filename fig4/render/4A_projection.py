"""fig4 render: Panel 4A - clean-deletion projection onto the reimaging landscape (open markers, 9 reps).

Renders FROM bundled tables: the reimaging landscape coordinates (background) and the projected
clean-deletion coordinates. Background keeps only WT + Biotin/Pyruvate/O-Antigen groups highlighted; the
rest is grey. Projected overlay = open shapes (no fill, black edges), thinned to 9 reps per mutant
(deterministic, random_state=42). Projected WT is omitted (represented by the reimaging WT).

Reads:  data/reimaging_landscape_coords.csv, data/cleanDeletions_projectedCoords.csv
Writes: figures/4A_projection.{png,svg}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig4/ for figlib
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import argparse
from figlib import config, plotting, CLEANDEL, CLEANDEL_DISPLAY

plotting.setStyle()
KEEP_FUNCS = ['Biotin Biosynthesis', 'Pyruvate Flux', 'O-Antigen Biosynthesis']
SEED = 42

ap = argparse.ArgumentParser()
ap.add_argument('--reps', type=int, default=0, help='replicates drawn per mutant (0 = all)')
ap.add_argument('--rep-alpha', type=float, default=0.38, help='opacity of the replicate cloud')
args = ap.parse_args()

# matplotlib sizes a marker by the area of its unit bounding box, so a square carries twice the ink of
# a triangle or diamond at the same `s` and reads as bigger. These factors even out the perceived size
# so the three shapes are comparable at a glance -- all are centred on their true coordinate.
SHAPE_SCALE = {'^': 1.45, 's': 1.0, 'D': 1.30, 'o': 1.15}
REP_S, CEN_S = 110, 1150

reim = pd.read_csv(config.TABLES / 'reimaging_landscape_coords.csv')
reim['geneLocus'] = reim['geneLocus'].fillna('')
isWT = reim['mutant'].astype(str) == 'WT'
keptLoci = {l: f for f in KEEP_FUNCS for l in plotting.HIGHLIGHT_SETS[f]}
funcOfRow = reim['geneLocus'].map(keptLoci)
isGrey = funcOfRow.isna() & ~isWT

proj = pd.read_csv(config.TABLES / 'cleanDeletions_projectedCoords.csv')
proj = proj[proj['mutant'].isin(CLEANDEL)]
if args.reps:
    proj = (proj.groupby('mutant', group_keys=False)
            .apply(lambda g: g.sample(n=min(args.reps, len(g)), random_state=SEED)))
proj = proj.reset_index(drop=True)
projLabs = [m for m in CLEANDEL if m in set(proj['mutant'].dropna())]

fig, ax = plt.subplots(figsize=(15, 14)); ax.set_box_aspect(1)
ax.scatter(reim.loc[isGrey, 'umap1'], reim.loc[isGrey, 'umap2'], c=plotting.BACKGROUND_COLOR,
           s=300, alpha=0.15, edgecolors='black', linewidth=0.5, zorder=1)
funcH = []
for f in KEEP_FUNCS:
    m = funcOfRow == f
    if m.any():
        ax.scatter(reim.loc[m, 'umap1'], reim.loc[m, 'umap2'], c=plotting.FUNCTION_COLORS[f], s=300, alpha=0.3,
                   edgecolors='black', linewidth=0.5, zorder=3)
        funcH.append(Line2D([0], [0], marker='o', linestyle='none', markersize=16,
                            markerfacecolor=plotting.FUNCTION_COLORS[f], markeredgecolor='black', markeredgewidth=0.6, label=f))
if isWT.any():
    ax.scatter(reim.loc[isWT, 'umap1'], reim.loc[isWT, 'umap2'], c='black', s=380, alpha=0.6,
               edgecolors='black', linewidth=0.5, zorder=5)
    funcH.append(Line2D([0], [0], marker='o', linestyle='none', markersize=16, markerfacecolor='black',
                        markeredgecolor='black', label='WT (reimaging)'))

# Each mutant: the full replicate cloud drawn faint and small, then its centroid as one large hollow
# black marker on top -- heavy stroke, no fill, so the landscape underneath stays visible. The centroid
# is what the eye should land on; the cloud shows the spread.
projH = []
for m in projLabs:
    sub = proj[proj['mutant'] == m]
    mk = CLEANDEL[m][0]
    scale = SHAPE_SCALE.get(mk, 1.0)
    ax.scatter(sub['umap1'], sub['umap2'], marker=mk, s=REP_S * scale, facecolors='none',
               edgecolors='black', linewidths=1.2, alpha=args.rep_alpha, zorder=7)
    cx, cy = sub['umap1'].mean(), sub['umap2'].mean()
    ax.scatter([cx], [cy], marker=mk, s=CEN_S * scale, facecolors='none', edgecolors='black',
               linewidths=4.0, alpha=1.0, zorder=9)
    projH.append(Line2D([0], [0], marker=mk, linestyle='none', markersize=20 * (scale ** 0.5),
                        markerfacecolor='none', markeredgecolor='black', markeredgewidth=3.0,
                        label=f'{CLEANDEL_DISPLAY[m]} (clean del.)'))
    print(f'  {m:6s} n={len(sub):3d}  centroid ({cx:.2f}, {cy:.2f})')

pad = 0.5
ax.set_xlim(reim['umap1'].min() - pad, reim['umap1'].max() + pad)
ax.set_ylim(reim['umap2'].min() - pad, reim['umap2'].max() + pad)
ax.set_xlabel('UMAP 1', fontsize=32); ax.set_ylabel('UMAP 2', fontsize=32); ax.tick_params(labelsize=28)
leg1 = ax.legend(handles=projH, title='Projected: clean deletions', frameon=False, fontsize=30,
                 title_fontsize=30, loc='upper left', bbox_to_anchor=(1.01, 1.0)); ax.add_artist(leg1)
ax.legend(handles=funcH, title='Reimaging function', frameon=False, fontsize=26, title_fontsize=28,
          loc='upper left', bbox_to_anchor=(1.01, 0.5))
fig.tight_layout()
out = config.ensure(config.FIGURES) / '4A_projection'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
fig.savefig(str(out) + '.svg', bbox_inches='tight')
print(f'Saved: {out}.png')
