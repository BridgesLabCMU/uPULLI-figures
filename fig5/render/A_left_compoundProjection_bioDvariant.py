"""fig5 render: alternative Panel 5A left - biotin view with the untreated ΔbioD clean deletion added.

Same landscape and same compound overlay as A_left_compoundProjection.py, but the biotin story is
carried by four projected conditions instead of three, and each gets its own color/shape:

  WT + DMSO                open black circle          (compound plate, cmpd_260522; all 16 wells)
  WT + 100 uM MAC13772     open purple diamond        (compound plate, anti-biotin compound; 9 wells)
  ΔbioD + biotin           open red down-triangle     (compound plate, bioD_biotin; all 16 wells)
  ΔbioD                    open grey up-triangle      (clean-deletion plates 260521/260522; 9 wells)

The two conditions that land in the crowded biotin cluster (anti-biotin, ΔbioD) are capped at 9
replicates each — a deterministic subsample, random_state=42, as in Fig 4A — and the purple diamonds
are drawn last so they read on top of the grey triangles they overlap.

Background is the reimaging landscape: grey atlas, the Tn bioA-D insertions in orange with black
edges, and reimaging WT as semi-transparent black. Note the two ΔbioD conditions come from different
experiments (treated = the compound plate, untreated = the clean-deletion plates) projected onto the
same manifold, not a treated/untreated pair from one plate.

Reads:  data/reimaging_landscape_coords.csv, data/compounds_projectedCoords.csv,
        data/cleanDeletions_projectedCoords.csv
Writes: figures/5A_left_compoundProjection_bioDvariant.{png,svg}

  python fig5/render/A_left_compoundProjection_bioDvariant.py            # all replicates
  python fig5/render/A_left_compoundProjection_bioDvariant.py --reps 9   # thin every condition to 9 (as Fig 4A)
"""
import sys
import argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig5/ for figlib
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from figlib import config, plotting

plotting.setStyle()
KEEP_FUNCS = ['Biotin Biosynthesis']          # the Tn bioA-D insertions (VC_1111/1113/1114/1115)
SEED = 42

# projected condition -> source table, label in that table, marker, edge color, legend text, cap, zorder.
# 'cmpd' = compounds_projectedCoords.csv, 'cleandel' = cleanDeletions_projectedCoords.csv.
# All four are drawn open (facecolors='none'); the color is the marker edge.
#   cap    = max replicates drawn (deterministic subsample, random_state=42; None = all).
#            The two conditions that pile into the biotin cluster are capped at 9 so they stay legible.
#   zorder = draw order within the overlay: purple (anti-biotin) must sit ON TOP of the grey ΔbioD
#            triangles it overlaps, so it is drawn last.
CONDITIONS = [
    ('cmpd',     'WT_DMSO',     'o', 'black',   'WT + DMSO',                          None, 8),
    ('cmpd',     'WT_antiBio',  'D', '#911eb4', r'WT + 100 $\mu$M MAC13772',             9, 11),
    ('cmpd',     'bioD_biotin', 'v', '#e6194B', r'$\Delta \mathit{bioD}$ + biotin',   None, 8),
    ('cleandel', 'BioD',        '^', '#787878', r'$\Delta \mathit{bioD}$',               9, 10),
]

ap = argparse.ArgumentParser()
ap.add_argument('--reps', type=int, default=0,
                help='further thin EVERY condition to N wells (0 = keep each condition\'s own cap)')
args = ap.parse_args()

# ── reimaging background ──────────────────────────────────────────────────────
reim = pd.read_csv(config.TABLES / 'reimaging_landscape_coords.csv')
reim['geneLocus'] = reim['geneLocus'].fillna('')
isWT = reim['mutant'].astype(str) == 'WT'
keptLoci = {l: f for f in KEEP_FUNCS for l in plotting.HIGHLIGHT_SETS[f]}
funcOfRow = reim['geneLocus'].map(keptLoci)
isGrey = funcOfRow.isna() & ~isWT

# ── projected conditions (already computed upstream; nothing is refit here) ────
tables = {'cmpd': pd.read_csv(config.TABLES / 'compounds_projectedCoords.csv'),
          'cleandel': pd.read_csv(config.TABLES / 'cleanDeletions_projectedCoords.csv')}


def wells(src, label, cap):
    sub = tables[src]
    sub = sub[sub['mutant'].astype(str) == label]
    limits = [n for n in (cap, args.reps or None) if n]
    if limits and len(sub) > min(limits):
        sub = sub.sample(n=min(limits), random_state=SEED)
    return sub


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
                            markerfacecolor=plotting.FUNCTION_COLORS[f], markeredgecolor='black',
                            markeredgewidth=0.6, label=r'Tn $\mathit{bioA}$-$\mathit{D}$'))
if isWT.any():
    ax.scatter(reim.loc[isWT, 'umap1'], reim.loc[isWT, 'umap2'], c='black', s=380, alpha=0.6,
               edgecolors='black', linewidth=0.5, zorder=5)
    funcH.append(Line2D([0], [0], marker='o', linestyle='none', markersize=16, markerfacecolor='black',
                        markeredgecolor='black', label='WT (reimaging)'))

projH = []
for src, label, mk, color, pretty, cap, z in CONDITIONS:
    sub = wells(src, label, cap)
    if sub.empty:
        print(f'[WARN] no wells for {label} in {src} table - skipped')
        continue
    ax.scatter(sub['umap1'], sub['umap2'], marker=mk, s=260, facecolors='none', edgecolors=color,
               linewidths=2.2, alpha=1.0, zorder=z)
    projH.append(Line2D([0], [0], marker=mk, linestyle='none', markersize=16, markerfacecolor='none',
                        markeredgecolor=color, markeredgewidth=2.2, label=pretty))
    print(f'{pretty}: {len(sub)} wells')

pad = 0.5
ax.set_xlim(reim['umap1'].min() - pad, reim['umap1'].max() + pad)
ax.set_ylim(reim['umap2'].min() - pad, reim['umap2'].max() + pad)
ax.set_xlabel('UMAP 1', fontsize=32); ax.set_ylabel('UMAP 2', fontsize=32); ax.tick_params(labelsize=28)
leg1 = ax.legend(handles=projH, title='Projected', frameon=False, fontsize=28, title_fontsize=30,
                 loc='upper left', bbox_to_anchor=(1.01, 1.0)); ax.add_artist(leg1)
leg2 = ax.legend(handles=funcH, title='Reimaging', frameon=False, fontsize=26, title_fontsize=28,
                 loc='upper left', bbox_to_anchor=(1.01, 0.45))
fig.tight_layout()
stem = '5A_left_compoundProjection_bioDvariant' + (f'_{args.reps}reps' if args.reps else '')
out = config.ensure(config.FIGURES) / stem
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight', bbox_extra_artists=(leg1, leg2))
fig.savefig(str(out) + '.svg', bbox_inches='tight', bbox_extra_artists=(leg1, leg2))
print(f'Saved: {out}.png')
