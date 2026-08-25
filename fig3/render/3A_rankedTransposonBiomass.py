#!/usr/bin/env python3
"""fig3 render: Panel 3A - ranked peak-biomass classification of the transposon screen.

Renders FROM the bundled per-well phenotype summary + thresholds. Each transposon mutant is a point,
ranked by normalized peak biofilm biomass (descending) and colored by its WT-anchored phenotypic class.
The B_max / B_min classification boundaries (from the WT replicate distribution) are drawn as dashed
lines.

The legend is written as its OWN file rather than placed on the axes, so the panel stays a clean square
plot and the key can be laid out independently at assembly.

Reads:  data/tn_phenotype_summary.csv, data/tn_thresholds.csv
Writes: figures/3A_rankedTransposonBiomass.{png,svg}  +  figures/3A_rankedTransposonBiomass_legend.{png,svg}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig3/ for figlib
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from figlib import config, plotting, PHENO_COLORS, PHENO_ORDER

plotting.setStyle(extra={'font.size': 44, 'axes.labelsize': 52, 'axes.titlesize': 52,
                         'xtick.labelsize': 44, 'ytick.labelsize': 44, 'legend.fontsize': 44,
                         'axes.linewidth': 3, 'xtick.major.size': 12, 'ytick.major.size': 12,
                         'xtick.major.width': 3, 'ytick.major.width': 3})

s = pd.read_csv(config.TABLES / 'tn_phenotype_summary.csv')
thr = pd.read_csv(config.TABLES / 'tn_thresholds.csv').iloc[0]
s = s.sort_values('Peak', ascending=False).reset_index(drop=True)
counts = s['Phenotype'].value_counts()

upperColor, lowerColor = PHENO_COLORS['High Biofilm'], PHENO_COLORS['Low Biofilm']

fig, ax = plt.subplots(figsize=(16, 16))
ax.set_box_aspect(1)
handles = []
for p in PHENO_ORDER:
    sub = s[s['Phenotype'] == p]
    if sub.empty:
        continue
    handles.append(ax.scatter(sub.index, sub['Peak'], s=170, color=PHENO_COLORS[p],
                              edgecolor='none', rasterized=True,
                              label=f'{p} (n = {int(counts.get(p, 0))})'))
ax.axhline(thr['Bmax'], ls='--', lw=4.5, color=upperColor)
ax.axhline(thr['Bmin'], ls='--', lw=4.5, color=lowerColor)
handles.append(Line2D([0], [0], ls='--', lw=4.5, color=upperColor, label=r'$B_{max}$ (High Biofilm threshold)'))
handles.append(Line2D([0], [0], ls='--', lw=4.5, color=lowerColor, label=r'$B_{min}$ (Low Biofilm threshold)'))
# wrapped: at this font size a single-line label is wider than the canvas and gets clipped
ax.set_xlabel('Transposon Ranking\n(by Peak Biofilm Biomass)')
ax.set_ylabel('Peak Biofilm Biomass (a.u.)')
ax.set_ylim(0, min(5, float(s['Peak'].max()) * 1.05))
ax.set_xlim(-len(s) * 0.01, len(s) * 1.02)
plt.tight_layout()

out = config.ensure(config.FIGURES) / '3A_rankedTransposonBiomass'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
fig.savefig(str(out) + '.svg', dpi=400, bbox_inches='tight')   # dpi = resolution of the rasterized point cloud
plt.close(fig)

# ── legend as a standalone file (same handles, same style; no axes) ──
figL = plt.figure(figsize=(14, 5))
figL.legend(handles=handles, loc='center', frameon=False, markerscale=2.2, handletextpad=0.9,
            labelspacing=0.8)
figL.savefig(str(out) + '_legend.png', dpi=300, bbox_inches='tight')
figL.savefig(str(out) + '_legend.svg', bbox_inches='tight')
plt.close(figL)
print(f'Saved: {out}.png  +  {out}_legend.png')
