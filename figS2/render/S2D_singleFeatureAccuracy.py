#!/usr/bin/env python3
"""figS2 render: Panel S2D - single-feature timecourse RF accuracy (ranked bar chart).

Renders FROM the bundled per-feature accuracy table: each feature base's balanced RF accuracy when it is
the ONLY feature (its timecourse alone), sorted best-first, bars colored by feature family.

Reads:  data/singleFeatureAccuracy.csv
Writes: figures/S2D_singleFeatureAccuracy.{png,svg}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS2/ for figlib
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
from figlib import config, plotting, FAMILY_COLORS

plotting.setStyle(extra={'font.size': 36, 'axes.titlesize': 40, 'axes.labelsize': 36,
                         'xtick.labelsize': 30, 'ytick.labelsize': 30, 'axes.linewidth': 2})

df = pd.read_csv(config.TABLES / 'singleFeatureAccuracy.csv')
df = df.sort_values('balancedAccuracy', ascending=False).reset_index(drop=True)

fig, ax = plt.subplots(figsize=(24, 24))
names = df['prettyName'].tolist()
colors = [FAMILY_COLORS.get(f, FAMILY_COLORS['other']) for f in df['family']]
ax.barh(names[::-1], df['balancedAccuracy'][::-1], color=colors[::-1])
ax.set_xlabel('Balanced Classification Accuracy')
ax.set_title('Single-Feature Timecourse Accuracy')
ax.legend(handles=[mpatches.Patch(color=FAMILY_COLORS['biomass'], label='Biofilm Biomass'),
                   mpatches.Patch(color=FAMILY_COLORS['whole'], label='Whole-Image Features'),
                   mpatches.Patch(color=FAMILY_COLORS['colony'], label='Colony-level Features')],
          frameon=False, loc='center left', bbox_to_anchor=(1.02, 0.5))
plt.tight_layout()
out = config.ensure(config.FIGURES) / 'S2D_singleFeatureAccuracy'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
fig.savefig(str(out) + '.svg', bbox_inches='tight')
print(f'Saved: {out}.png')
