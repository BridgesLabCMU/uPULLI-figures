"""fig4 render: Panel 4C left - normalized biofilm-biomass over time.

Renders FROM the bundled trace table: mean +- SD biomass(t) for reimaging WT and the three clean
deletions, normalized to the reimaging WT median peak. Square axes, no title, functional-group colors.

Reads:  data/biomassOverTime_normWTpeak.csv
Writes: figures/4C_left_biomassOverTime.{png,svg}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig4/ for figlib
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from figlib import config, plotting, CLEANDEL, CLEANDEL_DISPLAY

plotting.setStyle()
df = pd.read_csv(config.TABLES / 'biomassOverTime_normWTpeak.csv')
# (group, display, marker, color) in plot order
SERIES = [('WT', 'WT (reimaging)', 'o', 'black')] + \
         [(g, CLEANDEL_DISPLAY[g], CLEANDEL[g][0], CLEANDEL[g][1]) for g in ['BioD', 'ManA', 'PdhE2']]

fig, ax = plt.subplots(figsize=(10, 10)); ax.set_box_aspect(1)
handles = []
for grp, disp, mk, color in SERIES:
    s = df[df['group'] == grp].sort_values('frame')
    x, mean, sd = s['frame'].to_numpy(), s['mean'].to_numpy(), s['sd'].to_numpy()
    ax.fill_between(x, mean - sd, mean + sd, color=color, alpha=0.2, linewidth=0, zorder=2)
    ax.plot(x, mean, color=color, linewidth=2.5, zorder=3)
    ax.scatter(x, mean, marker=mk, s=260, facecolors=color, edgecolors='black', linewidths=1.0, zorder=4)
    handles.append(Line2D([0], [0], marker=mk, linestyle='-', color=color, markerfacecolor=color,
                          markeredgecolor='black', markeredgewidth=1.0, markersize=16, linewidth=2.5, label=disp))

ax.set_xlim(0, 30); ax.set_ylim(bottom=0); ax.set_xticks(range(0, 31, 5))
ax.set_xlabel('Time (h)', fontsize=32); ax.set_ylabel('Biofilm Biomass (a.u.)', fontsize=32)
ax.tick_params(labelsize=28)
ax.legend(handles=handles, frameon=False, fontsize=26, loc='upper left')
fig.tight_layout()
out = config.ensure(config.FIGURES) / '4C_left_biomassOverTime'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
fig.savefig(str(out) + '.svg', bbox_inches='tight')
print(f'Saved: {out}.png')
