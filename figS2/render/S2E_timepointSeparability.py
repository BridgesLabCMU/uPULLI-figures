#!/usr/bin/env python3
"""figS2 render: Panel S2E - mutant separability across time.

Renders FROM the bundled tables: the single-timepoint RF balanced accuracy per frame (dotted, with a
+/-1 std band) and the full-timecourse RF accuracy (dashed baseline, with its band). Built on the
quantitative features (biomass + whole-image + colony).

Reads:  data/timepointSalience_groupkfold.csv, data/fullTimecourse_accuracy.csv
Writes: figures/S2E_timepointSeparability.{png,svg}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS2/ for figlib
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from figlib import config, plotting

plotting.setStyle(extra={'font.size': 36, 'axes.titlesize': 40, 'axes.labelsize': 36,
                         'xtick.labelsize': 30, 'ytick.labelsize': 30, 'axes.linewidth': 2})

sal = pd.read_csv(config.TABLES / 'timepointSalience_groupkfold.csv')
full = pd.read_csv(config.TABLES / 'fullTimecourse_accuracy.csv')
frames = sal['frame'].values
rfMean = sal['rfMean'].values
rfStd = sal['rfStd'].values
fullMean = float(full['rfMean'].iloc[0])
fullStd = float(full['rfStd'].iloc[0])

fig, ax = plt.subplots(figsize=(16, 8))
ax.plot(frames, rfMean, linestyle=':', color='black', marker='o', markerfacecolor='white',
        markeredgecolor='black', linewidth=2, label='Random Forest (Single Timepoint)')
ax.fill_between(frames, rfMean - rfStd, rfMean + rfStd, color='black', alpha=0.15)
ax.axhline(fullMean, linestyle='--', linewidth=3, color='black', label='Random Forest (Full Timecourse)')
ax.fill_between(frames, fullMean - fullStd, fullMean + fullStd, color='black', alpha=0.05)
ax.set_xlabel('Time (h)'); ax.set_ylabel('Balanced Classification Accuracy')
ax.set_title('Mutant Separability Across Time')
ax.legend(frameon=False)
plt.tight_layout()
out = config.ensure(config.FIGURES) / 'S2E_timepointSeparability'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
fig.savefig(str(out) + '.svg', bbox_inches='tight')
print(f'Saved: {out}.png')
