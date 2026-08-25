#!/usr/bin/env python3
"""figS2 render: Panel S2B - colony-level (segmented) feature Random-Forest confusion matrix.

Renders FROM the bundled mean row-normalized confusion matrix. Balanced accuracy (title) = mean of the
diagonal. Same style as Fig 2E — this is the segmented-colony-feature slice.

Reads:  data/colony_confusion_cv.csv
Writes: figures/S2B_colonyConfusion.{png,svg}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS2/ for figlib
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from figlib import config, plotting, STRAIN_ORDER, DISPLAY_NAMES

CSV, TITLE, STEM = 'colony_confusion_cv.csv', 'Colony-level Features', 'S2B_colonyConfusion'
plotting.setStyle(extra={'font.size': 40, 'axes.titlesize': 46, 'axes.labelsize': 44,
                         'xtick.labelsize': 36, 'ytick.labelsize': 36, 'axes.linewidth': 2})

cm = pd.read_csv(config.TABLES / CSV, index_col=0)
labels = [s for s in STRAIN_ORDER if s in cm.index]
M = cm.loc[labels, labels].values
meanAcc = float(np.mean(np.diag(M)))

fig, ax = plt.subplots(figsize=(20, 16))
im = ax.imshow(M, cmap='viridis', vmin=0, vmax=1)
ax.set_xticks(range(len(labels))); ax.set_yticks(range(len(labels)))
ax.set_xticklabels([DISPLAY_NAMES[l] for l in labels], rotation=45, ha='right')
ax.set_yticklabels([DISPLAY_NAMES[l] for l in labels])
ax.set_xlabel('Predicted Label', labelpad=20); ax.set_ylabel('True Label', labelpad=20)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        ax.text(j, i, f'{M[i, j]:.2f}', ha='center', va='center', fontsize=34, fontweight='bold',
                color='black' if M[i, j] > 0.5 else 'white')
plt.title(f'{TITLE}\nBalanced Accuracy (CV) = {meanAcc:.3f}', pad=30)
cb = plt.colorbar(im); cb.ax.tick_params(labelsize=32)
plt.tight_layout()
out = config.ensure(config.FIGURES) / STEM
fig.savefig(str(out) + '.png', dpi=300)
fig.savefig(str(out) + '.svg')
print(f'Saved: {out}.png  (balanced accuracy = {meanAcc:.3f})')
