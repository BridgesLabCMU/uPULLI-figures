"""fig5 render: Panel 5C right - multispecies embedding RF confusion matrix (100% LB, 10X).

Balanced accuracy (title) = mean of the matrix diagonal. Reads the bundled mean row-normalized matrix.

Reads:  data/multispecies_100pctLB_10X_confusion_cv.csv   Writes: figures/5C_right_multispeciesConfusion.{png,svg}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig5/ for figlib
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from figlib import config, plotting, SPECIES, SPECIES_DISPLAY

plotting.setStyle(extra={'font.size': 40, 'axes.titlesize': 46, 'axes.labelsize': 44,
                         'xtick.labelsize': 36, 'ytick.labelsize': 36, 'axes.linewidth': 2})
cm = pd.read_csv(config.TABLES / 'multispecies_100pctLB_10X_confusion_cv.csv', index_col=0)
labels = [s for s in SPECIES if s in cm.index]
cm = cm.loc[labels, labels]; M = cm.values
meanAcc = float(np.mean(np.diag(M)))
fig, ax = plt.subplots(figsize=(20, 16))
im = ax.imshow(M, cmap='viridis', vmin=0, vmax=1)
ax.set_xticks(range(len(labels))); ax.set_yticks(range(len(labels)))
ax.set_xticklabels([SPECIES_DISPLAY[l] for l in labels], rotation=45, ha='right', fontsize=36)
ax.set_yticklabels([SPECIES_DISPLAY[l] for l in labels], fontsize=36)
ax.set_xlabel('Predicted Label', labelpad=20); ax.set_ylabel('True Label', labelpad=20)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        ax.text(j, i, f'{M[i, j]:.2f}', ha='center', va='center', fontsize=30, fontweight='bold',
                color='black' if M[i, j] > 0.5 else 'white')
plt.title(f'Multispecies (Embeddings), 100% LB, 10X\nBalanced Accuracy (CV) = {meanAcc:.3f}', pad=30)
cb = plt.colorbar(im); cb.ax.tick_params(labelsize=32)
plt.tight_layout()
out = config.ensure(config.FIGURES) / '5C_right_multispeciesConfusion'
fig.savefig(str(out) + '.png', dpi=300); fig.savefig(str(out) + '.svg')
print(f'Saved: {out}.png  (balanced accuracy = {meanAcc:.3f})')
