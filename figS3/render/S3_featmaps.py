#!/usr/bin/env python3
"""figS3 render: temporal mutant x frame heatmaps for ALL feature classes (one PNG + SVG each).

Renders FROM the bundled per-feature matrices + meta, in the same style as Fig 2F (per-family colormap,
masked low-biomass colony cells black, 2-line title: feature name (+unit) and single-feature RF
accuracy). Output files are named by PANEL LETTER (S3A..S3AA) in the figure's laid-out order (colony
features alphabetical, then Haralick 0-12, then global image entropy, then biofilm biomass), matching
the repo's <fig><panel> naming (e.g. 1D, 2E). Every heatmap shares the same figure size and frame axis.

Reads:  data/featmaps_meta.csv + data/featmap_<feat>.csv (x N)
Writes: figures/S3<letter>.{png,svg}  (+ prints a letter->feature manifest)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS3/ for figlib
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from figlib import config, plotting, STRAIN_ORDER, DISPLAY_NAMES

plotting.setStyle(extra={'font.size': 36, 'axes.titlesize': 40, 'axes.labelsize': 36,
                         'xtick.labelsize': 30, 'ytick.labelsize': 30, 'axes.linewidth': 2})

# Panel order as laid out in the figure (colony alphabetical -> Haralick 0-12 -> global image entropy
# -> biofilm biomass). Output filename = S3<panel letter>.
PANEL_ORDER = [
    'colony_area_um2_mean', 'colony_area_um2_std', 'colony_bgCV', 'colony_centroidOffset_um_mean',
    'colony_eccentricity_mean', 'colony_majorAxisLength_um_mean', 'colony_meanIntensity_kurtosis',
    'colony_meanIntensity_mean', 'colony_mstEdgeMax_um_mean', 'colony_nnDistance1_um_mean',
    'colony_nnDistance1_um_std', 'nColonies',
    'whole_haralick_0', 'whole_haralick_1', 'whole_haralick_2', 'whole_haralick_3', 'whole_haralick_4',
    'whole_haralick_5', 'whole_haralick_6', 'whole_haralick_7', 'whole_haralick_8', 'whole_haralick_9',
    'whole_haralick_10', 'whole_haralick_11', 'whole_haralick_12',
    'whole_entropy', 'biomass',
]


def panelLetter(i):
    return chr(65 + i) if i < 26 else 'A' + chr(65 + i - 26)


def panelStem(letter, feat, label):
    """S3<letter>_<English feature name>[_Haralick<N>] — panel letter + descriptive name."""
    name = label.replace('(', '').replace(')', '').replace(' ', '_')
    if feat.startswith('whole_haralick_'):
        name = f'{name}_Haralick{feat.split("_")[-1]}'
    return f'S3{letter}_{name}'


# per-family colormap: biomass = red-blue diverging, whole-image = viridis, colony = plasma
CMAP_BY_FAMILY = {'biomass': 'coolwarm', 'whole': 'viridis', 'colony': 'plasma'}

meta = pd.read_csv(config.TABLES / 'featmaps_meta.csv').set_index('feature')
outDir = config.ensure(config.FIGURES)
manifest = []

for i, feat in enumerate(PANEL_ORDER):
    if feat not in meta.index:
        print(f'[skip] {feat} not in meta'); continue
    r = meta.loc[feat]
    letter = panelLetter(i)
    label, unit = r['label'], ('' if pd.isna(r['unit']) else str(r['unit']))
    mat = pd.read_csv(config.TABLES / f'featmap_{feat}.csv', index_col=0)
    mat = mat.loc[[m for m in STRAIN_ORDER if m in mat.index]]
    frames = [int(c) for c in mat.columns]
    heatmap = np.ma.masked_invalid(mat.values.astype(float))
    vmin, vmax = float(np.nanmin(heatmap)), float(np.nanmax(heatmap))
    cmap = plt.get_cmap(CMAP_BY_FAMILY.get(r['family'], 'plasma')).copy(); cmap.set_bad('black')

    fig = plt.figure(figsize=(9, 8))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 14], left=0.25, right=0.95, bottom=0.10, top=0.86, hspace=0.15)
    cax = fig.add_subplot(gs[0]); ax = fig.add_subplot(gs[1])
    im = ax.imshow(heatmap, aspect='auto', cmap=cmap, interpolation='nearest', vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
    cbar.set_label(f'({unit})' if unit else '(unitless)', fontsize=24, labelpad=8)
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.tick_params(labelsize=20)
    unitTxt = f' ({unit})' if unit else ''
    rfTxt = '' if pd.isna(r['rfAccuracy']) else f'\nClassification Accuracy = {float(r["rfAccuracy"]):.3f}'
    cax.set_title(f'{label}{unitTxt}{rfTxt}', fontsize=26, pad=14)

    tickFrames = [t for t in frames if t >= 10 and t % 5 == 0]
    ax.set_xticks([frames.index(t) for t in tickFrames]); ax.set_xticklabels(tickFrames)
    ax.set_xlabel('Time (h)', fontsize=26)
    ax.set_yticks(range(len(mat.index))); ax.set_yticklabels([DISPLAY_NAMES[m] for m in mat.index])

    stem = panelStem(letter, feat, label)
    out = outDir / stem
    fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
    fig.savefig(str(out) + '.svg', bbox_inches='tight')
    plt.close(fig)
    manifest.append((letter, stem, label))
    print(f'Saved: {stem}.{{png,svg}}')

print(f'\n{len(manifest)} panels. Letter -> filename:')
for letter, stem, label in manifest:
    print(f'  {letter:>2s}: {stem}')
