#!/usr/bin/env python3
"""figS8 render: Panel S8E - peak PvpsL-lux activity per replicate, normalized to WT.

Renders FROM the bundled peak table: one dot per well (deterministic jitter, seed 42) with a wide
horizontal bar at each condition's MEDIAN. Colors are the same per-strain colors as every other Fig-4
panel (Fig 4, Fig 5, S8B-D): WT black, Δ*bioD* orange, Δ*manA* blue, Δ*pdhE2* green.

House conventions: Gillius ADF + stixsans, square axes, axis-label fs 32 / tick fs 28, PNG + SVG
@300 dpi, no title. A dotted line marks WT = 1 (the normalization reference).

Reads:  data/luxPeak_normWT.csv
Writes: figures/S8E_luxActivity.{png,svg}

Usage:
  python figS8/render/S8E_luxActivity.py
  python figS8/render/S8E_luxActivity.py --conditions WT,BioD,ManA,PdhE2,PdhE1   # include pdhE1
  python figS8/render/S8E_luxActivity.py --center mean
"""
import sys
import argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS8/ for figlib
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from figlib import config, plotting, LUX_COLORS, LUX_DISPLAY

ap = argparse.ArgumentParser()
ap.add_argument('--conditions', default='WT,BioD,ManA,PdhE2',
                help='comma-separated, in x-axis order (PdhE1 is in the data but off by default)')
ap.add_argument('--center', choices=['median', 'mean'], default='median')
ap.add_argument('--jitter', type=float, default=0.13)
args = ap.parse_args()

plotting.setStyle()
SEED = 42
order = [c.strip() for c in args.conditions.split(',') if c.strip()]

d = pd.read_csv(config.TABLES / 'luxPeak_normWT.csv')
missing = [c for c in order if c not in set(d['condition'])]
if missing:
    raise SystemExit(f'conditions absent from luxPeak_normWT.csv: {missing}')

rng = np.random.default_rng(SEED)
fig, ax = plt.subplots(figsize=(10, 10))
ax.set_box_aspect(1)
ax.axhline(1.0, ls=':', lw=2, color='#999999', zorder=1)     # WT normalization reference

for i, cond in enumerate(order):
    y = d.loc[d['condition'] == cond, 'peakNormWT'].to_numpy(dtype=float)
    color = LUX_COLORS[cond]
    ax.scatter(i + rng.uniform(-args.jitter, args.jitter, size=len(y)), y, s=200, facecolors=color,
               edgecolors='black', linewidths=1.2, alpha=0.9, zorder=3)
    center = np.median(y) if args.center == 'median' else np.mean(y)
    ax.plot([i - 0.32, i + 0.32], [center, center], lw=5, color='black', solid_capstyle='butt', zorder=4)
    print(f'{cond:6s} n={len(y):2d}  {args.center}={center:.2f}  range=[{y.min():.2f}, {y.max():.2f}]')

ax.set_xticks(range(len(order)))
ax.set_xticklabels([LUX_DISPLAY[c] for c in order], fontsize=30)
ax.set_xlim(-0.6, len(order) - 0.4)
ax.set_ylim(bottom=0)
ax.set_ylabel('Normalized P$_{\\mathit{vpsL}}$-lux Activity', fontsize=32)
ax.tick_params(axis='y', labelsize=28)
ax.tick_params(axis='x', length=0)
fig.tight_layout()
out = config.ensure(config.FIGURES) / 'S8E_luxActivity'
fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
fig.savefig(str(out) + '.svg', bbox_inches='tight')
plt.close(fig)
print(f'Saved: {out}.png')
