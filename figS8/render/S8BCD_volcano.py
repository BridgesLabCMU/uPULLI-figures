#!/usr/bin/env python3
"""figS8 render: Panels S8B/C/D - RNA-seq volcano plots for the three clean deletions vs WT.

Renders FROM the bundled per-mutant volcano tables. One square panel per mutant, house conventions
(Gillius ADF + stixsans, box_aspect 1, axis-label fs 32 / tick fs 28, PNG + SVG @300 dpi, no title).
Panels: S8B = Delta-bioD, S8C = Delta-pdhE2, S8D = Delta-manA.

Points past BOTH thresholds (|log2FC| > 2 and q < 0.05) are drawn in that mutant's own color -- the
same reimaging functional-group color it carries in Fig 4 and Fig 5 (Δ*bioD* orange,
Δ*manA* blue, Δ*pdhE2* green). Everything else is grey. Dashed guides mark both thresholds.

Reads:  data/rnaseq_volcano_<mutant>.csv
Writes: figures/S8{B,C,D}_volcano_<mutant>.{png,svg}

Usage:
  python figS8/render/S8BCD_volcano.py                      # all three panels
  python figS8/render/S8BCD_volcano.py --mutant BioD
  python figS8/render/S8BCD_volcano.py --fc 1 --q 0.01      # different thresholds
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
from figlib import config, plotting, CLEANDEL_COLORS, RNASEQ_PANELS

ap = argparse.ArgumentParser()
ap.add_argument('--mutant', default='all', choices=['all'] + list(RNASEQ_PANELS))
ap.add_argument('--fc', type=float, default=2.0, help='|log2FC| cutoff for coloring (default 2)')
ap.add_argument('--q', type=float, default=0.05, help='adjusted-p cutoff for coloring (default 0.05)')
args = ap.parse_args()

plotting.setStyle()
QFLOOR = 1e-12          # keeps -log10(0) off the plot if a q ever underflows


def yOf(df):
    return -np.log10(df['qvalue'].clip(lower=QFLOOR))


# Shared y-axis across S8B/C/D: computed from ALL THREE tables regardless of --mutant, so a single-panel
# re-render lands on the same scale as a full run and the three panels stay comparable side by side.
YMAX = max(yOf(pd.read_csv(config.TABLES / f'rnaseq_volcano_{m}.csv')).max() for m in RNASEQ_PANELS)
YMAX = float(np.ceil(YMAX * 1.05))
print(f'shared y-axis: 0 to {YMAX:g} (-log10 adjusted p)')


def render(mutant):
    panel, _ = RNASEQ_PANELS[mutant]
    color = CLEANDEL_COLORS[mutant]
    d = pd.read_csv(config.TABLES / f'rnaseq_volcano_{mutant}.csv')
    d['y'] = yOf(d)
    hit = (d['qvalue'] < args.q) & (d['logFC'].abs() > args.fc)

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_box_aspect(1)
    ax.scatter(d.loc[~hit, 'logFC'], d.loc[~hit, 'y'], s=55, c=plotting.BACKGROUND_COLOR,
               edgecolors='none', alpha=0.55, rasterized=True, zorder=2)
    ax.scatter(d.loc[hit, 'logFC'], d.loc[hit, 'y'], s=110, c=color, edgecolors='black',
               linewidths=0.6, alpha=0.95, zorder=3)
    for v in (-args.fc, args.fc):
        ax.axvline(v, ls='--', lw=2, color='#888888', zorder=1)
    ax.axhline(-np.log10(args.q), ls='--', lw=2, color='#888888', zorder=1)

    lim = np.ceil(np.abs(d['logFC']).max()) + 0.5      # symmetric x so up/down read comparably
    ax.set_xlim(-lim, lim)
    ax.set_ylim(0, YMAX)          # shared across the three panels

    ax.set_xlabel('log$_2$FC', fontsize=32)
    ax.set_ylabel('$-$log$_{10}$ adjusted $p$', fontsize=32)
    ax.tick_params(labelsize=28)
    fig.tight_layout()
    out = config.ensure(config.FIGURES) / f'{panel}_volcano_{mutant}'
    fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
    fig.savefig(str(out) + '.svg', dpi=400, bbox_inches='tight')   # dpi = rasterized grey cloud
    plt.close(fig)
    print(f'{panel} {mutant:6s} {int(hit.sum()):4d}/{len(d)} colored '
          f"(|log2FC|>{args.fc:g}, q<{args.q:g})  -> {out.name}.png/.svg")


for m in (list(RNASEQ_PANELS) if args.mutant == 'all' else [args.mutant]):
    render(m)
