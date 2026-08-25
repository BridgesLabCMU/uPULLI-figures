#!/usr/bin/env python3
"""figS5 render: per-chromosome biomass heatmaps of the LOW-BIOFILM transposon mutants.

Panels A/B — the Low-Biofilm class from the genome-wide transposon screen (Fig S4), split by chromosome:
(A) Chromosome I, (B) Chromosome II. One row per mutant, one column per hour (8-30 h), biomass normalized
to the WT peak mean. Low-biofilm mutants were not reimaged, so none carry a gene name (rows are labeled
by locus number) or a functional annotation (no markers). Biomass color scale is 0-3 (as in the other
figures); low-biofilm mutants sit below the low-biofilm threshold (~0.19 a.u.) and appear near the dark end.

Reads:  data/tn_biomass_matrix.csv, data/tn_locus_meta.csv  (transposon-screen tables, from figS4 build)
Writes: figures/S5A_ChromosomeI.{png,svg}, figures/S5B_ChromosomeII.{png,svg}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS5/ for figlib
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from figlib import config, plotting

plotting.setStyle(extra={'axes.linewidth': 1.2})

mat = pd.read_csv(config.TABLES / 'tn_biomass_matrix.csv', index_col=0)
meta = pd.read_csv(config.TABLES / 'tn_locus_meta.csv')
allFrames = [int(c) for c in mat.columns]
START_H = 8
keepIdx = [i for i, f in enumerate(allFrames) if f >= START_H]
frames = [allFrames[i] for i in keepIdx]
nFrames = len(frames)
cmap = plt.cm.RdYlBu_r.copy(); cmap.set_bad('white')

low = meta[meta['phenotype'] == 'Low Biofilm'].reset_index(drop=True)
VMAX = 3.0   # biomass color scale 0-3, matching the other figures (low-biofilm mutants sit near the bottom)
print(f'Low-biofilm color scale: vmin=0, vmax={VMAX} (n={len(low)})')

# layout (inches)
CELL_W, CELL_H = 0.22, 0.36
FONT = 20
ROWS_PER_COL = 46
LEFT_LABEL_W, COL_GAP, RIGHT_MARGIN = 1.9, 0.8, 0.4
TOP, BOTTOM = 2.4, 1.2


def renderChromosome(chrom, stem, title):
    sub = low[low['chromosome'] == chrom].reset_index(drop=True)
    loci = sub['geneLocus'].tolist()
    n = len(loci)
    if n == 0:
        print(f'[skip] {chrom}: no Low-biofilm mutants'); return
    nCols = max(1, -(-n // ROWS_PER_COL))
    rowsPerCol = -(-n // nCols)
    blockW, maxBlockH = nFrames * CELL_W, rowsPerCol * CELL_H
    colFoot = LEFT_LABEL_W + blockW
    figW = nCols * colFoot + (nCols - 1) * COL_GAP + RIGHT_MARGIN
    figH = TOP + maxBlockH + BOTTOM
    fig = plt.figure(figsize=(figW, figH))
    tickFrames = [f for f in frames if f % 5 == 0]
    tickPos = [frames.index(f) for f in tickFrames]
    im = None

    for ci in range(nCols):
        cl = loci[ci * rowsPerCol:(ci + 1) * rowsPerCol]
        nr = len(cl)
        if nr == 0:
            continue
        H = np.ma.masked_invalid(mat.loc[cl].values[:, keepIdx].astype(float))
        colOrigin = ci * (colFoot + COL_GAP)
        ax = fig.add_axes([(colOrigin + LEFT_LABEL_W) / figW, BOTTOM / figH, blockW / figW, (nr * CELL_H) / figH])
        im = ax.imshow(H, aspect='auto', cmap=cmap, vmin=0, vmax=VMAX, interpolation='nearest')
        ax.set_xticks(tickPos); ax.set_xticklabels(tickFrames, fontsize=14)
        ax.set_xlabel('Time (h)', fontsize=17)
        ax.set_yticks(range(nr)); ax.set_yticklabels(cl, fontsize=FONT)
        ax.tick_params(axis='y', length=2)

    cbW = min(3.2, figW * 0.4)
    cax = fig.add_axes([(figW - cbW) / 2 / figW, (BOTTOM + maxBlockH + 0.75) / figH, cbW / figW, 0.16 / figH])
    cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
    cbar.set_label('Biofilm Biomass (a.u.)', fontsize=15, labelpad=4)
    cbar.ax.xaxis.set_label_position('top'); cbar.ax.xaxis.set_ticks_position('top'); cbar.ax.tick_params(labelsize=12)
    fig.text(0.5, (BOTTOM + maxBlockH + 1.55) / figH, title, ha='center', fontsize=26, fontweight='bold')

    out = config.ensure(config.FIGURES) / stem
    fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
    fig.savefig(str(out) + '.svg', bbox_inches='tight')
    plt.close(fig)
    print(f'Saved: {stem}.{{png,svg}}  ({n} Low-biofilm mutants, {nCols} cols)')


renderChromosome('I', 'S5A_ChromosomeI', 'Chromosome I — Low-Biofilm Mutants')
renderChromosome('II', 'S5B_ChromosomeII', 'Chromosome II — Low-Biofilm Mutants')
