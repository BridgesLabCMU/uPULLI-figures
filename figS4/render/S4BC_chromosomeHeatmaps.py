#!/usr/bin/env python3
"""figS4 render: Panels B / C - per-chromosome biomass heatmaps of the reimaging-selected mutants.

The 157 reimaging-selected mutants (the Fig-3 set, non-WT), split by chromosome: (B) Chromosome I,
(C) Chromosome II. Trajectories are the reimaging feature data (per-gene mean, normalized
to the reimaging WT peak); phenotype is the binary reimaging classification (Dispersal Defect, else High
Biofilm). One row per mutant, one column per hour (8-30 h), wrapped across side-by-side columns. The
biomass color scale is 0-3 (as in the other figures; its WT normalization differs from S4A/S4D).

Row labels (large, black): the gene NAME if the locus has one, otherwise the locus number. Beside each
label is a marker encoding phenotype + functional annotation:
  * shape = phenotype: square = High Biofilm, circle = Dispersal Defect.
  * fill  = functional annotation color (Motility / O-Antigen / Polyamine / Biotin / Pyruvate /
    Vibriobactin, from the shared FUNCTION_COLORS); OPEN (black edge, no fill) if unannotated.

Build with build/S4BC_reimagingHits.py.
Reads:  data/reimHits_biomass_matrix.csv, data/reimHits_locus_meta.csv, data/reimagingGeneNames.csv
Writes: figures/S4B_ChromosomeI.{png,svg}, figures/S4C_ChromosomeII.{png,svg}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS4/ for figlib
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from figlib import config, plotting

plotting.setStyle(extra={'axes.linewidth': 1.2})

mat = pd.read_csv(config.TABLES / 'reimHits_biomass_matrix.csv', index_col=0)
meta = pd.read_csv(config.TABLES / 'reimHits_locus_meta.csv')   # 157 Fig-3 reimaging mutants; phenotype = High/Dispersal
allFrames = [int(c) for c in mat.columns]
START_H = 8                                                   # crop the empty early frames
keepIdx = [i for i, f in enumerate(allFrames) if f >= START_H]
frames = [allFrames[i] for i in keepIdx]                      # frames shown (8..30)
nFrames = len(frames)
cmap = plt.cm.RdYlBu_r.copy(); cmap.set_bad('white')
# biomass color scale 0-3 (matching the other figures); values above 3x WT saturate
VMAX = 3.0
print(f'reimaging-hit color scale: vmin=0, vmax={VMAX}')

# gene names (locus -> name) + functional-group membership (locus -> color)
FUNC_COLORS = plotting.FUNCTION_COLORS
loc2group = {loc: grp for grp, loci in plotting.HIGHLIGHT_SETS.items() for loc in loci}
ri = pd.read_csv(config.TABLES / 'reimagingGeneNames.csv').dropna(subset=['geneLocus'])
nameMap = dict(zip(ri['geneLocus'], ri['geneName'].astype(str)))


def dispName(locus):
    n = nameMap.get(locus, '')
    return n.strip() if (isinstance(n, str) and n.strip() and n.strip() != str(locus)) else str(locus)


# layout (inches). Cells are rectangular (narrower than tall) to compress width.
CELL_W = 0.22          # cell width (per timepoint)
CELL_H = 0.36          # cell height / row height (sized for the gene-name font)
FONT = 20              # gene-name label size
MARKER_S = 170         # phenotype/function marker size
ROWS_PER_COL = 46
LEFT_LABEL_W = 2.7     # room for gene names + marker
COL_GAP = 0.9
RIGHT_MARGIN = 0.5
TOP = 4.7              # title + colorbar + phenotype/function legends (stacked)
BOTTOM = 1.2           # x-axis
MARK_AX_W = 0.42       # inch: dedicated axis to the RIGHT of each heatmap holding the marker
MARK_GAP = 0.06        # inch: gap between the heatmap and its marker axis
TEXT_GAP = 0.16        # inch: gene-name right edge sits this far left of the heatmap


def renderChromosome(chrom, stem, title):
    sub = meta[meta['chromosome'] == chrom].reset_index(drop=True)   # all 157 are reimaging-selected hits
    loci = sub['geneLocus'].tolist()
    phenos = sub['phenotype'].tolist()
    n = len(loci)
    if n == 0:
        print(f'[skip] {chrom}: no hits'); return
    nCols = max(1, -(-n // ROWS_PER_COL))
    rowsPerCol = -(-n // nCols)
    blockW = nFrames * CELL_W
    maxBlockH = rowsPerCol * CELL_H
    colFoot = LEFT_LABEL_W + blockW + MARK_AX_W          # label | heatmap | marker axis
    figW = nCols * colFoot + (nCols - 1) * COL_GAP + RIGHT_MARGIN
    figH = TOP + maxBlockH + BOTTOM
    fig = plt.figure(figsize=(figW, figH))
    xText = -TEXT_GAP / blockW                            # gene name on the left (axes fraction)
    tickFrames = [f for f in frames if f % 5 == 0]
    tickPos = [frames.index(f) for f in tickFrames]
    im = None

    for ci in range(nCols):
        cl = loci[ci * rowsPerCol:(ci + 1) * rowsPerCol]
        cp = phenos[ci * rowsPerCol:(ci + 1) * rowsPerCol]
        nr = len(cl)
        if nr == 0:
            continue
        H = np.ma.masked_invalid(mat.loc[cl].values[:, keepIdx].astype(float))
        colOrigin = ci * (colFoot + COL_GAP)
        y0, h = BOTTOM / figH, (nr * CELL_H) / figH
        ax = fig.add_axes([(colOrigin + LEFT_LABEL_W) / figW, y0, blockW / figW, h])
        im = ax.imshow(H, aspect='auto', cmap=cmap, vmin=0, vmax=VMAX, interpolation='nearest')
        ax.set_xticks(tickPos); ax.set_xticklabels(tickFrames, fontsize=14)
        ax.set_xlabel('Time (h)', fontsize=17)
        ax.set_yticks([])
        # dedicated marker axis to the right (kept inside the layout so it is never clipped)
        axM = fig.add_axes([(colOrigin + LEFT_LABEL_W + blockW + MARK_GAP) / figW, y0,
                            (MARK_AX_W - 2 * MARK_GAP) / figW, h])
        axM.set_xlim(0, 1); axM.set_ylim(nr - 0.5, -0.5); axM.axis('off')
        tr = ax.get_yaxis_transform()
        for j, (locus, ph) in enumerate(zip(cl, cp)):
            grp = loc2group.get(locus)
            face = FUNC_COLORS[grp] if grp else 'none'
            shape = 's' if ph == 'High Biofilm' else 'o'
            ax.text(xText, j, dispName(locus), transform=tr, ha='right', va='center',
                    fontsize=FONT, color='black', clip_on=False)
            axM.scatter([0.5], [j], marker=shape, s=MARKER_S,
                        facecolors=face, edgecolors='black', linewidths=1.4, zorder=5)

    # ── top region (stacked, centered): colorbar, functional legend, phenotype legend, title ──
    def yfrac(inchesAbove):
        return (BOTTOM + maxBlockH + inchesAbove) / figH

    cbW = min(3.5, figW * 0.4)
    cax = fig.add_axes([(figW - cbW) / 2 / figW, yfrac(0.85), cbW / figW, 0.16 / figH])
    cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
    cbar.set_label('Biofilm Biomass (a.u.)', fontsize=15, labelpad=4)
    cbar.ax.xaxis.set_label_position('top'); cbar.ax.xaxis.set_ticks_position('top'); cbar.ax.tick_params(labelsize=12)

    phenoH = [Line2D([0], [0], marker='s', linestyle='none', markerfacecolor='none', markeredgecolor='black',
                     markersize=14, markeredgewidth=1.4, label='High Biofilm'),
              Line2D([0], [0], marker='o', linestyle='none', markerfacecolor='none', markeredgecolor='black',
                     markersize=14, markeredgewidth=1.4, label='Dispersal Defect')]
    funcH = [Patch(facecolor=FUNC_COLORS[g], edgecolor='black', label=g) for g in FUNC_COLORS]
    legF = fig.legend(handles=funcH, loc='center', bbox_to_anchor=(0.5, yfrac(2.05)), ncol=3,
                      title='Functional annotation (marker fill; open = none)', frameon=False,
                      fontsize=13, title_fontsize=14)
    fig.add_artist(legF)
    fig.legend(handles=phenoH, loc='center', bbox_to_anchor=(0.5, yfrac(3.1)), ncol=2,
               title='Phenotype (marker shape)', frameon=False, fontsize=13, title_fontsize=14)
    fig.text(0.5, yfrac(3.95), title, ha='center', fontsize=26, fontweight='bold')

    out = config.ensure(config.FIGURES) / stem
    fig.savefig(str(out) + '.png', dpi=300, bbox_inches='tight')
    fig.savefig(str(out) + '.svg', bbox_inches='tight')
    plt.close(fig)
    print(f'Saved: {stem}.{{png,svg}}  ({n} reimaging hits, {nCols} cols x ~{rowsPerCol} rows)')


renderChromosome('I', 'S4B_ChromosomeI', 'Chromosome I — Reimaging-Selected Mutants')
renderChromosome('II', 'S4C_ChromosomeII', 'Chromosome II — Reimaging-Selected Mutants')
