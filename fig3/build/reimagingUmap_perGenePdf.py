#!/usr/bin/env python3
"""fig3 copy — v2 reimaging per-gene UMAP PDF (one page per mutant) on the v2 manifold.

Copy of v2/reimaging/umap/reimagingUmap_perGenePdf.py for the figure package, with the import path
fixed for fig3/ and a per-replicate COORDINATES CSV added (for publication: every plotted dot's
umap position + gene/locus/function + the per-gene highlight color).

Reads:  data/v2/reimaging/reimaging_umapEmbeddings.parquet   (--embeddings to override)
Writes: results/v2/reimaging/umaps/reimagingUmap_nn{nn}_md{md}_perGene.pdf
        results/v2/reimaging/umaps/reimagingUmap_nn{nn}_md{md}_perGene_coords.csv
"""
import sys
import argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig3/ for figlib
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
from figlib import config, plotting

ap = argparse.ArgumentParser()
ap.add_argument('--embeddings', default=str(config.input('reimaging/umapEmbeddings.parquet')),
                help='UMAP embeddings parquet (standard schema). Defaults to the v2 numerical manifold.')
ap.add_argument('--outdir', default=str(config.FIGURES))
ap.add_argument('--nn', type=int, default=10)
ap.add_argument('--md', type=float, default=0.1)
args = ap.parse_args()

plotting.setStyle(extra={'font.size': 28, 'axes.linewidth': 2})

WT_LABEL = 'WT'
UNCLASSIFIED_DOT = '#ff00c3'
BACKGROUND_ALPHA = 0.15
locusToFunc = {l: f for f, loci in plotting.HIGHLIGHT_SETS.items() for l in loci}

outDir = Path(config.ensure(args.outdir))
pdfPath = outDir / f'reimagingUmap_nn{args.nn}_md{args.md:.2f}_perGene.pdf'
coordsCsv = config.TABLES / f'reimagingUmap_nn{args.nn}_md{args.md:.2f}_perGene_coords.csv'

emb = pd.read_parquet(args.embeddings)
emb = emb[(emb['n_neighbors'] == args.nn) & (emb['min_dist'] == args.md)].reset_index(drop=True)
for c in ('mutant', 'geneLocus', 'function'):
    emb[c] = emb[c].fillna('') if c in emb.columns else ''
print(f'Filtered to nn={args.nn}, md={args.md}: {len(emb)} rows')

wtMask = emb['mutant'].astype(str) == WT_LABEL
pad = 0.5
xlim = (emb['umap1'].min() - pad, emb['umap1'].max() + pad)
ylim = (emb['umap2'].min() - pad, emb['umap2'].max() + pad)


def dotColorFor(locus):
    return plotting.FUNCTION_COLORS.get(locusToFunc.get(locus, ''), UNCLASSIFIED_DOT)


# CSV: every plotted dot's coordinates + gene/locus/function + highlight color
coordsOut = emb[[c for c in ['plateId', 'wellId', 'mutant', 'geneLocus', 'function', 'umap1', 'umap2'] if c in emb.columns]].copy()
coordsOut['functionalGroup'] = coordsOut['geneLocus'].map(lambda l: locusToFunc.get(l, 'Unclassified'))
coordsOut.loc[wtMask.values, 'functionalGroup'] = 'WT'
coordsOut['dotColor'] = coordsOut['geneLocus'].map(dotColorFor)
coordsOut.loc[wtMask.values, 'dotColor'] = '#000000'
coordsOut.to_csv(coordsCsv, index=False)
print(f'Saved coordinates CSV: {coordsCsv} ({len(coordsOut)} replicates)')

uniqueGenes = sorted(g for g in emb.loc[~wtMask, 'mutant'].dropna().unique() if g)
print(f'Generating per-gene PDF: {len(uniqueGenes)} genes -> {pdfPath}')

with PdfPages(pdfPath) as pdf:
    for gene in uniqueGenes:
        geneMask = emb['mutant'] == gene
        otherMask = ~geneMask & ~wtMask

        fig, ax = plt.subplots(figsize=(15, 15))
        ax.set_box_aspect(1)
        ax.scatter(emb.loc[otherMask, 'umap1'], emb.loc[otherMask, 'umap2'],
                   c=plotting.BACKGROUND_COLOR, s=100, alpha=BACKGROUND_ALPHA,
                   edgecolors='none', linewidth=0, rasterized=True, zorder=1)
        if wtMask.any():
            ax.scatter(emb.loc[wtMask, 'umap1'], emb.loc[wtMask, 'umap2'],
                       c='black', s=100, alpha=0.25, linewidth=0, rasterized=True, zorder=2)

        locus = emb.loc[geneMask, 'geneLocus'].iloc[0] if geneMask.any() else ''
        dotColor = dotColorFor(locus)

        if geneMask.any():
            ax.scatter(emb.loc[geneMask, 'umap1'], emb.loc[geneMask, 'umap2'],
                       c=dotColor, s=160, alpha=0.9, edgecolors='black', linewidth=1.0,
                       zorder=5, label=gene)

        ax.set_xlim(xlim); ax.set_ylim(ylim)
        ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
        title = f'{gene} ({locus})' if locus and locus != gene else gene
        ax.set_title(title, fontsize=24)
        ax.legend(frameon=False, loc='upper right')
        plt.tight_layout()
        pdf.savefig(fig, dpi=150)
        plt.close(fig)

print(f'Saved: {pdfPath} ({len(uniqueGenes)} pages)')
