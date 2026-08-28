#!/usr/bin/env python3
"""fig3 render — Supplement: per-gene UMAP PDF (one landscape page per mutant).

Renders FROM the bundled coordinates table (which already carries each replicate's highlight color).
Each page = the full landscape (grey + faint WT) with one gene's replicates highlighted.

Reads:  data/reimagingUmap_nn10_md0.10_perGene_coords.csv
Writes: figures/reimagingUmap_perGene.pdf
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # -> fig3/ for figlib (archive/ is one level deeper)
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
from figlib import config, plotting

plotting.setStyle(extra={'font.size': 28, 'axes.linewidth': 2})

df = pd.read_csv(config.TABLES / 'reimagingUmap_nn10_md0.10_perGene_coords.csv')
for c in ('mutant', 'geneLocus', 'function', 'dotColor'):
    df[c] = df[c].fillna('')
wt = df['mutant'].astype(str) == 'WT'
pad = 0.5
xlim = (df.umap1.min() - pad, df.umap1.max() + pad)
ylim = (df.umap2.min() - pad, df.umap2.max() + pad)

genes = sorted(g for g in df.loc[~wt, 'mutant'].dropna().unique() if g)
out = config.ensure(config.FIGURES) / 'reimagingUmap_perGene.pdf'
print(f'Per-gene PDF: {len(genes)} genes -> {out}')

with PdfPages(out) as pdf:
    for gene in genes:
        gmask = df['mutant'] == gene
        other = ~gmask & ~wt
        fig, ax = plt.subplots(figsize=(15, 15)); ax.set_box_aspect(1)
        ax.scatter(df.loc[other, 'umap1'], df.loc[other, 'umap2'], c=plotting.BACKGROUND_COLOR,
                   s=100, alpha=0.15, edgecolors='none', rasterized=True, zorder=1)
        if wt.any():
            ax.scatter(df.loc[wt, 'umap1'], df.loc[wt, 'umap2'], c='black', s=100, alpha=0.25,
                       linewidth=0, rasterized=True, zorder=2)
        sub = df[gmask]
        locus = sub['geneLocus'].iloc[0] if len(sub) else ''
        dotColor = sub['dotColor'].iloc[0] if len(sub) else '#ff00c3'
        ax.scatter(sub.umap1, sub.umap2, c=dotColor, s=160, alpha=0.9, edgecolors='black',
                   linewidth=1.0, zorder=5, label=gene)
        ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
        ax.set_title(f'{gene} ({locus})' if locus and locus != gene else gene, fontsize=24)
        ax.legend(frameon=False, loc='upper right')
        plt.tight_layout()
        pdf.savefig(fig, dpi=150)
        plt.close(fig)
print(f'Saved: {out} ({len(genes)} pages)')
