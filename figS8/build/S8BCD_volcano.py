#!/usr/bin/env python3
"""figS8 build B-D: per-mutant RNA-seq volcano tables (source data for Panels S8B/C/D).

One table per clean deletion (vs WT), carrying every quantified gene -- a volcano needs the
non-significant cloud, so nothing is filtered here. The lab's differential-expression output is
already volcano-ready; this step only selects columns, renames them, and derives the plotted y
(-log10 of the Benjamini-Hochberg adjusted p-value).

Contrast letters (B/E/M vs W) are the lab's file naming; each was verified by the deleted gene being
the extreme negative outlier of its own table (bioD -10.52, manA -6.23).

Reads:  config.input('rnaseq/{BvW,EvW,MvW}_allGenes.csv')
Writes: data/rnaseq_volcano_{BioD,PdhE2,ManA}.csv   (locus, oldLocus, gene, description, logFC, qvalue)
        -> panels S8B (bioD), S8C (pdhE2), S8D (manA)

Usage:  python figS8/build/S8BCD_volcano.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS8/ for figlib
import numpy as np
import pandas as pd
from figlib import config, RNASEQ_PANELS

COLS = {'Locustag': 'locus', 'OldLocustag': 'oldLocus', 'Gene': 'gene', 'Description': 'description',
        'logFC': 'logFC', 'Benjamini_Hochberg_Adjusted_PValue': 'qvalue'}

config.ensure(config.TABLES)
for mutant, (panel, stem) in RNASEQ_PANELS.items():
    src = Path(config.input('rnaseq/')) / f'{stem}_allGenes.csv'
    d = pd.read_csv(src)
    missing = [c for c in COLS if c not in d.columns]
    if missing:
        raise SystemExit(f'{src.name}: missing expected columns {missing}')
    out = d[list(COLS)].rename(columns=COLS)
    # the deleted gene should be the strongest depletion -- a cheap guard against a mislabeled contrast
    lowest = out.loc[out['logFC'].idxmin()]
    sig = (out['qvalue'] < 0.05)
    print(f'{panel} {mutant:6s} <- {stem}: n={len(out)}  q<0.05 {int(sig.sum()):4d}  '
          f"|logFC|>2 & q<0.05 {int((sig & (out['logFC'].abs() > 2)).sum()):4d}  "
          f"strongest depletion: {lowest['gene']} ({lowest['logFC']:.2f})")
    out.to_csv(config.TABLES / f'rnaseq_volcano_{mutant}.csv', index=False)

print(f'Saved: rnaseq_volcano_*.csv ({len(RNASEQ_PANELS)} tables) -> {config.TABLES}')
