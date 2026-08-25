"""fig4 build 4E: RNA-seq log2 fold-change matrix for the three clean-deletion mutants.

Port of the lab RNA-seq heatmap (RNAseq_Heatmap_V1.py + TSV_CSV_addColumn.py). Reads the per-mutant
"all genes" differential-expression tables (mutant vs WT), pulls the curated biofilm-gene set, and
builds a mutant x gene log2FC matrix ordered ΔbioD (top), ΔpdhE2, ΔmanA (bottom) — matching the
reordered Fig-4C feature heatmap.

Reads:  rnaseq/{BvW,EvW,MvW}_allGenes.csv  (build input; see inputs.json)
Writes: data/rnaseq_logFC_matrix.csv   (mutant x gene, log2FC)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig4/ for figlib
import pandas as pd
from figlib import config

RNASEQ_DIR = config.input('rnaseq/')      # directory of <strain>_allGenes.csv differential-expression tables

# curated biofilm-gene set (order = x-axis order), from RNAseq_Heatmap_V1.py
GENES = ['vpsU', 'vpsA', 'vpsB', 'vpsC', 'vpsD', 'vpsE', 'vpsF', 'vpsG', 'vpsH', 'vpsI', 'vpsJ', 'vpsK',
         'vpsL', 'vpsM', 'vpsN', 'vpsO', 'vpsP', 'vpsQ',
         'rbmA', 'rbmB', 'rbmC', 'rbmD', 'rbmEF',
         'mshI', 'mshJ', 'mshK', 'mshL', 'mshM', 'mshN', 'mshE', 'mshG', 'mshF',
         'mshB', 'mshA', 'mshC', 'mshD', 'mshO', 'mshP', 'mshQ',
         'lapD', 'lapG', 'lapB', 'lapC', 'craA', 'frhA', 'bap1']

# genotype file stem -> mutant key (order: bioD top, pdhE2 middle, manA bottom)
GENOTYPES = [('BvW', 'bioD'), ('EvW', 'pdhE2'), ('MvW', 'manA')]

rows = {}
for stem, mut in GENOTYPES:
    df = pd.read_csv(RNASEQ_DIR / f'{stem}_allGenes.csv')
    df['Display_Name'] = df['Gene'].fillna(df['OldLocustag'])
    s = (df[df['Display_Name'].isin(GENES)].groupby('Display_Name')['logFC'].mean())
    rows[mut] = s
    print(f'{mut}: {s.notna().sum()} of {len(GENES)} genes')

mat = pd.DataFrame(rows).T                                  # mutant (rows) x gene (cols)
mat = mat.reindex(index=[m for _, m in GENOTYPES], columns=[g for g in GENES if g in mat.columns])
config.ensure(config.TABLES)
mat.to_csv(config.TABLES / 'rnaseq_logFC_matrix.csv')
print(f'Saved: rnaseq_logFC_matrix.csv  ({mat.shape[0]} mutants x {mat.shape[1]} genes)')
