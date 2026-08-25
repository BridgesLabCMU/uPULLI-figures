#!/usr/bin/env python3
"""figS4 build: Data S1 — full results of the initial genome-wide transposon screen.

The supplementary dataset behind the screen: every quantified feature for every transposon mutant,
carrying the phenotype call that selected the reimaging library.

Two files, because one table cannot serve both readers:

  DataS1_transposonScreen_perWell.csv        one row per well (~3,000). Identity + the screen's
                                             summary measures (Peak / Final / Early) + the phenotype
                                             call + every feature at that well's peak-biomass frame.
                                             This is the table to open and read.

  DataS1_transposonScreen_allFrames.csv.gz   one row per well x frame (~95,000). The same identity
                                             and label columns, plus every feature at every
                                             timepoint — the complete record, gzipped.

Phenotype calls are read from the screen analysis of record (`transposons/results_10x.csv`) and are
NOT recomputed: they are the criterion as it was applied when the reimaging library was chosen.

Written to `deposit/` rather than `data/`, because these are deposit artifacts (tens of MB) rather
than per-panel source data bundled with the code.

Reads:  config.input('transposons/master_frame_features.csv'), config.input('transposons/results_10x.csv')
        data/tn_geneIndex.csv, data/reimagingGeneNames.csv  (bundled)
Writes: ../deposit/DataS1_transposonScreen_perWell.csv
        ../deposit/DataS1_transposonScreen_allFrames.csv.gz

Usage:  python figS4/build/DataS1_transposonScreen.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS4/ for figlib
import numpy as np
import pandas as pd
from figlib import config

OUT = Path(__file__).resolve().parents[2] / 'deposit'
HIT_CLASSES = ['High Biofilm', 'Dispersal Defect']
IDS = ['plateLabel', 'plateID', 'wellId', 'geneLocus', 'geneName', 'chromosome']
CALLS = ['Peak', 'Final', 'Early', 'Phenotype', 'isHit']


def stripMag(s):
    return s.astype(str).str.replace(r'_\d+$', '', regex=True)


def normLocus(s):
    return s.astype(str).str.replace(r'\s*;\s*', '/', regex=True).str.strip()


def chromosomeOf(locus):
    return 'II' if str(locus).startswith(('VC_A', 'VCA')) else 'I'


# ── the calls of record ───────────────────────────────────────────────────────
res = pd.read_csv(config.input('transposons/results_10x.csv'))
calls = (res[['Plate', 'Well', 'Gene.Locus', 'Peak', 'Final', 'Early', 'Phenotype']]
         .drop_duplicates().rename(columns={'Well': 'wellId', 'Gene.Locus': 'geneLocus'}))
calls['geneLocus'] = normLocus(calls['geneLocus'])
calls['plateLabel'] = calls['Plate'].apply(lambda n: f'TN-Plate{int(n):02d}')
calls = calls.drop(columns='Plate')
calls['isHit'] = calls['Phenotype'].isin(HIT_CLASSES)
print(f'calls of record: {len(calls)} wells\n' + calls['Phenotype'].value_counts().to_string())

# ── features, joined to locus via the screen index ────────────────────────────
gene = pd.read_csv(config.TABLES / 'tn_geneIndex.csv')          # plateID, plateLabel, wellId, geneLocus
gene['geneLocus'] = normLocus(gene['geneLocus'])
feat = pd.read_csv(config.input('transposons/master_frame_features.csv'))
feat['wellId'] = stripMag(feat['wellID'])
featCols = [c for c in feat.columns if c not in ('drawerID', 'plateID', 'wellID', 'wellId', 'mag', 'frame')]
print(f'features: {len(feat)} rows x {len(featCols)} feature columns')

df = feat.merge(gene[['plateID', 'plateLabel', 'wellId', 'geneLocus']], on=['plateID', 'wellId'], how='inner')
df = df.merge(calls[['plateLabel', 'wellId', 'geneLocus'] + CALLS], on=['plateLabel', 'wellId', 'geneLocus'],
              how='left')
names = pd.read_csv(config.TABLES / 'reimagingGeneNames.csv')   # geneLocus -> geneName (named genes only)
df = df.merge(names, on='geneLocus', how='left')
df['chromosome'] = df['geneLocus'].map(chromosomeOf)
unlabelled = int(df['Phenotype'].isna().sum())
if unlabelled:
    print(f'[note] {unlabelled} well-frames have features but no call of record (kept, Phenotype empty)')

OUT.mkdir(parents=True, exist_ok=True)

# ── per well: everything at that well's peak-biomass frame ────────────────────
peakIdx = df.groupby(['plateID', 'wellId'])['biomass'].idxmax()
perWell = df.loc[peakIdx].rename(columns={'frame': 'peakFrame'})
perWell = perWell[IDS + CALLS + ['peakFrame'] + featCols].sort_values(['plateLabel', 'wellId'])
perWell.to_csv(OUT / 'DataS1_transposonScreen_perWell.csv', index=False)
print(f'perWell: {len(perWell)} wells x {perWell.shape[1]} cols  '
      f'({(OUT / "DataS1_transposonScreen_perWell.csv").stat().st_size / 1e6:.1f} MB)')

# ── all frames: the complete record ───────────────────────────────────────────
allFrames = df[IDS + CALLS + ['frame'] + featCols].sort_values(['plateLabel', 'wellId', 'frame'])
allFrames.to_csv(OUT / 'DataS1_transposonScreen_allFrames.csv.gz', index=False, compression='gzip')
print(f'allFrames: {len(allFrames)} rows x {allFrames.shape[1]} cols  '
      f'({(OUT / "DataS1_transposonScreen_allFrames.csv.gz").stat().st_size / 1e6:.1f} MB gzipped)')
