"""figS4 build: genome-wide transposon biomass screen (source data for Fig S4 + Fig S5).

Two sources, deliberately:

  * **Phenotype calls + thresholds** come from `TransposonResults_10x.csv`, the screen analysis of
    record. These calls are what selected the reimaging library, so they are a fixed historical
    result and are NOT recomputed here; they color the trajectories in S4A / S5. The per-well call
    table and the recovered B_max / B_min thresholds now belong to Fig 3A and are written by
    fig3/build/3A_rankedTransposonBiomass.py from the same source.
  * **Biomass trajectories** come from the screen's `master_frame_features.csv`, normalized to the
    WT peak mean of the training set. Panels S4A and S5 render these, colored by the calls above —
    the same split panels S4B/S4C already use (calls from the reimaging analysis, trajectories from
    the feature data), so their color scale differs from S4D's axis.

Reads:  transposons/results_10x.csv, transposons/master_frame_features.csv,
        training/master_frame_features.csv, training/layout.csv  (build inputs; see inputs.json)
        data/tn_geneIndex.csv  (bundled)
Writes: data/tn_biomass_matrix.csv      (geneLocus x frame, mean biomassNorm, genome-ordered)
        data/tn_locus_meta.csv          (geneLocus, chromosome, peak, phenotype, isHit; matrix row order)
        data/tn_normalization.csv       (trajectory normalizer + the screen's class counts)
"""
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS4/ for figlib
import numpy as np
import pandas as pd
from figlib import config

HIT_CLASSES = ['High Biofilm', 'Dispersal Defect']
RANK = {'High Biofilm': 3, 'Dispersal Defect': 2, 'Low Biofilm': 1, 'Normal': 0}


def stripMag(s):
    return s.astype(str).str.replace(r'_\d+$', '', regex=True)


def plateKey(s):
    return s.astype(str).str.replace(' ', '_')


def normLocus(s):
    """Multi-locus insertions are written 'VC_0255; VC_A0373' in the results table and
    'VC_0255/VC_A0373' in the screen index. Normalize both to the slash form."""
    return s.astype(str).str.replace(r'\s*[;/]\s*', '/', regex=True).str.strip()


def locusSortKey(locus):
    first = str(locus).split('/')[0]
    m = re.match(r'VC_A?(\d+)', first)
    return (1 if 'VC_A' in first else 0, int(m.group(1)) if m else 0)


def chromosome(locus):
    return 'II' if 'VC_A' in str(locus).split('/')[0] else 'I'


# ── calls of record: the screen analysis that selected the reimaging library ──
res = pd.read_csv(config.input('transposons/results_10x.csv'))
summary = (res[['Plate', 'Well', 'Gene.Locus', 'Peak', 'Final', 'Early', 'Phenotype']]
           .drop_duplicates()
           .rename(columns={'Well': 'well', 'Gene.Locus': 'geneLocus'}))
summary['geneLocus'] = normLocus(summary['geneLocus'])
summary['plateLabel'] = summary['Plate'].apply(lambda n: f'TN-Plate{int(n):02d}')
summary = summary.drop(columns='Plate').rename(columns={'plateLabel': 'plate'})
counts = summary['Phenotype'].value_counts()
print(f'Calls of record: {len(summary)} wells, {summary["geneLocus"].nunique()} loci')
print(counts.to_string())

# ── trajectories: screen feature table, normalized to the training-set WT peak mean ──
layout = pd.read_csv(config.input('training/layout.csv'))
layout = layout[layout['mutant'] == 'WT'][['plateId', 'wellId']].drop_duplicates()
layout['plateKey'] = plateKey(layout['plateId'])

train = pd.read_csv(config.input('training/master_frame_features.csv'),
                    usecols=['plateID', 'wellID', 'frame', 'biomass'])
train['plateKey'] = plateKey(train['plateID'])
train['well'] = stripMag(train['wellID'])
wt = train.merge(layout[['plateKey', 'wellId']], left_on=['plateKey', 'well'],
                 right_on=['plateKey', 'wellId'], how='inner')
if wt.empty:
    raise SystemExit('no WT wells matched between the training layout and the training master')
wtPeakMean = float(wt.groupby(['plateKey', 'well'])['biomass'].max().mean())
print(f'\nTrajectory normalizer: WT peak mean {wtPeakMean:.6f}')

gene = pd.read_csv(config.TABLES / 'tn_geneIndex.csv')      # plateID, plateLabel, wellId, geneLocus
gene['geneLocus'] = normLocus(gene['geneLocus'])
tn = pd.read_csv(config.input('transposons/master_frame_features.csv'),
                 usecols=['plateID', 'wellID', 'frame', 'biomass'])
tn['well'] = stripMag(tn['wellID'])
df = tn.merge(gene[['plateID', 'wellId', 'geneLocus']], left_on=['plateID', 'well'],
              right_on=['plateID', 'wellId'], how='inner')
df['biomassNorm'] = df['biomass'] / wtPeakMean

# keep the loci that carry a call — meta and matrix must stay row-aligned
called = set(summary['geneLocus'])
missing = called - set(df['geneLocus'])
extra = set(df['geneLocus']) - called
if missing:
    print(f'[warn] {len(missing)} called loci have no trajectory; dropped from the matrix')
if extra:
    print(f'[note] {len(extra)} loci have trajectories but no call; not shown in S4A/S5')
loci = sorted(called - missing, key=locusSortKey)

timepoints = np.sort(df['frame'].unique())
mat = (df.groupby(['geneLocus', 'frame'])['biomassNorm'].mean().unstack('frame')
       .reindex(index=loci, columns=timepoints))
print(f'Matrix: {len(loci)} loci x {len(timepoints)} frames')

# ── per-locus meta: calls of record, genome-ordered, aligned to the matrix ──
locusPheno = (summary.assign(_r=summary['Phenotype'].map(RANK))
              .sort_values('_r').drop_duplicates('geneLocus', keep='last')
              .set_index('geneLocus')['Phenotype'])
locusMeta = pd.DataFrame({'geneLocus': loci})
locusMeta['chromosome'] = locusMeta['geneLocus'].map(chromosome)
locusMeta['peak'] = mat.max(axis=1).values
locusMeta['phenotype'] = locusMeta['geneLocus'].map(locusPheno).fillna('Normal').values
locusMeta['isHit'] = locusMeta['phenotype'].isin(HIT_CLASSES).values
print('Per-locus calls:\n' + locusMeta['phenotype'].value_counts().to_string())

config.ensure(config.TABLES)
mat.to_csv(config.TABLES / 'tn_biomass_matrix.csv')
locusMeta.to_csv(config.TABLES / 'tn_locus_meta.csv', index=False)
# The per-well call table and the recovered B_max/B_min thresholds moved to fig3 with Panel 3A
# (fig3/build/3A_rankedTransposonBiomass.py) -- no S4 panel draws them. What stays here is the
# normalizer S4A's and S5's trajectories are expressed in, plus the screen's class counts.
pd.DataFrame([{'trajectoryNormalizer_wtPeakMean': wtPeakMean,
               'nWells': int(len(summary)), 'nLoci': int(len(loci)),
               'nHigh': int(counts.get('High Biofilm', 0)), 'nLow': int(counts.get('Low Biofilm', 0)),
               'nDispersal': int(counts.get('Dispersal Defect', 0)),
               'nNormal': int(counts.get('Normal', 0))}]).to_csv(
    config.TABLES / 'tn_normalization.csv', index=False)
print(f'Saved 3 tables to {config.TABLES}')
