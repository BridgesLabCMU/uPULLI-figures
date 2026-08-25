"""figS3 build: temporal (mutant x frame) matrices for ALL feature classes (source data for Fig S3).

The all-feature-class extension of fig2/build/2F_featmaps.py. For every quantitative feature class
(biomass + 14 whole-image + 12 colony bases), builds a mutant (row) x frame (col) matrix of the mean
value across replicates. Identical conventions to Fig 2F so the four features shared with Fig 2F match:
colony-feature timepoints with < minRepsPerTp replicates above the biomass threshold are masked NaN;
frames start at timeStart=8. Also writes a meta table (display label, unit, family, single-feature RF
accuracy from ACC_PATH) and a family-ordered manifest for single-page layout.

Reads:  config.input('training/wide.parquet') (training_wide.parquet); ACC_PATH (singleFeatureAccuracy.csv, optional)
Writes: data/featmap_<feat>.csv (x N)  +  data/featmaps_meta.csv
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS3/ for figlib
import numpy as np
import pandas as pd
from figlib import config, features, STRAIN_ORDER

ACC_PATH = ACC_PATH   # bundled Fig S2D output, not a deposit input

timeStart, biomassThresh, minRepsPerTp = 8, 5e-3, 10

# All feature classes, family-ordered for single-page layout: biomass, then whole-image, then colony.
FEATURES = ['biomass'] + list(features.WHOLE_FEATURE_BASES) + list(features.ALLOWED_COLONY_BASES)


def familyOf(feat):
    return 'biomass' if feat == 'biomass' else ('whole' if feat.startswith('whole_') else 'colony')


df = pd.read_parquet(config.input('training/wide.parquet'))
df = df[df['mutant'].isin(STRAIN_ORDER)].reset_index(drop=True)


def frameOf(c):
    return int(c.split('_t')[-1]) if '_t' in c and c.split('_t')[-1].isdigit() else None


allFrames = sorted({frameOf(c) for c in df.columns if frameOf(c) is not None})

# per-mutant low-biomass timepoint exclusion (colony features only)
excludedTp = {}
for m in STRAIN_ORDER:
    dfm = df[df['mutant'] == m]
    excl = set()
    for t in allFrames:
        col = f'biomass_t{t}'
        if col not in dfm.columns or int(np.sum(dfm[col].astype(float).values > biomassThresh)) < minRepsPerTp:
            excl.add(t)
    excludedTp[m] = excl

acc = {}
if ACC_PATH.exists():
    a = pd.read_csv(ACC_PATH)
    acc = dict(zip(a['featureBase'], a['balancedAccuracy']))
else:
    print(f'[WARN] {ACC_PATH} not found — RF accuracy omitted from meta')

config.ensure(config.TABLES)
metaRows = []
for feat in FEATURES:
    family = familyOf(feat)
    validFrames = [t for t in allFrames if t >= timeStart and f'{feat}_t{t}' in df.columns]
    if not validFrames:
        print(f'[skip] {feat}: no columns in wide table')
        continue
    rows = []
    for m in STRAIN_ORDER:
        dfm = df[df['mutant'] == m]
        traj = np.full(len(validFrames), np.nan)
        for i, t in enumerate(validFrames):
            if family == 'colony' and t in excludedTp[m]:
                continue
            vals = dfm[f'{feat}_t{t}'].astype(float).values
            if np.isfinite(vals).any():
                traj[i] = np.nanmean(vals)
        rows.append(traj)
    mat = pd.DataFrame(rows, index=STRAIN_ORDER, columns=validFrames)
    mat.to_csv(config.TABLES / f'featmap_{feat}.csv')
    metaRows.append({'feature': feat, 'label': features.pretty_name(feat), 'unit': features.feature_unit(feat),
                     'family': family, 'rfAccuracy': acc.get(feat, np.nan)})
    accTxt = f'{acc.get(feat):.3f}' if feat in acc else 'n/a'
    print(f'Saved: featmap_{feat}.csv  ({mat.shape[0]}x{mat.shape[1]}, RF acc {accTxt})')

pd.DataFrame(metaRows).to_csv(config.TABLES / 'featmaps_meta.csv', index=False)
print(f'Saved: featmaps_meta.csv  ({len(metaRows)} feature classes)')
