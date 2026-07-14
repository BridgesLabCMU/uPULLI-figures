"""fig2 build 2F: single-feature mutant x frame matrices (source data for Panel 2F).

For each of the four featured single features, builds a mutant (row) x frame (col) matrix of the median
value across replicates (colony features: timepoints with < minRepsPerTp replicates above the biomass
threshold are masked NaN). Also writes a meta table with the display label, unit, group, and the
single-feature RF balanced accuracy (from config.ACC) used in each panel's title.

Reads:  config.WIDE (training_wide.parquet); config.ACC (singleFeatureAccuracy.csv, optional)
Writes: data/featmap_<feat>.csv (x4)  +  data/featmaps_meta.csv
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig2/ for figlib
import numpy as np
import pandas as pd
from figlib import config, features, STRAIN_ORDER

FEATURES = ['colony_meanIntensity_mean', 'nColonies', 'colony_eccentricity_mean', 'whole_haralick_8']
timeStart, biomassThresh, minRepsPerTp = 8, 5e-3, 10

df = pd.read_parquet(config.WIDE)
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
if config.ACC.exists():
    a = pd.read_csv(config.ACC)
    acc = dict(zip(a['featureBase'], a['balancedAccuracy']))
else:
    print(f'[WARN] {config.ACC} not found — RF accuracy omitted from meta')

config.ensure(config.TABLES)
metaRows = []
for feat in FEATURES:
    group = 'whole' if feat.startswith('whole_') else 'colony'
    validFrames = [t for t in allFrames if t >= timeStart and f'{feat}_t{t}' in df.columns]
    rows = []
    for m in STRAIN_ORDER:
        dfm = df[df['mutant'] == m]
        traj = np.full(len(validFrames), np.nan)
        for i, t in enumerate(validFrames):
            if group == 'colony' and t in excludedTp[m]:
                continue
            vals = dfm[f'{feat}_t{t}'].astype(float).values
            if np.isfinite(vals).any():
                traj[i] = np.nanmean(vals)
        rows.append(traj)
    mat = pd.DataFrame(rows, index=STRAIN_ORDER, columns=validFrames)
    mat.to_csv(config.TABLES / f'featmap_{feat}.csv')
    metaRows.append({'feature': feat, 'label': features.pretty_name(feat), 'unit': features.feature_unit(feat),
                     'group': group, 'rfAccuracy': acc.get(feat, np.nan)})
    print(f'Saved: featmap_{feat}.csv  ({mat.shape[0]}x{mat.shape[1]}, RF acc {acc.get(feat, float("nan")):.3f})')

pd.DataFrame(metaRows).to_csv(config.TABLES / 'featmaps_meta.csv', index=False)
print('Saved: featmaps_meta.csv')
