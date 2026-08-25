#!/usr/bin/env python3
"""fig2 build step 0: training collapsed-wide table from the ORIGINAL full feature set.

Turns the per-(well, frame) training master feature set (`master_frame_features.csv`, µPULLI output,
a KiltHub deposit) into the wide table the other build scripts consume, and recovers the
(plateId, wellId) -> mutant labels from the legacy cleaned table.

Upstream of this (NOT in this repo): raw images -> µPULLI image pipeline -> processed images/masks
(BioImage Archive); processed images -> µPULLI feature extraction -> master_frame_features.csv.

Inputs (see inputs.json):
  training/master_frame_features.csv   per-(plate, well, frame) features
  training/layout.csv                  plateId, wellId -> mutant
Output:
  training/wide.parquet                consumed by build/2*.py
  data/trainingLayout.csv              bundled copy of the layout, for the record
"""
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig2/ for figlib
import pandas as pd
from figlib import config

parquetOut = Path(config.ensure(config.input_path('training/wide.parquet')))
layoutOut = config.TABLES / 'trainingLayout.csv'
config.ensure(config.TABLES)


def renameCol(c):
    if c == 'colAgg_nColonies':
        return 'nColonies'
    if c.startswith('colAgg_'):
        c = 'colony_' + c[len('colAgg_'):]
    if c.endswith('_skewness'):
        c = c[:-len('_skewness')] + '_skew'
    return c


def stripMag(w):
    m = re.match(r'^([A-Z]\d+)(?:_(\d+))?$', str(w))
    return m.group(1) if m else str(w)


masterFrame = config.input('training/master_frame_features.csv')
print(f'Loading {masterFrame}')
frame = pd.read_csv(masterFrame)
frame['plateId'] = frame['plateID'].astype(str).str.replace(' ', '_', regex=False)  # 'Plate 1' -> 'Plate_1'
frame['wellId'] = frame['wellID'].apply(stripMag)
frame = frame[frame['frame'].between(0, 30)].copy()

idCols = ['drawerID', 'plateID', 'wellID', 'mag', 'plateId', 'wellId', 'frame']
featCols = [c for c in frame.columns if c not in idCols]
long = frame[['plateId', 'wellId', 'frame'] + featCols].rename(columns={c: renameCol(c) for c in featCols})
long = long.melt(id_vars=['plateId', 'wellId', 'frame'], var_name='feature', value_name='value')
long['feature'] = long['feature'] + '_t' + long['frame'].astype(int).astype(str)
wide = long.pivot_table(index=['plateId', 'wellId'], columns='feature', values='value', aggfunc='first').reset_index()
wide.columns.name = None

lab = pd.read_csv(config.input('training/layout.csv'))[['plateId', 'wellId', 'mutant']].drop_duplicates(['plateId', 'wellId'])
lab.to_csv(layoutOut, index=False)
merged = wide.merge(lab, on=['plateId', 'wellId'], how='left')
front = ['plateId', 'wellId', 'mutant']
merged = merged[front + [c for c in merged.columns if c not in front]]
merged.to_parquet(parquetOut, index=False)
print(f'Saved: {parquetOut}  shape={merged.shape} | labeled {merged.mutant.notna().sum()}/{len(merged)}')
print(f'Saved: {layoutOut}')
