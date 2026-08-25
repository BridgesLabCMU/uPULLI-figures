#!/usr/bin/env python3
"""fig3 build step 0: collapsed-wide table from the ORIGINAL full feature set.

First step of the build layer: turns the per-(well, frame) master feature set produced by the µPULLI
feature-extraction pipeline (`master_frame_features.csv`, a KiltHub deposit) into the analysis-ready
wide table the other `build/` scripts consume. Logic (unchanged from the analysis pipeline):
  * strip the `_<mag>` suffix from wellID (`A10_02` -> `A10`)
  * canonical renames: colAgg_* -> colony_*, *_skewness -> *_skew, colAgg_nColonies -> nColonies
  * keep frames 0-30, pivot to `<feature>_t<frame>`, one row per (plateId, wellId)
  * join the reimaging gene index for geneLocus/geneName/function; derive `mutant`

Upstream of this (NOT in this repo): raw images -> µPULLI image pipeline -> processed images/masks
(BioImage Archive); processed images -> µPULLI feature extraction (multiWellAnalysis + DINOv2) ->
master_frame_features.csv. See the top-level README / the pipeline repository.

Inputs (see inputs.json):
  reimaging/master_frame_features.csv  per-(plate, well, frame) features
  reimaging/geneIndex.csv              well -> geneLocus / geneName / function
Output:
  reimaging/collapsedWide.parquet      consumed by build/3*.py
"""
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig3/ for figlib
import pandas as pd
from figlib import config

masterFramePath = config.input('reimaging/master_frame_features.csv')
indexPath = config.input('reimaging/geneIndex.csv')
parquetOut = Path(config.ensure(config.input_path('reimaging/collapsedWide.parquet')))
framesExpected = list(range(31))


def renameMasterColumn(col):
    if col == 'colAgg_nColonies':
        return 'nColonies'
    if col.startswith('colAgg_'):
        col = 'colony_' + col[len('colAgg_'):]
    if col.endswith('_skewness'):
        col = col[:-len('_skewness')] + '_skew'
    return col


def stripMagSuffix(wellToken):
    m = re.match(r'^([A-Z]\d+)(?:_(\d+))?$', str(wellToken))
    if m is None:
        return wellToken, None
    return m.group(1), m.group(2)


print(f'Loading master frame features: {masterFramePath}')
frame = pd.read_csv(masterFramePath)
print(f'master_frame_features shape: {frame.shape}')

frame['plateId'] = frame['plateID']
parsed = frame['wellID'].apply(stripMagSuffix)
frame['wellId'] = parsed.str[0]
fieldRepFromMag = frame['mag'].astype(str).str.lstrip('_').replace('', '1').astype(int)
frame['fieldRep'] = parsed.str[1].fillna(fieldRepFromMag).astype(int)

frame = frame[frame['frame'].between(min(framesExpected), max(framesExpected))].copy()
idCols = ['drawerID', 'plateID', 'wellID', 'mag', 'plateId', 'wellId', 'fieldRep', 'frame']
masterFeatureCols = [c for c in frame.columns if c not in idCols]
renameMap = {c: renameMasterColumn(c) for c in masterFeatureCols}
renamedFeatures = [renameMap[c] for c in masterFeatureCols]

print('Pivoting wide ...')
long = frame[['plateId', 'wellId', 'frame'] + masterFeatureCols].rename(columns=renameMap)
long = long.melt(id_vars=['plateId', 'wellId', 'frame'], value_vars=renamedFeatures,
                 var_name='feature', value_name='value')
long['feature'] = long['feature'] + '_t' + long['frame'].astype(int).astype(str)
wide = long.pivot_table(index=['plateId', 'wellId'], columns='feature', values='value',
                        aggfunc='first').reset_index()
wide.columns.name = None

print(f'Joining reimaging index: {indexPath}')
idx = pd.read_csv(indexPath, usecols=['plateId', 'wellId', 'geneLocus', 'geneName', 'function', 'srcPlate', 'srcWell'])
merged = wide.merge(idx.drop_duplicates(subset=['plateId', 'wellId']), on=['plateId', 'wellId'], how='left')
nMissing = merged['geneLocus'].isna().sum()
if nMissing:
    print(f'[WARN] {nMissing} wells have no index entry — dropping.')
    merged = merged[merged['geneLocus'].notna()].reset_index(drop=True)

merged['mutant'] = merged['geneName'].where(merged['geneName'].notna(), merged['geneLocus'])
frontCols = ['plateId', 'wellId', 'mutant', 'geneName', 'geneLocus', 'function', 'srcPlate', 'srcWell']
merged = merged[frontCols + [c for c in merged.columns if c not in frontCols]]

merged.to_parquet(parquetOut, index=False)
print(f'Saved: {parquetOut}  {merged.shape}')
