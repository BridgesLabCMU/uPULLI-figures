#!/usr/bin/env python3
"""fig2 build 2C: peak-biomass-frame feature matrix (source data for Panel 2C).

Per mutant, samples every feature (biomass + whole entropy/haralick + the 12 colony bases) at that
mutant's peak-biomass frame, medianed across replicates, then z-scored across the 8 mutants per feature
row. Fixed mutant order (STRAIN_ORDER); rows sorted into the canonical feature order.

Reads:  config.WIDE (training_wide.parquet)
Writes: data/trainingHeatmap_peakBiomass_featureMatrix.csv  +  _peakFrames.csv
"""
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig2/ for figlib
import numpy as np
import pandas as pd
from figlib import config, features, STRAIN_ORDER

canonHaralick = ['energy', 'contrast', 'correlation', 'variance', 'inverse_difference_moment', 'sum_average',
                 'sum_variance', 'sum_entropy', 'entropy', 'difference_variance', 'difference_entropy', 'imc1', 'imc2']


def featGroup(b):
    return 'biomass' if b == 'biomass' else 'haralick' if 'haralick' in b else 'entropy' if 'entropy' in b else 'colony'


def sortKey(b):
    if b == 'biomass':
        return (0, 0)
    if featGroup(b) == 'colony':
        return (1, b)
    if 'entropy' in b:
        return (2, 0)
    k = re.sub(r'_(mean|std|var)$', '', b.replace('whole_haralick_', ''))
    idx = int(k) if k.isdigit() else (canonHaralick.index(k) if k in canonHaralick else 99)
    return (3, idx)


df = pd.read_parquet(config.WIDE)
df = df[df['mutant'].isin(STRAIN_ORDER)].reset_index(drop=True)

bases = set()
for c in df.columns:
    m = re.search(r'_t(\d+)$', c)
    if not m:
        continue
    b = c[:m.start()]
    if b == 'biomass' or b.startswith('whole_haralick') or b.startswith('whole_entropy') or b in features.ALLOWED_COLONY_BASES:
        bases.add(b)
featureBases = sorted(bases, key=sortKey)

biomassCols = {int(re.search(r'_t(\d+)$', c).group(1)): c for c in df.columns if re.match(r'^biomass_t\d+$', c)}
peakFrames = {}
for m in STRAIN_ORDER:
    g = df[df['mutant'] == m]
    traj = pd.Series({fr: g[col].median() for fr, col in biomassCols.items()}).sort_index()
    peakFrames[m] = int(traj.idxmax()) if not traj.dropna().empty else np.nan

rows = []
for base in featureBases:
    vals = []
    for m in STRAIN_ORDER:
        g = df[df['mutant'] == m]
        pf = peakFrames[m]
        col = f'{base}_t{int(pf)}' if not pd.isna(pf) else None
        vals.append(g[col].median() if col and col in g.columns else np.nan)
    rows.append(vals)
mat = pd.DataFrame(rows, index=featureBases, columns=STRAIN_ORDER)
z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, np.nan), axis=0).fillna(0)

config.ensure(config.TABLES)
z.to_csv(config.TABLES / 'trainingHeatmap_peakBiomass_featureMatrix.csv')
pd.DataFrame({'mutant': list(peakFrames), 'peakFrame': list(peakFrames.values())}).to_csv(
    config.TABLES / 'trainingHeatmap_peakBiomass_peakFrames.csv', index=False)
print(f'Saved featureMatrix ({z.shape[0]} features x {z.shape[1]} mutants) + peakFrames')
