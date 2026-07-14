#!/usr/bin/env python3
"""fig2 build 2D: training UMAP coordinates for the all-features view (source data for Panel 2D).

Growth-filters the training wide table, standardizes the three feature modalities over frames 9-27
(biomass log1p, whole-image haralick/entropy, the 12 colony bases), concatenates all three, and fits
UMAP (n_neighbors=25, min_dist=0.25, random_state=0). Emits per-well coordinates + mutant label.

Reads:  config.WIDE (training_wide.parquet)
Writes: data/trainingUmap_all_three_coords.csv
"""
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig2/ for figlib
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
import umap.umap_ as umap
from figlib import config, features, STRAIN_ORDER

FRAMES, NN, MD, RS = features.KEEP_FRAMES, 25, 0.25, 0
df = pd.read_parquet(config.WIDE)

bcols = [c for c in df.columns if re.match(r'^biomass_t\d+$', c)]
maxBio = df[bcols].max(axis=1)
wtMed = (pd.DataFrame({'p': df['plateId'], 'b': maxBio, 'm': df['mutant']})
         .query("m=='WT'").groupby('p')['b'].median())
keep = (df['mutant'] == 'vpsL') | (maxBio >= 0.15 * df['plateId'].map(wtMed))
df = df[keep & df['mutant'].isin(STRAIN_ORDER)].reset_index(drop=True)
print(f'{len(df)} wells after growth filter')


def cols_for(bases):
    return [f'{b}_t{t}' for t in FRAMES for b in bases if f'{b}_t{t}' in df.columns]


def std(cols, log=False):
    X = df[cols].astype(float)
    X = X.fillna(X.median())
    if log:
        X = np.log1p(X)
    return StandardScaler().fit_transform(X.values)


Xbio = std(cols_for(['biomass']), log=True)
Xwhole = std(cols_for([f'whole_haralick_{i}' for i in range(13)] + ['whole_entropy']))
Xcol = std(cols_for(features.ALLOWED_COLONY_BASES))
emb = umap.UMAP(n_neighbors=NN, min_dist=MD, random_state=RS).fit_transform(np.hstack([Xbio, Xwhole, Xcol]))

out = df[['plateId', 'wellId', 'mutant']].copy()
out['umap1'], out['umap2'] = emb[:, 0], emb[:, 1]
config.ensure(config.TABLES)
out.to_csv(config.TABLES / 'trainingUmap_all_three_coords.csv', index=False)
print(f'Saved: trainingUmap_all_three_coords.csv ({len(out)} wells)')
