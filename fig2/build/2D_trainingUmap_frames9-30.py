#!/usr/bin/env python3
"""fig2 build 2D (frames 9-30 variant): training UMAP over a single consistent frame window (9-30)
for all three feature modalities, colored by mutant.

This is the standard 2D UMAP (9-27) extended to the last frame: biomass(log1p) + whole-image
haralick/entropy + the 12 colony bases, all over frames 9-30. Start-frame 9 matches the colony
window, the canonical UMAP filter (9-27), and multispecies; end-frame 30 matches the 2E confusion
matrix. UMAP params (nn=25, md=0.25, rs=0) and growth filter identical to 2D_trainingUmap.py, which
this does NOT overwrite.

Reads:  config.input('training/wide.parquet') (training_wide.parquet)
Writes: data/trainingUmap_all_three_frames9-30_coords.csv
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

FRAMES, NN, MD, RS = list(range(9, 31)), 25, 0.25, 0   # single 9-30 window for all modalities
df = pd.read_parquet(config.input('training/wide.parquet'))

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
print(f'features (frames 9-30): biomass={Xbio.shape[1]}, whole={Xwhole.shape[1]}, colony={Xcol.shape[1]}')
emb = umap.UMAP(n_neighbors=NN, min_dist=MD, random_state=RS).fit_transform(np.hstack([Xbio, Xwhole, Xcol]))

out = df[['plateId', 'wellId', 'mutant']].copy()
out['umap1'], out['umap2'] = emb[:, 0], emb[:, 1]
config.ensure(config.TABLES)
out.to_csv(config.TABLES / 'trainingUmap_all_three_frames9-30_coords.csv', index=False)
print(f'Saved: trainingUmap_all_three_frames9-30_coords.csv ({len(out)} wells)')
