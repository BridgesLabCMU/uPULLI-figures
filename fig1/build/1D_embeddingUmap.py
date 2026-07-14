"""fig1 build 1D: training embedding-UMAP coordinates (source data for Panel 1D).

Stacks DINOv2 CLS embeddings over frames 9-27 (19x768 = 14592 dims) for the growth-filtered 8-mutant
training wells, then fits UMAP (n_neighbors=25, min_dist=0.25, random_state=0) under two representations:
cosine (L2-normalize + cosine metric) and euclid (StandardScaler + euclidean). Emits per-well coords for
both; the paper panel uses the euclid view.

Reads:  config.CLS (training_cls.npy), config.EMBIDX (training_embIndex.csv), config.WIDE (labels + growth filter)
Writes: data/trainingEmbeddingUmap_coords.csv   (metric, plateId, wellId, mutant, umap1, umap2)
"""
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig1/ for figlib
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler, normalize
import umap.umap_ as umap
from figlib import config, features, STRAIN_ORDER

FRAMES, NN, MD, RS = features.KEEP_FRAMES, 25, 0.25, 0

cls = np.load(config.CLS)
idx = pd.read_csv(config.EMBIDX)
wide = pd.read_parquet(config.WIDE)
bcols = [c for c in wide.columns if re.match(r'^biomass_t\d+$', c)]
maxBio = wide[bcols].max(axis=1)
wtMed = pd.DataFrame({'p': wide.plateId, 'b': maxBio, 'm': wide.mutant}).query("m=='WT'").groupby('p')['b'].median()
wide = wide.assign(_maxBio=maxBio.values)
wide['_keep'] = (wide.mutant == 'vpsL') | (wide._maxBio >= 0.15 * wide.plateId.map(wtMed))
lbl = wide[['plateId', 'wellId', 'mutant', '_keep']].drop_duplicates(['plateId', 'wellId'])

idx = idx.merge(lbl, on=['plateId', 'wellId'], how='left')
mask = idx['mutant'].isin(STRAIN_ORDER).to_numpy() & idx['_keep'].fillna(False).to_numpy()
sub = idx.loc[mask].reset_index(drop=True)
rep = cls[sub['row'].to_numpy()][:, FRAMES, :].reshape(len(sub), -1).astype(np.float32)
print(f'{mask.sum()} wells matched; stacked rep {rep.shape}')

meta = sub[['plateId', 'wellId', 'mutant']].reset_index(drop=True)
parts = []
for metric, X in [('cosine', normalize(rep)), ('euclid', StandardScaler().fit_transform(rep))]:
    emb = umap.UMAP(n_neighbors=NN, min_dist=MD, random_state=RS,
                    metric='cosine' if metric == 'cosine' else 'euclidean').fit_transform(X)
    d = meta.copy(); d['metric'] = metric; d['umap1'] = emb[:, 0]; d['umap2'] = emb[:, 1]
    parts.append(d)
    print(f'  {metric} UMAP done')

res = pd.concat(parts, ignore_index=True)
config.ensure(config.TABLES)
res.to_csv(config.TABLES / 'trainingEmbeddingUmap_coords.csv', index=False)
print(f'Saved: trainingEmbeddingUmap_coords.csv  {res.shape}')
