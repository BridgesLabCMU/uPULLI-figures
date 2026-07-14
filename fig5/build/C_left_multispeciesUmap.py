"""fig5 build C-left: multispecies embedding-UMAP coordinates, 100% LB / 10X (source data for Panel 5C left).

Stacks DINOv2 CLS embeddings over frames 9-23 for the 100% LB wells at 10X, StandardScaler's them, and
fits UMAP (nn=10, md=0.1, random_state=42, euclidean) - the direct (non-PCA) embedding manifold used
for the panel. Species + LB condition come from the embedding index.

Reads:  config.MULTI_CLS (multispecies_10X_cls.npy), config.MULTI_EMBIDX
Writes: data/multispecies_100pctLB_10X_umap_coords.csv   (species, plateId, wellId, umap1, umap2)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig5/ for figlib
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
import umap.umap_ as umap
from figlib import config, MULTI_FRAMES

NN, MD, SEED = 10, 0.1, 42
COND = '100% LB'

cls = np.load(config.MULTI_CLS)
idx = pd.read_csv(config.MULTI_EMBIDX)
m = (idx['LB_condition'] == COND).to_numpy()
sub = idx[m].reset_index(drop=True)
rep = cls[m][:, MULTI_FRAMES, :].reshape(int(m.sum()), -1).astype(np.float32)
print(f'{COND} 10X: {rep.shape[0]} wells; stacked {rep.shape}')

Xs = StandardScaler().fit_transform(rep).astype(np.float32)
emb = umap.UMAP(n_neighbors=NN, min_dist=MD, n_components=2, metric='euclidean',
                random_state=SEED, low_memory=True).fit_transform(Xs)

cols = [c for c in ['species', 'plateId', 'wellId', 'LB_condition'] if c in sub.columns]
out = sub[cols].copy()
out['umap1'], out['umap2'] = emb[:, 0], emb[:, 1]
config.ensure(config.TABLES)
out.to_csv(config.TABLES / 'multispecies_100pctLB_10X_umap_coords.csv', index=False)
print(f'Saved: multispecies_100pctLB_10X_umap_coords.csv  ({len(out)} wells)')
