"""fig5 build B-left: K. pneumoniae embedding-UMAP coordinates (source data for Panel 5B left).

Stacks DINOv2 CLS embeddings over frames 9-24 for the growth-filtered 5-strain kleb set (waaL dropped),
then fits UMAP (nn=10, md=0.1, random_state=0) under cosine (L2-norm) and euclid (StandardScaler)
representations; the panel uses euclid. Kleb has one mutant per plate, so labels come from
plate-timestamp -> NV id -> gene; growth filter = max biomass >= 0.005 (mrkA exempt, expected non-former).

Reads:  config.KLEB_CLS, config.KLEB_EMBIDX, config.KLEB_FRAME (biomass for growth filter)
Writes: data/kleb_embeddingUmap_coords.csv   (metric, plateId, wellId, mutant, umap1, umap2)
"""
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig5/ for figlib
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler, normalize
import umap.umap_ as umap
from figlib import config, KLEB_FRAMES, KLEB_PLATE_NV, KLEB_MUTANT_MAP, KLEB_NO_GROWTH_FLOOR, KLEB_ORDER

NN, MD, RS = 10, 0.1, 0
cls = np.load(config.KLEB_CLS)
idx = pd.read_csv(config.KLEB_EMBIDX)


def _ts(p):
    m = re.match(r'(\d{6}_\d{6})', str(p)); return m.group(1) if m else None


idx['mutant'] = idx['plateId'].map(lambda p: KLEB_PLATE_NV.get(_ts(p))).map(KLEB_MUTANT_MAP)
fr = pd.read_csv(config.KLEB_FRAME, usecols=['plateID', 'wellID', 'biomass'])
fr['plateId'] = fr['plateID'].astype(str).str.replace(' ', '_', regex=False)
fr['wellId'] = fr['wellID'].astype(str).str.replace(r'_\d+$', '', regex=True)
maxBio = fr.groupby(['plateId', 'wellId'])['biomass'].max().rename('_maxBio')
idx = idx.merge(maxBio, on=['plateId', 'wellId'], how='left')
idx['_keep'] = (idx['mutant'] == 'mrkA') | (idx['_maxBio'].fillna(0.0) >= KLEB_NO_GROWTH_FLOOR)

mask = idx['mutant'].isin(KLEB_ORDER).to_numpy() & idx['_keep'].fillna(False).to_numpy()
sub = idx.loc[mask].reset_index(drop=True)
rep = cls[sub['row'].to_numpy()][:, KLEB_FRAMES, :].reshape(len(sub), -1).astype(np.float32)
print(f'{mask.sum()} kleb wells; stacked rep {rep.shape}')

meta = sub[['plateId', 'wellId', 'mutant']].reset_index(drop=True)
parts = []
for metric, X in [('cosine', normalize(rep)), ('euclid', StandardScaler().fit_transform(rep))]:
    emb = umap.UMAP(n_neighbors=NN, min_dist=MD, random_state=RS,
                    metric='cosine' if metric == 'cosine' else 'euclidean').fit_transform(X)
    d = meta.copy(); d['metric'] = metric; d['umap1'] = emb[:, 0]; d['umap2'] = emb[:, 1]
    parts.append(d); print(f'  {metric} done')

res = pd.concat(parts, ignore_index=True)
config.ensure(config.TABLES)
res.to_csv(config.TABLES / 'kleb_embeddingUmap_coords.csv', index=False)
print(f'Saved: kleb_embeddingUmap_coords.csv  {res.shape}')
