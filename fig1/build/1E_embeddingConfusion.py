"""fig1 build 1E: GroupKFold Random-Forest confusion matrix from DINOv2 embeddings (Panel 1E).

Same CV protocol as the numerical confusion matrix (GroupKFold by plateId, 5 folds x 5 repeats,
RandomForest 200 trees balanced, per-fold StandardScaler), but the feature matrix is the stacked CLS
embedding (frames 9-27 -> 14592 dims). Emits the mean row-normalized confusion matrix.

Reads:  config.CLS, config.EMBIDX, config.WIDE (labels + growth filter)
Writes: data/embeddings_confusion_cv.csv
"""
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig1/ for figlib
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import confusion_matrix, balanced_accuracy_score
from figlib import config, features, STRAIN_ORDER

randomState, nRepeats, nSplits = 0, 5, 5
FRAMES = features.KEEP_FRAMES

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
sub = idx.loc[mask]
X = cls[sub['row'].to_numpy()][:, FRAMES, :].reshape(int(mask.sum()), -1).astype(np.float32)
y = sub['mutant'].to_numpy()
groupsOrig = sub['plateId'].to_numpy()
labels = [s for s in STRAIN_ORDER if s in y]
print(f'{mask.sum()} wells; X={X.shape}')

uniquePlates = np.unique(groupsOrig)
rng = np.random.default_rng(randomState)
scores, allCm = [], np.zeros((len(labels), len(labels)))
for _ in range(nRepeats):
    shuffled = uniquePlates.copy(); rng.shuffle(shuffled)
    gmap = dict(zip(uniquePlates, shuffled))
    groups = np.array([gmap[g] for g in groupsOrig])
    for tr, te in GroupKFold(n_splits=nSplits).split(X, y, groups):
        sc = StandardScaler()
        clf = RandomForestClassifier(n_estimators=200, min_samples_leaf=2, class_weight='balanced',
                                     random_state=randomState, n_jobs=4)
        clf.fit(sc.fit_transform(X[tr]), y[tr])
        pred = clf.predict(sc.transform(X[te]))
        scores.append(balanced_accuracy_score(y[te], pred))
        cm = confusion_matrix(y[te], pred, labels=labels).astype(float)
        allCm += cm / cm.sum(axis=1, keepdims=True)
meanCm = allCm / (nRepeats * nSplits)
print(f'embeddings: balanced accuracy (CV) = {np.mean(scores):.4f}')
config.ensure(config.TABLES)
pd.DataFrame(meanCm, index=labels, columns=labels).to_csv(config.TABLES / 'embeddings_confusion_cv.csv')
print('Saved: embeddings_confusion_cv.csv')
