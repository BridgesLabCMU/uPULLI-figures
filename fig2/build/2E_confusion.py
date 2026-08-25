"""fig2 build 2E: GroupKFold Random-Forest confusion matrix, all-features (source data for Panel 2E).

Port of the 'all' combination from the training confusion-matrix analysis: GroupKFold by plateId
(5 folds x 5 repeats), RandomForest (200 trees, balanced), predicting the 8-mutant genotype from
biomass(log1p) + whole haralick/entropy + the 12 colony bases (frames 9-27; colony 9-28). Emits the
mean row-normalized confusion matrix (the panel's balanced accuracy = mean of its diagonal).

Reads:  config.input('training/wide.parquet') (training_wide.parquet)
Writes: data/all_confusion_cv.csv
"""
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig2/ for figlib
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import confusion_matrix, balanced_accuracy_score
from figlib import config, features, STRAIN_ORDER

randomState, nRepeats, nSplits = 0, 5, 5
minFrame, maxFrame = 9, 30   # colony window: exclude pre-colony frames (colonies segmentable from t9); start-frame consistent with the 9-27 UMAP filter + multispecies. biomass/whole use full 0-30.
allowed = features.ALLOWED_COLONY_BASES


def selBio(cols):
    return [c for c in cols if re.match(r'^biomass_t\d+$', c)]


def selWhole(cols):
    return [c for c in cols if c.startswith('whole_') and c.rsplit('_t', 1)[0].startswith(('whole_haralick', 'whole_entropy'))]


def selColony(cols):
    out = []
    for c in cols:
        m = re.search(r'_t(\d+)$', c)
        if m and minFrame <= int(m.group(1)) <= maxFrame and c.rsplit('_t', 1)[0] in allowed:
            out.append(c)
    return out


df = pd.read_parquet(config.input('training/wide.parquet'))
bcols = selBio(df.columns)
maxBio = df[bcols].max(axis=1)
wtMed = (pd.DataFrame({'plateId': df['plateId'], 'maxBio': maxBio, 'mutant': df['mutant']})
         .query("mutant=='WT'").groupby('plateId')['maxBio'].median())
keep = (df['mutant'] == 'vpsL') | (maxBio >= 0.15 * df['plateId'].map(wtMed))
df = df[keep & df['mutant'].isin(STRAIN_ORDER)].reset_index(drop=True)
print(f'{len(df)} wells after growth filter')

Xbio = np.log1p(df[selBio(df.columns)].fillna(df[selBio(df.columns)].median())).values
whole, col = selWhole(df.columns), selColony(df.columns)
Xwhole = df[whole].fillna(df[whole].median()).values
Xcol = df[col].fillna(df[col].median()).values
X = np.hstack([Xbio, Xwhole, Xcol])
y = df['mutant'].values
groupsOrig = df['plateId'].values
labels = [s for s in STRAIN_ORDER if s in y]

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
print(f'all: balanced accuracy (CV) = {np.mean(scores):.4f}  (== mean diag = {np.mean(np.diag(meanCm)):.4f})')
config.ensure(config.TABLES)
pd.DataFrame(meanCm, index=labels, columns=labels).to_csv(config.TABLES / 'all_confusion_cv.csv')
print('Saved: all_confusion_cv.csv')
