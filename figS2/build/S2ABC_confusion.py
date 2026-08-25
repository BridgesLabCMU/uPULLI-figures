"""figS2 build S2A/B/C: per-feature-family RF confusion matrices (source data for Panels S2A-C).

Port of v2/training/classification/confusionMatrices.py, restricted to the three single-family feature
sets (the all-feature matrix is Fig 2E). Same CV protocol as Fig 2: GroupKFold by plateId (5 folds x 5
repeats), RandomForest (200 trees, balanced), per-fold StandardScaler, balanced accuracy, predicting the
8-mutant genotype. Feature windows match the numerical pipeline exactly:
  * biomass  (Panel S2A) : biomass_t0..t30, log1p
  * colony   (Panel S2B) : the 12 ALLOWED_COLONY_BASES over frames 9-30 (colonies segmentable from t9)
  * whole    (Panel S2C) : whole_haralick_* + whole_entropy_* over frames 0-30
Growth filter (shared with Fig 2): keep a well if mutant==vpsL OR max biomass >= 0.15x per-plate WT median.

Reads:  config.input('training/wide.parquet') (training_wide.parquet)
Writes: data/{bio,colony,whole}_confusion_cv.csv   (8x8 mean row-normalized confusion matrices)
"""
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS2/ for figlib
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import confusion_matrix, balanced_accuracy_score
from figlib import config, features, STRAIN_ORDER

randomState, nRepeats, nSplits = 0, 5, 5
minFrame, maxFrame = 9, 30   # colony window (matches confusionMatrices.py); biomass/whole use full 0-30
allowedColonyBases = features.ALLOWED_COLONY_BASES


def selectBiomass(cols):
    return [c for c in cols if re.match(r'^biomass_t\d+$', c)]


def selectWhole(cols):
    out = []
    for c in cols:
        if not c.startswith('whole_'):
            continue
        base = c.rsplit('_t', 1)[0]
        if base.startswith('whole_haralick') or base.startswith('whole_entropy'):
            out.append(c)
    return out


def selectColony(cols):
    out = []
    for c in cols:
        m = re.search(r'_t(\d+)$', c)
        if not m:
            continue
        if not (minFrame <= int(m.group(1)) <= maxFrame):
            continue
        if c.rsplit('_t', 1)[0] in allowedColonyBases:
            out.append(c)
    return out


def buildMatrix(X, y, groupsOriginal, labels, prefix):
    uniquePlates = np.unique(groupsOriginal)
    rng = np.random.default_rng(randomState)
    allScores, allCm = [], np.zeros((len(labels), len(labels)))
    for _ in range(nRepeats):
        shuffled = uniquePlates.copy(); rng.shuffle(shuffled)
        plateMap = dict(zip(uniquePlates, shuffled))
        groups = np.array([plateMap[g] for g in groupsOriginal])
        for tr, te in GroupKFold(n_splits=nSplits).split(X, y, groups):
            sc = StandardScaler()
            clf = RandomForestClassifier(n_estimators=200, min_samples_leaf=2, class_weight='balanced',
                                         random_state=randomState, n_jobs=4)
            clf.fit(sc.fit_transform(X[tr]), y[tr])
            pred = clf.predict(sc.transform(X[te]))
            allScores.append(balanced_accuracy_score(y[te], pred))
            cm = confusion_matrix(y[te], pred, labels=labels).astype(float)
            allCm += cm / cm.sum(axis=1, keepdims=True)
    meanAcc = float(np.mean(allScores))
    meanCm = allCm / (nRepeats * nSplits)
    print(f'{prefix}: balanced accuracy (CV) = {meanAcc:.4f}')
    config.ensure(config.TABLES)
    pd.DataFrame(meanCm, index=labels, columns=labels).to_csv(config.TABLES / f'{prefix}_confusion_cv.csv')
    print(f'Saved: {prefix}_confusion_cv.csv')


wideDf = pd.read_parquet(config.input('training/wide.parquet'))
bcols = selectBiomass(wideDf.columns)
maxBio = wideDf[bcols].max(axis=1)
wtMedByPlate = (pd.DataFrame({'plateId': wideDf['plateId'], 'maxBio': maxBio, 'mutant': wideDf['mutant']})
                .query("mutant == 'WT'").groupby('plateId')['maxBio'].median())
wtRef = wideDf['plateId'].map(wtMedByPlate)
keep = (wideDf['mutant'] == 'vpsL') | (maxBio >= 0.15 * wtRef)
wideDf = wideDf[keep & wideDf['mutant'].isin(STRAIN_ORDER)].reset_index(drop=True)
print(f'After growth + strain filter: {len(wideDf)} wells')

bio = selectBiomass(wideDf.columns); whole = selectWhole(wideDf.columns); col = selectColony(wideDf.columns)
print(f'features: biomass={len(bio)}, whole={len(whole)}, colony={len(col)}')
Xbio = np.log1p(wideDf[bio].fillna(wideDf[bio].median())).values
Xwhole = wideDf[whole].fillna(wideDf[whole].median()).values
Xcol = wideDf[col].fillna(wideDf[col].median()).values
y = wideDf['mutant'].values
groups = wideDf['plateId'].values
labels = [s for s in STRAIN_ORDER if s in y]

buildMatrix(Xbio, y, groups, labels, 'bio')       # Panel S2A - Biofilm Biomass
buildMatrix(Xcol, y, groups, labels, 'colony')    # Panel S2B - Colony-level Features
buildMatrix(Xwhole, y, groups, labels, 'whole')   # Panel S2C - Whole-image Features
print('Done — S2A/B/C confusion matrices.')
