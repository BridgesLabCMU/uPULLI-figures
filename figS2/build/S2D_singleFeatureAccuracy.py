"""figS2 build S2D: single-feature timecourse RF accuracy (source data for Panel S2D).

Port of v2/training/classification/singleFeatureAccuracy.py. For each feature base, train a GroupKFold
RandomForest on that feature's timecourse ALONE and record balanced accuracy. Feature bases: biomass +
whole haralick/entropy + the 12 ALLOWED_COLONY_BASES. Colony features use frames 9-30 (segmentable from
t9); biomass/whole use 0-30. Same CV protocol as the confusion panels (GroupKFold by plateId, 5x5,
RandomForest 200 trees balanced). Growth filter shared with Fig 2.

Reads:  config.input('training/wide.parquet') (training_wide.parquet)
Writes: data/singleFeatureAccuracy.csv   (featureBase, prettyName, family, balancedAccuracy)
"""
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS2/ for figlib
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import balanced_accuracy_score
from figlib import config, features, STRAIN_ORDER

RANDOM_STATE, nRepeats, nSplits, nTrees = 0, 5, 5, 200
minFrameAll, maxFrameAll = 0, 30
minFrameColony, maxFrameColony = 9, 30
colonyBases = list(features.ALLOWED_COLONY_BASES)


def isColony(base):
    return base.startswith('colony_') or base == 'nColonies'


def family(base):
    return 'biomass' if base == 'biomass' else ('whole' if base.startswith('whole') else ('colony' if isColony(base) else 'other'))


def featureBases(cols):
    fams = set()
    for c in cols:
        m = re.search(r'_t(\d+)$', c)
        if not m:
            continue
        base = c.rsplit('_t', 1)[0]
        if base == 'biomass' or base.startswith('whole_haralick') or base.startswith('whole_entropy'):
            fams.add(base)
        elif base in colonyBases:
            fams.add(base)
    return sorted(fams)


def evaluateSingleFeature(base, wideDf, meta):
    frames = range(minFrameColony, maxFrameColony + 1) if isColony(base) else range(minFrameAll, maxFrameAll + 1)
    cols = [f'{base}_t{t}' for t in frames if f'{base}_t{t}' in wideDf.columns]
    if not cols:
        return None
    X = wideDf[cols].replace([np.inf, -np.inf], np.nan)
    X = X.loc[:, ~X.isna().all(axis=0)]
    if X.shape[1] == 0:
        return None
    X = X.fillna(X.median())
    X = X.loc[:, X.var(axis=0) > 0]
    if X.shape[1] == 0:
        return None
    if base == 'biomass':
        X = np.log1p(X)
    X = np.nan_to_num(X.values)
    y = meta['mutant'].values
    groupsOrig = meta['plateId'].values
    uniq = np.unique(groupsOrig)
    rng = np.random.default_rng(RANDOM_STATE)
    scores = []
    for _ in range(nRepeats):
        shuffled = uniq.copy(); rng.shuffle(shuffled)
        pm = dict(zip(uniq, shuffled))
        groups = np.array([pm[g] for g in groupsOrig])
        for tr, te in GroupKFold(n_splits=nSplits).split(X, y, groups):
            sc = StandardScaler()
            rf = RandomForestClassifier(n_estimators=nTrees, min_samples_leaf=2, class_weight='balanced',
                                        random_state=RANDOM_STATE, n_jobs=1)
            rf.fit(sc.fit_transform(X[tr]), y[tr])
            scores.append(balanced_accuracy_score(y[te], rf.predict(sc.transform(X[te]))))
    return base, float(np.mean(scores))


wideDf = pd.read_parquet(config.input('training/wide.parquet'))
bcols = [c for c in wideDf.columns if re.match(r'^biomass_t\d+$', c)]
maxBio = wideDf[bcols].max(axis=1)
wtMed = (pd.DataFrame({'p': wideDf['plateId'], 'b': maxBio, 'm': wideDf['mutant']})
         .query("m=='WT'").groupby('p')['b'].median())
keep = (wideDf['mutant'] == 'vpsL') | (maxBio >= 0.15 * wideDf['plateId'].map(wtMed))
wideDf = wideDf[keep & wideDf['mutant'].isin(STRAIN_ORDER)].reset_index(drop=True)
meta = wideDf[['plateId', 'wellId', 'mutant']].copy()
bases = featureBases(wideDf.columns)
print(f'{len(wideDf)} wells; evaluating {len(bases)} feature bases...')

results = Parallel(n_jobs=-1)(delayed(evaluateSingleFeature)(f, wideDf, meta) for f in bases)
results = [r for r in results if r is not None]
results.sort(key=lambda x: x[1], reverse=True)

resDf = pd.DataFrame({'featureBase': [r[0] for r in results],
                      'prettyName': [features.pretty_name(r[0]) for r in results],
                      'family': [family(r[0]) for r in results],
                      'balancedAccuracy': [r[1] for r in results]})
config.ensure(config.TABLES)
resDf.to_csv(config.TABLES / 'singleFeatureAccuracy.csv', index=False)
print(resDf.to_string(index=False))
print('Saved: singleFeatureAccuracy.csv')
