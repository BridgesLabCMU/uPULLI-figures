"""figS2 build S2E: mutant separability across time (source data for Panel S2E).

Port of v2/training/classification/timepointSeparability.py.
  * single-timepoint: for each frame t, train a GroupKFold RF on the quantitative features AT THAT FRAME
    ONLY (biomass + whole haralick/entropy + colony), balanced accuracy over 5 repeats x 5 folds;
  * full-timecourse: one RF on all those features across all frames (the dashed baseline).
biomass is log1p'd. GroupKFold by plateId with the standard plate-shuffle, RandomForest 200 trees.
Growth filter shared with Fig 2.

Reads:  config.input('training/wide.parquet') (training_wide.parquet)
Writes: data/timepointSalience_groupkfold.csv   (frame, rfMean, rfStd)
        data/fullTimecourse_accuracy.csv         (rfMean, rfStd)
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
from figlib import config, STRAIN_ORDER

RANDOM_STATE, nRepeats, nSplits, nTrees = 0, 5, 5, 200


def isFeature(base):
    return (base == 'biomass' or base.startswith('whole_haralick')
            or base.startswith('whole_entropy') or base.startswith('colony_') or base == 'nColonies')


def cvAccuracy(X, y, groupsOrig):
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
    return float(np.mean(scores)), float(np.std(scores))


def cleanMatrix(sub, logBiomass=True):
    X = sub.replace([np.inf, -np.inf], np.nan)
    X = X.loc[:, ~X.isna().all(axis=0)]
    X = X.fillna(X.median())
    X = X.loc[:, X.var(axis=0) > 0]
    bcols = [c for c in X.columns if c.startswith('biomass_')]
    if logBiomass and bcols:
        X[bcols] = np.log1p(X[bcols])
    return np.nan_to_num(X.values)


wideDf = pd.read_parquet(config.input('training/wide.parquet'))
bcols = [c for c in wideDf.columns if re.match(r'^biomass_t\d+$', c)]
maxBio = wideDf[bcols].max(axis=1)
wtMed = (pd.DataFrame({'p': wideDf['plateId'], 'b': maxBio, 'm': wideDf['mutant']})
         .query("m=='WT'").groupby('p')['b'].median())
keep = (wideDf['mutant'] == 'vpsL') | (maxBio >= 0.15 * wideDf['plateId'].map(wtMed))
wideDf = wideDf[keep & wideDf['mutant'].isin(STRAIN_ORDER)].reset_index(drop=True)
y = wideDf['mutant'].values
groups = wideDf['plateId'].values
print(f'{len(wideDf)} wells after growth + strain filter')

frameCols = {}
for c in wideDf.columns:
    m = re.search(r'_t(\d+)$', c)
    if not m:
        continue
    if isFeature(c[:m.start()]) and pd.api.types.is_numeric_dtype(wideDf[c]):
        frameCols.setdefault(int(m.group(1)), []).append(c)
frames = sorted(frameCols)


def perFrame(t):
    X = cleanMatrix(wideDf[frameCols[t]])
    mean, std = cvAccuracy(X, y, groups)
    return t, mean, std


print(f'Single-timepoint RF over {len(frames)} frames ...')
res = Parallel(n_jobs=-1)(delayed(perFrame)(t) for t in frames)
res.sort()
validFrames = np.array([r[0] for r in res])
rfMean = np.array([r[1] for r in res])
rfStd = np.array([r[2] for r in res])
config.ensure(config.TABLES)
pd.DataFrame({'frame': validFrames, 'rfMean': rfMean, 'rfStd': rfStd}).to_csv(
    config.TABLES / 'timepointSalience_groupkfold.csv', index=False)

print('Full-timecourse RF (all frames) ...')
allCols = sorted({c for cols in frameCols.values() for c in cols})
fullMean, fullStd = cvAccuracy(cleanMatrix(wideDf[allCols]), y, groups)
pd.DataFrame({'rfMean': [fullMean], 'rfStd': [fullStd]}).to_csv(config.TABLES / 'fullTimecourse_accuracy.csv', index=False)
print(f'full-timecourse balanced accuracy = {fullMean:.4f} +/- {fullStd:.4f}')
print('Saved: timepointSalience_groupkfold.csv, fullTimecourse_accuracy.csv')
