#!/usr/bin/env python3
"""figS6 build — replicate-consistency source tables for panels B and C.

Regenerates the four tables the S6B/C/D/E render scripts read, from the deposited reimaging wide
table. Self-contained: it reproduces the manifold's own preprocessing chain (the same one behind
Fig 3), then does every statistic on one pairwise-distance matrix.

Preprocessing, identical to the reimaging manifold build:
  drop rows with no mutant -> drop EXCLUDE_LOCI -> select_umap_feature_columns (frames 9-27) ->
  fillna(0) -> drop zero-variance columns -> growth filter -> min-5-replicates filter ->
  StandardScaler.  Result: 3669 wells x 285 standardized features, 158 mutants, 48 plates.

Statistics (all on Euclidean distances in that standardized space):
  * within- vs between-mutant pairwise distance distributions; AUC = normalized Mann-Whitney U.
    Pairwise distances are NOT independent (each well is in 3668 pairs), so every p-value here is
    obtained by permuting the mutant labels -- never from an analytic Mann-Whitney/KS test.
  * per-mutant mean within-replicate distance vs a per-mutant permutation null. This is 158
    separate tests, so it carries a Benjamini-Hochberg FDR column; note a permutation p cannot fall
    below 1/(nPerms+1), which is why the default is 10,000 permutations (at 1000 the floor sits
    ABOVE the Bonferroni threshold 0.05/158 and no mutant could pass regardless of effect size).
  * paired within-vs-between per mutant: Wilcoxon signed-rank with the MUTANT as the unit of
    analysis (n = 158), which is the one place a classical paired test is legitimate here.
  * two null models throughout: labels shuffled freely, and shuffled only within a plate (the
    latter holds plate composition fixed and is the control against a batch explanation).

Inputs (logical names, resolved via config.input -> see ../../INPUTS.md):
  reimaging/collapsedWide.parquet
Writes -> ../data/replicate_{distanceHistogram,perMutantNullHistogram,perMutant,summary}.csv

Usage:  python build/S6BC_replicateConsistency.py [--perms 10000]
        (~5 min at 10,000 permutations; --perms 1000 is ~40 s but disables the Bonferroni column)
"""
import sys
import argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS6/ for figlib
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import ks_2samp, rankdata, wilcoxon
from sklearn.preprocessing import StandardScaler
from figlib import config, features

NBINS, PM_NBINS = 80, 40          # panel B bins (millions of pairs) / panel C bins (158 values)
SEED = 42

ap = argparse.ArgumentParser()
ap.add_argument('--perms', type=int, default=10000)
args = ap.parse_args()

rng = np.random.default_rng(SEED)
out = config.TABLES
summary = {}


def put(key, value, note=''):
    summary[key] = value
    print(f'  {key:<32s} {value:>14.6g}  {note}')


# ── 1. Wide table -> the standardized feature matrix the manifold is built on ──
print('Loading the reimaging wide table ...')
wide = pd.read_parquet(config.input('reimaging/collapsedWide.parquet'))
wide = wide[wide['mutant'].notna()].reset_index(drop=True)
if 'geneLocus' in wide.columns:
    wide = wide[~wide['geneLocus'].isin(features.EXCLUDE_LOCI)].reset_index(drop=True)

featureCols = features.select_umap_feature_columns(wide)
X = wide[featureCols].copy().fillna(0)
X = X.drop(columns=X.columns[X.nunique(dropna=False) <= 1].tolist())

gmask = features.growth_mask(X)
wide, X = wide[gmask].reset_index(drop=True), X[gmask].reset_index(drop=True)
low = features.low_replicate_genes(wide['mutant'])
if low:
    keep = ~wide['mutant'].isin(low)
    wide, X = wide[keep].reset_index(drop=True), X[keep].reset_index(drop=True)

Xs = StandardScaler().fit_transform(X).astype(np.float64)
lab = wide['mutant'].astype(str).to_numpy()
plate = wide['plateId'].astype(str).to_numpy()
N = len(wide)
mutants = np.array(sorted(pd.unique(lab)))
groupIdx = [np.where(lab == m)[0] for m in mutants]
sizes = np.array([len(g) for g in groupIdx])
a = len(mutants)
print(f'  {N} wells x {Xs.shape[1]} features | {a} mutants | {len(np.unique(plate))} plates')

# ── 2. Pairwise distances ─────────────────────────────────────────────────────
print('Pairwise Euclidean distances ...')
D = squareform(pdist(Xs, 'euclidean')).astype(np.float32)
iu = np.triu_indices(N, 1)
dv = D[iu]
sameMut = (lab[:, None] == lab[None, :])[iu]
samePlate = (plate[:, None] == plate[None, :])[iu]
w, b = dv[sameMut], dv[~sameMut]
nW, nB = w.size, b.size

print('\nPanel B — distance distributions')
put('nWells', float(N)); put('nMutants', float(a))
put('nWithinPairs', float(nW)); put('nBetweenPairs', float(nB))
put('meanWithin', float(w.mean())); put('meanBetween', float(b.mean()))
put('cohensD', float((b.mean() - w.mean()) / np.sqrt((w.var(ddof=1) + b.var(ddof=1)) / 2)))
put('meanBetweenSamePlate', float(dv[samePlate & ~sameMut].mean()))
put('meanBetweenDiffPlate', float(dv[~samePlate & ~sameMut].mean()))
put('fracWithinPairsCrossPlate', float((~samePlate[sameMut]).mean()),
    '1.0 => replicate tightness cannot be a batch effect')

r = rankdata(dv.astype(np.float64))
R = squareform(r.astype(np.float32))


def aucFromRankSum(s):
    return float(1 - (s - nW * (nW + 1) / 2) / (nW * nB))


obsAuc = aucFromRankSum(r[sameMut].sum())
put('aucWithinCloser', obsAuc, '= normalized Mann-Whitney U')
put('cliffsDelta', float(2 * obsAuc - 1))
put('mannWhitneyU', float(obsAuc * nW * nB))
put('ksStatistic', float(ks_2samp(w, b).statistic), 'descriptive only -- analytic p NOT used')

edges = np.linspace(float(dv.min()), float(dv.max()), NBINS + 1)
hw, _ = np.histogram(w, bins=edges)
hb, _ = np.histogram(b, bins=edges)
width = np.diff(edges)
pd.DataFrame({'binLeft': edges[:-1], 'binRight': edges[1:], 'binCenter': (edges[:-1] + edges[1:]) / 2,
              'withinCount': hw, 'betweenCount': hb,
              'withinDensity': hw / hw.sum() / width,
              'betweenDensity': hb / hb.sum() / width}).to_csv(out / 'replicate_distanceHistogram.csv',
                                                              index=False)

# ── 3. Per-mutant observed values ─────────────────────────────────────────────
triCache = {n: np.triu_indices(n, 1) for n in np.unique(sizes)}


def blockMeans(blocks):
    return np.array([D[np.ix_(i, i)][triCache[len(i)]].mean() for i in blocks])


def blockRankSum(blocks):
    return float(sum(R[np.ix_(i, i)][triCache[len(i)]].sum() for i in blocks))


obsWithin = blockMeans(groupIdx)
allIdx = np.arange(N)
obsToOthers = np.array([D[np.ix_(i, np.setdiff1d(allIdx, i))].mean() for i in groupIdx])

# ── 4. Permutation nulls (two models) ─────────────────────────────────────────
# Permuting well indices into blocks of the observed group sizes is exactly equivalent to shuffling
# the mutant labels, and much cheaper.
splitAt = np.cumsum(sizes)[:-1]
print(f'\n{args.perms} permutations x 2 null models ...')
nullPerMut = np.empty((args.perms, a))
nullAuc = np.empty(args.perms)
for i in range(args.perms):
    blocks = np.split(rng.permutation(N), splitAt)
    nullPerMut[i] = blockMeans(blocks)
    nullAuc[i] = aucFromRankSum(blockRankSum(blocks))
    if (i + 1) % 2000 == 0:
        print(f'  shuffled-label null {i + 1}/{args.perms}')

plateWells = [np.where(plate == p)[0] for p in np.unique(plate)]
nullPerMutStrat = np.empty((args.perms, a))
for i in range(args.perms):
    permLab = lab.copy()
    for pw in plateWells:
        permLab[pw] = rng.permutation(lab[pw])
    nullPerMutStrat[i] = blockMeans([np.where(permLab == m)[0] for m in mutants])
    if (i + 1) % 2000 == 0:
        print(f'  within-plate null {i + 1}/{args.perms}')

pFloor = 1 / (args.perms + 1)
put('nPerms', float(args.perms)); put('permP_floor', pFloor)
put('auc_z', float((obsAuc - nullAuc.mean()) / nullAuc.std(ddof=1)))
put('auc_p', float((1 + int((nullAuc >= obsAuc).sum())) * pFloor), 'permutation test, panel B')
put('nullAuc_mean', float(nullAuc.mean())); put('nullAuc_sd', float(nullAuc.std(ddof=1)))

print('\nPanel C — per-mutant distributions')
obsMean = float(obsWithin.mean())
put('obsMeanWithinPerMutant', obsMean)
for name, null in (('global', nullPerMut), ('withinPlate', nullPerMutStrat)):
    m = null.mean(axis=1)
    put(f'null_{name}_mean', float(m.mean())); put(f'null_{name}_sd', float(m.std(ddof=1)))
    put(f'null_{name}_z', float((obsMean - m.mean()) / m.std(ddof=1)))
    put(f'null_{name}_p', float((1 + int((m <= obsMean).sum())) * pFloor))

pool = np.sort(nullPerMut.ravel())


def aucVsPool(v):
    return float((pool.size - np.searchsorted(pool, v, side='right')).sum() / (v.size * pool.size))


obsPmAuc = aucVsPool(obsWithin)
nullPmAuc = np.array([aucVsPool(nullPerMut[i]) for i in range(args.perms)])
put('perMutantAuc', obsPmAuc, 'P(observed per-mutant value < null per-mutant value)')
put('perMutantAuc_z', float((obsPmAuc - nullPmAuc.mean()) / nullPmAuc.std(ddof=1)))
put('perMutantAuc_p', float((1 + int((nullPmAuc >= obsPmAuc).sum())) * pFloor))
put('nMutantsBelowNullMean', float(int((obsWithin < nullPerMut.mean()).sum())))

pmEdges = np.linspace(min(obsWithin.min(), nullPerMut.min(), nullPerMutStrat.min()) - 0.2,
                      max(obsWithin.max(), nullPerMut.max(), nullPerMutStrat.max()) + 0.2, PM_NBINS + 1)
pmWidth = np.diff(pmEdges)
hObs, _ = np.histogram(obsWithin, bins=pmEdges)
hNull, _ = np.histogram(nullPerMut.ravel(), bins=pmEdges)
hNullS, _ = np.histogram(nullPerMutStrat.ravel(), bins=pmEdges)
pd.DataFrame({'binLeft': pmEdges[:-1], 'binRight': pmEdges[1:],
              'binCenter': (pmEdges[:-1] + pmEdges[1:]) / 2,
              'observedCount': hObs, 'nullGlobalCount': hNull, 'nullWithinPlateCount': hNullS,
              'observedDensity': hObs / hObs.sum() / pmWidth,
              'nullGlobalDensity': hNull / hNull.sum() / pmWidth,
              'nullWithinPlateDensity': hNullS / hNullS.sum() / pmWidth,
              }).to_csv(out / 'replicate_perMutantNullHistogram.csv', index=False)


# ── 5. Per-mutant table with FDR (panels D, E) ────────────────────────────────
def benjaminiHochberg(p):
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    res = np.empty(n)
    res[order] = np.clip(q, 0, 1)
    return res


perMut = pd.DataFrame({
    'mutant': mutants,
    'geneLocus': wide.groupby('mutant')['geneLocus'].first().reindex(mutants).fillna('').to_numpy(),
    'nReplicates': sizes, 'meanWithin': obsWithin, 'meanToOtherMutants': obsToOthers,
    'nullMeanWithin': nullPerMut.mean(axis=0), 'nullSdWithin': nullPerMut.std(axis=0, ddof=1),
})
perMut['p'] = [(1 + int((nullPerMut[:, i] <= obsWithin[i]).sum())) * pFloor for i in range(a)]
perMut['qBH'] = benjaminiHochberg(perMut['p'])
perMut['pBonferroni'] = np.clip(perMut['p'] * a, 0, 1)
perMut['withinMinusOthers'] = perMut.meanWithin - perMut.meanToOtherMutants
perMut.sort_values('meanWithin').to_csv(out / 'replicate_perMutant.csv', index=False)

print('\nPanels D/E — per-mutant tests (158 tests -> FDR)')
put('nMutantsBelowNull', float(int((perMut.meanWithin < perMut.nullMeanWithin).sum())))
put('nMutantsPermP_lt_0.05', float(int((perMut.p < 0.05).sum())), 'UNCORRECTED')
put('nMutantsFDR_lt_0.05', float(int((perMut.qBH < 0.05).sum())), 'Benjamini-Hochberg')
put('nMutantsBonferroni_lt_0.05', float(int((perMut.pBonferroni < 0.05).sum())))
put('bonferroniThreshold', float(0.05 / a), 'unreachable unless permP_floor is below it')

diffs = perMut.withinMinusOthers.to_numpy()
stat, pw = wilcoxon(perMut.meanWithin, perMut.meanToOtherMutants, alternative='less')
ranks = rankdata(np.abs(diffs))
boot = np.array([np.median(rng.choice(diffs, diffs.size, replace=True)) for _ in range(10000)])
put('nMutantsTighterThanOthers', float(int((diffs < 0).sum())))
put('wilcoxon_stat', float(stat))
put('wilcoxon_p', float(pw), 'normal approximation; see README before quoting at face value')
put('rankBiserial', float((ranks[diffs > 0].sum() - ranks[diffs < 0].sum()) / ranks.sum()))
put('medianPairedDiff', float(np.median(diffs)))
put('medianPairedDiff_ci_lo', float(np.percentile(boot, 2.5)))
put('medianPairedDiff_ci_hi', float(np.percentile(boot, 97.5)))
put('smallestPairedGap', float(np.abs(diffs).min()))

pd.DataFrame({'statistic': list(summary), 'value': list(summary.values())}).to_csv(
    out / 'replicate_summary.csv', index=False)
print(f'\nWrote 4 tables -> {out}/replicate_*.csv')
