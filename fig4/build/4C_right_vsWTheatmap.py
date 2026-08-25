"""fig4 build 4C-right: WT-normalized peak-feature matrix (source data for Panel 4C right).

For each clean-deletion condition, samples every feature (biomass + whole entropy/haralick + the 12
colony bases) at that condition's peak-biomass frame, then reports (condition - plate WT) / reimaging-
atlas sigma. WT is the reference and is dropped from the matrix. Rows = conditions, columns = features
(canonical order). NOTE: fixes the `whole_entropy_` filter typo in the original so the Global Entropy
row IS included.

Reads:  [config.input('cluster/cleanDel_260521_collapsedWide.parquet'),
                       config.input('cluster/cleanDel_260522_collapsedWide.parquet')] (conditions + WT reference), config.input('reimaging/collapsedWide.parquet') (atlas sigma)
Writes: data/cleanDel_vsWT_horizontal_matrix.csv   (conditions x features)
"""
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig4/ for figlib
import numpy as np
import pandas as pd
from figlib import config, features

PEAK_MIN, PEAK_MAX = 0, 30
allowed = features.ALLOWED_COLONY_BASES
canonHaralick = ['energy', 'contrast', 'correlation', 'variance', 'inverse_difference_moment', 'sum_average',
                 'sum_variance', 'sum_entropy', 'entropy', 'difference_variance', 'difference_entropy', 'imc1', 'imc2']


def splitFF(c):
    m = re.match(r'(.+)_t(\d+)$', c); return (m.group(1), int(m.group(2))) if m else (c, None)


def inRange(c):
    m = re.search(r'_t(\d+)$', c); return m and PEAK_MIN <= int(m.group(1)) <= PEAK_MAX


def featGroup(b):
    return 'biomass' if b == 'biomass' else 'haralick' if 'haralick' in b else 'entropy' if 'entropy' in b else 'colony'


def sortKey(b):
    if b == 'biomass':
        return (0, 0)
    if featGroup(b) == 'colony':
        return (1, b)
    if 'entropy' in b:
        return (2, 0)
    if 'haralick' in b:
        k = re.sub(r'_(mean|std|var)$', '', b.replace('whole_haralick_', ''))
        return (3, int(k) if k.isdigit() else (canonHaralick.index(k) if k in canonHaralick else 99))
    return (4, b)


df = pd.concat([pd.read_parquet(p) for p in [config.input('cluster/cleanDel_260521_collapsedWide.parquet'),
                       config.input('cluster/cleanDel_260522_collapsedWide.parquet')]], ignore_index=True)
cols = []
for c in df.columns:
    if not inRange(c):
        continue
    b = splitFF(c)[0]
    # FIX: whole_entropy (not whole_entropy_) so the Global Entropy row is kept
    if (b == 'biomass' or b.startswith('whole_haralick_') or b.startswith('whole_entropy') or b in allowed) \
            and pd.api.types.is_numeric_dtype(df[c]):
        cols.append(c)
bases = sorted(set(splitFF(c)[0] for c in cols), key=sortKey)

biomassCols = {splitFF(c)[1]: c for c in df.columns if c.startswith('biomass_') and inRange(c)}
peakFrame = {}
for m, g in df.groupby('mutant'):
    traj = pd.Series({fr: g[col].median() for fr, col in biomassCols.items()}).sort_index()
    peakFrame[m] = int(traj.idxmax()) if traj.notna().any() else PEAK_MAX

mutants = list(df['mutant'].dropna().unique())
mutants = (['WT'] if 'WT' in mutants else []) + sorted(m for m in mutants if m != 'WT')
mat = pd.DataFrame(index=bases, columns=mutants, dtype=float)
for m in mutants:
    g = df[df['mutant'] == m]; pf = peakFrame[m]
    for b in bases:
        col = f'{b}_t{pf}'
        mat.loc[b, m] = g[col].median() if col in g.columns else np.nan


def atlasStats(bases):
    reim = pd.read_parquet(config.input('reimaging/collapsedWide.parquet'))
    bc = {int(re.search(r'_t(\d+)$', c).group(1)): c for c in reim.columns if re.match(r'^biomass_t\d+$', c)}
    frames = sorted(bc)
    pf = np.array(frames)[np.nanargmax(reim[[bc[f] for f in frames]].to_numpy(), axis=1)]
    out = {}
    for b in bases:
        cm = {int(re.search(r'_t(\d+)$', c).group(1)): c for c in reim.columns if re.match(rf'^{re.escape(b)}_t\d+$', c)}
        if not cm:
            out[b] = (np.nan, np.nan); continue
        vals = np.array([reim[cm[f]].to_numpy()[i] if f in cm else np.nan for i, f in enumerate(pf)])
        out[b] = (np.nanmean(vals), np.nanstd(vals))
    return pd.DataFrame(out, index=['mean', 'std']).T


st = atlasStats(bases)
wtCol = 'WT' if 'WT' in mat.columns else next((c for c in mat.columns if str(c).startswith('WT')), None)
if wtCol is None:
    raise SystemExit('vsWT needs a WT/WT_* control column')
matz = mat.sub(mat[wtCol], axis=0).div(st['std'].replace(0, np.nan), axis=0).fillna(0).drop(columns=[wtCol])
condOrder = [m for m in mutants if m != wtCol]
H = matz[condOrder].T                                # rows = conditions, cols = features
config.ensure(config.TABLES)
H.to_csv(config.TABLES / 'cleanDel_vsWT_horizontal_matrix.csv')
print(f'Saved: cleanDel_vsWT_horizontal_matrix.csv  ({H.shape[0]} conditions x {H.shape[1]} features; WT={wtCol} dropped)')
