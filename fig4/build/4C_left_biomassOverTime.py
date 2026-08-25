"""fig4 build 4C-left: normalized biomass-over-time traces (source data for Panel 4C left).

Per group (WT from the reimaging atlas; Δ*bioD*/Δ*manA*/Δ*pdhE2* from the two clean-deletion plates),
computes mean +- SD biofilm biomass at each frame, normalized to the median peak biomass of the
reimaging WT wells (median over WT wells of each well's max-over-time biomass). No-growth wells
(max biomass <= 0.005) are dropped first.

Reads:  [config.input('cluster/cleanDel_260521_collapsedWide.parquet'),
                       config.input('cluster/cleanDel_260522_collapsedWide.parquet')] (2 plates), config.input('reimaging/collapsedWide.parquet') (WT wells)
Writes: data/biomassOverTime_normWTpeak.csv   (group, frame, mean, sd, n)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig4/ for figlib
import numpy as np
import pandas as pd
from figlib import config, features

FRAMES = list(range(31))
BIO = [f'biomass_t{t}' for t in FRAMES]


def growth(df):
    return df[df[BIO].max(axis=1) > features.NO_GROWTH_FLOOR]


cd = growth(pd.concat([pd.read_parquet(p) for p in [config.input('cluster/cleanDel_260521_collapsedWide.parquet'),
                       config.input('cluster/cleanDel_260522_collapsedWide.parquet')]], ignore_index=True))
reim = pd.read_parquet(config.input('reimaging/collapsedWide.parquet'), columns=['mutant'] + BIO)
wt = growth(reim[reim['mutant'].astype(str) == 'WT'])
normConst = float(wt[BIO].max(axis=1).median())
print(f'Reimaging WT wells {len(wt)} | normalization (median peak WT biomass) = {normConst:.5f}')

groups = {'WT': wt, 'BioD': cd[cd['mutant'] == 'BioD'], 'ManA': cd[cd['mutant'] == 'ManA'], 'PdhE2': cd[cd['mutant'] == 'PdhE2']}
rows = []
for grp, g in groups.items():
    vals = g[BIO].to_numpy(dtype=float) / normConst
    n = vals.shape[0]
    mean = np.nanmean(vals, axis=0)
    # SD, not SEM: the panel's claim is that replicate trajectories differ reproducibly between
    # genotypes, which is a statement about replicate spread. SEM (= SD/sqrt(n), n = 24-32 here) is
    # ~5-6x narrower and reads as far more precision than the biology has. It is also not comparable
    # across these series: reimaging WT's wells come from 24 different plates, while each clean
    # deletion's 32 wells come from only 2, so their effective n is much smaller than 32.
    sd = np.nanstd(vals, axis=0, ddof=1)
    for i, t in enumerate(FRAMES):
        rows.append({'group': grp, 'frame': t, 'mean': mean[i], 'sd': sd[i], 'n': n})
    print(f'  {grp}: n={n}, peak mean={np.nanmax(mean):.3f}')

out = pd.DataFrame(rows)
config.ensure(config.TABLES)
out.to_csv(config.TABLES / 'biomassOverTime_normWTpeak.csv', index=False)
print(f'Saved: biomassOverTime_normWTpeak.csv ({len(out)} rows)')
