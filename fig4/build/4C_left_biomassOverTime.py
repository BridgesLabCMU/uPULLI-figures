"""fig4 build 4C-left: normalized biomass-over-time traces (source data for Panel 4C left).

Per group (WT from the reimaging atlas; Δ*bioD*/Δ*manA*/Δ*pdhE2* from the two clean-deletion plates),
computes mean +- SEM biofilm biomass at each frame, normalized to the median peak biomass of the
reimaging WT wells (median over WT wells of each well's max-over-time biomass). No-growth wells
(max biomass <= 0.005) are dropped first.

Reads:  config.CLEANDEL_WIDES (2 plates), config.REIM_WIDE (WT wells)
Writes: data/biomassOverTime_normWTpeak.csv   (group, frame, mean, sem, n)
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


cd = growth(pd.concat([pd.read_parquet(p) for p in config.CLEANDEL_WIDES], ignore_index=True))
reim = pd.read_parquet(config.REIM_WIDE, columns=['mutant'] + BIO)
wt = growth(reim[reim['mutant'].astype(str) == 'WT'])
normConst = float(wt[BIO].max(axis=1).median())
print(f'Reimaging WT wells {len(wt)} | normalization (median peak WT biomass) = {normConst:.5f}')

groups = {'WT': wt, 'BioD': cd[cd['mutant'] == 'BioD'], 'ManA': cd[cd['mutant'] == 'ManA'], 'PdhE2': cd[cd['mutant'] == 'PdhE2']}
rows = []
for grp, g in groups.items():
    vals = g[BIO].to_numpy(dtype=float) / normConst
    n = vals.shape[0]
    mean = np.nanmean(vals, axis=0)
    sem = np.nanstd(vals, axis=0, ddof=1) / np.sqrt(n)
    for i, t in enumerate(FRAMES):
        rows.append({'group': grp, 'frame': t, 'mean': mean[i], 'sem': sem[i], 'n': n})
    print(f'  {grp}: n={n}, peak mean={np.nanmax(mean):.3f}')

out = pd.DataFrame(rows)
config.ensure(config.TABLES)
out.to_csv(config.TABLES / 'biomassOverTime_normWTpeak.csv', index=False)
print(f'Saved: biomassOverTime_normWTpeak.csv ({len(out)} rows)')
