#!/usr/bin/env python3
"""interactive build step 0 — per-timepoint heatmap matrices for Interactive Plot 3.

Interactive Plot 3 animates the atlas heatmap across the timecourse, so it needs one z-scored
mutant x feature matrix per hour, not just the peak-biomass one bundled for Fig. S6A. This writes
them once as a small array; IP3 then needs no access to the 103 MB wide table.

Convention matches the published heatmap animations exactly:
  per-mutant MEDIAN of each feature at frame t  ->  z-score ACROSS MUTANTS WITHIN THAT FRAME
  -> clip to [-3, 3].
Each frame is standardized independently, so the color scale is comparable across mutants within a
timepoint (which is the comparison the figure is making) but not across timepoints; absolute drift
over the timecourse is deliberately removed. The peak-biomass frame is not recomputed here — IP3
appends the Fig. S6A matrix as the final slider position.

Reads:  reimaging/collapsedWide.parquet (config.input),
        data/fullAtlas_pcaLinkage_cluster_order.csv, data/fullAtlas_peakBiomass_featureMatrix.csv
Writes: data/fullAtlas_frameMatrices.npz  (frames 0-30, float16, ~250 kB)

Usage:  python build/0_buildFrameMatrices.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> interactive/ for figlib
import numpy as np
import pandas as pd
from figlib import config

FRAMES = list(range(31))
VMIN, VMAX = -3, 3

order = pd.read_csv(config.TABLES / 'fullAtlas_pcaLinkage_cluster_order.csv')
mutants = list(order['mutant'])                      # dendrogram leaf order, 158
features = list(pd.read_csv(config.TABLES / 'fullAtlas_peakBiomass_featureMatrix.csv')
                .rename(columns={'Unnamed: 0': 'feature'})['feature'])

import pyarrow.parquet as pq
have = set(pq.ParquetFile(config.input('reimaging/collapsedWide.parquet')).schema_arrow.names)
need = [f'{b}_t{t}' for b in features for t in FRAMES if f'{b}_t{t}' in have]
wide = pd.read_parquet(config.input('reimaging/collapsedWide.parquet'), columns=['mutant'] + need)
wide = wide[wide['mutant'].isin(mutants)]
print(f'{len(wide)} wells, {wide.mutant.nunique()}/{len(mutants)} mutants, '
      f'{len(features)} features x {len(FRAMES)} frames')

# NOT every feature exists at every frame: the 11 colony-segmentation features start at t5, since
# colonies are not segmentable before then. Those cells are left NaN and drawn as an explicit
# "no data" color, rather than filled with 0 -- a z of 0 would read as "average", which is wrong.
out = np.full((len(FRAMES), len(mutants), len(features)), np.nan, dtype=np.float32)
for ti, t in enumerate(FRAMES):
    cols = [(j, f'{b}_t{t}') for j, b in enumerate(features) if f'{b}_t{t}' in have]
    if not cols:
        continue
    js = [j for j, _ in cols]
    med = wide.groupby('mutant')[[c for _, c in cols]].median().reindex(mutants).to_numpy(dtype=float)
    sd = np.nanstd(med, axis=0)
    sd[sd == 0] = np.nan
    out[ti][:, js] = ((med - np.nanmean(med, axis=0)) / sd).clip(VMIN, VMAX)
nMissing = int(np.isnan(out).sum())
print(f'{nMissing} of {out.size} cells have no data ({100 * nMissing / out.size:.1f}%: '
      f'colony features before t5)')

np.savez_compressed(config.TABLES / 'fullAtlas_frameMatrices.npz',
                    matrices=out.astype(np.float16),
                    frames=np.array(FRAMES), mutants=np.array(mutants), features=np.array(features))
p = config.TABLES / 'fullAtlas_frameMatrices.npz'
print(f'Wrote {p}  ({p.stat().st_size / 1e3:.0f} kB)  shape {out.shape}')
