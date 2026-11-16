#!/usr/bin/env python3
"""interactive build step 0 — per-timepoint heatmap matrices for Interactive Plot 3.

Interactive Plot 3 animates the atlas heatmap across the timecourse, so it needs one z-scored
mutant x feature matrix per hour, not just the peak-biomass one bundled for Fig. S6A. This writes
them once as a small array; IP3 then needs no access to the 103 MB wide table.

Written for BOTH of the plot's views — the 158-mutant full atlas (Fig. S6A) and the 50-leaf
functional-annotation subset (Fig. 3D), which the page toggles between.

Convention matches the published heatmap animations exactly:
  per-mutant MEDIAN of each feature at frame t  ->  z-score ACROSS MUTANTS WITHIN THAT FRAME
  -> clip to [-3, 3].
Each frame is standardized independently, so the color scale is comparable across mutants within a
timepoint (which is the comparison the figure is making) but not across timepoints; absolute drift
over the timecourse is deliberately removed. The peak-biomass frame is not recomputed here — IP3
appends the Fig. S6A / Fig. 3D matrix as the final slider position.

The z-score is always taken over the WHOLE atlas and only then subset to a view's leaves, which is
what `generateDendogramHeatmaps-withPCA.py` does for the functional panel (it z-scores the full
`workerDf` and selects columns afterwards). Standardizing within the 50 leaves instead would give
the same mutant a different color in the two views.

Reads:  reimaging/collapsedWide.parquet (config.input),
        data/<set>_pcaLinkage_cluster_order.csv, data/<set>_peakBiomass_featureMatrix.csv
Writes: data/<set>_frameMatrices.npz  (frames 0-30, float16)

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
SETS = ['fullAtlas', 'functionalSubset']       # fullAtlas first: it defines the z-score population

leafOrder = {s: list(pd.read_csv(config.TABLES / f'{s}_pcaLinkage_cluster_order.csv')['mutant'])
             for s in SETS}
featureNames = {s: list(pd.read_csv(config.TABLES / f'{s}_peakBiomass_featureMatrix.csv')
                        .rename(columns={'Unnamed: 0': 'feature'})['feature']) for s in SETS}
if featureNames['functionalSubset'] != featureNames['fullAtlas']:
    raise SystemExit('the two views must share one feature set — they share one heatmap axis')
features = featureNames['fullAtlas']

# the z-score population: every mutant the atlas clusters, in full-atlas leaf order
population = leafOrder['fullAtlas']

import pyarrow.parquet as pq
have = set(pq.ParquetFile(config.input('reimaging/collapsedWide.parquet')).schema_arrow.names)
need = [f'{b}_t{t}' for b in features for t in FRAMES if f'{b}_t{t}' in have]
wide = pd.read_parquet(config.input('reimaging/collapsedWide.parquet'), columns=['mutant'] + need)
wide = wide[wide['mutant'].isin(population)]
print(f'{len(wide)} wells, {wide.mutant.nunique()}/{len(population)} mutants, '
      f'{len(features)} features x {len(FRAMES)} frames')

# NOT every feature exists at every frame: the 11 colony-segmentation features start at t5, since
# colonies are not segmentable before then. Those cells are left NaN and drawn as an explicit
# "no data" color, rather than filled with 0 -- a z of 0 would read as "average", which is wrong.
out = np.full((len(FRAMES), len(population), len(features)), np.nan, dtype=np.float32)
for ti, t in enumerate(FRAMES):
    cols = [(j, f'{b}_t{t}') for j, b in enumerate(features) if f'{b}_t{t}' in have]
    if not cols:
        continue
    js = [j for j, _ in cols]
    med = wide.groupby('mutant')[[c for _, c in cols]].median().reindex(population).to_numpy(dtype=float)
    sd = np.nanstd(med, axis=0)
    sd[sd == 0] = np.nan
    out[ti][:, js] = ((med - np.nanmean(med, axis=0)) / sd).clip(VMIN, VMAX)
nMissing = int(np.isnan(out).sum())
print(f'{nMissing} of {out.size} cells have no data ({100 * nMissing / out.size:.1f}%: '
      f'colony features before t5)')

rowOf = {m: i for i, m in enumerate(population)}
for s in SETS:
    rows = [rowOf[m] for m in leafOrder[s]]          # reorder/subset to this view's leaves
    mats = out[:, rows, :]
    p = config.TABLES / f'{s}_frameMatrices.npz'
    np.savez_compressed(p, matrices=mats.astype(np.float16), frames=np.array(FRAMES),
                        mutants=np.array(leafOrder[s]), features=np.array(features))
    print(f'Wrote {p}  ({p.stat().st_size / 1e3:.0f} kB)  shape {mats.shape}')
