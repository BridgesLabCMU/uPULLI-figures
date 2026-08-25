#!/usr/bin/env python3
"""fig2 build 2B: per-replicate biomass-trace matrix (source data for Panel 2B).

Emits the WT-peak-normalized biofilm-biomass trajectory of every replicate of the 8 training mutants:
one row per well (mutant, plateId, wellId, biomass_t0..t30 / WT-peak-mean). The render script orders
rows by strain then (plateId, wellId), forward/back-fills trailing NaNs, and draws the heatmap.

Reads:  config.input('training/wide.parquet') (training_wide.parquet)
Writes: data/biomassTraces_normWTpeak.csv
"""
import sys
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig2/ for figlib
import pandas as pd
from figlib import config, STRAIN_ORDER

df = pd.read_parquet(config.input('training/wide.parquet'))
df = df[df['mutant'].isin(STRAIN_ORDER)].reset_index(drop=True)
bcols = sorted([c for c in df.columns if re.match(r'^biomass_t\d+$', c)], key=lambda c: int(c.split('_t')[1]))

wtPeakMean = float(df[df['mutant'] == 'WT'][bcols].astype(float).max(axis=1).mean())
out = df[['mutant', 'plateId', 'wellId']].copy()
out[bcols] = df[bcols].astype(float) / wtPeakMean

outPath = config.ensure(config.TABLES) / 'biomassTraces_normWTpeak.csv'
out.to_csv(config.TABLES / 'biomassTraces_normWTpeak.csv', index=False)
print(f'Saved: {outPath}  ({len(out)} replicates; WT peak mean = {wtPeakMean:.4g})')
