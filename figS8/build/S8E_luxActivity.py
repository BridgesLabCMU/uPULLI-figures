#!/usr/bin/env python3
"""figS8 build E: per-replicate peak PvpsL-lux activity (source data for Panel S8E).

Ports the lab's peak calculation (Continuous_Peak_Plotting.R, 251212_lux_37C_Shaking) so the figure
matches the analysis of record:

  RLU(well, t) = Lum / OD, with inf / 0 / negative set to NaN
  peak(well)   = max over the 41 hourly timepoints
  normalized   = peak(well) / MEAN of the WT wells' peaks     <- WT reference is the mean, not median

Well -> condition comes from that script's layout. Row B (pdhR) is deliberately excluded: it is
commented out of the lab analysis and is not part of this figure.

Reads:  config.input('lux/Lum.csv'), config.input('lux/OD.csv')
Writes: data/luxPeak_normWT.csv   (condition, well, peakRLU, peakNormWT)

Usage:  python figS8/build/S8E_luxActivity.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS8/ for figlib
import numpy as np
import pandas as pd
from figlib import config, LUX_WELLS

WT = 'WT'


def readPlate(path):
    d = pd.read_csv(path)
    d = d.rename(columns={d.columns[0]: 'Time'})
    wells = [c for c in d.columns if c != 'Time']
    return d[wells].apply(pd.to_numeric, errors='coerce')


lum, od = readPlate(config.input('lux/Lum.csv')), readPlate(config.input('lux/OD.csv'))
if len(lum) != len(od):
    raise SystemExit(f'Lum.csv ({len(lum)} rows) and OD.csv ({len(od)} rows) disagree on timepoints')
shared = [w for w in lum.columns if w in od.columns]
rlu = (lum[shared] / od[shared]).replace([np.inf, -np.inf], np.nan)
rlu = rlu.mask(rlu <= 0)                      # zero / negative RLU is not a measurement
print(f'{len(rlu)} timepoints x {len(shared)} wells')

rows = []
for cond, wells in LUX_WELLS.items():
    for w in wells:
        if w not in rlu.columns:
            print(f'[WARN] {cond}: well {w} absent from the plate tables - skipped')
            continue
        peak = rlu[w].max(skipna=True)
        if pd.isna(peak):
            print(f'[WARN] {cond}: well {w} has no usable RLU - skipped')
            continue
        rows.append({'condition': cond, 'well': w, 'peakRLU': float(peak)})

out = pd.DataFrame(rows)
wtRef = out.loc[out['condition'] == WT, 'peakRLU'].mean()
if not np.isfinite(wtRef) or wtRef <= 0:
    raise SystemExit('WT reference peak is not usable')
out['peakNormWT'] = out['peakRLU'] / wtRef
print(f'WT reference (mean peak RLU over {int((out.condition == WT).sum())} wells) = {wtRef:.1f}')
for c, g in out.groupby('condition', sort=False):
    print(f'  {c:6s} n={len(g):2d}  median={g.peakNormWT.median():.3f}  '
          f'mean={g.peakNormWT.mean():.3f}  range=[{g.peakNormWT.min():.3f}, {g.peakNormWT.max():.3f}]')

config.ensure(config.TABLES)
out.to_csv(config.TABLES / 'luxPeak_normWT.csv', index=False)
print(f'Saved: luxPeak_normWT.csv ({len(out)} wells)')
