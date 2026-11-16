#!/usr/bin/env python3
"""fig3 build - three-dimensional version of the reimaging UMAP landscape (nn=10, md=0.1).

Same manifold recipe as `3B_reimagingUmap_functionalAnnotations.py`, with `n_components=3` instead of
2, so the extra axis is the ONLY difference: identical row filters, identical feature selection,
identical scaler, identical UMAP hyperparameters and random_state. It is nonetheless a SEPARATE fit -
UMAP does not nest, so umap1/umap2 here are not the published 2-D coordinates with a third column
bolted on, and the two panels will not superimpose. Treat this as an alternative view of the same
data, never as a re-labelling of Figure 3B.

One table serves both colorings the atlas is drawn with, so the fit is done once:
  - functional annotation (the Fig 3B look), via `geneLocus` + HIGHLIGHT_SETS;
  - WT-normalized peak biofilm biomass (the retired 3A look), via `peakBiomassNorm`.

NOTE ON THE NAME: this is *not* Panel 3D (the dendrogram + heatmap). Everything here is spelled
`...3d...` in lower case and carries no panel letter, because it is an exploratory view rather than a
manuscript panel.

Reads:  reimaging/collapsedWide.parquet   (build input; see ../../INPUTS.md)
Writes: data/reimagingLandscape3d_nn10_md0.10_coords.csv
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig3/ for figlib
import re
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
import umap.umap_ as umap
from figlib import config, features

NN, MD, RANDOM_STATE = 10, 0.1, 42
stem = f'reimagingLandscape3d_nn{NN}_md{MD:.2f}'

# ── rows: exactly the filter chain of the 2-D build ───────────────────────────
wide = pd.read_parquet(config.input('reimaging/collapsedWide.parquet'))
wide = wide[wide['mutant'].notna()].reset_index(drop=True)
if 'geneLocus' in wide.columns:
    wide = wide[~wide['geneLocus'].isin(features.EXCLUDE_LOCI)].reset_index(drop=True)

featureCols = features.select_umap_feature_columns(wide)
X = wide[featureCols].copy().fillna(0)
zv = X.columns[X.nunique(dropna=False) <= 1].tolist()
if zv:
    X = X.drop(columns=zv)

gmask = features.growth_mask(X)
wide, X = wide[gmask].reset_index(drop=True), X[gmask].reset_index(drop=True)

low = features.low_replicate_genes(wide['mutant'])
if low:
    keep = ~wide['mutant'].isin(low)
    wide, X = wide[keep].reset_index(drop=True), X[keep].reset_index(drop=True)

print(f'Feature matrix: {X.shape}  (WT replicates: {int((wide["mutant"].astype(str) == "WT").sum())})')

# ── the fit: the one line that differs from the 2-D build ─────────────────────
XScaled = StandardScaler().fit_transform(X).astype(np.float32)
reducer = umap.UMAP(n_neighbors=NN, min_dist=MD, n_components=3, metric='euclidean',
                    random_state=RANDOM_STATE, low_memory=True)
embedding = reducer.fit_transform(XScaled)

metaCols = [c for c in ['plateId', 'wellId', 'mutant', 'geneLocus', 'function'] if c in wide.columns]
coords = wide[metaCols].copy()
coords['n_neighbors'], coords['min_dist'] = NN, MD
coords['umap1'], coords['umap2'], coords['umap3'] = embedding[:, 0], embedding[:, 1], embedding[:, 2]

# ── the biomass color channel (same definition as the retired 3A panel) ───────
# peak biomass = max over frames, normalized to the MEAN peak of the WT replicates that survived the
# same filters. Computed here so the render layer never has to touch the wide table.
biomassCols = [c for c in wide.columns if re.match(r'^biomass_t\d+$', c)]
coords['peakBiomass'] = wide[biomassCols].max(axis=1)
wtMask = coords['mutant'].astype(str) == 'WT'
wtMeanPeak = coords.loc[wtMask, 'peakBiomass'].mean()
coords['peakBiomassNorm'] = coords['peakBiomass'] / wtMeanPeak
print(f'WT reps: {int(wtMask.sum())} | WT mean peak biomass: {wtMeanPeak:.4f}')

coordsCsv = config.TABLES / f'{stem}_coords.csv'
coords.to_csv(coordsCsv, index=False)
print(f'Saved coordinates: {coordsCsv}  ({len(coords)} replicates)')
