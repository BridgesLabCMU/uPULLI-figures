# build/: regenerate the source-data tables

These scripts recompute the tables in `../data/` from the **original full feature set** (not needed to
render the figures — use `../render/` for that). They are the analysis that produced the deposited
tables: reshaping the feature set into the wide table (step 0), then the UMAP fit, per-mutant centroids,
peak-biomass feature matrices, and the PCA→Ward clustering. Upstream of step 0 (raw images → processed
images → the feature set) is the µPULLI pipeline in a separate repository — see `../README.md`.

## Inputs (obtain from KiltHub, then point figlib at them)

- `master_frame_features.csv`: the original full per-(well, frame) feature set (µPULLI output);
  input to step 0.
- `reimagingIndex.csv`: gene locus/name/function annotations.
- `reimaging_umapEmbeddings.parquet`: precomputed UMAP embeddings (bundled in `../data/`, used by the
  centroid/per-gene builders so they need only the small file).
- `reimaging_collapsedWide.parquet`: the wide table; produced by step 0, then used by 3A/3B/3C.

Set the paths for a standalone checkout:

```bash
export FIG3_MASTER_FRAME=/path/to/master_frame_features.csv
export FIG3_REIMAGING_INDEX=/path/to/reimagingIndex.csv
export FIG3_WIDE_TABLE=/path/to/reimaging_collapsedWide.parquet   # where step 0 writes / 3* read
```

(In the original monorepo these default to the in-repo/`/mnt` paths, so the scripts run without env vars.)

## Run

```bash
python build/0_buildCollapsedWide.py                       # master_frame_features.csv -> reimaging_collapsedWide.parquet (WIDE)
python build/3B_reimagingUmap_functionalAnnotations.py     # -> data/reimagingLandscape_..._coords.csv (refits UMAP)
python build/3A_reimagingUmap_coloredByBiomass.py          # -> data/coloredByBiomass..._coords.csv
python build/3Ctop_reimagingUmap_centroidsByFunction.py    # -> data/centroidsByFunction..._centroids.csv
python build/3Cbottom_reimagingUmap_centroidsByLocus.py    # -> data/centroidsByLocus..._centroids.csv
python build/3C_generateDendogramHeatmaps-withPCA.py       # -> data/functional_*.csv + linkage.npy
python build/reimagingUmap_perGenePdf.py                   # -> data/reimagingUmap_..._perGene_coords.csv
```

CSVs are written to `../data/`; any figures these emit go to `../figures/`.

**Note:** UMAP output can drift across `umap-learn`/`numba`/`numpy` versions even with
`random_state=42`, so a rebuild may not be bit-identical to the deposited coordinates. Use the pinned
`environment.yml` for the closest match; the published figures are pinned to the coordinates in
`../data/` (see the full chain in `../README.md`).
