# build/: regenerate the source-data tables

These scripts recompute the tables in `../data/` from the **original full feature set** (not needed to
render the figures — use `../render/` for that). They are the analysis that produced the deposited
tables: reshaping the feature set into the wide table (step 0), then the UMAP fit, per-mutant centroids,
peak-biomass feature matrices, and the PCA→Ward clustering. Upstream of step 0 (raw images → processed
images → the feature set) is the µPULLI pipeline in a separate repository — see `../README.md`.

## Inputs

Fetch them once, from the repository root:

```bash
python fetch_data.py fig3
```

That populates `build/inputs/` with the logical names these scripts resolve through
`config.input(...)` — see [`../../INPUTS.md`](../../INPUTS.md) for sizes and contents:

- `reimaging/master_frame_features.csv` — per-(well, frame) feature set; input to step 0.
- `reimaging/geneIndex.csv` — gene locus / name / function annotations.
- `reimaging/collapsedWide.parquet` — the wide table; written by step 0, read by 3A/3B/3D.
- `reimaging/umapEmbeddings.parquet` — UMAP coordinates, used by the centroid/per-gene builders.

Already have the deposit unpacked elsewhere? `export UPULLI_DATA_ROOT=/path/to/deposit` instead.

## Run

```bash
python build/0_buildCollapsedWide.py                       # master_frame_features.csv -> reimaging_collapsedWide.parquet (WIDE)
python build/3B_reimagingUmap_functionalAnnotations.py     # -> data/reimagingLandscape_..._coords.csv (refits UMAP)
python build/3A_reimagingUmap_coloredByBiomass.py          # -> data/coloredByBiomass..._coords.csv
python build/3Btop_reimagingUmap_centroidsByFunction.py    # -> data/centroidsByFunction..._centroids.csv
python build/3Bbottom_reimagingUmap_centroidsByLocus.py    # -> data/centroidsByLocus..._centroids.csv
python build/3D_generateDendogramHeatmaps-withPCA.py       # -> data/functional_*.csv + linkage.npy
python build/reimagingUmap_perGenePdf.py                   # -> data/reimagingUmap_..._perGene_coords.csv
```

CSVs are written to `../data/`; any figures these emit go to `../figures/`.

**Note:** UMAP output can drift across `umap-learn`/`numba`/`numpy` versions even with
`random_state=42`, so a rebuild may not be bit-identical to the deposited coordinates. Use the pinned
`environment.yml` for the closest match; the published figures are pinned to the coordinates in
`../data/` (see the full chain in `../README.md`).
