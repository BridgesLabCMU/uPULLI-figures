# Build inputs

The `render/` scripts that redraw each panel read the source-data tables bundled in `fig*/data/`, so a
fresh clone reproduces every figure with **no downloads and no configuration**. The `build/` scripts
that regenerate those tables from the underlying feature data need the inputs below.

```bash
python fetch_data.py            # download everything (~996 MB) into fig*/build/inputs/
python fetch_data.py figS4      # or just what one figure needs
python fetch_data.py --list     # show this table from the manifest
python fetch_data.py --verify   # re-check checksums
```

Already have the deposit unpacked? Point at it instead and skip the download:

```bash
export UPULLI_DATA_ROOT=/path/to/deposit
```

Deposit: **uPULLI figure build inputs** — DOI `[pending]`.

`inputs.json` is the machine-readable manifest (logical name, size, SHA-256, consuming figures) and is
the single source of truth: this table is generated from it, and every build script resolves its inputs
through `config.input("<logical name>")`. No script contains an absolute path.

| Logical name | Size | Used by | Contents |
|---|---|---|---|
| `cluster/cleanDel_260521_collapsedWide.parquet` | 5 MB | fig4, fig5 | Clean-deletion plate 260521, pivoted wide. |
| `cluster/cleanDel_260522_collapsedWide.parquet` | 5 MB | fig4, fig5 | Clean-deletion plate 260522, pivoted wide. |
| `cluster/cmpd_260522_collapsedWide.parquet` | 6 MB | fig5 | Compound-treatment plate 260522, pivoted wide. |
| `kleb/embeddings/cls.npy` | 44 MB | fig5 | DINOv2 CLS embeddings for the K. pneumoniae set. |
| `kleb/embeddings/index.csv` | 0 MB | fig5 | Row index for the K. pneumoniae embedding array. |
| `kleb/master_frame_features.csv` | 20 MB | fig5 | K. pneumoniae mutant set: per-(plate, well, frame) features (growth filter). |
| `multispecies/embeddings/cls_10X.npy` | 31 MB | fig5 | DINOv2 CLS embeddings for the 8-species panel, 10X wells. |
| `multispecies/embeddings/index_10X.csv` | 0 MB | fig5 | Row index for the multispecies 10X embedding array. |
| `reimaging/collapsedWide.parquet` | 103 MB | fig3, fig4, fig5, figS4, figS6 | Reimaging atlas pivoted wide to <feature>_t<frame>, joined to the gene index. |
| `reimaging/embeddings/dendroClusterOrder.csv` | 0 MB | figS7 | Leaf order + cluster assignment for the embedding-PC dendrogram. |
| `reimaging/embeddings/dendroPcaCentroids.csv` | 0 MB | figS7 | Per-mutant centroids in embedding-PC space, input to the S7 dendrograms. |
| `reimaging/geneIndex.csv` | 4 MB | fig3 | Reimaging well -> geneLocus / geneName / functional annotation. |
| `reimaging/master_frame_features.csv` | 223 MB | fig3 | Reimaging atlas: per-(plate, well, frame) features. |
| `reimaging/reimagingResults_10x.csv` | 22 MB | figS4 | Reimaging analysis of record: set membership + phenotype per mutant. |
| `reimaging/thumbnails112/` | 12 MB | interactive | Peak-biofilm-biomass thumbnail pack: one 112px grayscale JPEG per reimaging well, keyed `<plate-well>_t<peakFrame>.jpg`, plus manifest.csv. Inlined as base64 by the interactive plots. |
| `reimaging/umapEmbeddings.parquet` | 0 MB | fig3 | Reimaging UMAP coordinates over the (n_neighbors, min_dist) grid. |
| `rnaseq/` | — | fig4 | Directory of per-strain differential-expression tables (<strain>_allGenes.csv) for the Fig 4E biofilm-gene heatmap. Sequencing data are deposited separately under their own accession. |
| `training/embeddings/cls.npy` | 155 MB | fig1 | DINOv2 CLS embeddings for the training set, [nWells, 31, 768] float32. |
| `training/embeddings/index.csv` | 0 MB | fig1, figS1 | Row index for the training embedding arrays (row -> plate, well). |
| `training/embeddings/patchmean.npy` | 155 MB | figS1 | Mean-pooled DINOv2 patch tokens for the training set, [nWells, 31, 768] float32. |
| `training/layout.csv` | 0 MB | fig2, figS4 | Training set well -> mutant map (plateId, wellId, mutant). |
| `training/master_frame_features.csv` | 89 MB | fig2, figS4 | Training set (8 V. cholerae strains): per-(plate, well, frame) features. |
| `training/wide.parquet` | 42 MB | fig2, figS1, figS2, figS3 | Training set pivoted wide to <feature>_t<frame>; labelled via training/layout.csv. |
| `transposons/master_frame_features.csv` | 80 MB | figS4 | Transposon screen: per-(plate, well, frame) features (trajectories for S4A / S5). |
| `transposons/results_10x.csv` | 11 MB | figS4 | Genome-wide screen analysis of record: per-well Peak/Final/Early/Phenotype. These calls selected the reimaging library and are read as given, never recomputed. |

24 inputs, 996 MB total.

Two inputs are *produced* by a step-0 build rather than downloaded — `training/wide.parquet`
(`fig2/build/0_buildTrainingWide.py`) and `reimaging/collapsedWide.parquet`
(`fig3/build/0_buildCollapsedWide.py`). They are deposited as well so the later steps can run without
repeating step 0.

## Also in the deposit (not build inputs)

The KiltHub item carries three things `fetch_data.py` does **not** download, because no build step
reads them:

* **Data S1** — `DataS1_transposonScreen_perWell.csv` (2,939 wells × 203 columns) and
  `DataS1_transposonScreen_allFrames.csv.gz` (91,097 well-frames): every quantified feature for every
  screened transposon mutant, with its phenotype call. Built by `figS4/build/DataS1_transposonScreen.py`.
* **`source-data-tables.zip`** — the bundled `fig*/data` tables plus each figure's README (their column
  dictionary), so the deposit is readable without cloning.
* **`opentrons-protocols.zip`** — the liquid-handling protocols, with a README mapping each to the
  figures it produced.
