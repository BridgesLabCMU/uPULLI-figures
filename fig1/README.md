# Figure 1: Training Set (DINOv2 Embeddings)

The 8-mutant *V. cholerae* training set (WT, Δ*vpsL*, Δ*rbmB*, Δ*hapR*, Δ*potD1*, Δ*flaA*, *luxO*^D47E,
*vpvC*^W240R) characterized from **DINOv2 CLS embeddings** (the µPULLI deep-representation route): an
embedding UMAP and an RF genotype-classification confusion matrix. Compare to Figure 2, which does the
same from hand-engineered numerical features.

## Layout

```
figlib.py            thin shim: fig1 paths + strain order/labels/colors + re-export of ../figlib_shared.py
data/                source-data tables (CSVs); also on KiltHub
render/              draw each panel from data/  -> figures/   (what you run)
build/               regenerate data/ from the CLS embeddings (KiltHub; heavy, optional)
figures/             rendered PNG/SVG output
```
(Shared `requirements.txt` / `environment.yml` / `LICENSE` live at the repository root.)

## From raw data to panels (full chain)

1. **Raw brightfield timelapses** -> *(µPULLI image pipeline, separate repo `[pipeline repo]`)* ->
   **processed images + masks** (BioImage Archive `[S-BIAD#####]`).
2. Processed images -> *(µPULLI feature extraction, DINOv2 route)* -> **`training_cls.npy`** (per-well
   CLS embeddings [nWells, 31, 768]) + `training_embIndex.csv` (KiltHub `[DOI]`).
3. `training_cls.npy` -> **`build/1D`** (stack frames 9-27 -> UMAP) / **`build/1E`** (stack -> GroupKFold
   RF) -> **source-data CSVs in `data/`**. Labels + the growth filter are joined from `training_wide.parquet`.
4. `data/` CSVs -> **`render/1*.py`** -> **panels** in `figures/`.

There is no wide-table reshaping step here (the embeddings are the feature representation); `build/`
consumes the CLS array directly. Steps 1-2 are the µPULLI pipeline (separate repo). Fill bracketed IDs
at submission.

### Panels -> scripts -> tables

| Panel | `render/` script | `data/` table | `build/` script |
|---|---|---|---|
| **1D** embedding UMAP | `1D_embeddingUmap.py` | `trainingEmbeddingUmap_coords.csv` | `1D_embeddingUmap.py` |
| **1E** embedding RF confusion matrix | `1E_embeddingConfusion.py` | `embeddings_confusion_cv.csv` | `1E_embeddingConfusion.py` |

Growth filter (shared with Fig 2): keep a well if `mutant==vpsL` OR max biomass >= 0.15x the per-plate
WT median max. UMAP (nn=25/md=0.25/rs=0) and the RF CV are library-version sensitive, so the figures are
pinned to the tables in `data/`.

## Column dictionaries

- **`trainingEmbeddingUmap_coords.csv`** - `metric` (`cosine` or `euclid`), `plateId`, `wellId`,
  `mutant`, `umap1`, `umap2`; one row per growth-filtered well per representation. The panel uses the
  `euclid` rows (StandardScaler + euclidean UMAP on the stacked CLS embedding, frames 9-27).
- **`embeddings_confusion_cv.csv`** - 8x8 mean row-normalized confusion matrix (rows = true, cols =
  predicted) from GroupKFold RF on the stacked CLS embedding; the panel's balanced accuracy = mean of
  the diagonal.
