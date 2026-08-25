# Supplemental Figure S1: Training Set (DINOv2 Patch Tokens)

Supplement to Figure 1. The same 8-mutant *V. cholerae* training set (WT, Δ*vpsL*, Δ*rbmB*, Δ*hapR*,
Δ*potD1*, Δ*flaA*, *luxO*^D47E, *vpvC*^W240R), characterized from **DINOv2 patch tokens** instead of the
CLS token used in Figure 1: an embedding UMAP and an RF genotype-classification confusion matrix. Shows
whether the patch-token representation recovers the same genotype structure the CLS token does.

## Patch-token representation

DINOv2 emits, per (well, frame), one global **CLS** token plus a grid of **patch** tokens (9 here). Fig 1
uses the CLS token; this figure mean-pools the 9 patch tokens per (well, frame) to a single 768-d
"patch-mean" descriptor, then runs the **identical** pipeline as Fig 1 — stack the descriptor over the
growth-phase window frames 9–30 (22×768 = 16896 dims), UMAP (nn=25, md=0.25, rs=0), and GroupKFold RF.
Only CLS vs mean-patch differs, so S1 is a like-for-like comparison to Fig 1.

## Layout

```
figlib.py            thin shim: figS1 paths + strain order/labels/colors + re-export of ../figlib_shared.py
data/                source-data tables (CSVs); also on KiltHub
render/              draw each panel from data/  -> figures/   (what you run)
build/               regenerate data/ from the patch-mean embeddings (KiltHub; heavy, optional)
figures/             rendered PNG/SVG output
```
(Shared `requirements.txt` / `environment.yml` / `LICENSE` live at the `paper-figures/` root.)

## From raw data to panels (full chain)

1. **Raw brightfield timelapses** -> *(µPULLI image pipeline, separate repo `[pipeline repo]`)* ->
   **processed images + masks** (BioImage Archive [`S-BIAD3830`](https://doi.org/10.6019/S-BIAD3830)).
2. Processed images -> *(µPULLI feature extraction, DINOv2 route)* -> **`training_patchmean.npy`** (per-well
   mean-pooled patch descriptors [nWells, 31, 768]) + `training_embIndex.csv` (KiltHub `[DOI]`). The
   patch-mean tensor averages the 9 DINOv2 patch tokens per (well, frame); it is deposited as
   `training/embeddings/patchmean.npy` (see `../INPUTS.md`).
3. `training_patchmean.npy` -> **`build/S1B`** (stack frames 9-30 -> UMAP) / **`build/S1C`** (stack ->
   GroupKFold RF) -> **source-data CSVs in `data/`**. Labels + the growth filter are joined from
   `training_wide.parquet`.
4. `data/` CSVs -> **`render/S1*.py`** -> **panels** in `figures/`.

### Panels -> scripts -> tables

| Panel | `render/` script | `data/` table | `build/` script |
|---|---|---|---|
| **S1B** patch-token UMAP | `S1B_patchUmap.py` | `trainingPatchUmap_coords.csv` | `S1B_patchUmap.py` |
| **S1C** patch-token RF confusion matrix | `S1C_patchConfusion.py` | `patch_confusion_cv.csv` | `S1C_patchConfusion.py` |

Panel **S1A** is not generated here.

Growth filter (shared with Fig 1/2): keep a well if `mutant==vpsL` OR max biomass >= 0.15x the per-plate
WT median max. UMAP (nn=25/md=0.25/rs=0) and the RF CV are library-version sensitive, so the figures are
pinned to the tables in `data/`.

## Column dictionaries

- **`trainingPatchUmap_coords.csv`** - `metric` (`cosine` or `euclid`), `plateId`, `wellId`, `mutant`,
  `umap1`, `umap2`; one row per growth-filtered well per representation. The panel uses the `euclid` rows
  (StandardScaler + euclidean UMAP on the stacked mean-patch descriptor, frames 9-30).
- **`patch_confusion_cv.csv`** - 8x8 mean row-normalized confusion matrix (rows = true, cols = predicted)
  from GroupKFold RF on the stacked mean-patch descriptor; the panel's balanced accuracy = mean of the
  diagonal.
