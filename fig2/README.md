# Figure 2: Training Set (Numerical Phenotyping)

The 8-mutant *V. cholerae* training set (WT, Δ*vpsL*, Δ*rbmB*, Δ*hapR*, Δ*potD1*, Δ*flaA*, *luxO*^D47E,
*vpvC*^W240R) characterized from hand-engineered numerical features: biomass traces, a peak-biomass
feature heatmap, a multimodal UMAP, an RF genotype-classification confusion matrix, and single-feature
timecourse heatmaps.

## Layout

```
figlib.py            thin shim — fig2 paths + strain order/labels/colors + re-export of ../figlib_shared.py
data/                source-data tables (CSVs); also on KiltHub
render/              draw each panel from data/  → figures/   (what you run)
build/               regenerate data/ from the original feature set (KiltHub; optional)
figures/             rendered PNG/SVG output
```
(Shared `requirements.txt` / `environment.yml` / `LICENSE` live at the `paper-figures/` root.)

## From raw data to panels (full chain)

1. **Raw brightfield timelapses** → *(µPULLI image pipeline — separate repo `[pipeline repo]`)* →
   **processed images + masks** (BioImage Archive [`S-BIAD3830`](https://doi.org/10.6019/S-BIAD3830)).
2. Processed images → *(µPULLI feature extraction — `multiWellAnalysis` + DINOv2)* →
   **`master_frame_features.csv`**: the original full per-(well, frame) training feature set (KiltHub `[DOI]`).
3. `master_frame_features.csv` → **`build/0_buildTrainingWide.py`** → **`training_wide.parquet`**
   (canonical renames, pivot to `<feature>_t<frame>`, recover `(plateId,wellId)→mutant` labels).
4. Wide table → **`build/2*.py`** (aggregation, UMAP fit, GroupKFold RF) → **source-data CSVs in `data/`**
   (KiltHub `[DOI]` + bundled here).
5. `data/` CSVs → **`render/2*.py`** → **panels** in `figures/`.

Steps 1–2 are the µPULLI pipeline (separate repo); steps 3–5 are here. `render/` needs only `data/`;
`build/` needs the KiltHub feature set (`build/README.md`). Fill bracketed IDs at submission.

### Panels → scripts → tables

| Panel | `render/` script | `data/` table(s) | `build/` script |
|---|---|---|---|
| **2B** biomass traces | `2B_biomassTraces.py` | `biomassTraces_normWTpeak.csv` | `2B_biomassTraces.py` |
| **2C** peak-biomass feature heatmap, features clustered | `2C_trainingHeatmap.py` | `trainingHeatmap_peakBiomass_featureMatrix.csv` (+ `_peakFrames.csv`) | `2C_trainingHeatmap.py` |
| 2C variant, features in family order (no dendrogram) | `2C_trainingHeatmap_familyOrder.py` | same table | — |
| **2D** multimodal UMAP | `2D_trainingUmap.py` | `trainingUmap_all_three_coords.csv` | `2D_trainingUmap.py` |
| **2E** RF confusion matrix | `2E_confusion.py` | `all_confusion_cv.csv` | `2E_confusion.py` |
| **2F** single-feature heatmaps (×4) | `2F_featmaps.py` | `featmap_<feat>.csv` ×4 + `featmaps_meta.csv` | `2F_featmaps.py` |

**2C clustering.** Rows (features) are ordered by hierarchical clustering — Ward linkage on Euclidean
distances between each feature's z-score profile across the eight mutants; the bundled matrix is already
standardized, so no further scaling is applied. Columns keep the fixed strain order and are *not*
clustered. Because clustering breaks up the feature families, family membership moves from contiguous
brackets to the color strip on the right (same three classes/colors as Fig. S2D). vpsL's 12 colony
features are undefined (it forms no microcolonies): they are drawn grey and treated as the row mean
(0 in z-space) for the linkage only. The pre-clustering family-ordered layout is still renderable via
`2C_trainingHeatmap_familyOrder.py`.

Growth filter: 2D and 2E keep a well if `mutant==vpsL` OR its max biomass ≥ 0.15× the per-plate WT
median max (1149 wells). 2B/2C/2F use all replicates of the 8 mutants. UMAP (2D, nn=25/md=0.25/rs=0) and
the RF CV (2E) are library-version sensitive, so the figures are pinned to the tables in `data/`.

## Column dictionaries

Shared IDs: **`mutant`** (one of the 8 genotypes), **`plateId`**, **`wellId`** (mag suffix stripped).

- **`biomassTraces_normWTpeak.csv`** (one row per replicate well) — `mutant, plateId, wellId,
  biomass_t0…biomass_t30`; biomass is divided by the WT peak mean (mean over WT wells of each well's max
  biomass), so WT peaks ≈ 1.
- **`trainingHeatmap_peakBiomass_featureMatrix.csv`** — row index = 27 feature bases (biomass, 12 colony,
  global entropy, 13 haralick); columns = 8 mutants; cell = z-score (across the 8 mutants) of the median
  value at that mutant's peak-biomass frame. **`_peakFrames.csv`**: `mutant, peakFrame` (hour of peak).
- **`trainingUmap_all_three_coords.csv`** — `plateId, wellId, mutant, umap1, umap2` for each
  growth-filtered well (UMAP on standardized biomass[log1p] + whole-image + colony features, frames 9–27).
- **`all_confusion_cv.csv`** — 8×8 mean row-normalized confusion matrix (rows = true, cols = predicted);
  the panel's balanced accuracy = mean of the diagonal.
- **`featmap_<feat>.csv`** (×4: `colony_meanIntensity_mean`, `nColonies`, `colony_eccentricity_mean`,
  `whole_haralick_8`) — row index = 8 mutants, columns = frame (h), cell = median across replicates
  (colony features: timepoints with < 10 replicates above the biomass threshold are blank/masked).
  **`featmaps_meta.csv`**: `feature, label, unit, group, rfAccuracy` (single-feature RF balanced accuracy
  shown in each panel title).
