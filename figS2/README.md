# Supplemental Figure S2: Training-Set Classification, Dissected (Numerical Features)

Supplement to Figure 2. Where Fig 2E reports the all-feature RF confusion matrix, S2 breaks the
hand-engineered numerical classification of the 8-mutant *V. cholerae* training set into its parts:
which feature families carry the genotype signal, which individual features are most discriminative,
and when in the timecourse the mutants become separable.

## Layout

```
figlib.py            thin shim: figS2 paths + strain order/labels/colors + re-export of ../figlib_shared.py
data/                source-data tables (CSVs); also on KiltHub
render/              draw each panel from data/  -> figures/   (what you run)
build/               regenerate data/ from the training wide table (KiltHub; heavy, optional)
figures/             rendered PNG/SVG output
```
(Shared `requirements.txt` / `environment.yml` / `LICENSE` live at the `paper-figures/` root.)

## From raw data to panels (full chain)

1. **Raw brightfield timelapses** -> *(µPULLI image + feature pipeline, separate repo `[pipeline repo]`)* ->
   **`training_wide.parquet`** (per-well `<feature>_t<frame>` matrix; KiltHub `[DOI]`).
2. `training_wide.parquet` -> **`build/S2*.py`** (GroupKFold RF classification) -> **source-data CSVs in
   `data/`**.
3. `data/` CSVs -> **`render/S2*.py`** -> **panels** in `figures/`.

### Panels -> scripts -> tables

| Panel | `render/` script | `data/` table | `build/` script |
|---|---|---|---|
| **S2A** biomass-only RF confusion | `S2A_bioConfusion.py` | `bio_confusion_cv.csv` | `S2ABC_confusion.py` |
| **S2B** colony-level (segmented) RF confusion | `S2B_colonyConfusion.py` | `colony_confusion_cv.csv` | `S2ABC_confusion.py` |
| **S2C** whole-image RF confusion | `S2C_wholeConfusion.py` | `whole_confusion_cv.csv` | `S2ABC_confusion.py` |
| **S2D** single-feature timecourse accuracy | `S2D_singleFeatureAccuracy.py` | `singleFeatureAccuracy.csv` | `S2D_singleFeatureAccuracy.py` |
| **S2E** mutant separability across time | `S2E_timepointSeparability.py` | `timepointSalience_groupkfold.csv` + `fullTimecourse_accuracy.csv` | `S2E_timepointSeparability.py` |

**Shared CV protocol** (all panels): GroupKFold by `plateId` (5 folds × 5 repeats, plate-shuffled),
RandomForest (200 trees, `min_samples_leaf=2`, `class_weight='balanced'`), per-fold StandardScaler,
**balanced accuracy**, predicting the 8-mutant genotype. Growth filter (shared with Fig 2): keep a well
if `mutant==vpsL` OR max biomass ≥ 0.15× the per-plate WT median max. Biomass is log1p'd. Colony
features use frames 9–30 (colonies segmentable from t9); biomass and whole-image features use 0–30.
The RF CV is library-version sensitive, so the figures are pinned to the tables in `data/`.

## Column dictionaries

- **`{bio,colony,whole}_confusion_cv.csv`** — 8×8 mean row-normalized confusion matrix (rows = true,
  cols = predicted) for that single feature family; each panel's balanced accuracy = mean of the diagonal.
- **`singleFeatureAccuracy.csv`** — `featureBase`, `prettyName`, `family` (`biomass`/`whole`/`colony`),
  `balancedAccuracy`; one row per feature base, each trained as the sole predictor across its timecourse.
- **`timepointSalience_groupkfold.csv`** — `frame`, `rfMean`, `rfStd`: RF balanced accuracy using only
  the features at that single frame.
- **`fullTimecourse_accuracy.csv`** — `rfMean`, `rfStd`: the full-timecourse RF baseline (all frames).
