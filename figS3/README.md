# Supplemental Figure S3: Temporal Heatmaps for All Feature Classes (Training Set)

The all-feature-class extension of Figure 2F. Where Fig 2F shows four curated feature classes, S3 renders
a temporal (mutant × frame) heatmap for **every** quantitative feature class — biofilm biomass, the 14
whole-image texture features (global entropy + 13 Haralick descriptors), and the 12 colony-segmentation
features (27 total) — in the identical Fig-2F style, so the four features shared with Fig 2F render
identically. Each heatmap is emitted as its own PNG **and** SVG for single-page layout.

## Layout

```
figlib.py            thin shim: figS3 paths + strain order/labels + re-export of ../figlib_shared.py
data/                source-data tables (CSVs): featmap_<feat>.csv (x27) + featmaps_meta.csv
render/              draw each heatmap from data/  -> figures/   (what you run)
build/               regenerate data/ from the training wide table (KiltHub; optional)
figures/             rendered PNG + SVG output, named panel letter + feature (S3A_Average_Colony_Area, …)
```
(Shared `requirements.txt` / `environment.yml` / `LICENSE` live at the `paper-figures/` root.)

## From raw data to panels

1. `training_wide.parquet` (µPULLI feature pipeline; KiltHub `[DOI]`) -> **`build/S3_featmaps.py`** ->
   per-feature `data/featmap_<feat>.csv` matrices + `featmaps_meta.csv`.
2. `data/` CSVs -> **`render/S3_featmaps.py`** -> one `figures/S3<letter>_<FeatureName>.{png,svg}` per
   panel (Haralick panels also carry the Haralick number, e.g. `S3U_Global_Texture_Entropy_Haralick8`).
   The panel letter follows the figure's laid-out order (`PANEL_ORDER` in the render): colony features
   alphabetical (A–L), Haralick 0–12 (M–Y), global image entropy (Z), biofilm biomass (AA). The render
   prints the letter→filename manifest.

`build/S3_featmaps.py` also reads the single-feature RF accuracies from Fig S2D
(`figS2/data/singleFeatureAccuracy.csv`, override with the bundled `figS2/data/singleFeatureAccuracy.csv`) for each panel's title;
if absent, titles fall back to the feature name.

## Conventions (identical to Fig 2F)

- Rows = the 8 *V. cholerae* strains; columns = imaging time, frames 8–30 h.
- Cell = mean feature value across replicate wells (up to n = 144 per strain; 141 for ΔpotD1).
- Colormap = plasma; each heatmap self-scaled to its own min/max (colorbar in the feature's units).
- Colony-segmentation features are masked **black** at timepoints with < 10 replicates above the biomass
  threshold (colonies are not reliably segmentable before ~t8, and low-biomass wells have no colonies).
- Title = feature name (+ unit) and the single-feature RF balanced accuracy (Fig S2D).

## Single-page layout

Panels are lettered A–AA in the laid-out order: colony features alphabetical (A–L), Haralick 0–12
(M–Y), global image entropy (Z), biofilm biomass (AA). Files are named
`S3<letter>_<FeatureName>[_Haralick<N>].{png,svg}`. All heatmaps share the same figure size and frame
axis so they tile cleanly onto one page.

## Column dictionary

- **`featmap_<feat>.csv`** — 8 (strain) × 23 (frame, t8–t30) matrix of mean feature value; NaN = masked.
- **`featmaps_meta.csv`** — `feature`, `label` (pretty name), `unit`, `family`
  (`biomass`/`whole`/`colony`), `rfAccuracy` (single-feature RF balanced accuracy).
