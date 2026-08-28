# Figure 4: Clean-Deletion Projection

In-frame clean deletions (Δ*bioD*, Δ*manA*, Δ*pdhE2*) validated against the reimaging transposon
landscape: each clean deletion is projected onto the Figure-3 UMAP (4A), its biofilm-biomass trajectory
is compared to reimaging WT (4C left), and its peak-feature profile is shown relative to WT (4C right).

## Layout

```
figlib.py            thin shim: fig4 paths + clean-deletion markers/colors + re-export of ../figlib_shared.py
data/                source-data tables (CSVs); also on KiltHub
render/              draw each panel from data/  -> figures/   (what you run)
build/               regenerate the 4C tables from the clean-deletion + reimaging feature tables (KiltHub)
figures/             rendered PNG/SVG output
```
(Shared `requirements.txt` / `environment.yml` / `LICENSE` live at the `paper-figures/` root.)

## From raw data to panels (full chain)

1. **Raw brightfield timelapses** -> *([µPULLI-I](https://github.com/BridgesLabCMU/uPULLI-I))* ->
   **processed images + masks** (BioImage Archive [`S-BIAD3830`](https://doi.org/10.6019/S-BIAD3830)).
2. Processed images -> *(µPULLI feature extraction)* -> the clean-deletion plate feature tables +
   the reimaging atlas feature table (CMU KiltHub `[DOI]`).
3a. **4A:** clean-deletion wells are transformed onto the reimaging UMAP by the projection step
   (`projectOnReimaging`, scaler + fitted UMAP `transform`) -> **`cleanDeletions_projectedCoords.csv`**.
   The 4A background is the Figure-3 reimaging landscape (`reimaging_landscape_coords.csv`).
3b. **4C:** clean-deletion + reimaging feature tables -> **`build/4C_*`** -> the 4C source tables.
4. `data/` CSVs -> **`render/4*.py`** -> **panels** in `figures/`.

Steps 1-2 (and the 4A projection in 3a) are upstream analysis; `render/` needs only `data/`. Fill
bracketed IDs at submission.

### Panels -> scripts -> tables

| Panel | `render/` script | `data/` table(s) | `build/` script |
|---|---|---|---|
| **4A** projection onto reimaging landscape | `4A_projection.py` | `cleanDeletions_projectedCoords.csv` + `reimaging_landscape_coords.csv` | upstream projection (see note) |
| **4C left** biomass over time | `4C_left_biomassOverTime.py` | `biomassOverTime_normWTpeak.csv` | `4C_left_biomassOverTime.py` |
| **4C right** WT-normalized feature heatmap | `4C_right_vsWTheatmap.py` | `cleanDel_vsWT_horizontal_matrix.csv` | `4C_right_vsWTheatmap.py` |
| **4C right (alt)** same matrix, conditions reordered Δ*bioD* / Δ*pdhE2* / Δ*manA* | `4C_right_vsWTheatmap_ordered.py` | same table | same build |
| **4E** RNA-seq biofilm-gene heatmap | `4E_rnaseqHeatmap.py` | `rnaseq_logFC_matrix.csv` | `4E_rnaseqHeatmap.py` |

**Note on 4A:** `cleanDeletions_projectedCoords.csv` is produced by the projection step of the analysis
pipeline (clean-deletion wells passed through the reimaging canonical scaler + fitted UMAP `transform`,
nn=10/md=0.1). That fitted reducer is an upstream artifact, so the table is provided as source data
rather than rebuilt here. The three clean deletions land on their transposon counterparts (BioD→bioD,
ManA→manA, PdhE2→pdhE2).

## Column dictionaries

- **`cleanDeletions_projectedCoords.csv`** - `srcPlate, plateId, wellId, mutant, umap1, umap2`; each
  clean-deletion replicate's coordinates in the reimaging UMAP. 4A subsamples 9 reps/mutant
  (random_state=42) and draws open markers (Δ*bioD* triangle, Δ*manA* square, Δ*pdhE2* diamond).
- **`reimaging_landscape_coords.csv`** - the Figure-3 reimaging background (`plateId, wellId, mutant,
  geneLocus, function, n_neighbors, min_dist, umap1, umap2`); grey + WT/Biotin/Pyruvate/O-Antigen shown.
- **`biomassOverTime_normWTpeak.csv`** - `group` (WT / BioD / ManA / PdhE2), `frame`, `mean`, `sd`, `n`;
  biofilm biomass normalized to the reimaging WT median peak (mean ± **SD** across wells per frame).
  SD rather than SEM: the panel's claim is about replicate spread, and SEM (= SD/√n, n = 24–32) is
  ~5–6× narrower — it would imply precision the biology does not have. SEM would also not mean the same
  thing across series here: reimaging WT's 24 wells come from 24 different plates, whereas each clean
  deletion's 32 wells come from only 2 (16 per plate), so their effective n is far below 32.
- **`cleanDel_vsWT_horizontal_matrix.csv`** - rows = clean-deletion conditions, columns = feature bases
  (27, incl. Global Entropy); cell = (condition − WT) / reimaging-atlas σ at the condition's peak-biomass
  frame (WT dropped). Corrects the original `whole_entropy_` filter typo, so Global Entropy is included.
