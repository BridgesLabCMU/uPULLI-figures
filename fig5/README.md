# Figure 5: Generalization Across Systems

Three panel pairs, each a UMAP/projection (left) + a classification/heatmap (right), showing the
approach generalizes beyond the *V. cholerae* transposon atlas:
- **A** compounds - chemical treatments projected onto the reimaging landscape, and a biotin-pathway
  feature barcode.
- **B** *K. pneumoniae* - 5 kleb transposon mutants from DINOv2 embeddings (UMAP + RF confusion).
- **C** multispecies - an 8-species panel from DINOv2 embeddings (UMAP + RF confusion), 100% LB / 10X.

## Layout

```
figlib.py            shim: fig5 paths + the three datasets' display constants + re-export of ../figlib_shared.py
data/                source-data tables (CSVs); also on KiltHub
render/              draw each panel from data/  -> figures/   (what you run)
build/               regenerate the two embedding-UMAP coord tables from the CLS caches (KiltHub)
figures/             rendered PNG/SVG output
```
(Shared `requirements.txt` / `environment.yml` / `LICENSE` live at the `paper-figures/` root.)

## Panels -> scripts -> tables

| Panel | `render/` script | `data/` table(s) | source of the table |
|---|---|---|---|
| **5A left** compound projection (biotin view): WT, WT+DMSO, WT + 100 µM MAC13772, Tn *bioA–D*, Δ*bioD*, Δ*bioD*+biotin | `A_left_compoundProjection_bioDvariant.py` | `compounds_projectedCoords.csv` + `cleanDeletions_projectedCoords.csv` + `reimaging_landscape_coords.csv` | upstream projections (fitted reducer) |
| 5A-left earlier version, three compound conditions only (no untreated Δ*bioD*) — **superseded** | `A_left_compoundProjection.py` | `compounds_projectedCoords.csv` + `reimaging_landscape_coords.csv` | upstream projection (fitted reducer) |
| **5A right** biotin feature barcode (vsWT) | `A_right_biotinBarcode.py` | `compounds_biotinBarcode_vsWT_matrix.csv` | compounds vsWT-barcode analysis |
| **5B left** *K. pneumoniae* embedding UMAP | `B_left_klebUmap.py` | `kleb_embeddingUmap_coords.csv` | **`build/B_left_klebUmap.py`** |
| **5B right** *K. pneumoniae* RF confusion | `B_right_klebConfusion.py` | `kleb_embeddings_confusion_cv.csv` | kleb embedding RF (RepeatedStratifiedKFold) |
| **5C left** multispecies embedding UMAP | `C_left_multispeciesUmap.py` | `multispecies_100pctLB_10X_umap_coords.csv` | **`build/C_left_multispeciesUmap.py`** |
| **5C right** multispecies RF confusion | `C_right_multispeciesConfusion.py` | `multispecies_100pctLB_10X_confusion_cv.csv` | multispecies embedding RF (GroupKFold) |

`build/` regenerates the two embedding-UMAP coordinate tables (5B/5C left) from the DINOv2 CLS caches
(`kleb_cls.npy`, `multispecies_10X_cls.npy`; KiltHub). The three classification/barcode tables (5A right,
5B/5C right) are produced by their datasets' analysis pipelines and are provided here as source data
(the RF cross-validations are slow; the deposited CSVs are the tables). 5A-left projected coords, like
Fig 4's, come from the reimaging projection step (a fitted UMAP `transform`, an upstream artifact).

## Datasets & conventions

- **B (kleb):** 5 *K. pneumoniae* transposon mutants (waaL dropped); one mutant per plate, so labels
  come from plate-timestamp → NV id → gene, and the RF CV is **RepeatedStratifiedKFold over wells** (not
  GroupKFold-by-plate — accepted caveat, no plate-level batch control). Growth floor: max biomass ≥
  0.005, mrkA exempt (expected non-former). Embedding window: CLS frames 9–24. UMAP uses the euclid rep.
- **C (multispecies):** 8-species panel, analyzed **within** 100% LB at 10X; embedding window frames
  9–23; RF CV is GroupKFold by plate; balanced accuracy = mean of the confusion diagonal.
- **5A left** (`A_left_compoundProjection_bioDvariant.py`) is the published panel. It carries the
  **untreated Δ*bioD* clean deletion** alongside the compound conditions, so the chemical block, the deletion and the rescue appear on one
  landscape — WT+DMSO open black ○, WT + 100 µM MAC13772 open purple ◇, Δ*bioD*+biotin open red ▽,
  Δ*bioD* open grey (#787878) △; background = grey atlas, Tn *bioA–D* orange with black edges, reimaging
  WT semi-transparent black. The two Δ*bioD* conditions come from **different experiments** (treated from
  the compound plate, untreated from the clean-deletion plates) projected onto the same manifold — not a
  treated/untreated pair from one plate. The two conditions that land in the crowded biotin cluster —
  anti-biotin and Δ*bioD* — are **capped at 9 wells** each (deterministic subsample, `random_state=42`,
  Fig-4A's declutter convention); WT+DMSO and Δ*bioD*+biotin keep all 16. Purple is drawn last so it sits
  above the grey triangles it overlaps. `--reps N` thins every condition further and suffixes the output
  `_Nreps`.

## Column dictionaries

- **`compounds_projectedCoords.csv`** - `srcPlate, plateId, wellId, mutant, umap1, umap2`; compound wells
  projected onto the reimaging UMAP. 5A-left overlays WT+DMSO, WT+anti-biotin, Δ*bioD*+biotin (open markers).
- **`cleanDeletions_projectedCoords.csv`** - same schema; the Fig-4 clean-deletion wells (plates 260521 /
  260522) on the same manifold. Bundled here too so Fig 5 stands alone (a copy of `fig4/data/`'s table).
- **`reimaging_landscape_coords.csv`** - the Fig-3 reimaging background (Biotin Biosynthesis + WT shown).
- **`compounds_biotinBarcode_vsWT_matrix.csv`** - features (rows) × conditions (columns, display labels);
  cell = (condition − WT+DMSO) / reimaging-atlas σ at the condition's peak-biomass frame (WT+DMSO = zero
  baseline column). 27 features incl. Global Entropy.
- **`kleb_embeddingUmap_coords.csv`** - `metric` (cosine/euclid), `plateId, wellId, mutant, umap1, umap2`;
  the panel uses euclid.
- **`kleb_embeddings_confusion_cv.csv`** - 5×5 mean row-normalized confusion matrix.
- **`multispecies_100pctLB_10X_umap_coords.csv`** - `species, plateId, wellId, LB_condition, umap1, umap2`.
- **`multispecies_100pctLB_10X_confusion_cv.csv`** - 8×8 mean row-normalized confusion matrix.
