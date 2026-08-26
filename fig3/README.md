# Figure 3: Reimaging Set UMAP Landscape (paper figure package)

Self-contained package for the *V. cholerae* reimaging-atlas UMAP figure. Every panel is built from
**one shared manifold**: the reimaging atlas embedded with UMAP at **n_neighbors = 10, min_dist = 0.1**
(3669 replicate wells, 158 mutants incl. WT). All panels share the same per-replicate `umap1/umap2`
coordinates; they differ only in how points are colored/aggregated. Every plotted quantity has a
companion source-data table — open a figure and its `*_coords.csv` / `*_centroids.csv` and each dot/cell
is a row.

## Layout

```
figlib.py            thin shim — fig3 paths + re-export of ../figlib_shared.py
data/                source-data tables (CSVs + linkage.npy + embeddings parquet); also on KiltHub
render/              draw each panel from data/  → figures/   (what you run)
build/               regenerate data/ from the original feature set (KiltHub; optional)
figures/             rendered PNG/SVG/PDF output
```
(Shared `requirements.txt` / `environment.yml` / `LICENSE` live at the `paper-figures/` root.)

## Reproduce the figures from processed numerical data (no raw data)

```bash
conda env create -f environment.yml && conda activate paper-figures   # or: pip install -r requirements.txt
python render/3A_rankedTransposonBiomass.py                   # Panel 3A (+ its separate legend)
python render/3B_reimagingUmap_functionalAnnotations.py      # Panel 3B (main landscape)
python render/3Btop_reimagingUmap_centroidsByFunction.py     # Panel 3B top inset
python render/3Bbottom_reimagingUmap_centroidsByLocus.py     # Panel 3B bottom inset
python render/3D_dendrogramHeatmap.py                        # Panel 3D (dendrogram + heatmap)
python render/reimagingUmap_perGenePdf.py                    # Supplement (per-gene PDF)
```

Each reads only the bundled `data/` tables and writes to `figures/`. To regenerate the tables
themselves from the raw feature matrix, see `build/README.md`. The figures assume the **Gillius ADF**
font (matplotlib falls back silently if it's missing; data/coordinates are unaffected).

---

## From raw data to panels (full chain)

Each layer is a deposited artifact; each arrow is a script. You can enter at any level.

1. **Raw brightfield timelapses** → *([µPULLI-I](https://github.com/BridgesLabCMU/uPULLI-I))* →
   **processed images + segmentation masks** (BioImage Archive [`S-BIAD3830`](https://doi.org/10.6019/S-BIAD3830)).
2. Processed images → *(µPULLI feature extraction — `multiWellAnalysis` + DINOv2, same pipeline repo)* →
   **`master_frame_features.csv`** — the *original full* per-(well, frame) numerical feature set
   (+ per-colony table + embeddings) (CMU KiltHub `[DOI]`).
3. `master_frame_features.csv` → **`build/0_buildCollapsedWide.py`** →
   **`reimaging_collapsedWide.parquet`** (wide table: strip mag suffix, canonical renames, pivot to
   `<feature>_t<frame>`, join the gene index).
4. Wide table → **`build/3*.py`** (UMAP fit, per-mutant centroids, PCA→Ward clustering) → the
   **source-data CSVs in `data/`** (KiltHub `[DOI]` + bundled in this repo).
5. `data/` CSVs → **`render/3*.py`** → the **panels** in `figures/`.

Steps 1–2 are the µPULLI pipeline and live in separate repositories ([µPULLI-I](https://github.com/BridgesLabCMU/uPULLI-I) for imaging + features, [µPULLI-DL](https://github.com/BridgesLabCMU/uPULLI-DL) for the DINOv2 embeddings); **steps
3–5 are in this repo.** `render/` (step 5) needs only the bundled `data/`; `build/` (steps 3–4) needs
the KiltHub feature set — see `build/README.md`. Fill the bracketed deposit IDs at submission.

### Panels → scripts → tables

**3A** was previously supplemental (Fig S4D) and now opens Figure 3, replacing the UMAP-colored-by-peak-
biomass panel that used to be 3A. That retired panel now lives in **`archive/`**
(`archive/{render,build}/reimagingUmap_coloredByBiomass.py`) — still runnable for reference, but not a
panel; its source table `coloredByBiomassNormWT_plasma_nn10_md0.10_coords.csv` stays in `data/` so it
can render. **3C** is
representative microscopy and has no script here — those stills come from the image tree via
`v2/reimaging/representativeImages/pullRepresentativeFrames.py`, outside this package.

3A's two tables are produced by **fig3's own** build step from `transposons/results_10x.csv`, the screen
analysis of record: the calls are read as given (they selected the reimaging library, so they are a fixed
historical result), and `B_max`/`B_min`/`maxFinal` are recovered as the midpoints of the gaps that bracket
each class boundary. figS4's build parses the same table for its trajectory panels — a deliberate
duplication, so neither package reaches into the other's `data/`. 3A writes its legend as a **separate**
pair of files (`3A_rankedTransposonBiomass_legend.{png,svg}`).


| Panel | `render/` script | `data/` table(s) | `build/` script |
|---|---|---|---|
| **3A** ranked transposon peak biomass | `3A_rankedTransposonBiomass.py` | `tn_phenotype_summary.csv`, `tn_thresholds.csv` | `3A_rankedTransposonBiomass.py` |
| **3B** functional landscape | `3B_reimagingUmap_functionalAnnotations.py` | `reimagingLandscape_nn10_md0.10_coords.csv` | `3B_reimagingUmap_functionalAnnotations.py` |
| **3B** top inset (by function) | `3Btop_reimagingUmap_centroidsByFunction.py` | `centroidsByFunction_nn10_md0.10_centroids.csv` (+ landscape coords) | `3Btop_reimagingUmap_centroidsByFunction.py` |
| **3B** bottom inset (by locus) | `3Bbottom_reimagingUmap_centroidsByLocus.py` | `centroidsByLocus_nn10_md0.10_centroids.csv` (+ landscape coords) | `3Bbottom_reimagingUmap_centroidsByLocus.py` |
| **3C** representative images | — | — | not in this repo (image panels) |
| **3D** dendrogram + heatmap | `3D_dendrogramHeatmap.py` | `functional_peakBiomass_featureMatrix.csv`, `functional_pcaLinkage_cluster_order.csv`, `functional_pcaLinkage_linkage.npy` | `3D_generateDendogramHeatmaps-withPCA.py` |
| **Supp.** per-gene PDF | `reimagingUmap_perGenePdf.py` | `reimagingUmap_nn10_md0.10_perGene_coords.csv` | `reimagingUmap_perGenePdf.py` |

**The manifold fit** (step 4, applied identically for every panel): feature filter frames 9–27 (biomass
+ whole-image entropy/haralick; drop skew/kurtosis/cv) → drop 4 mislocalized loci → `fillna(0)` → drop
zero-variance → growth filter (max biomass > 0.005) → min-replicate filter (≥5 reps) → `StandardScaler`
→ `UMAP(nn=10, md=0.1, metric='euclidean', random_state=42)`. UMAP is **not** bit-stable across
`umap-learn`/`numba`/`numpy` versions, so the figures are pinned to the coordinates in `data/`; `build/`
is a best-effort re-derivation (run it in the pinned `environment.yml`).

---

## Column dictionaries

Common identifier columns (shared by most tables):
- **`plateId`** — reimaging plate (drawer + acquisition timestamp), e.g. `Plate1_Drawer1 02-Oct-2025 11-48-52`.
- **`wellId`** — 96-well position, magnification suffix stripped (`A10_02` → `A10`).
- **`mutant`** — gene name if annotated, else the locus tag; `WT` for wild type. One transposon insertion per non-WT mutant.
- **`geneLocus`** — VC locus tag (`VC_1115`), `''` if unknown. `WT` rows may be blank.
- **`function`** — free-text gene annotation (from `reimagingIndex.csv`); may contain commas (CSV is quoted).
- **`umap1`, `umap2`** — UMAP coordinates (nn=10, md=0.1). Per-replicate in `*_coords.csv`; per-mutant mean/median in `*_centroids.csv`.

### `reimagingLandscape_nn10_md0.10_coords.csv`:  main UMAP landscape (3669 rows, one per replicate)
| column | meaning |
|---|---|
| plateId, wellId, mutant, geneLocus, function | identifiers (above) |
| n_neighbors, min_dist | UMAP params (10, 0.1) — constant, kept for provenance |
| umap1, umap2 | replicate coordinates |

Coloring in the figure: grey = unclassified; six functional-group colors (see below); black = WT.

### `coloredByBiomassNormWT_plasma_nn10_md0.10_coords.csv`: biomass overlay (3669 rows)
| column | meaning |
|---|---|
| plateId, wellId, mutant, geneLocus, function, umap1, umap2 | as above |
| **peakBiomass** | max biofilm biomass over the timecourse for that well (a.u.) |
| **peakBiomassNorm** | `peakBiomass` ÷ (mean WT peak biomass) — the value mapped to the plasma colormap (2nd–98th pct range) |

### `centroidsByFunction_nn10_md0.10_centroids.csv`: per-mutant centroids (158 rows, one per mutant)
| column | meaning |
|---|---|
| mutant, geneLocus | identifiers |
| **nReps** | number of surviving replicate wells averaged into the centroid |
| umap1, umap2 | **mean** UMAP position of the mutant's replicates (the plotted centroid) |
| **functionalGroup** | one of the six highlighted pathways, `WT`, or `Unclassified` |
| **color** | hex color used for the centroid dot (matches `functionalGroup`) |

### `centroidsByLocus_nn10_md0.10_centroids.csv`: centroids colored by genome position (158 rows)
| column | meaning |
|---|---|
| mutant, geneLocus, nReps, umap1, umap2 | as in centroidsByFunction |
| **isChrII** | `True` if the locus is on chromosome II (`VC_Axxxx`), else chromosome I |
| **locusNumRaw** | numeric part of the locus tag (per-chromosome), drives the split Chr I / Chr II colormaps (`_byChr` plot) |
| **locusNum** | Chr I locus number as-is; Chr II shifted up by max Chr I number → single continuous scale for the viridis plot |

### `reimagingUmap_nn10_md0.10_perGene_coords.csv`: per-gene supplement PDF (3669 rows)
| column | meaning |
|---|---|
| plateId, wellId, mutant, geneLocus, function, umap1, umap2 | as above |
| **functionalGroup** | pathway / `WT` / `Unclassified` |
| **dotColor** | hex highlight color used on that gene's PDF page (functional-group color, or `#ff00c3` magenta if unclassified, black for WT) |

Each PDF page = the full landscape (grey + faint WT) with one gene's replicates highlighted.

### Dendrogram / heatmap tables (`results/v2/reimaging/dendogram/functionalAnnotations_pcaLinkage/`)
Clusters the **functional-subset** = the 6 highlighted pathways' mutants + WT (50 leaves), via
PCA-50 → per-mutant median centroid → Ward linkage. The heatmap shows, per mutant, its features at the
mutant's own peak-biomass frame, z-scored per feature.

- **`functional_peakBiomass_featureMatrix.csv`** : the heatmap values. Row index = 27 feature bases
  (biomass, 12 colony features, global entropy, 13 Haralick); columns = all 158 mutants; each cell =
  z-score (across the 158 mutants) of that mutant's median value at its peak-biomass frame. The rendered
  heatmap shows the 50-mutant functional subset in dendrogram-leaf order (see cluster_order).
- **`functional_pcaLinkage_cluster_order.csv`** (50 rows, dendrogram left→right leaf order) —
  `mutant`, `display` (gene name shown), `gene` (locus), `peakFrame` (hour of peak biomass),
  `annotation` (pathway / `WT` / `Compound` / `Other`), `color` (leaf/strip hex).
- **`functional_pcaLinkage_centroids.csv`** (50 rows) — `mutant` + `PC1…PC50`: the PCA centroids that
  were clustered.
- **`functional_pcaLinkage_explainedVariance.csv`** (50 rows): `pc`, `explainedVarianceRatio`,
  `cumulativeExplainedVariance` for the PCA.
- **`functional_peakBiomass_frames.csv`** (158 rows): `mutant`, `peakFrame` (peak-biomass hour used
  to sample features).
- **`functional_pcaLinkage_linkage.npy`** — the SciPy Ward linkage matrix (binary, not CSV).

---

## Functional groups & colors (shared across panels)

| group | color | loci |
|---|---|---|
| Motility | `#ff0004` | VC_2059…VC_2208 (26 loci) |
| O-Antigen Biosynthesis | `#0096ff` | VC_0212…VC_0269 (12) |
| Polyamine Import | `#14f7f0` | VC_1424, VC_1426, VC_1427, VC_1428 |
| Biotin Biosynthesis | `#ff9f1c` | VC_1111, VC_1113, VC_1114, VC_1115 |
| Pyruvate Flux | `#39ff14` | VC_2413, VC_0943 |
| Vibriobactin Biosynthesis | `#ba17f6` (`#a200ff` in the dendrogram) | VC_0771, VC_0772 |
| WT | `#000000` | — |
| Unclassified | grey (`#939090` / `#d0d0d0`) | everything else |

(The exact locus lists are duplicated at the top of each script and match `v2/common/plotting.HIGHLIGHT_SETS`.)
