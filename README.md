# Paper figures — reproduction code + source data

Self-contained code and source-data tables that regenerate the figures for the µPULLI manuscript
(*V. cholerae* biofilm phenotyping). One package per figure — `fig1/`…`fig5/` and `figS1/`…`figS8/`,
plus `interactive/` for the HTML supplements; a shared library (`figlib_shared.py`) holds the
functional-group colors, the UMAP feature filter, and the plot style so they are defined once.

## Layout

```
paper-figures/
  figlib_shared.py     shared: colors / feature filter / plot style / per-figure config factory
  inputs.json          manifest of build-layer inputs (logical name → size, SHA-256, figures)
  INPUTS.md            that manifest as a table, plus how to get the inputs
  fetch_data.py        download the build inputs from the data deposit
  reproduce/           deposit_to_cytation.py: deposited raw stacks → µPULLI-I input (see Reproducing the paper)
  requirements.txt  environment.yml  LICENSE (MIT)  .gitignore
  fig1/ … fig5/  figS1/ … figS8/     one package per figure:
    figlib.py            thin shim — declares this figure's build inputs, re-exports figlib_shared
    data/                source-data tables (CSVs + small arrays; also on KiltHub)
    render/              draw each panel FROM data/  → figures/   (fast, deterministic; what you run)
    build/               regenerate data/ FROM the deposited feature data (optional)
    figures/             rendered PNG/SVG output
    README.md            per-figure column dictionary + panel details
  interactive/           the six interactive HTML supplements (build/ only — the page IS the output)
  archive/movies/        retired Movie S3; builds, but nothing in the paper depends on it
```

## System requirements

- **Operating system:** any system that runs conda or Python 3.9 (Linux, macOS, Windows). Tested on
  Ubuntu 24.04.3 LTS; the code is pure Python with no compiled parts of its own and no absolute paths.
- **Software:** Python 3.9 and the versions pinned in [`environment.yml`](environment.yml) /
  [`requirements.txt`](requirements.txt). Tested with Python 3.9.23, numpy 2.0.2, pandas 2.3.1,
  matplotlib 3.9.4, scipy 1.13.1, and, for `build/` only, scikit-learn 1.6.1, umap-learn 0.5.9.post2,
  numba 0.60.0 and pyarrow 21.0.0. Optional: the **Gillius ADF** font (text spacing only).
- **Hardware:** a standard desktop computer; no GPU. Most panels need under 1 GB of RAM. The full-atlas
  dendrogram heatmaps need more: Fig 3D peaks at 4.5 GB and S6A at 12 GB, so rendering every panel
  needs **16 GB of RAM**.
- **Install time:** about 2–5 minutes (`conda env create`, measured at 2 min on a fresh package cache).

## Reproduce a figure

```bash
conda env create -f environment.yml && conda activate paper-figures   # or: pip install -r requirements.txt
python fig3/render/3B_reimagingUmap_functionalAnnotations.py          # a Fig-3 panel; see each fig*/README.md
```

That is the whole setup: render scripts read only the bundled `fig*/data/` tables, so there is nothing to
download, nothing to configure, and no network access. Output is deterministic.

To regenerate those tables from the underlying feature data, fetch the build inputs first — see
[`INPUTS.md`](INPUTS.md):

```bash
python fetch_data.py            # ~996 MB into fig*/build/inputs/  (or set UPULLI_DATA_ROOT)
python fig3/build/0_buildCollapsedWide.py
```

No script contains an absolute path: every build input is resolved from `inputs.json` by logical name.

## Demo

Run times below are for a desktop computer with an 8-core CPU (measured with the pinned environment,
restricted to 8 cores).

**1. Render one panel from the bundled data** (about 3 seconds, 0.3 GB RAM):

```bash
python fig3/render/3B_reimagingUmap_functionalAnnotations.py
```

Expected output: `fig3/figures/3B_landscape.png` and `.svg`, the reimaging phenotype landscape of
Fig. 3B. Render scripts are deterministic, so the panel matches the published one; text spacing differs
slightly if the Gillius ADF font is not installed.

**2. Rebuild a figure's data from the deposit, then render it** (about 15 seconds plus a 3 MB download):

```bash
python fetch_data.py figS8                 # Fig. S8 inputs from KiltHub (SHA-256 checked on download)
python figS8/build/S8BCD_volcano.py && python figS8/build/S8E_luxActivity.py
python figS8/render/S8BCD_volcano.py && python figS8/render/S8E_luxActivity.py
```

Expected output: the four tables in `figS8/data/` are rewritten **byte-identical** to the copies
shipped in the repository (`git status` shows no change), and the RNA-seq volcanoes (S8B–D) and
P*vpsL*-lux panel (S8E) appear in `figS8/figures/`.

**3. Every panel** (all 51 render scripts):

```bash
for f in fig*/render/*.py; do python "$f"; done
```

About **35 minutes** in total. The 45 ordinary panels take 2–45 seconds each (under 5 minutes
together); the six full-atlas dendrogram heatmaps (Fig. 3D, S6A ×2, S7A ×2, S7C) take 1–9 minutes each
and account for the rest.
Figures assume the **Gillius ADF** font (matplotlib falls back silently if it is missing; the data and
coordinates are unaffected).

## Reproducing the paper

The paper's numbers come from three stages. Each can be rerun from public material, and each was
checked against the published results:

| Stage | Tool | Check |
|---|---|---|
| A. Raw images → per-frame features | [µPULLI-I](https://github.com/BridgesLabCMU/uPULLI-I) | Two wells each from the training, screen, reimaging, *K. pneumoniae* and multispecies (10×) sets, reprocessed with the commands below (the screen and reimaging wells starting from the deposited raw stacks): all 191 features **identical** to the deposited tables |
| B. Processed images → DINOv2 embeddings | [µPULLI-DL](https://github.com/BridgesLabCMU/uPULLI-DL) | Reimaging wells re-embedded: every frame's vector matches the deposited one to within 4 × 10⁻⁶ (relative), cosine similarity 1.0000000; the residual is floating-point difference between PyTorch versions |
| C. Features and embeddings → tables → figures | this repository | Training and reimaging wide tables rebuilt **identical**; figure tables rebuilt **byte-identical** (Fig. S8 tested end to end) |

The paper's data were processed with development builds of µPULLI-I (0.5.0 and 1.0.0). Use the current
release: it gives identical features on the same images. For exact numerics, install µPULLI-I from its
pinned `environment.yml` (Python 3.9.23).

### A. Raw images → features (µPULLI-I)

The raw brightfield timelapses are in the BioImage Archive entry
[`S-BIAD3830`](https://doi.org/10.6019/S-BIAD3830) under `raw-images/`, one OME-TIFF stack per well.
µPULLI-I reads per-frame Cytation 5 files instead, so convert first. `reproduce/deposit_to_cytation.py`
writes every frame with its original pixel values and the pixel size from the deposit. It needs only
numpy and tifffile; run it in the µPULLI-I environment.

```bash
# 1. Convert the plates of one dataset (example: the transposon reimaging set)
python reproduce/deposit_to_cytation.py <deposit>/raw-images/Vcholerae/transposonReimaging \
    --out raw_frames/reimaging --mag _02

# 2. Process them: the flags common to every dataset, plus that dataset's row from the table below
biofilm-processing-run --plates raw_frames/reimaging/* --output-dir features/reimaging --mag _02 \
    --block-diam 101 --fft-stride 1 --shift-thresh 250 \
    --fixed-thresh 0.025 --downsample 1 --min-colony-area 200 --prop-radius 50 \
    --whole-image --colony-tracking --colony-feats --workers 8
```

The run writes `features/<dataset>/master_frame_features.csv`, the table deposited on KiltHub.

| Dataset | Deposit folder (`raw-images/…`) | `--mag` | `--fixed-thresh` | `--downsample` | `--min-colony-area` | `--prop-radius` | Deposited output | Figures |
|---|---|---|---|---|---|---|---|---|
| Training set (8 strains) | `Vcholerae/knownMutant` | `_03` | 0.025 | 1 | 200 | 50 | `training/master_frame_features.csv` | 1, 2, S1–S5 |
| Genome-wide transposon screen | `Vcholerae/transposonScreen` | `_03` | 0.025 | 1 | 200 | 50 | `transposons/master_frame_features.csv` | S4, S5, Data S1 |
| Transposon reimaging | `Vcholerae/transposonReimaging` | `_02` | 0.025 | 1 | 200 | 50 | `reimaging/master_frame_features.csv` | 3, 4, 5, S4, S6, S7, Interactive 1–3 |
| Clean deletions, compounds | `Vcholerae/cleanDeletion`, `Vcholerae/compound` | `_02` | 0.025 | 1 | 200 | 50 | input to `cluster/*` (see "Deposited as is") | 4, 5A |
| *K. pneumoniae* mutants | `Kpneumoniae` | `_02` | 0.04 | 4 | 200 | 25 | `kleb/master_frame_features.csv` | 5B |
| Multispecies panel, 10× | `Multispecies` | `_04` | 0.03 | 4 | 400 | 25 | input to stage B | 5C |

The multispecies deposit is filed by species. Point the converter at the whole `Multispecies` folder and
it reassembles each plate. Processing takes about 5 minutes per 10× well (1992 × 1992 px, 31 frames)
per worker and about 3 GB of RAM per worker; the sets run with `--downsample 4` are faster.

### B. Processed images → embeddings (µPULLI-DL)

Run µPULLI-DL on each stage-A output folder. Its defaults are the paper's settings: `facebook/dinov2-base`
(ViT-B/14), 518 px input, 3 × 3 patch grid.

```bash
biofilm-embeddings-run features/training          # likewise features/reimaging and features/kleb
biofilm-embeddings-run features/cluster --n-frames 31
biofilm-embeddings-run features/multispecies --mag _04 --n-frames 24 --cache-dir features/multispecies/10X
```

Each run writes `embeddings/cls_cache.pt`. The deposited arrays are a format conversion of that cache:
`cls.npy` is `cache['cls']` as float32 (wells × frames × 768), `patchmean.npy` is `cache['patches']`
averaged over its 9 patch cells, and `index.csv` lists `cache['wells']` in row order. Read the cache with
`biofilm_embeddings.embeddings.extractor.loadCache`. A GPU is strongly recommended for the full sets.

### C. Features → tables → figures (this repository)

Put the stage-A and stage-B outputs where the build layer expects them: under the logical names in
[`INPUTS.md`](INPUTS.md), in `fig*/build/inputs/` or under `UPULLI_DATA_ROOT`. Alternatively run
`python fetch_data.py` to take them from KiltHub. Then:

```bash
python fig2/build/0_buildTrainingWide.py     # training/wide.parquet    (verified identical to the deposit)
python fig3/build/0_buildCollapsedWide.py    # reimaging/collapsedWide.parquet (verified identical)
```

Then run each figure's `build/` scripts in the order its `fig*/README.md` gives, and the `render/`
scripts as in [Demo](#demo). Analysis choices (growth filters, replicate minimums, excluded loci,
frame windows) are in the build scripts and documented in each `fig*/README.md`.

### Deposited as is

These KiltHub inputs are used as published. The code in this repository and the two tools does not
regenerate them.

| Input | What it is | Origin |
|---|---|---|
| `transposons/results_10x.csv`, `reimaging/reimagingResults_10x.csv` | Phenotype calls from the original screen analysis; these selected the reimaging library | Screen analysis of record |
| `training/layout.csv`, `reimaging/geneIndex.csv` | Well → strain and well → gene/function annotations | Plate layouts and genome annotation |
| `reimaging/umapEmbeddings.parquet`, `reimaging/embeddings/dendro*.csv`, `cluster/*_collapsedWide.parquet`, `reimaging/thumbnails*/` | UMAP grid, embedding-dendrogram inputs, clean-deletion and compound wide tables, image thumbnails | Lab analysis scripts not included in this repository |
| `rnaseq/` | Differential-expression tables (Fig. 4E, S8B–D) | RNA-seq pipeline in the Methods; reads on NCBI SRA, BioProject `PRJNA1513060` |
| `lux/Lum.csv`, `lux/OD.csv` | Plate-reader P*vpsL*-lux timecourse (Fig. S8E) | Raw instrument export |

## Figures

| Figure | Package | Status |
|---|---|---|
| 1 — training set, DINOv2 embeddings | `fig1/` | ✅ complete |
| 2 — training set, numerical | `fig2/` | ✅ complete |
| 3 — reimaging phenotype landscape | `fig3/` | ✅ complete |
| 4 — clean-deletion projection | `fig4/` | ✅ complete (4E RNA-seq pending) |
| 5 — generalization (compounds / kleb / multispecies) | `fig5/` | ✅ complete |

| Supplemental | Package | Content |
|---|---|---|
| S1 — training patch-embedding UMAP + confusion | `figS1/` | DINOv2 patch tokens |
| S2 — training classification dissected | `figS2/` | per-family confusion, single-feature accuracy, timepoint separability |
| S3 — training temporal heatmaps, all feature classes | `figS3/` | one page of per-feature time×mutant heatmaps |
| S4 — genome-wide transposon biomass screen | `figS4/` | genome heatmap, Chr I/II, ranking scatter |
| S5 — low-biofilm transposon mutants by chromosome | `figS5/` | Chr I / Chr II biomass trajectories |
| S6 — full reimaging atlas + replicate reproducibility | `figS6/` | A dendrogram/heatmap over all 158 mutants; B–E replicates are closer than chance and than other mutants |
| S7 — full reimaging atlas, DINOv2-embedding views | `figS7/` | A embedding-PC dendrogram, B embedding UMAP, C functional subset |
| S8 — clean-deletion transcriptional + reporter phenotypes | `figS8/` | B–D RNA-seq volcanoes (Δ*bioD*, Δ*pdhE2*, Δ*manA* vs WT), E peak P*vpsL*-lux activity |

| Interactive | Package | Content |
|---|---|---|
| 1 — UMAP, quantitative features (µPULLI-I) | `interactive/` | click a replicate → its peak biofilm biomass image |
| 2 — UMAP, DINOv2 CLS embeddings PCA-50 (µPULLI-DL) | `interactive/` | same explorer, embedding manifold |
| 1T, 2T — the same two UMAPs, timelapse | `interactive/` | click a representative replicate → **plays its whole timecourse**; other replicates keep their peak still |
| 3 — dendrogram + heatmap, all timepoints | `interactive/` | **toggle full atlas ↔ functional-annotation subset**; scrub or play through 0–30 h + peak; click a mutant → its replicate images, which follow the slider |
| 4 — RNA-seq volcano, Δ*bioD* vs WT | `interactive/` | hover a gene → identifiers, annotation, log₂FC, adjusted *p*; pathway highlights; gene/locus/description search |
| 5 — RNA-seq volcano, Δ*pdhE2* vs WT | `interactive/` | same explorer |
| 6 — RNA-seq volcano, Δ*manA* vs WT | `interactive/` | same explorer |

The interactive supplements are self-contained single-file HTML — for Plots 1–3 every thumbnail is
inlined as a base64 data URI, so each opens offline by double-clicking with no server and no asset
folder. They are built rather than bundled (see [`interactive/README.md`](interactive/README.md));
Plots 4–6 need no build input at all, since their tables ship in `interactive/data/`.

Plots 1T/2T are alternates of 1/2, not extra figures — ship either the still or the timelapse
version of each UMAP, and renumber to 7/8 if both go in the same supplement.

`IP12_interactiveUmap.py --nas` additionally builds lab-only variants that **link** the image trees
instead of inlining them: 12.5 MB, all 3,669 replicates playing at full resolution, but they open
only where the NAS is mounted. Useful internally, **never depositable** — see
[`interactive/README.md`](interactive/README.md).

> The dendrogram-timecourse movie is **not part of the manuscript** and has been retired to
> [`archive/movies/`](archive/movies/README.md); it still builds (it needs `ffmpeg`), but nothing in the
> paper depends on it.

Each `fig*/README.md` carries that figure's full lineage — the panel→script→table map, the column
dictionary for every table, and the processing chain from the original feature set down to the panels.

## Data / processing chain (repo-wide)

| Layer | Artifact | Deposited at | Produced by |
|---|---|---|---|
| Raw | brightfield timelapses | BioImage Archive [`S-BIAD3830`](https://doi.org/10.6019/S-BIAD3830) | — |
| Images | processed images + segmentation masks | BioImage Archive [`S-BIAD3830`](https://doi.org/10.6019/S-BIAD3830) | [µPULLI-I](https://github.com/BridgesLabCMU/uPULLI-I) |
| Features | full feature sets + gene index | CMU KiltHub [`10.1184/R1/33333711`](https://doi.org/10.1184/R1/33333711.v1) | [µPULLI-I](https://github.com/BridgesLabCMU/uPULLI-I) |
| Embeddings | DINOv2 CLS + patch descriptors | CMU KiltHub [`10.1184/R1/33333711`](https://doi.org/10.1184/R1/33333711.v1) | [µPULLI-DL](https://github.com/BridgesLabCMU/uPULLI-DL) |
| Tables | figure source-data CSVs (`fig*/data/`) | KiltHub [`10.1184/R1/33333711`](https://doi.org/10.1184/R1/33333711.v1) + this repo | `fig*/build/` |
| Datasets | **Data S1** — full transposon-screen results (every feature + phenotype call, per well and per frame) | KiltHub [`10.1184/R1/33333711`](https://doi.org/10.1184/R1/33333711.v1) | `figS4/build/DataS1_transposonScreen.py` |
| Figures | panels (`fig*/figures/`) | this repo | `fig*/render/` |

Steps above the "Tables" row are the µPULLI pipeline (a separate repository); this repo covers
Tables → Figures (and the feature-set → table reshaping in each `build/`).

## Data & Code Availability

Derived/tabular data (feature vectors, UMAP coordinates, per-microcolony measurements, biomass
trajectories) → CMU KiltHub [`10.1184/R1/33333711`](https://doi.org/10.1184/R1/33333711.v1). Processed image data + segmentation masks → BioImage Archive
[`S-BIAD3830`](https://doi.org/10.6019/S-BIAD3830); raw timelapses in the same BioImage Archive deposit. Raw RNA sequencing reads → NCBI SRA, BioProject `PRJNA1513060`. Code (this repository) →
Bridges Lab GitHub, MIT license.
