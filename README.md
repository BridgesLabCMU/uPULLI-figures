# Paper figures — reproduction code + source data

Self-contained code and source-data tables that regenerate the figures for the µPULLI manuscript
(*V. cholerae* biofilm phenotyping). One package per figure — `fig1/`…`fig5/` and `figS1/`…`figS7/`; a
shared library (`figlib_shared.py`) holds the functional-group colors, the UMAP feature filter, and the
plot style so they are defined once.

## Layout

```
paper-figures/
  figlib_shared.py     shared: colors / feature filter / plot style / per-figure config factory
  inputs.json          manifest of build-layer inputs (logical name → size, SHA-256, figures)
  INPUTS.md            that manifest as a table, plus how to get the inputs
  fetch_data.py        download the build inputs from the data deposit
  requirements.txt  environment.yml  LICENSE (MIT)  .gitignore
  fig1/ … fig5/  figS1/ … figS7/     one package per figure:
    figlib.py            thin shim — declares this figure's build inputs, re-exports figlib_shared
    data/                source-data tables (CSVs + small arrays; also on KiltHub)
    render/              draw each panel FROM data/  → figures/   (fast, deterministic; what you run)
    build/               regenerate data/ FROM the deposited feature data (optional)
    figures/             rendered PNG/SVG output
    README.md            per-figure column dictionary + panel details
```

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
Figures assume the **Gillius ADF** font (matplotlib falls back silently if it is missing; the data and
coordinates are unaffected).

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

| Interactive | Package | Content |
|---|---|---|
| 1 — UMAP, quantitative features (µPULLI-I) | `interactive/` | click a replicate → its peak biofilm biomass image |
| 2 — UMAP, DINOv2 CLS embeddings PCA-50 (µPULLI-DL) | `interactive/` | same explorer, embedding manifold |
| 3 — dendrogram + heatmap, all timepoints | `interactive/` | scrub or play through 0–30 h + peak; click a mutant → its replicate images, which follow the slider |
| 4 — RNA-seq volcano, Δ*bioD* vs WT | `interactive/` | hover a gene → identifiers, annotation, log₂FC, adjusted *p*; pathway highlights; gene/locus/description search |
| 5 — RNA-seq volcano, Δ*pdhE2* vs WT | `interactive/` | same explorer |
| 6 — RNA-seq volcano, Δ*manA* vs WT | `interactive/` | same explorer |

The interactive supplements are self-contained single-file HTML — for Plots 1–3 every thumbnail is
inlined as a base64 data URI, so each opens offline by double-clicking with no server and no asset
folder. They are built rather than bundled (see [`interactive/README.md`](interactive/README.md));
Plots 4–6 need no build input at all, since their tables ship in `interactive/data/`.

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
| Features | full feature sets + gene index | CMU KiltHub `[DOI]` | [µPULLI-I](https://github.com/BridgesLabCMU/uPULLI-I) |
| Embeddings | DINOv2 CLS + patch descriptors | CMU KiltHub `[DOI]` | [µPULLI-DL](https://github.com/BridgesLabCMU/uPULLI-DL) |
| Tables | figure source-data CSVs (`fig*/data/`) | KiltHub `[DOI]` + this repo | `fig*/build/` |
| Datasets | **Data S1** — full transposon-screen results (every feature + phenotype call, per well and per frame) | KiltHub `[DOI]` | `figS4/build/DataS1_transposonScreen.py` |
| Figures | panels (`fig*/figures/`) | this repo | `fig*/render/` |

Steps above the "Tables" row are the µPULLI pipeline (a separate repository); this repo covers
Tables → Figures (and the feature-set → table reshaping in each `build/`).

## Data & Code Availability (fill bracketed IDs at submission)

Derived/tabular data (feature vectors, UMAP coordinates, per-microcolony measurements, biomass
trajectories) → CMU KiltHub `[DOI]`. Processed image data + segmentation masks → BioImage Archive
[`S-BIAD3830`](https://doi.org/10.6019/S-BIAD3830); raw timelapses in the same BioImage Archive deposit. Raw RNA sequencing reads → NCBI SRA, BioProject `PRJNA1513060`. Code (this repository) →
Bridges Lab GitHub, MIT license.
