# uPULLI-figures: reproduction code + source data

Self-contained code and source-data tables that regenerate the figures for the µPULLI manuscript
(*V. cholerae* biofilm phenotyping). One package per figure (`fig1/`…`fig5/`); a shared library
(`figlib_shared.py`) holds the functional-group colors, the UMAP feature filter, and the plot style so
they are defined once.

## Layout

```
uPULLI-figures/
  figlib_shared.py     shared: colors / feature filter / plot style / per-figure path factory
  requirements.txt  environment.yml  LICENSE (MIT)  .gitignore
  fig3/                one package per figure:
    figlib.py            thin shim — per-figure paths + re-export of figlib_shared
    data/                source-data tables (CSVs + small arrays; also deposited on KiltHub)
    render/              draw each panel FROM data/  → figures/   (fast, deterministic; what you run)
    build/               regenerate data/ FROM the raw feature matrix (KiltHub; heavy, optional)
    figures/             rendered PNG/SVG/PDF output
    README.md            per-figure column dictionary + panel details
```

## Reproduce a figure

```bash
conda env create -f environment.yml && conda activate upulli-figures   # or: pip install -r requirements.txt
python fig3/render/3B_reimagingUmap_functionalAnnotations.py          # a Fig-3 panel; see each fig*/README.md
```

Render scripts read only the bundled `fig*/data/` tables — no heavy data, no network, deterministic
output. To regenerate the tables themselves from the raw feature matrix, see `fig*/build/README.md`.
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

Each `fig*/README.md` carries that figure's full lineage — the panel→script→table map, the column
dictionary for every table, and the processing chain from the original feature set down to the panels.

## Data / processing chain (repo-wide)

| Layer | Artifact | Deposited at | Produced by |
|---|---|---|---|
| Raw | brightfield timelapses | institutional NAS (on request) | — |
| Images | processed images + segmentation masks | BioImage Archive `[S-BIAD#####]` | µPULLI image pipeline `[pipeline repo]` |
| Features | full feature sets + DINOv2 embeddings + gene index | CMU KiltHub `[DOI]` | µPULLI feature extraction `[pipeline repo]` |
| Tables | figure source-data CSVs (`fig*/data/`) | KiltHub `[DOI]` + this repo | `fig*/build/` |
| Figures | panels (`fig*/figures/`) | this repo | `fig*/render/` |

Steps above the "Tables" row are the µPULLI pipeline (a separate repository); this repo covers
Tables → Figures (and the feature-set → table reshaping in each `build/`).

## Data & Code Availability (fill bracketed IDs at submission)

Derived/tabular data (feature vectors, UMAP coordinates, per-microcolony measurements, biomass
trajectories) → CMU KiltHub `[DOI]`. Processed image data + segmentation masks → BioImage Archive
`[S-BIAD#####]`; raw timelapses on institutional storage, available on request. Code (this repository) →
Bridges Lab GitHub, MIT license.
