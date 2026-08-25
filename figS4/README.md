# Supplemental Figure S4: Genome-Wide Transposon Biomass Screen

The label-free biofilm-biomass screen of an ordered *V. cholerae* transposon library (one insertion per
non-essential gene, ~2,900 mutants). Each mutant is classified against the WT replicate distribution
(Normal / Low Biofilm / High Biofilm / Dispersal Defect); the High-Biofilm + Dispersal-Defect hits are
the strains taken forward to the reimaging atlas (Fig 3).

**The phenotype calls are a fixed result, not recomputed here.** They are the screen analysis that
selected the reimaging library, so they are read from `TransposonResults_10x.csv` as given. The biomass
*trajectories* rendered in panels A and S5 come from the screen's feature table — the same split panels
B/C use, so their color scale is independent of the Fig-3A ranked-scatter axis.

## Layout

```
figlib.py            thin shim: figS4 paths + phenotype colors + re-export of ../figlib_shared.py
data/                source-data tables (CSVs); also on KiltHub
render/              draw each panel from data/  -> figures/   (what you run)
build/               regenerate data/ from the screen master feature table (optional)
figures/             rendered PNG/SVG output
```
(Shared `requirements.txt` / `environment.yml` / `LICENSE` live at the `paper-figures/` root.)

## From raw data to panels

1. **Raw brightfield timelapses** -> *(µPULLI pipeline)* -> two inputs (both KiltHub `[DOI]`):
   `TransposonResults_10x.csv`, the screen analysis of record (per-well `Peak`/`Final`/`Early`/
   `Phenotype` for 2,849 wells / 2,847 loci), and the screen's `master_frame_features.csv`
   (per plate/well/frame, 3,070 wells across 32 plates).
2. `build/S4_screen.py` takes the calls and thresholds from the results table as given, builds
   trajectories from the feature table (joined to the bundled `data/tn_geneIndex.csv` for
   well -> geneLocus, normalized to the training-set WT peak mean), and writes the source-data
   tables in `data/`.
3. `data/` CSVs -> `render/S4*.py` -> panels in `figures/`.

Multi-locus insertions are written `VC_0255; VC_A0373` in the results table and `VC_0255/VC_A0373` in
the screen index; the build normalizes both to the slash form before joining.

### Panels -> scripts -> tables

| Panel | `render/` script | `data/` table(s) | `build/` script |
|---|---|---|---|
| **A** (`S4A`) genome-wide biomass heatmap, horizontal (genome on x, time on y; + reimaging strip) | `S4A_genomeHeatmap.py` | `tn_biomass_matrix.csv`, `tn_locus_meta.csv` | `S4_screen.py` |
| **B** (`S4B`) Chromosome I, reimaging trajectories (wrapped columns) | `S4BC_chromosomeHeatmaps.py` | `reimHits_biomass_matrix.csv`, `reimHits_locus_meta.csv` | `S4BC_reimagingHits.py` |
| **C** (`S4C`) Chromosome II, reimaging trajectories | `S4BC_chromosomeHeatmaps.py` | `reimHits_biomass_matrix.csv`, `reimHits_locus_meta.csv` | `S4BC_reimagingHits.py` |

> **Panel D was promoted to Figure 3A.** The ranked peak-biomass classification scatter is now
> `fig3/render/3A_rankedTransposonBiomass.py`, and it took its data with it: `tn_phenotype_summary.csv`
> and `tn_thresholds.csv` are produced by **fig3's own** build (`fig3/build/3A_rankedTransposonBiomass.py`)
> from the same `transposons/results_10x.csv`, and no longer live here. `build/S4_screen.py` still parses
> the calls of record — S4A and S5 color their trajectories by them — but now writes only what S4 draws,
> plus `tn_normalization.csv` for provenance.

Panels **B/C** show the **reimaging-selected mutants** — the 157 non-WT members of the Fig-3 reimaging
set (134 on Chr I, 23 on Chr II). Membership + phenotype come from the authoritative reimaging analysis
(`ReimagingResults_10x.csv`; binary: Dispersal Defect, else High Biofilm), and the **trajectories are the
reimaging feature data** (`data/v2/reimaging/reimaging_collapsedWide`, per-gene mean over
replicates, normalized to the reimaging WT peak — so B/C use a different color scale from A/D). One cell
per timepoint (8–30 h), one row per mutant, wrapped across side-by-side columns. Each row is labeled
(black, large) with its **gene name** (from `REIMAGING_INDEX`) or **locus number**. A marker at the right
encodes **phenotype by shape** (square = High Biofilm, circle = Dispersal Defect) and **functional
annotation by fill color** (the six `FUNCTION_COLORS` groups; open = unannotated). Panel **A** is the
full-genome transposon-screen (all ~2,900 loci) overview and comes from the transposon screen; **B/C**
come from the reimaging analysis. Fig 3A ranks the same screen A summarizes, on its own axis.

## Data S1 (supplementary dataset, not a panel)

`build/DataS1_transposonScreen.py` writes the paper's **Data S1** — every quantified feature for every
screened well, plus its phenotype call — as `DataS1_transposonScreen_perWell.csv` (2,939 wells x 203
columns, features at each well's peak-biomass frame) and `DataS1_transposonScreen_allFrames.csv.gz`
(91,097 well-frames). It writes to `../deposit/` rather than `data/`: at 8.4 MB + 27.4 MB these are
deposit artifacts, not per-panel source data to bundle with the code. Same inputs and same calls of
record as the panels above.

## Normalization + classification (manuscript Methods)

As applied in the screen analysis of record:

- `biomassNorm = biomass / mean(peak biomass of the WT control wells)`.
- Wells with `biomassNorm > 0.15` at t = 5 h are excluded as early-onset anomalies.
- `B_max = max(WT peak)/mean(WT peak)`; `B_min = min(WT peak)/mean(WT peak) - 0.25` (a conservative
  margin below the observed WT floor).
- **High Biofilm** if peak > B_max; **Low Biofilm** if peak < B_min; **Dispersal Defect** if late-phase
  (t = 25 h) biomass exceeds the WT maximum while peak is within [B_min, B_max]; **Normal** otherwise.

The WT reference set that produced those boundaries is not archived with the results table, so
`B_max`/`B_min`/`maxFinal` in `tn_thresholds.csv` are **recovered from the class boundaries** — each
threshold is bracketed by the gap between adjacent classes and reported as its midpoint
(`B_max` 1.3347, `B_min` 0.3015, `maxFinal` 0.1271). Figure 3A's dashed lines use those values, and its
build recovers them the same way from the same table — see `fig3/data/tn_thresholds.csv`.

The trajectories in panels A and S5 are separately normalized to the training-set WT peak mean
(`trajectoryNormalizer_wtPeakMean` in `tn_normalization.csv`), which is why their color scale is not on
the Fig-3A ranked-scatter axis.

## Column dictionaries

- **`tn_biomass_matrix.csv`** — geneLocus (rows, genome-ordered) × frame (cols); mean `biomassNorm`.
- **`tn_locus_meta.csv`** — `geneLocus`, `chromosome` (I/II), `peak`, `phenotype`, `isHit`
  (High Biofilm or Dispersal Defect); row-aligned to the matrix.
- **`tn_normalization.csv`** — `trajectoryNormalizer_wtPeakMean` (the training-set WT peak mean that
  A's and S5's trajectories are expressed in) plus the screen's per-class counts. The per-well call table
  and the recovered `Bmax`/`Bmin`/`maxFinal` moved to `fig3/data/` with Panel 3A.
- **`tn_geneIndex.csv`** — screen annotation: `plateID`, `plateLabel` (TN-Plate01…), `wellId`,
  `geneLocus`. One row per annotated well (2,939); the join key into the screen master.
- **`reimagingGeneNames.csv`** — `geneLocus` -> `geneName` for row labeling in panels A and B/C.

## Note on counts

The per-class counts come from the screen analysis of record and match the manuscript body: 2,849 wells
over 2,847 loci — 53 Low Biofilm, 120 High Biofilm, 46 Dispersal Defect, 2,630 Normal. Because these
calls selected the reimaging library, they are held fixed rather than recomputed; reclassifying on
re-measured features would move mutants across the thresholds and no longer describe the set that was
actually reimaged.
