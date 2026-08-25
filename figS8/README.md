# Figure S8: transcriptional and reporter phenotypes of the clean deletions

**B–D** RNA-seq volcano plots for the three in-frame clean deletions vs WT, and **E** peak
P*vpsL*-lux reporter activity in the same strains. Both address the same question from different
angles: what each deletion does to gene expression, and what it does to *vps* transcription
specifically.

## Layout

```
figlib.py            shim: figS8 paths + clean-deletion colors/labels + lux layout + re-export of ../figlib_shared.py
data/                source-data tables (CSVs); also on KiltHub
render/              draw each panel from data/  -> figures/   (what you run)
build/               regenerate the tables from the deposited RNA-seq / plate-reader inputs
figures/             rendered PNG/SVG output
```
(Shared `requirements.txt` / `environment.yml` / `LICENSE` live at the `paper-figures/` root.)

## Panels -> scripts -> tables

| Panel | `render/` script | `data/` table | `build/` script |
|---|---|---|---|
| **S8B** Δ*bioD* volcano | `S8BCD_volcano.py --mutant BioD` | `rnaseq_volcano_BioD.csv` | `S8BCD_volcano.py` |
| **S8C** Δ*pdhE2* volcano | `S8BCD_volcano.py --mutant PdhE2` | `rnaseq_volcano_PdhE2.csv` | `S8BCD_volcano.py` |
| **S8D** Δ*manA* volcano | `S8BCD_volcano.py --mutant ManA` | `rnaseq_volcano_ManA.csv` | `S8BCD_volcano.py` |
| **S8E** peak P*vpsL*-lux activity | `S8E_luxActivity.py` | `luxPeak_normWT.csv` | `S8E_luxActivity.py` |

Running `S8BCD_volcano.py` with no `--mutant` renders all three panels; the panel letter is attached
per mutant in `figlib.RNASEQ_PANELS`, so outputs are `S8B_volcano_BioD`, `S8C_volcano_PdhE2`,
`S8D_volcano_ManA`.

> Also in this package: `render/S8_trainingFeatureDendrogram.py`, the feature-clustered training
> heatmap. It was **promoted to Figure 2C** (`fig2/render/2C_trainingHeatmap.py`) and is kept here only
> for reference — it is not an S8 panel. ⚠ Its layout has a known bug (`yaxis.tick_right()` collides
> with the family strip and colorbar); the fixed version is the fig2 copy.

## Panel details

**S8B–D (volcanoes).** One square panel per mutant, drawn from the full quantified-gene table (all
3,572 CDS — a volcano needs the non-significant cloud, so `build/` filters nothing). Points past
**both** thresholds (|log2FC| > 2 **and** q < 0.05) take that mutant's own reimaging functional-group
color, the same one it carries in Fig 4 and Fig 5 (Δ*bioD* orange, Δ*manA* blue, Δ*pdhE2* green);
everything else is grey. Dashed guides mark both cutoffs — `--fc` / `--q` to change them. Up to
`--labels` gene names are drawn per side, de-duplicated by gene name (names recur across loci, e.g.
two *groL*) and skipped where they would overlap a label already placed. The mutant label sits
bottom-left, which is empty in every volcano by construction — the top-left corner is exactly where
the deleted gene lands.

Colored counts at the defaults: **Δ*bioD* 189**, **Δ*pdhE2* 658**, **Δ*manA* 1**. Δ*manA*'s only
significant gene is *manA* itself — that deletion produces essentially no transcriptional response, so
S8D is deliberately near-empty. In Δ*bioD*, *bioD* is at −10.5 while the rest of the operon is strongly
up (*bioA* +3.9, *bioB* +3.8, *bioF* +3.4, *bioC* +2.3) — derepression of the biotin pathway.

Contrast letters (B/E/M vs W) are the lab's file naming; `build/` re-checks each by asserting the
deleted gene is the strongest depletion (bioD −10.52; manA −6.23; for Δ*pdhE2* it prints **aceF**
−7.26, which is the same gene — pyruvate dehydrogenase E2 / dihydrolipoamide acetyltransferase).

**S8E (lux).** Peak P*vpsL*-lux activity per replicate well, normalized to WT. Ported from the lab's
`Continuous_Peak_Plotting.R` so the panel matches the analysis of record: RLU = Lum/OD per well per
timepoint (inf / 0 / negative → NaN), peak = max over the 41 hourly timepoints, then divided by the
**mean** of the 9 WT wells' peaks. Note the reference is a mean while the plotted center line is a
**median** (`--center mean` to switch). One dot per well with deterministic jitter (seed 42), a wide
bar at the median, dotted line at WT = 1.

Row B of the plate is Δ*pdhR* and is **excluded** — it is commented out of the lab analysis. Δ*pdhE1*
is present in `data/` but off by default (`--conditions WT,BioD,ManA,PdhE2,PdhE1` to include it), since
it has no counterpart in the imaging clean-deletion set.

## Column dictionaries

- **`rnaseq_volcano_<mutant>.csv`** (×3, one per clean deletion vs WT) — `locus` (VC_RS-style),
  `oldLocus` (the `VC_0002`-style ID that matches the reimaging atlas), `gene`, `description`, `logFC`
  (log2 vs WT), `qvalue` (Benjamini–Hochberg adjusted p). All 3,572 quantified CDS per contrast,
  unfiltered, no missing q-values. Only 1,404 rows carry a gene name.
- **`luxPeak_normWT.csv`** — `condition` (WT / PdhE1 / BioD / ManA / PdhE2), `well`, `peakRLU` (max
  Lum/OD over 41 timepoints), `peakNormWT` (peak ÷ mean WT peak). 45 wells, 9 per condition; WT
  reference 735.2 RLU. Medians: WT 1.00, Δ*bioD* 7.19, Δ*pdhE2* 3.98, Δ*manA* 1.66 (Δ*pdhE1* 4.90).
- **`trainingHeatmap_peakBiomass_featureMatrix.csv`** — the Fig 2C matrix, used only by the retired
  training-dendrogram script above.

## Build inputs (deposited)

- **`rnaseq/`** — the per-contrast `<strain>_allGenes.csv` differential-expression tables (shared with
  Fig 4E).
- **`lux/Lum.csv`**, **`lux/OD.csv`** — matched luminescence and optical-density plate-reader
  timecourses (41 hourly timepoints × wells A1–H12) for the P*vpsL*-lux reporter.

```bash
python stage_inputs.py figS8            # maintainer: symlink the inputs
python figS8/build/S8BCD_volcano.py     # -> data/rnaseq_volcano_*.csv
python figS8/build/S8E_luxActivity.py   # -> data/luxPeak_normWT.csv
python figS8/render/S8BCD_volcano.py    # -> figures/S8{B,C,D}_volcano_*.{png,svg}
python figS8/render/S8E_luxActivity.py  # -> figures/S8E_luxActivity.{png,svg}
```
