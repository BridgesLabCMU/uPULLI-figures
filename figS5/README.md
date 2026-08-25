# Supplemental Figure S5: Low-Biofilm Transposon Mutants, by Chromosome

Per-gene biomass trajectories of the **Low-Biofilm** class from the genome-wide transposon screen
(Fig S4), split by chromosome. **(A)** Chromosome I, **(B)** Chromosome II.

## Layout

```
figlib.py   thin shim (paths + shared lib)
data/       tn_biomass_matrix.csv + tn_locus_meta.csv + tn_geneIndex.csv (bundled)
build/      S5_screenTrajectories.py — regenerates those tables from the deposited screen data
render/     S5_lowBiofilmChromosome.py -> figures/
figures/    S5A_ChromosomeI.{png,svg}, S5B_ChromosomeII.{png,svg}
```

## Notes

- **Membership** = the Low-Biofilm class (peak biomass below the WT-anchored low threshold) from the
  screen classification in `tn_locus_meta.csv` (53 mutants: 48 on Chr I, 5 on Chr II).
- Low-biofilm mutants were **not reimaged**, so they have no gene names or functional annotations —
  rows are labeled by locus number, with no phenotype/function markers.
- **Color scale 0–3** (matching the other figures). Because low-biofilm mutants sit below the low-biofilm
  threshold (see `B_min` in `figS4/data/tn_thresholds.csv`), they appear near the bottom (dark) end of
  the scale.
- One cell per timepoint (8–30 h); the locus list is wrapped across side-by-side columns.

The `data/` tables are the transposon-screen source tables shared with Fig S4 (regenerate them with
`build/S5_screenTrajectories.py`, the same parse figS4 runs for its own panels — duplicated so
neither package reaches into the other's `data/`).
