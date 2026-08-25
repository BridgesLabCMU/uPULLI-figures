# Supplemental Figure S6: Full Reimaging Atlas — Quantitative Features

Two things about the atlas built from the hand-engineered **quantitative** features: how all 158
mutants relate to one another (**A**), and that those measurements are reproducible (**B, C**).

B and C form a 2×2 — **rows are the claim, columns are the granularity**:

| | pooled (left) | per mutant (right) |
|---|---|---|
| **B** — replicates are closer to each other than to *other mutants* | all 6.7M pairwise distances | paired, one line per mutant |
| **C** — replicates are closer than *chance* | 158 per-mutant values vs the permutation null | each mutant vs its own null, FDR-corrected |

Read left→right, top→bottom, each panel needs only concepts the previous ones introduced: B-left
defines "distance" and the within/between split; B-right introduces the per-mutant unit; C-left
distributes exactly the quantity B-right's "Same Mutant" column plots; C-right resolves it per mutant
and adds significance.

## Layout

```
figlib.py   thin shim (paths + shared lib)
data/       fullAtlas_* (panel A) + replicate_* (panels B, C) source-data tables
render/     one script per panel -> figures/   (what you run; reads only data/)
build/      S6BC_replicateConsistency.py — regenerates the replicate_* tables (optional)
figures/    rendered PNG + SVG
```

| Panel | `render/` script | `data/` table(s) | `build/` |
|---|---|---|---|
| **A** dendrogram + heatmap, 158 mutants | `S6A_fullDendrogramHeatmap_vertical.py` (+ `_horizontal.py`) | `fullAtlas_pcaLinkage_{linkage.npy, cluster_order.csv}`, `fullAtlas_peakBiomass_featureMatrix.csv` | upstream (v2 full-atlas dendrogram) |
| **B left** within- vs between-mutant pairwise distances | `S6B_left_distanceDistributions.py` | `replicate_distanceHistogram.csv` | `S6BC_replicateConsistency.py` |
| **B right** paired within vs between, per mutant | `S6B_right_pairedWithinVsBetween.py` | `replicate_perMutant.csv`, `replicate_summary.csv` | ″ |
| **C left** observed vs null per-mutant distance | `S6C_left_perMutantNullOverlay.py` | `replicate_perMutantNullHistogram.csv` | ″ |
| **C right** each mutant vs its own null | `S6C_right_perMutantVsNull.py` | `replicate_perMutant.csv` | ″ |

- **A, horizontal vs vertical:** horizontal = dendrogram on top, features as rows (grouped: biofilm
  biomass, colony-segmentation, whole-image Haralick — brackets), mutants as columns in leaf order,
  functional-annotation color strip under the labels. Vertical = the transpose (mutants on y with the
  dendrogram and strip on the left, large italic labels on the right, feature labels rotated 30°).
  RdBu_r, ±3 z-score; heatmap/strip rasterized, SVG at dpi 200 (raster-size safe at this extent).
- **B (right)** takes `--sig {value,stars,none}` (default `value`) for the significance bracket.

## Rebuilding the replicate_* tables

```bash
python fetch_data.py figS6                    # from the repo root -> reimaging/collapsedWide.parquet
python build/S6BC_replicateConsistency.py   # ~4 min at the default 10,000 permutations
```

The build is self-contained: it reproduces the manifold's own preprocessing (drop unlabelled rows →
drop excluded loci → frames 9–27 feature filter → fillna(0) → drop zero-variance → growth filter →
≥5 replicates → StandardScaler, giving 3669 wells × 285 features) and derives every statistic from
one pairwise-distance matrix. Panel A's `fullAtlas_*` tables come from the upstream v2 full-atlas
dendrogram (`results/v2/reimaging/dendogram/fullReimaging_pcaLinkage/`) and are provided as source
data.

**Do not lower `--perms` below ~3200 if you care about the Bonferroni column.** A permutation *p*
cannot go below 1/(nPerms+1); at 1000 permutations that floor (1.0 × 10⁻³) sits *above* the
Bonferroni threshold 0.05/158 = 3.2 × 10⁻⁴, so zero mutants can pass regardless of the true effect.
`replicate_summary.csv` carries `permP_floor` and `bonferroniThreshold` so the two can be compared.

## Column dictionaries

- **`replicate_distanceHistogram.csv`** — shared-bin histogram of all pairwise distances.
  `binLeft, binRight, binCenter`; `withinCount / betweenCount` (raw counts, same-mutant pairs vs
  different-mutant pairs); `withinDensity / betweenDensity` (each normalized to integrate to 1, which
  is what the panel draws — the two classes differ 160-fold in count).
- **`replicate_perMutantNullHistogram.csv`** — shared-bin histogram of *one mutant's* mean
  within-replicate distance. `binLeft, binRight, binCenter`; `observedCount` (the 158 real mutants),
  `nullGlobalCount` / `nullWithinPlateCount` (158 × nPerms values pooled over permutations, labels
  shuffled freely / shuffled only within a plate); matching `*Density` columns.
- **`replicate_perMutant.csv`** — one row per mutant, sorted by `meanWithin`. `mutant, geneLocus,
  nReplicates`; `meanWithin` (mean pairwise distance among its own replicates), `meanToOtherMutants`
  (mean distance from its wells to every other mutant's wells), `withinMinusOthers` (their
  difference; negative = tighter within); `nullMeanWithin, nullSdWithin` (that mutant's own
  permutation null — the SD depends on replicate count); `p` (per-mutant permutation p, floored at
  1/(nPerms+1)), **`qBH`** (Benjamini–Hochberg FDR across the 158 tests — the one to report),
  `pBonferroni`.
- **`replicate_summary.csv`** — `statistic, value`, one row each; every number quoted in the figure
  legend below is here, so captions can be regenerated rather than retyped.

## Figure legend — the panels are spare, the numbers live here

B and C carry no stats blocks and no axis parentheticals **by design**; all of it belongs in the
manuscript legend. Shared preamble:

> Distances are Euclidean in the 285-dimensional standardized feature space used to build the
> reimaging manifold (19 timepoints, frames 9–27, of biomass + whole-image entropy/Haralick features,
> z-scored per feature) and are therefore **dimensionless** — not "a.u." and not standard deviations.
> *n* = 3669 wells, 158 mutants (median 24 replicates each), 48 plates, after the atlas filters.
> All *p* values are from 10,000 permutations of the mutant labels and are floored at 1/10,001;
> report them as *p* < 0.0001.

- **(B, left)** *Pairwise distances between wells of the same mutant (blue, n = 41,110 pairs) and of
  different mutants (grey, n = 6,687,836). Replicate pairs are closer: mean 14.0 vs 22.1
  (Cohen's d = 1.13). A random replicate pair is the closer one 81% of the time (AUC = 0.806, the
  normalized Mann–Whitney U; Cliff's δ = 0.61; KS D = 0.459), exceeding all 10,000 permutations
  (null AUC 0.500 ± 0.001, z = 285, p < 0.0001). Significance is by permutation, not from the
  analytic Mann–Whitney or KS p-value: each well appears in 3668 pairs, so the distances are not
  independent.*
- **(B, right)** *For each mutant, the mean distance from its wells to its own replicates vs to all other
  mutants' wells; one line per mutant. All 158/158 are tighter within than between (Wilcoxon
  signed-rank W = 0, n = 158, p < 10⁻²⁰; rank-biserial r = −1.0; median difference −7.65,
  95% CI [−8.0, −6.9] by 10,000 bootstrap resamples over mutants; smallest single-mutant gap 1.08).
  The mutant is the unit of analysis, so unlike pairwise distances these observations are
  approximately independent and a paired test is valid.*
- **(C, left)** *Distribution of one mutant's mean within-replicate distance — observed across the 158
  mutants (blue; median 13.6, range 9.7–24.8) and under the null pooled over 10,000 permutations
  (grey, labels shuffled; dark grey, shuffled within plate). 157 of 158 observed values fall below
  the null mean of 22.09; a random observed value is smaller than a random null value 98.7% of the
  time (AUC = 0.987; z = 68.3, p < 0.0001). The within-plate null holds plate composition fixed and
  so controls for batch.*
- **(C, right)** *Each mutant's mean within-replicate distance, sorted, colored by functional group. Dashed
  line = the null expectation (22.09); grey band = ± 2 SD of a single mutant's null, drawn at the
  across-mutant mean SD (2.13) as a visual guide while each mutant is tested against its own null
  (SD 2.08 at 24 replicates to 4.37 at 6). 157 of 158 mutants fall below chance and **150 of 158 are
  significant at Benjamini–Hochberg q < 0.05** (129 under Bonferroni; 151 uncorrected). Mutants near
  the chance line are phenotypically WT-like, not poorly measured.*

Two statements worth adding to the text, both from `replicate_summary.csv`: the design makes a batch
explanation impossible — each plate carries ~76 mutants at one well each and each mutant recurs on
~24 plates, so **100% of within-mutant pairs are cross-plate pairs** — and consistently, cross-mutant
pairs on the same plate average 21.5 vs 22.2 on different plates, i.e. plate identity moves the
distance by 0.6 where mutant identity moves it by 8.1.

## Cautions

- **The reported Wilcoxon p (5.6 × 10⁻²⁸) is correct but over-precise.** It is the normal
  approximation with W = 0, n = 158 (z = −10.90), and is conservative against the exact signed-rank
  floor of 2⁻¹⁵⁸. But it assumes the 158 paired differences are independent, and every mutant's
  "distance to other mutants" is computed against nearly the same pool of wells. Panel B (right) therefore
  prints an inequality; lead with 158/158 and the median gap.
- **Only panel C (right) needs multiplicity control** — it is 158 tests. B (left), B (right) and C (left) are each one
  pre-specified test; no correction is applied across panels, which are facets of one hypothesis
  computed from a single distance matrix.
- **The feature set is 87% Haralick texture** (247 of 285 columns; biomass and entropy 19 each), so
  "distance" here is dominated by texture. The direction of every result is robust; magnitudes would
  shift under per-family weighting.
- Distances are computed in the standardized feature space, **not** in UMAP coordinates — UMAP
  optimizes local neighborhood preservation, so measuring replicate tightness there would be close
  to circular.
