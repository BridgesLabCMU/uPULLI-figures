# Interactive Plots 1–6: the reimaging atlas and the clean-deletion transcriptomes

Six self-contained single-file HTML supplements. Every thumbnail is inlined as a base64 data URI,
so each file works **offline by double-clicking** — no server, no sibling asset folder, no network.

| Plot | `build/` script | What it shows |
|---|---|---|
| **1** | `IP12_interactiveUmap.py --track numerical` | UMAP of the atlas from the quantitative µPULLI-I features; click a replicate → its peak biofilm biomass image |
| **2** | `IP12_interactiveUmap.py --track embedding` | the same explorer over the PCA-50 DINOv2 CLS embedding manifold (µPULLI-DL) |
| **3** | `IP3_interactiveDendrogramTimelapse.py` | dendrogram + heatmap of the atlas across the whole timecourse (slider + play, adjustable fps); click a mutant → its replicate images, which **follow the time slider** |
| **4** | `IP456_interactiveVolcano.py --mutant BioD` | RNA-seq volcano for Δ*bioD* vs WT; hover a gene → name, locus, description, log₂FC, adjusted *p*; sidebar toggles pathway highlights; search matches gene/locus/description |
| **5** | `IP456_interactiveVolcano.py --mutant PdhE2` | the same explorer for Δ*pdhE2* vs WT |
| **6** | `IP456_interactiveVolcano.py --mutant ManA` | the same explorer for Δ*manA* vs WT |

```
figlib.py   thin shim (paths + shared lib); declares the one build input
data/       coordinates, well labels and dendrogram tables (all small, all bundled)
build/      the three generators
archive/    retired builders, kept for reference (not part of the package's output)
figures/    the rendered .html (gitignored — rebuild, or take from the deposit)
```

> **Retired:** the static Plot 3 — same viewer, but every replicate shown at its own peak-biomass frame
> and no image animation — is no longer built. Its dual-mode generator is at
> `archive/IP3_interactiveDendrogram.py` (no flag reproduces that old build, `--timelapse` the current
> one, though it still writes the superseded `interactivePlot3b_` filename). The timelapse viewer is
> **Interactive Plot 3** from here on. Note the archived copy predates
> the axis-label and scroll fixes below.

## Build

```bash
python fetch_data.py interactive          # from the repo root -> thumbnail pack + wide table
python build/0_buildFrameMatrices.py                      # ~1 s  -> data/fullAtlas_frameMatrices.npz
python build/IP12_interactiveUmap.py --track numerical    # ~5 s  -> 27 MB
python build/IP12_interactiveUmap.py --track embedding    # ~5 s  -> 27 MB
python build/IP3_interactiveDendrogramTimelapse.py         # ~9 s  -> 30 MB
python build/IP456_interactiveVolcano.py                  # ~5 s  -> 3 x 0.8 MB
```

Plots 4–6 need **no build input at all** — their tables are bundled in `data/`, so they rebuild on a
fresh clone with nothing downloaded.

Step 0 is only needed for Plot 3 and only if `data/fullAtlas_frameMatrices.npz` is absent — it is
bundled, so a normal rebuild starts at `IP3`.

### Plot 3 playback

The heatmap carries **32 timepoints**: hourly frames 0–30 h plus the peak-biomass frame as the final
slider position, which is where the page opens so it matches Fig. S6A on load. A time slider scrubs,
**Play** animates, and a second slider sets the rate from 1–24 fps (default 8, changeable mid-play).
Arrow keys step; space toggles play. The dendrogram, leaf order and replicate galleries are fixed —
only the heatmap cells change with time.

Frames follow the published animation convention: **per-mutant median** of each feature at that
frame, **z-scored across mutants within that frame**, clipped to ±3. Each frame is standardized
independently, so colors compare mutants *within* a timepoint but not *across* timepoints.

### Why three replicates per mutant

Making the images move with the slider cannot be had for free — they are inlined, and **every well at
every frame is ~330 MB of JPEG (~440 MB base64)**, far past what a single file can carry. Measured at
112 px:

| replicates/mutant | × 31 frames | inlined HTML |
|---|---|---|
| all (~23) | 334 MB | ~445 MB — not viable |
| 6 | 83 MB | ~110 MB |
| 4 | 56 MB | ~75 MB |
| **3 @ 96 px** | **~29 MB** | **~40 MB** ← the shipped build |

So the viewer keeps every timepoint for three representative replicates rather than every replicate at
one timepoint (which is what the retired Plot 3 did). The representatives are the wells whose peak
biomass is closest to the mutant's median peak biomass (the "median" selection the representativeImages
pullers use), so they are typical rather than extreme. Regenerate the pack at a different trade-off with
`v2/reimaging/umap/exportTimelapsePack.py --reps N --px P --quality Q`.

The caption under each image reports the frame being shown and marks that well's own peak
(`t=17h (peak)`); at the slider's final "Peak biomass" position each well shows its own peak frame,
which is why the images there differ from one another in time.

Only the thumbnail pack is a build input; everything else is in `data/`. The builders need **no
image-processing dependency** — the pack is already encoded at final size, so they only base64 it,
and the pinned core environment is enough (`IP3` additionally uses scipy + matplotlib, both pinned).


### Plots 4–6: the volcano explorers

Ported from the lab's `Interactive_Volcano_Explorer_light_V2.py` (RNAseq/Scripts on the lab share),
with the package's rules applied: no absolute paths, tables bundled in `data/`, nothing to download.

Every quantified gene is a point (3,572 CDS per contrast); the non-significant cloud is kept, since a
volcano needs it. Dashed guides mark **|log₂FC| = 2** and adjusted *p* = 0.05 — the same cutoffs as
Fig S8B–D, so the interactive and printed volcanoes call exactly the same genes significant. (The
lab's original explorer used 1, which disagreed with the panels.) Points past **both** are drawn in
**that mutant's own color** — orange for Δ*bioD*, green for Δ*pdhE2*, blue for Δ*manA*, matching
Figs 3–5 and S8B–D — with a thin dark edge, as the printed panels use. Search hits are pink rather
than the lab's green, which would have collided with the Δ*pdhE2* page. The sidebar toggles pathway highlights; the search box matches gene name, either locus form
(`VC_RS00005` or `VC_0002`) or description text, and highlights hits in green with a live count.

The six locus-keyed functional groups carry **the same loci and the same colors** they have in Fig 3
and Fig 4, so a pathway reads identically across the paper. Five further sets, keyed by gene name, are
specific to this transcriptional view: *vps I*, *vps II*, Matrix Proteins, Type IV Pilus, Adhesins.
Matching normalizes ids (strips `_`/spaces, uppercases) and tries the locus, the gene name and every
`;`-separated old locus, so all eleven sets resolve fully against the tables.

Counts past both guides: **Δ*bioD* 189**, **Δ*pdhE2* 658**, **Δ*manA* 1** — Δ*manA*'s single hit is
*manA* itself, which is the point of that panel.

**All three share one set of axes** (log₂FC −11…10, −log₁₀ *p* 0…11, computed over all three
contrasts), so the panels are directly comparable and "reset" lands in the same place in each.

Highlighted genes are drawn in a pass *after* the main cloud and search hits after those, so a
highlighted gene is never hidden under a neighbour; hover resolves ties the same way, so the point you
see on top is the one the card describes. **Clicking a gene pins its card** (click empty space to
dismiss); the card follows its gene as you pan and is clamped inside the plot so edge genes are not
cut off. Each sidebar row reports `genes · significant` for that pathway, and **all** / **none**
toggle every set at once.

Navigation matches Plots 1–3: scroll pans, ctrl/⌘ + scroll zooms at the cursor, drag pans, shift-drag
box-zooms, double-click resets.

## Bundled tables

- **`wells.csv`** — one row per reimaging well (3,669): `plateId, wellId, peakFrame, peakBiomass,
  file` (the thumbnail's filename in the pack), `mutant, geneLocus, function`. Shared by all three
  plots; it is what links a point or a dendrogram leaf to an image.
- **`umapCoords_numerical.parquet`, `umapCoords_embedding.parquet`** — `plateId, wellId, mutant,
  geneLocus, function, n_neighbors, min_dist, umap1, umap2`, over the full 3 × 3 (n_neighbors,
  min_dist) grid, which is what the in-page selector switches between. Plots 1 and 2 differ *only*
  in which of these two they read.
- **`fullAtlas_pcaLinkage_{linkage.npy, cluster_order.csv}`, `fullAtlas_peakBiomass_featureMatrix.csv`**
  — the Fig. S6A clustering, reused so the interactive and printed dendrograms are the same object.
  ⚠ The feature matrix is keyed by **`mutant`** (gene name where one exists, else the locus), not by
  `gene`/locus — those agree for only 68 of 158 columns.
- **`rnaseq_volcano_<mutant>.csv`** (×3, plots 4–6) — `locus`, `oldLocus`, `gene`, `description`,
  `logFC` (log₂ vs WT), `qvalue` (Benjamini–Hochberg adjusted *p*). All 3,572 quantified CDS per
  contrast, unfiltered. Copies of the Fig S8B–D tables, bundled so this package stands alone.
- **`fullAtlas_frameMatrices.npz`** (208 kB) — `matrices` (31 × 158 × 27 float16: frame × mutant in
  leaf order × feature), plus `frames`, `mutants`, `features`. Written by `build/0_buildFrameMatrices.py`;
  drives Plot 3's animation. NaN marks a cell with no data — see below.

### "No data" cells in Plot 3

The 11 colony-segmentation features do not exist before **t5** (colonies are not segmentable that
early), so ~16% of the animated cells have no value. Those are drawn in a distinct dark grey
(palette index 255) rather than filled with 0. This is a deliberate departure from the legacy
animation, which substituted 0 — a z of 0 renders mid-scale and reads as *average*, which is
precisely wrong for a missing measurement. A handful of additional cells at later frames are
genuinely NaN in the source and are shown the same way.

Cells ship as **base64 palette indices**, one byte each, against a 256-entry RdBu_r ramp: all 32
timepoints cost ~180 kB instead of ~1.4 MB of hex strings, and the colors still match the static
figure exactly.

### Navigating the plot

**Scroll pans** the view (two-finger scroll pans in both axes); **ctrl/cmd + scroll zooms** about the
cursor, as does a trackpad pinch. Dragging still pans, double-click or **Fit** resets. Feature names are
printed on **both** x axes — above the heatmap rising to the right, below it falling to the left — so a
column can be read without scrolling back to the top. The header and footer are sized at runtime from
the measured label widths (`ctx.measureText`), so no name is ever clipped or overrun by the cells.

## Build inputs (deposited)

`reimaging/thumbnails112/` (Plots 1–2) — one 112 px grayscale JPEG per well at its peak-biomass frame,
named `<plate-well>_t<peakFrame>.jpg`, plus `manifest.csv`. ~12 MB for 3,867 wells. Produced upstream by
`v2/reimaging/umap/exportThumbnailPack.py` from the reprocessed `processed.tif` stack using the
standard 0.5–99.5 percentile stretch, i.e. the same recipe as every other still in the paper.

`reimaging/thumbnailsTimelapse/` (Plot 3) — every frame of three representative wells per mutant at
96 px, same naming and same stretch, ~65 MB; from `v2/reimaging/umap/exportTimelapsePack.py`.

Pre-encoding at final size is deliberate: the alternative is shipping the 200 GB image tree or the
200 MB 256 px cache, and re-deriving the stretch at build time.

## Suggested captions

> **Interactive Plot 1.** Interactive UMAP of the transposon reimaging landscape based on
> quantitative features extracted using µPULLI-I. Clicking any replicate displays its peak biofilm
> biomass image.
>
> **Interactive Plot 2.** Interactive UMAP of the transposon reimaging landscape based on PCA-50
> dimensionality reduction of DINOv2 CLS embeddings extracted using µPULLI-DL. Clicking any
> replicate displays its peak biofilm biomass image.
>
> **Interactive Plot 3.** Interactive dendrogram and heatmap of the transposon reimaging atlas
> across the imaging timecourse, linking every mutant to its genotype and to the biofilm images its
> features were measured from. The heatmap and the images can be scrubbed or played through all
> hourly timepoints at an adjustable frame rate.

> **Interactive Plots 4–6.** Interactive RNA-seq volcano plots of Δ*bioD* (4), Δ*pdhE2* (5) and
> Δ*manA* (6) relative to wild type. Every quantified gene is shown; hovering reports its identifiers,
> annotation, log₂ fold change and adjusted *p*-value, the sidebar highlights curated pathways, and the
> search box locates genes by name, locus or description.

The Plot 3 caption is deliberately not the "…and timelapse video" wording: the viewer carries stills, not
video, and its leaves are the 158 **mutants** (clicking one opens its representative replicates), so
the images are reached through the genotype rather than being individually clickable.

## Notes

- **Sizes:** 27 MB (plots 1, 2), 30 MB (plot 3 — 14,694 thumbnails, every frame of three
  representative wells per mutant), 0.8 MB each (plots 4–6, which carry no images). Each thumbnail is
  stored once and referenced by index.
- **Where they go:** these are journal supplementary files and/or the GitHub Pages viewer. They are
  gitignored here like `fig*/figures/` — regenerate them, or take them from the deposit.
- **Verification:** after building, the JavaScript can be checked with
  `node --check` on the extracted `<script>` block, and every thumbnail index should be < the length
  of `THUMBS`. Both hold for the current build.
