# Interactive Plots 1–6: the reimaging atlas and the clean-deletion transcriptomes

Six self-contained single-file HTML supplements (plus two alternate builds of Plots 1 and 2). Every
thumbnail is inlined as a base64 data URI, so each file works **offline by double-clicking** — no
server, no sibling asset folder, no network.

| Plot | `build/` script | What it shows |
|---|---|---|
| **1** | `IP12_interactiveUmap.py --track numerical` | UMAP of the atlas from the quantitative µPULLI-I features; click a replicate → its peak biofilm biomass image |
| **2** | `IP12_interactiveUmap.py --track embedding` | the same explorer over the PCA-50 DINOv2 CLS embedding manifold (µPULLI-DL) |
| **1T** | `IP12_interactiveUmap.py --track numerical --timelapse` | Plot 1, but clicking a **representative** replicate plays its whole timecourse (slider, play/pause, adjustable fps); every other replicate keeps its peak still |
| **2T** | `IP12_interactiveUmap.py --track embedding --timelapse` | the same, on the embedding manifold |
| **3** | `IP3_interactiveDendrogramTimelapse.py` | dendrogram + heatmap of the atlas across the whole timecourse (slider + play, adjustable fps), **switchable between the full 158-mutant atlas and the 50-leaf functional-annotation subset**; click a mutant → its replicate images, which **follow the time slider** |
| **4** | `IP456_interactiveVolcano.py --mutant BioD` | RNA-seq volcano for Δ*bioD* vs WT; hover a gene → name, locus, description, log₂FC, adjusted *p*; sidebar toggles pathway highlights; search matches gene/locus/description |
| **5** | `IP456_interactiveVolcano.py --mutant PdhE2` | the same explorer for Δ*pdhE2* vs WT |
| **6** | `IP456_interactiveVolcano.py --mutant ManA` | the same explorer for Δ*manA* vs WT |

**1T/2T are alternates of 1/2, not additional figures.** They come from the same script and the same
coordinates; the only difference is which thumbnail pack is inlined. Ship one version of each UMAP,
or renumber them 7/8 if both go into the supplement.

```
figlib.py   thin shim (paths + shared lib); declares the build inputs
data/       coordinates, well labels and dendrogram tables (all small, all bundled)
build/      the three generators (IP12 emits both the still and the timelapse UMAP builds)
archive/    retired builders, kept for reference (not part of the package's output)
figures/    the rendered .html (TRACKED, unlike every other figures/ — see "Where they go")
```

> **Retired:** the static Plot 3 — same viewer, but every replicate shown at its own peak-biomass frame
> and no image animation — is no longer built. Its dual-mode generator is at
> `archive/IP3_interactiveDendrogram.py` (no flag reproduces that old build, `--timelapse` the current
> one, though it still writes the superseded `interactivePlot3b_` filename). The timelapse viewer is
> **Interactive Plot 3** from here on. Note the archived copy predates
> the axis-label and scroll fixes below.

## Build

```bash
python fetch_data.py interactive          # from the repo root -> thumbnail packs + wide table
python build/0_buildFrameMatrices.py                      # ~2 s  -> data/*_frameMatrices.npz (both views)
python build/IP12_interactiveUmap.py --track numerical    # ~5 s  -> 27 MB
python build/IP12_interactiveUmap.py --track embedding    # ~5 s  -> 27 MB
python build/IP12_interactiveUmap.py --track numerical --timelapse   # ~25 s -> 56 MB
python build/IP12_interactiveUmap.py --track embedding --timelapse   # ~25 s -> 56 MB
python build/IP3_interactiveDendrogramTimelapse.py         # ~9 s  -> 30 MB
python build/IP456_interactiveVolcano.py                  # ~5 s  -> 3 x 0.8 MB
```

Plots 4–6 need **no build input at all** — their tables are bundled in `data/`, so they rebuild on a
fresh clone with nothing downloaded.

Step 0 is only needed for Plot 3 and only if the `*_frameMatrices.npz` are absent — they are
bundled, so a normal rebuild starts at `IP3`. It writes one per view (`fullAtlas`,
`functionalSubset`), z-scoring across the whole atlas and only then subsetting to each view's leaves.

### Plot 3: the two views

**View** in the top-left switches the heatmap between the printed panels' two scopes:

| Button | Leaves | Linkage | Printed counterpart |
|---|---|---|---|
| **Full atlas** (opens here) | all 158 clustered mutants | `data/fullAtlas_pcaLinkage_*` | Fig. S6A |
| **Functional annotations** | 50 — the six curated pathways + WT | `data/functionalSubset_pcaLinkage_*` | Fig. 3D |

Each view carries its own Ward-on-PCA-50 linkage and leaf order, taken from the static panel's own
tables, so the interactive dendrogram is the same object the paper prints in both cases. What it does
*not* re-do is the standardization: cells are z-scored **across all 158 mutants** in both views, which
is exactly what `generateDendogramHeatmaps-withPCA.py` does for the functional panel (it z-scores the
whole atlas and selects columns afterwards). So a mutant's colors do not change when you switch — only
which rows are shown and how they are ordered. Standardizing within the 50 leaves instead would make
the same cell two different colors in the two views and break the comparison with Fig. S6A.

The feature axis, the time slider, the search box and the replicate galleries are shared; switching
keeps the current timepoint, re-runs the search against the new leaf set, and keeps the open genotype's
gallery open when that genotype is in the new view.

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

### Plots 1T/2T: the UMAP with timelapses

`--timelapse` swaps the single peak-biomass still in the sidebar for a player: **Play/Pause**, a frame
slider, and a 1–24 fps rate control (arrow keys step, space toggles, exactly as in Plot 3). It opens on
the well's own peak frame, so the first thing shown is what Plots 1/2 show.

The cost is coverage, not the manifold. All **3,669** replicates are still plotted, still searchable and
still clickable — but only the **471** wells in the timelapse pack (3 representative replicates for each
of **157 of 158** mutants) carry every frame. Those are drawn with a thin white ring, and
**Timelapse wells only** hides the rest when the rings get lost in a dense region. Clicking any other
replicate shows its own peak-biomass still from the 112 px pack, with a line saying why it does not
play — no point is worse off than it is in Plots 1/2.

| build | inlined | file |
|---|---|---|
| `--timelapse` (default `--stills all`) | 471 timecourses + 3,198 peak stills | **56 MB** |
| `--timelapse --stills none` | 471 timecourses only | 44 MB |
| no flag (Plots 1/2) | 3,669 peak stills | 27 MB |

`--stills none` is the option to reach for if a journal caps supplementary file size; it leaves the
other 3,198 replicates plotted and hoverable but imageless, so prefer the default unless the size is
actually blocking.

Frames are the pack's 96 px previews scaled up to fill the sidebar, with the browser's default
smoothing — at 4.8× nearest-neighbour turns microcolony texture into noise, while the interpolated
image still reads as a biofilm. The footer under the player states the source resolution so the frame
is not mistaken for full detail. Plots 1/2 are untouched and still show their 112 px stills at native
size.

### `--nas`: the lab builds that link images instead of inlining them

⚠ **Not depositable, and not part of the published package.** These open only where the image trees
are mounted, which is exactly what the deposited plots are designed not to require. They exist
because the inlining budget — not the science — is what limits Plots 1T/2T to 471 wells at 96 px.
Drop the budget and every one of the 3,669 replicates plays at full detail:

```bash
python build/IP12_interactiveUmap.py --track numerical --nas mp4     # -> ...1N_numericalUmapLinked_mp4.html
python build/IP12_interactiveUmap.py --track numerical --nas frames  # -> ...1N_numericalUmapLinked_frames.html
```

| mode | what a click plays | source | needs |
|---|---|---|---|
| `mp4` | the well's own `processed.mp4`, full resolution (~1980², 31 frames over ~15 s), with the browser's transport controls | the reimaging image tree | nothing — the videos already exist |
| `frames` | the same frame slider as Plots 1T/2T, but every well at 256 px and figure-grade contrast | the frame pack | `v2/reimaging/umap/exportFramePack.py` (119,877 images, 2.0 GB, ~1 h) |

Both land at **12.5 MB** — no images at all, just the 9-combo coordinate grid, which is now the
whole file. Coverage goes from 471 wells to all 3,669, so the ring and the **Timelapse wells only**
filter disappear: in these builds everything plays.

`mp4` costs nothing to produce but carries the flatter tonal range that
`v2/common/stills.py` exists to avoid — fine for browsing, not for reading texture. `frames` applies
the standard 0.5–99.5 percentile stretch, so it matches the figures.

### The page finds its own images — put it anywhere

Nothing about the asset location is baked in. The build ships a list of directory **tails**, and on
load the page hunts for them: its build-time guess first, then `./`, `../`, `../../` … up to
`--search-depth` (8) levels above wherever the file actually sits. The first candidate that yields a
real image or video wins, and a status chip in the toolbar reports which. So the same file works from
the mirror, from a local checkout, from a mounted lab share, from a USB stick — anywhere the
assets are at or above it. Tails searched:

| mode | tails |
|---|---|
| `frames` | `data/v2/reimaging/framesPack`, `v2/reimaging/framesPack`, `framesPack` |
| `mp4` | `Sehna/multiphenotype-data-06-25-26/reimagingData`, `multiphenotype-data-06-25-26/reimagingData`, `reimagingData` |

If the hunt fails, the toolbar offers a **folder box**: paste the path, it is probed before being
accepted, and it is remembered (`localStorage`) for next time. `--asset-tail` adds more tails at build
time for a layout you expect to reuse.

Verified by relocating a built file and letting it re-find its assets:

| where the file was put | resolved to |
|---|---|
| mirror, `frames` | `../../../data/v2/reimaging/framesPack` |
| mirror, `mp4` | `../../../../Sehna/multiphenotype-data-06-25-26/reimagingData` |
| beside the pack | `framesPack` |
| 3 unrelated dirs below the pack | `../../../framesPack` |
| 2 unrelated dirs below `Sehna/` | `../../Sehna/multiphenotype-data-06-25-26/reimagingData` |

Every candidate is **relative**, so no absolute path enters the package and nothing is pinned to one
machine's mount spelling — the Linux mount point on the acquisition box and
`/Volumes/phenotyper` on a Mac are both fine.

⚠ The v2-side `interactiveLandscape.py` does **not** do this search — it writes fixed relative URLs
via `config.nas_relative_url()`, so that one still expects the mirror layout.

**Linked frames are preloaded, and playback waits for them.** An inlined frame is a data URI and is
there the instant it is assigned; a linked one is a network fetch. Swapping `img.src` on a fixed timer
therefore outruns the source on the first pass: the slider advances (it is only a number) while the
visible image sits on the last frame that actually decoded — the peak frame the sidebar opened on,
which reads as "Play does nothing for a while, then starts". So opening a well now pulls its whole
timecourse into the browser cache (progress shows under the player), and the play timer holds instead
of advancing to a frame that has not arrived. A frame that fails to load counts as ready, so one bad
file cannot stall playback on that index forever. Verified over SMB from the mirror: all 31 frames
display with Play pressed immediately on a cold cache.

`v2/reimaging/umap/interactiveLandscape.py` is the v2-side viewer with the same idea; it derives its
relative URLs with `config.nas_relative_url()` rather than hardcoding them.

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
- **`functionalSubset_pcaLinkage_{linkage.npy, cluster_order.csv}`,
  `functionalSubset_peakBiomass_featureMatrix.csv`** — the same three tables for Fig. 3D's 50-leaf
  functional subset, copied from `fig3/data/functional_pcaLinkage_*`. The linkage (49 merges) and the
  leaf order are the subset's own; the feature matrix is byte-identical to the full atlas's, because
  both panels z-score across all 158 mutants and differ only in which columns they draw. It is
  bundled rather than aliased so the two views stay independent if that ever changes.
- **`rnaseq_volcano_<mutant>.csv`** (×3, plots 4–6) — `locus`, `oldLocus`, `gene`, `description`,
  `logFC` (log₂ vs WT), `qvalue` (Benjamini–Hochberg adjusted *p*). All 3,572 quantified CDS per
  contrast, unfiltered. Copies of the Fig S8B–D tables, bundled so this package stands alone.
- **`fullAtlas_frameMatrices.npz`** (208 kB), **`functionalSubset_frameMatrices.npz`** (69 kB) —
  `matrices` (31 × *nLeaves* × 27 float16: frame × mutant in that view's leaf order × feature), plus
  `frames`, `mutants`, `features`. Written together by `build/0_buildFrameMatrices.py`; they drive
  Plot 3's animation in the two views. The subset file is a row-subset of the full one, in its own
  leaf order. NaN marks a cell with no data — see below.

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

`reimaging/thumbnails112/` (Plots 1, 2, 1T, 2T) — one 112 px grayscale JPEG per well at its
peak-biomass frame, named `<plate-well>_t<peakFrame>.jpg`, plus `manifest.csv`. ~12 MB for 3,867 wells.
Produced upstream by `v2/reimaging/umap/exportThumbnailPack.py` from the reprocessed `processed.tif`
stack using the standard 0.5–99.5 percentile stretch, i.e. the same recipe as every other still in the
paper.

`reimaging/thumbnailsTimelapse/` (Plots 3, 1T, 2T) — every frame of three representative wells per
mutant at 96 px, same naming and same stretch, ~25 MB / 15,438 images; from
`v2/reimaging/umap/exportTimelapsePack.py`. 471 of its 498 wells survive the manifold filters and so
appear in 1T/2T; Plot 3 keys it by mutant instead and uses all of them whose mutant is a leaf.

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
> hourly timepoints at an adjustable frame rate, and the view can be switched between the full atlas
> and the subset of mutants carrying the annotated functional groups.
>
> **Interactive Plots 1T and 2T** (alternates of 1 and 2). Interactive UMAP of the transposon
> reimaging landscape based on quantitative features extracted using µPULLI-I (1T) or on PCA-50
> dimensionality reduction of DINOv2 CLS embeddings extracted using µPULLI-DL (2T). Clicking a
> representative replicate, marked with a ring, plays its biofilm timecourse through all hourly
> timepoints at an adjustable frame rate; clicking any other replicate displays its peak biofilm
> biomass image.

> **Interactive Plots 4–6.** Interactive RNA-seq volcano plots of Δ*bioD* (4), Δ*pdhE2* (5) and
> Δ*manA* (6) relative to wild type. Every quantified gene is shown; hovering reports its identifiers,
> annotation, log₂ fold change and adjusted *p*-value, the sidebar highlights curated pathways, and the
> search box locates genes by name, locus or description.

The Plot 3 caption is deliberately not the "…and timelapse video" wording: the viewer carries stills, not
video, and its leaves are the 158 **mutants** (clicking one opens its representative replicates), so
the images are reached through the genotype rather than being individually clickable.

## Notes

- **Sizes:** 27 MB (plots 1, 2), 56 MB (plots 1T, 2T — 17,799 thumbnails; 44 MB with
  `--stills none`), 30 MB (plot 3 — 14,694 thumbnails, every frame of three representative wells per
  mutant), 0.8 MB each (plots 4–6, which carry no images). Each thumbnail is stored once and
  referenced by index.
- **Where they go:** these are journal supplementary files and/or the GitHub Pages viewer. Unlike
  `fig*/figures/`, they are **tracked** — `paper-figures/.gitignore` un-ignores them with
  `!interactive/figures/*.html`, because an interactive page cannot be rebuilt without the deposited
  thumbnail packs and a reader who clones the repo should get it. Only Plots **1, 2, 3, 4, 5, 6**
  belong in the published repo. ⚠️ That un-ignore is a blanket `*.html`, so the `1T`/`2T` alternates
  and the lab-only `--nas` builds (`*1N_*Linked_*.html`, `*2N_*`) are also copied by
  `export_public.sh` and show up as untracked files there. Ship one version of each UMAP, and never
  commit an `N` build — it only opens where the NAS is mounted.
- **Verification:** after building, the JavaScript can be checked with
  `node --check` on the extracted `<script>` block, and every thumbnail index should be < the length
  of `THUMBS`. Both hold for the current build.
