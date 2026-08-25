# Supplemental movies

Frame-by-frame renders of the atlas figures across the imaging timecourse, encoded to **MP4**
(H.264) and **AVI** (MJPEG) — the same frames, two containers, for journals that require one or the
other.

| Movie | `build/` script | Content |
|---|---|---|
| **S3** | `S3_dendrogramTimecourse.py` | the Figure 3D functional-annotation dendrogram + heatmap animated over hourly frames 0–30 h, 8 fps |

```
figlib.py   thin shim; adds config.VIDEOS (videos/) to the usual TABLES (data/)
data/       clustering + per-frame matrices (all bundled; no build inputs)
build/      one script per movie
videos/      the encoded .mp4 / .avi (gitignored — rebuild, or take from the deposit)
```

## Build

```bash
python build/S3_dendrogramTimecourse.py                       # ~40 s -> 1.3 MB mp4 + 9.6 MB avi
python build/S3_dendrogramTimecourse.py --fps 4 --width 1280  # slower playback, smaller frame
python build/S3_dendrogramTimecourse.py --keep-frames         # also leave the PNGs in videos/
```

Needs `ffmpeg` on PATH (6.x tested); everything else is the pinned core environment. No build inputs
— the package is self-contained.

## Movie S3

The static Fig. 3D shows the 50-mutant functional-annotation subset at each mutant's *peak-biomass*
frame. Movie S3 renders that identical figure once per hour instead, so each feature class can be
watched developing. **Only the heatmap changes**: dendrogram, leaf order, labels and the annotation
strip are fixed, so any motion on screen is change in the features, not re-clustering. The only text
is the title, `Time = X h`.

- **31 frames, 0–30 h, 8 fps** → 3.9 s. The peak-biomass frame is *not* appended (Fig. 3D already
  shows it); to add it as a final held frame, append the `peak` matrix the way
  `interactive/build/IP3_interactiveDendrogram.py` does.
- **Normalization matches the printed panel**: per-mutant median at each frame, z-scored across the
  full 158-mutant atlas, clipped ±3, then subset to the 50 functional leaves. A color therefore
  means the same thing in the movie and in Fig. 3D. Note each frame is standardized independently,
  so colors compare mutants *within* a timepoint, not across timepoints.
- **Grey cells = no data**, drawn in `#dcdcdc` rather than filled with 0 (a z of 0 sits mid-scale and
  would read as "average"). Two causes, both genuine:
  1. the feature does not exist yet — the 11 colony-segmentation features begin at **t5**, since
     colonies are not segmentable before then;
  2. the feature has no variance across mutants at that frame, so a z-score is undefined —
     `nColonies` through t5 (no colonies detected anywhere) and `biomass` at t0.

  In practice: 13 fully-grey rows at t0, 12 at t3–t4, 1 at t5, none from t6 on.
- **Fixed frame geometry.** Fig. 3D crops with `bbox_inches='tight'`, but a per-frame crop would let
  the frame size drift with the title width ("Time = 0 h" vs "Time = 30 h") and every frame in a
  video must be identical. The builder measures the tight box on the shortest- and longest-titled
  frames, unions them, and saves every frame against that one box — cropped like the figure, yet
  pixel-identical. It asserts all frames match before calling ffmpeg.

### Encoding

| | codec | pixel format | notes |
|---|---|---|---|
| `.mp4` | H.264 (libx264), CRF 18 | yuv420p | `+faststart`; plays everywhere |
| `.avi` | MJPEG, q 3 | yuvj444p | intra-only, larger, maximally compatible |

Both are 1920 px wide (height follows the aspect, forced even via `scale=1920:-2`) at 8 fps.
Rendering happens at `--dpi 60` on the full 47.5 × 25.9 in canvas and is downscaled by ffmpeg with
Lanczos, which is sharper than rendering small.

## Bundled tables

- **`functional_pcaLinkage_{linkage.npy, cluster_order.csv}`** — the Fig. 3D clustering: Ward linkage
  and the 50 leaves in order with display name, locus, annotation and strip color. Copied from
  `fig3/data/` so this package stands alone.
- **`fullAtlas_frameMatrices.npz`** — 31 × 158 × 27 float16 (frame × mutant × feature), NaN where
  there is no data. Identical to the file in `interactive/data/`, produced by
  `interactive/build/0_buildFrameMatrices.py` from the deposited wide table.

> The Fig. 3D and full-atlas peak-biomass feature matrices are byte-for-byte the same 27 × 158 table
> — the functional panel is a 50-column selection from it — which is why the per-frame matrices
> computed for the full atlas can drive this movie directly.
