"""interactive figlib shim — per-figure paths + re-export of the shared library.

Build scripts import: `from figlib import config, plotting`.

The three interactive supplements, all of the reimaging transposon atlas and all self-contained
single-file HTML (every thumbnail inlined as a base64 data URI, so they work offline by
double-clicking, with no server and no sibling asset folder):
  1  UMAP of the atlas from the quantitative uPULLI-I features; click a replicate -> its peak
     biofilm biomass image,
  2  the same explorer over the PCA-50 DINOv2 CLS embedding manifold (uPULLI-DL),
  3  dendrogram + heatmap of the atlas; click a mutant -> the peak-biomass image of every one of
     its replicates,
  4-6  RNA-seq volcano explorers for the three clean deletions (bioD, pdhE2, manA vs WT), with
     pathway highlighting, gene search and per-gene tooltips.

Only the thumbnail pack is a build input; the coordinates, labels, dendrogram tables and volcano
tables are bundled in data/.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # -> paper-figures/ for figlib_shared
import figlib_shared as _S

features = _S.features
plotting = _S.plotting

HERE = Path(__file__).resolve().parent      # paper-figures/interactive/

config = _S.make_config(
    HERE,
    inputs=['reimaging/thumbnails112/', 'reimaging/thumbnailsTimelapse/', 'reimaging/collapsedWide.parquet'],
)

# ── Plots 4-6: RNA-seq volcano explorers ──────────────────────────────────────
# mutant key -> (plot number, header shown in the page). Tables come from figS8's build layer and are
# bundled here as data/rnaseq_volcano_<mutant>.csv so this package stands alone.
VOLCANO_PANELS = {
    'BioD':  (4, 'Δ<i>bioD</i> relative to WT'),
    'PdhE2': (5, 'Δ<i>pdhE2</i> relative to WT'),
    'ManA':  (6, 'Δ<i>manA</i> relative to WT'),
}

# Significant genes take the mutant's own color -- the matching reimaging functional group, the same
# color it carries in Figs 3-5 and in the printed volcanoes (Fig S8B-D).
VOLCANO_SIG_COLORS = {
    'BioD':  _S.FUNCTION_COLORS['Biotin Biosynthesis'],     # orange
    'PdhE2': _S.FUNCTION_COLORS['Pyruvate Flux'],           # green
    'ManA':  _S.FUNCTION_COLORS['O-Antigen Biosynthesis'],  # blue
}
# Search hits need to stay distinct from every one of those, so they are pink rather than the green
# the lab's explorer used (which would collide with the ΔpdhE2 page).
VOLCANO_SEARCH_COLOR = '#ff2d95'

# Pathways selectable in the sidebar. The six locus-keyed functional groups are the shared ones used
# by every other figure (identical loci AND identical colors, so a group reads the same everywhere);
# the biofilm sets below them are keyed by GENE NAME and are specific to the transcriptional view.
VOLCANO_HIGHLIGHTS = {
    **{_S.functionLabel(k): v for k, v in _S.HIGHLIGHT_SETS.items()},
    '<i>vps I</i>': ['vpsU', 'vpsA', 'vpsB', 'vpsC', 'vpsD', 'vpsE', 'vpsF', 'vpsG', 'vpsH', 'vpsI', 'vpsJ', 'vpsK'],
    '<i>vps II</i>': ['vpsL', 'vpsM', 'vpsN', 'vpsO', 'vpsP', 'vpsQ'],
    'Matrix Proteins': ['rbmA', 'rbmB', 'rbmC', 'rbmD', 'rbmEF'],
    'Type IV Pilus': ['mshI', 'mshJ', 'mshK', 'mshL', 'mshM', 'mshN', 'mshE', 'mshG', 'mshF', 'mshB',
                      'mshA', 'mshC', 'mshD', 'mshO', 'mshP', 'mshQ'],
    'Adhesins': ['lapD', 'lapG', 'lapB', 'lapC', 'craA', 'frhA', 'bap1'],
}
VOLCANO_COLORS = {
    **{_S.functionLabel(k): v for k, v in _S.FUNCTION_COLORS.items()},
    '<i>vps I</i>': '#27272A', '<i>vps II</i>': '#52525B', 'Matrix Proteins': '#71717A',
    'Type IV Pilus': '#A1A1AA', 'Adhesins': '#D4D4D8',
}
