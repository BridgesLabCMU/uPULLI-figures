"""fig3 figlib shim — per-figure paths + re-export of the shared library.

Panel scripts import unchanged: `from figlib import config, features, plotting`.
The shared, figure-agnostic code (colors, feature filter, plot style) lives in
paper-figures/figlib_shared.py; this file only defines fig3's paths.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # -> paper-figures/ for figlib_shared
import figlib_shared as _S

features = _S.features
plotting = _S.plotting

HERE = Path(__file__).resolve().parent      # paper-figures/fig3/

config = _S.make_config(
    HERE,
    inputs=[
        'reimaging/master_frame_features.csv',
        'reimaging/geneIndex.csv',
        'reimaging/collapsedWide.parquet',
        'reimaging/umapEmbeddings.parquet',
        # Panel 3A: transposon screen calls of record (shared with figS4, read as given)
        'transposons/results_10x.csv',
    ],
)

# ── Panel 3A: transposon-screen phenotype classes (mirrors figS4/figlib.py; keep the two in step) ──
PHENO_ORDER = ['Low Biofilm', 'Normal', 'Dispersal Defect', 'High Biofilm']
PHENO_COLORS = {'Low Biofilm': '#2c6fbb', 'Normal': '#b0b0b0',
                'Dispersal Defect': '#f39c12', 'High Biofilm': '#c0392b'}
