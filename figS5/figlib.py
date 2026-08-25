"""figS5 figlib shim — per-figure paths + re-export of the shared library.

Panel scripts import: `from figlib import config, plotting`.

Supplemental Figure S5 = per-chromosome biomass heatmaps of the Low-Biofilm transposon mutants
(the low-biofilm class from the genome-wide screen, Fig S4). Source tables (`tn_biomass_matrix.csv`,
`tn_locus_meta.csv`) are the transposon-screen tables produced by figS4/build/S4_screen.py.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # -> paper-figures/ for figlib_shared
import figlib_shared as _S

features = _S.features
plotting = _S.plotting

HERE = Path(__file__).resolve().parent      # paper-figures/figS5/

config = _S.make_config(
    HERE,
    inputs=[
        # calls of record (phenotype classes) + the screen trajectories they color
        'transposons/results_10x.csv',
        'transposons/master_frame_features.csv',
        # WT reference for the trajectory normalizer
        'training/master_frame_features.csv',
        'training/layout.csv',
    ],
)
