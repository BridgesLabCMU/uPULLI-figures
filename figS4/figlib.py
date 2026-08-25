"""figS4 figlib shim — per-figure paths + re-export of the shared library.

Panel scripts import: `from figlib import config, plotting, PHENO_COLORS, PHENO_ORDER`.
Shared, figure-agnostic code lives in paper-figures/figlib_shared.py.

Supplemental Figure S4 is the genome-wide transposon biomass screen: all ~2,900 V. cholerae transposon
mutants' biofilm-biomass trajectories, normalized to the WT peak mean of the training
set, genome-ordered, with WT-anchored phenotypic classification (Normal / Low Biofilm / High Biofilm /
Dispersal Defect).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # -> paper-figures/ for figlib_shared
import figlib_shared as _S

features = _S.features
plotting = _S.plotting

HERE = Path(__file__).resolve().parent      # paper-figures/figS4/

config = _S.make_config(
    HERE,
    inputs=[
        # calls of record — these selected the reimaging library and are read as given, not recomputed
        'transposons/results_10x.csv',
        # trajectories for panels A / S5
        'transposons/master_frame_features.csv',
        # WT reference for the trajectory normalizer
        'training/master_frame_features.csv',
        'training/layout.csv',
        # panels B / C: reimaging set membership + phenotype, and their trajectories
        'reimaging/reimagingResults_10x.csv',
        'reimaging/collapsedWide.parquet',
    ],
)

# WT wells on the training plates that define the normalization baseline. The authoritative labels come
# from TRAINING_LAYOUT (mutant == 'WT'); this list is the documented cross-check the build asserts against.
WT_WELLS = ['A5', 'B5', 'C5', 'D5', 'E11', 'F11', 'G11', 'H11']

# phenotype classes: display order + colors (used by both panels).
PHENO_ORDER = ['Low Biofilm', 'Normal', 'Dispersal Defect', 'High Biofilm']
PHENO_COLORS = {'Low Biofilm': '#2c6fbb', 'Normal': '#b0b0b0',
                'Dispersal Defect': '#f39c12', 'High Biofilm': '#c0392b'}
