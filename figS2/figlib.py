"""figS2 figlib shim — per-figure paths + re-export of the shared library.

Panel scripts import: `from figlib import config, features, plotting, STRAIN_ORDER, DISPLAY_NAMES, STRAIN_COLORS`.
Shared, figure-agnostic code lives in paper-figures/figlib_shared.py.

Supplemental Figure S2 dissects the hand-engineered numerical classification of the 8-mutant training
set (the companion to Figure 2's all-feature confusion): per-feature-family RF confusion matrices
(biomass / colony / whole-image), the single-feature timecourse accuracy ranking, and mutant
separability across time. All panels are built from the training wide table.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # -> paper-figures/ for figlib_shared
import figlib_shared as _S

features = _S.features
plotting = _S.plotting

HERE = Path(__file__).resolve().parent      # paper-figures/figS2/

config = _S.make_config(
    HERE,
    inputs=[
        'training/wide.parquet',
    ],
)

# figS2 shares the 8-mutant training set display with fig1/fig2.
STRAIN_ORDER = ['WT', 'vpsL', 'rbmB', 'hapR', 'potD1', 'flaA', 'luxO_D47E', 'vpvC_W240R']
DISPLAY_NAMES = {
    'WT': r'WT', 'vpsL': r'$\Delta \mathit{vpsL}$', 'rbmB': r'$\Delta \mathit{rbmB}$',
    'hapR': r'$\Delta \mathit{hapR}$', 'potD1': r'$\Delta \mathit{potD1}$', 'flaA': r'$\Delta \mathit{flaA}$',
    'luxO_D47E': r'$\mathit{luxO}^{D47E}$', 'vpvC_W240R': r'$\mathit{vpvC}^{W240R}$',
}
STRAIN_COLORS = {
    'WT': '#808588', 'vpsL': '#ff6f00', 'rbmB': '#180ae0', 'hapR': '#00a2ff',
    'potD1': '#b5ff60', 'flaA': '#00ffb3', 'luxO_D47E': '#ffbc3e', 'vpvC_W240R': '#CF0000',
}

# feature-family colors for the single-feature accuracy bar chart (Panel S2D).
FAMILY_COLORS = _S.FAMILY_COLORS      # shared with Fig S8; defined in figlib_shared.py
