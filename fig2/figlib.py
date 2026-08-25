"""fig2 figlib shim — per-figure paths + re-export of the shared library.

Panel scripts import unchanged: `from figlib import config, features, plotting`.
Shared, figure-agnostic code lives in paper-figures/figlib_shared.py; this file only defines fig2's paths.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # -> paper-figures/ for figlib_shared
import figlib_shared as _S

features = _S.features
plotting = _S.plotting

HERE = Path(__file__).resolve().parent      # paper-figures/fig2/

config = _S.make_config(
    HERE,
    inputs=[
        'training/master_frame_features.csv',
        'training/layout.csv',
        'training/wide.parquet',
    ],
)

# fig2-specific: the 8-mutant training set — display order, math-italic labels, per-strain colors.
# feature families: shared with Fig S2D so the colors and names cannot drift (2C's row strip)
FAMILY_COLORS = _S.FAMILY_COLORS
FAMILY_LABELS = _S.FAMILY_LABELS
FAMILY_ORDER = _S.FAMILY_ORDER
featureFamily = _S.feature_family

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

