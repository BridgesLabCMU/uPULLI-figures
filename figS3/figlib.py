"""figS3 figlib shim — per-figure paths + re-export of the shared library.

Panel scripts import: `from figlib import config, features, plotting, STRAIN_ORDER, DISPLAY_NAMES`.
Shared, figure-agnostic code lives in paper-figures/figlib_shared.py.

Supplemental Figure S3 is the all-feature-class extension of Figure 2F: one temporal (mutant × frame)
heatmap per quantitative feature class, rendered in the same style as Fig 2F so the four features shared
with Fig 2F look identical. Each heatmap is emitted as its own PNG + SVG for single-page layout.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # -> paper-figures/ for figlib_shared
import figlib_shared as _S

features = _S.features
plotting = _S.plotting

HERE = Path(__file__).resolve().parent      # paper-figures/figS3/

config = _S.make_config(
    HERE,
    inputs=[
        'training/wide.parquet',
    ],
)

# figS3 shares the 8-mutant training set display with fig1/fig2.
STRAIN_ORDER = ['WT', 'vpsL', 'rbmB', 'hapR', 'potD1', 'flaA', 'luxO_D47E', 'vpvC_W240R']
DISPLAY_NAMES = {
    'WT': r'WT', 'vpsL': r'$\Delta \mathit{vpsL}$', 'rbmB': r'$\Delta \mathit{rbmB}$',
    'hapR': r'$\Delta \mathit{hapR}$', 'potD1': r'$\Delta \mathit{potD1}$', 'flaA': r'$\Delta \mathit{flaA}$',
    'luxO_D47E': r'$\mathit{luxO}^{D47E}$', 'vpvC_W240R': r'$\mathit{vpvC}^{W240R}$',
}
