"""fig2 figlib shim - per-figure paths + re-export of the shared library.

Panel scripts import: `from figlib import config, features, plotting, STRAIN_ORDER, DISPLAY_NAMES, STRAIN_COLORS`.
Shared, figure-agnostic code lives in ../figlib_shared.py; this file only defines fig2's paths.
"""
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # -> repo root for figlib_shared
import figlib_shared as _S

features = _S.features
plotting = _S.plotting

HERE = Path(__file__).resolve().parent      # <repo>/fig2/
INPUTS = HERE / 'build' / 'inputs'          # place the KiltHub feature tables here (or set the env vars)

config = _S.make_config(
    HERE,
    # build-layer inputs (from KiltHub). Override via env, e.g. FIG2_WIDE_TABLE=/path/to/training_wide.parquet
    WIDE=Path(os.environ.get('FIG2_WIDE_TABLE', INPUTS / 'training_wide.parquet')),
    # MASTER_FRAME = the original full per-(well,frame) training feature set (step-0 input).
    MASTER_FRAME=Path(os.environ.get('FIG2_MASTER_FRAME', INPUTS / 'training_master_frame_features.csv')),
    # OLD_CLEANED = legacy cleaned table used only by step 0 to recover (plateId,wellId)->mutant labels.
    OLD_CLEANED=Path(os.environ.get('FIG2_OLD_CLEANED', INPUTS / 'training_wide_cleaned.parquet')),
    # ACC = single-feature RF accuracies (for the Panel F title lines).
    ACC=Path(os.environ.get('FIG2_SINGLEFEAT_ACC', INPUTS / 'singleFeatureAccuracy.csv')),
)

# the 8-mutant training set: display order, math-italic labels, per-strain colors.
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
