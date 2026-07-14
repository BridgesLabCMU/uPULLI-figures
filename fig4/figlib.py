"""fig4 figlib shim - per-figure paths + re-export of the shared library.

Panel scripts import: `from figlib import config, features, plotting`.
Shared, figure-agnostic code lives in ../figlib_shared.py.
"""
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # -> repo root for figlib_shared
import figlib_shared as _S

features = _S.features
plotting = _S.plotting

HERE = Path(__file__).resolve().parent      # <repo>/fig4/
INPUTS = HERE / 'build' / 'inputs'          # place the KiltHub feature tables here (or set the env vars)

config = _S.make_config(
    HERE,
    # build-layer inputs (from KiltHub): the two clean-deletion plates + the reimaging atlas.
    CLEANDEL_WIDES=[Path(os.environ.get('FIG4_CLEANDEL_260521', INPUTS / 'cleanDel_260521_collapsedWide.parquet')),
                    Path(os.environ.get('FIG4_CLEANDEL_260522', INPUTS / 'cleanDel_260522_collapsedWide.parquet'))],
    REIM_WIDE=Path(os.environ.get('FIG4_REIM_WIDE', INPUTS / 'reimaging_collapsedWide.parquet')),
)

# clean-deletion display: label -> (marker, color). Colors = the matching reimaging functional group.
CLEANDEL = {
    'BioD':  ('^', plotting.FUNCTION_COLORS['Biotin Biosynthesis']),    # -> bioD region
    'ManA':  ('s', plotting.FUNCTION_COLORS['O-Antigen Biosynthesis']),  # -> manA region
    'PdhE2': ('D', plotting.FUNCTION_COLORS['Pyruvate Flux']),           # -> pdhE2 region
}
CLEANDEL_DISPLAY = {'BioD': r'$\Delta\mathit{bioD}$', 'ManA': r'$\Delta\mathit{manA}$', 'PdhE2': r'$\Delta\mathit{pdhE2}$'}
