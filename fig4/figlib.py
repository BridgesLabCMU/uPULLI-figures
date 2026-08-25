"""fig4 figlib shim — per-figure paths + re-export of the shared library.

Panel scripts import: `from figlib import config, features, plotting`.
Shared, figure-agnostic code lives in paper-figures/figlib_shared.py.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # -> paper-figures/ for figlib_shared
import figlib_shared as _S

features = _S.features
plotting = _S.plotting

HERE = Path(__file__).resolve().parent      # paper-figures/fig4/

config = _S.make_config(
    HERE,
    inputs=[
        'cluster/cleanDel_260521_collapsedWide.parquet',
        'cluster/cleanDel_260522_collapsedWide.parquet',
        'reimaging/collapsedWide.parquet',
        'rnaseq/',
    ],
)

# clean-deletion display: label -> (marker, color). Colors = the matching reimaging functional group.
CLEANDEL = {
    'BioD':  ('^', plotting.FUNCTION_COLORS['Biotin Biosynthesis']),    # -> bioD region
    'ManA':  ('s', plotting.FUNCTION_COLORS['O-Antigen Biosynthesis']),  # -> manA region
    'PdhE2': ('D', plotting.FUNCTION_COLORS['Pyruvate Flux']),           # -> pdhE2 region
}
CLEANDEL_DISPLAY = {'BioD': r'$\Delta\mathit{bioD}$', 'ManA': r'$\Delta\mathit{manA}$', 'PdhE2': r'$\Delta\mathit{pdhE2}$'}
