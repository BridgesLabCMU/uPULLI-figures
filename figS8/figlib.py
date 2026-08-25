"""figS8 figlib shim — per-figure paths + re-export of the shared library.

Panel scripts import: `from figlib import config, plotting, STRAIN_ORDER, DISPLAY_NAMES,
FAMILY_COLORS, FAMILY_LABELS, FAMILY_ORDER, featureFamily`.

Supplemental Figure S8 = the Figure 2C training-set heatmap with the FEATURES hierarchically
clustered instead of held in their fixed family order. Same matrix, same colors; the question is
which measurements behave alike across the eight mutants, rather than how the mutants relate.
No build layer -- the source table is the one bundled for Fig. 2C.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # -> paper-figures/ for figlib_shared
import figlib_shared as _S

features = _S.features
plotting = _S.plotting

# feature families: shared with Fig S2D so the colors and names cannot drift
FAMILY_COLORS = _S.FAMILY_COLORS
FAMILY_LABELS = _S.FAMILY_LABELS
FAMILY_ORDER = _S.FAMILY_ORDER
featureFamily = _S.feature_family

HERE = Path(__file__).resolve().parent      # paper-figures/figS8/

config = _S.make_config(
    HERE,
    inputs=[
        'rnaseq/',        # B-D: per-contrast differential-expression tables
        'lux/Lum.csv',    # E: PvpsL-lux plate reader
        'lux/OD.csv',
    ],
)

# ── B-E: the clean-deletion set (shared with fig4, whose CLEANDEL this mirrors) ──
# Colors are the matching reimaging functional group, so a mutant keeps one color across every figure.
CLEANDEL_COLORS = {'BioD': plotting.FUNCTION_COLORS['Biotin Biosynthesis'],
                   'ManA': plotting.FUNCTION_COLORS['O-Antigen Biosynthesis'],
                   'PdhE2': plotting.FUNCTION_COLORS['Pyruvate Flux']}
CLEANDEL_DISPLAY = {'BioD': r'$\Delta\mathit{bioD}$', 'ManA': r'$\Delta\mathit{manA}$',
                    'PdhE2': r'$\Delta\mathit{pdhE2}$'}

# B-D: mutant -> (panel letter, RNA-seq contrast stem). Letters B/E/M vs W are the lab's file naming;
# each was verified by the deleted gene being the extreme negative outlier of its table.
RNASEQ_PANELS = {'BioD': ('S8B', 'BvW'), 'PdhE2': ('S8C', 'EvW'), 'ManA': ('S8D', 'MvW')}

# E: lux plate layout, from Continuous_Peak_Plotting.R next to the raw CSVs. Row B (pdhR) is NOT used.
LUX_WELLS = {'WT': [f'A{i}' for i in range(1, 10)], 'PdhE1': [f'C{i}' for i in range(1, 10)],
             'BioD': [f'D{i}' for i in range(1, 10)], 'ManA': [f'E{i}' for i in range(1, 10)],
             'PdhE2': [f'F{i}' for i in range(1, 10)]}
LUX_COLORS = {'WT': 'black', **CLEANDEL_COLORS, 'PdhE1': plotting.FUNCTION_COLORS['Pyruvate Flux']}
LUX_DISPLAY = {'WT': 'WT', 'PdhE1': r'$\Delta\mathit{pdhE1}$', **CLEANDEL_DISPLAY}

# the 8-mutant training set — display order and math-italic labels (mirrors fig2/figlib.py)
STRAIN_ORDER = ['WT', 'vpsL', 'rbmB', 'hapR', 'potD1', 'flaA', 'luxO_D47E', 'vpvC_W240R']
DISPLAY_NAMES = {
    'WT': r'WT', 'vpsL': r'$\Delta \mathit{vpsL}$', 'rbmB': r'$\Delta \mathit{rbmB}$',
    'hapR': r'$\Delta \mathit{hapR}$', 'potD1': r'$\Delta \mathit{potD1}$', 'flaA': r'$\Delta \mathit{flaA}$',
    'luxO_D47E': r'$\mathit{luxO}^{D47E}$', 'vpvC_W240R': r'$\mathit{vpvC}^{W240R}$',
}
