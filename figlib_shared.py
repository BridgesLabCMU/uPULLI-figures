"""Shared support library for the uPULLI-figures packages (fig1..fig5).

Figure-agnostic: the functional-group colors, the UMAP feature filter, the house plot style, and a
`make_config()` factory that builds a per-figure paths namespace. Each figure has a thin `figN/figlib.py`
that imports this and defines its own `config`, so the panel scripts just `from figlib import config,
features, plotting`. No dependency on the analysis monorepo.
"""
import re
from pathlib import Path
from types import SimpleNamespace

import matplotlib as mpl
import numpy as np  # noqa: F401  (available to importers)
import pandas as pd


def ensure(*paths):
    for p in paths:
        p = Path(p)
        (p.parent if p.suffix else p).mkdir(parents=True, exist_ok=True)
    return paths[0] if len(paths) == 1 else paths


# ── features: UMAP feature filter + growth/replicate filters (single source of truth) ──
EXCLUDE_LOCI = ['VC_0185', 'VC_1797', 'VC_2111', 'VC_A1031']
KEEP_FRAMES = list(range(9, 28))
DROP_SUFFIXES = ['_skew', '_kurtosis', '_cv']
NO_GROWTH_FLOOR = 0.005
MIN_REPLICATES = 5
WT_LABEL = 'WT'
BIOMASS_BASE = 'biomass'
ALLOWED_COLONY_BASES = [
    'nColonies', 'colony_area_um2_mean', 'colony_area_um2_std', 'colony_bgCV',
    'colony_centroidOffset_um_mean', 'colony_majorAxisLength_um_mean',
    'colony_meanIntensity_mean', 'colony_meanIntensity_kurtosis', 'colony_mstEdgeMax_um_mean',
    'colony_nnDistance1_um_mean', 'colony_nnDistance1_um_std', 'colony_eccentricity_mean',
]
_TP = re.compile(r'_t(\d+)$')


def _is_allowed_whole(base_name):
    lower = base_name.lower()
    return 'entropy' in lower or 'haralick' in lower


def select_umap_feature_columns(df, include_colony=False):
    cols = []
    for col in df.columns:
        m = _TP.search(col)
        if m is None or int(m.group(1)) not in KEEP_FRAMES:
            continue
        if not pd.api.types.is_numeric_dtype(df[col]):
            continue
        base = _TP.sub('', col)
        if any(col.lower().endswith(s) for s in DROP_SUFFIXES) and not (include_colony and base in ALLOWED_COLONY_BASES):
            continue
        if base.lower().startswith('biomass'):
            cols.append(col)
        elif base.lower().startswith('whole_') and _is_allowed_whole(base):
            cols.append(col)
        elif include_colony and base in ALLOWED_COLONY_BASES:
            cols.append(col)
    return cols


def biomass_columns(columns):
    pat = re.compile(rf'^{re.escape(BIOMASS_BASE)}_t(\d+)$')
    return {int(m.group(1)): c for c in columns if (m := pat.match(c))}


def growth_mask(X):
    bcols = biomass_columns(X.columns)
    if not bcols:
        return pd.Series(True, index=X.index)
    allbio = [bcols[t] for t in sorted(bcols)]
    return X[allbio].max(axis=1) > NO_GROWTH_FLOOR


def low_replicate_genes(gene_series):
    counts = gene_series.value_counts()
    return counts[counts < MIN_REPLICATES].index.tolist()


# display names + units for figure labels (biomass + whole entropy/haralick + the 12 colony bases)
_UM, _UM2 = r'$\mu$m', r'$\mu$m$^2$'
FEATURE_PRETTY = {
    'biomass': 'Biofilm Biomass', 'whole_entropy': 'Global Image Entropy',
    'whole_haralick_0': 'Global Texture Energy', 'whole_haralick_1': 'Global Texture Contrast',
    'whole_haralick_2': 'Global Texture Correlation', 'whole_haralick_3': 'Global Texture Variance',
    'whole_haralick_4': 'Global Texture Homogeneity', 'whole_haralick_5': 'Global Texture Sum Average',
    'whole_haralick_6': 'Global Texture Sum Variance', 'whole_haralick_7': 'Global Texture Sum Entropy',
    'whole_haralick_8': 'Global Texture Entropy', 'whole_haralick_9': 'Global Texture Difference Variance',
    'whole_haralick_10': 'Global Texture Difference Entropy', 'whole_haralick_11': 'Global Texture IMC1',
    'whole_haralick_12': 'Global Texture IMC2', 'nColonies': 'Number of Colonies',
    'colony_area_um2_mean': 'Average Colony Area', 'colony_area_um2_std': 'Colony Area Variability',
    'colony_bgCV': 'Background Intensity Variability (CV)', 'colony_centroidOffset_um_mean': 'Average Colony Radial Offset',
    'colony_majorAxisLength_um_mean': 'Average Colony Major Axis Length', 'colony_meanIntensity_mean': 'Average Colony Intensity',
    'colony_meanIntensity_kurtosis': 'Colony Intensity Tailedness', 'colony_mstEdgeMax_um_mean': 'Largest Inter-Colony Distance',
    'colony_nnDistance1_um_mean': 'Mean Nearest-Neighbor Distance', 'colony_nnDistance1_um_std': 'Nearest-Neighbor Distance Variability',
    'colony_eccentricity_mean': 'Average Colony Eccentricity',
}
FEATURE_UNITS = {
    'biomass': 'a.u.', 'whole_entropy': 'a.u.', 'nColonies': 'count', 'colony_area_um2_mean': _UM2,
    'colony_area_um2_std': _UM2, 'colony_bgCV': '', 'colony_centroidOffset_um_mean': _UM,
    'colony_majorAxisLength_um_mean': _UM, 'colony_meanIntensity_mean': 'a.u.', 'colony_meanIntensity_kurtosis': '',
    'colony_mstEdgeMax_um_mean': _UM, 'colony_nnDistance1_um_mean': _UM, 'colony_nnDistance1_um_std': _UM,
    'colony_eccentricity_mean': '',
}
WHOLE_FEATURE_BASES = ['whole_entropy'] + [f'whole_haralick_{i}' for i in range(13)]
DENDROGRAM_FEATURE_BASES = WHOLE_FEATURE_BASES + list(ALLOWED_COLONY_BASES)


def pretty_name(base):
    return FEATURE_PRETTY.get(base, base)


def feature_unit(base):
    if base in FEATURE_UNITS:
        return FEATURE_UNITS[base]
    if base.startswith('whole_haralick') or base.startswith('whole_entropy'):
        return 'a.u.'
    return ''


features = SimpleNamespace(
    EXCLUDE_LOCI=EXCLUDE_LOCI, KEEP_FRAMES=KEEP_FRAMES, DROP_SUFFIXES=DROP_SUFFIXES,
    NO_GROWTH_FLOOR=NO_GROWTH_FLOOR, MIN_REPLICATES=MIN_REPLICATES, WT_LABEL=WT_LABEL,
    ALLOWED_COLONY_BASES=ALLOWED_COLONY_BASES, select_umap_feature_columns=select_umap_feature_columns,
    biomass_columns=biomass_columns, growth_mask=growth_mask, low_replicate_genes=low_replicate_genes,
    FEATURE_PRETTY=FEATURE_PRETTY, FEATURE_UNITS=FEATURE_UNITS, WHOLE_FEATURE_BASES=WHOLE_FEATURE_BASES,
    DENDROGRAM_FEATURE_BASES=DENDROGRAM_FEATURE_BASES, pretty_name=pretty_name, feature_unit=feature_unit,
)

# ── plotting: shared functional groups + house style ──
HIGHLIGHT_SETS = {
    'Motility': ['VC_2059','VC_2066','VC_2067','VC_2069','VC_2120','VC_2121','VC_2122','VC_2123',
                 'VC_2129','VC_2130','VC_2134','VC_2136','VC_2137','VC_2138','VC_2140','VC_2188','VC_2191',
                 'VC_2196','VC_2197','VC_2198','VC_2200','VC_2203','VC_2204','VC_2206','VC_2207','VC_2208'],
    'O-Antigen Biosynthesis': ['VC_0212','VC_0223','VC_0239','VC_0241','VC_0242','VC_0245','VC_0247','VC_0249',
                               'VC_0250','VC_0251','VC_0259','VC_0269'],
    'Polyamine Import': ['VC_1424','VC_1426','VC_1427','VC_1428'],
    'Biotin Biosynthesis': ['VC_1111','VC_1113','VC_1114','VC_1115'],
    'Pyruvate Flux': ['VC_2413','VC_0943'],
    'Vibriobactin Biosynthesis': ['VC_0771','VC_0772'],
}
FUNCTION_COLORS = {'Motility': '#ff0004', 'O-Antigen Biosynthesis': '#0096ff', 'Polyamine Import': '#14f7f0',
                   'Biotin Biosynthesis': '#ff9f1c', 'Pyruvate Flux': '#39ff14', 'Vibriobactin Biosynthesis': '#ba17f6'}
BACKGROUND_COLOR = '#d0d0d0'


def setStyle(extra=None):
    rc = {'font.family': 'Gillius ADF', 'mathtext.fontset': 'stixsans', 'axes.linewidth': 1.5, 'savefig.dpi': 300}
    if extra:
        rc.update(extra)
    mpl.rcParams.update(rc)
    return rc


plotting = SimpleNamespace(setStyle=setStyle, HIGHLIGHT_SETS=HIGHLIGHT_SETS,
                           FUNCTION_COLORS=FUNCTION_COLORS, BACKGROUND_COLOR=BACKGROUND_COLOR)


def make_config(figdir, **extra):
    """Per-figure paths namespace. figdir = the figN/ directory (contains data/, render/, build/, figures/).
    Standard paths: TABLES (bundled source-data CSVs), FIGURES (output). Extra build-layer inputs
    (WIDE feature matrix, INDEX, EMB, …) are passed as kwargs by the per-figure figlib shim."""
    figdir = Path(figdir)
    base = dict(FIG=figdir, TABLES=figdir / 'data', FIGURES=figdir / 'figures', ensure=ensure)
    base.update(extra)
    return SimpleNamespace(**base)
