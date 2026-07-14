"""fig5 figlib shim - per-figure paths + dataset display constants + re-export of the shared library.

Figure 5 spans three datasets: compounds (A), K. pneumoniae kleb (B), multispecies (C).
Panel scripts: `from figlib import config, features, plotting` plus the dataset constants they need.
"""
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # -> repo root for figlib_shared
import figlib_shared as _S

features = _S.features
plotting = _S.plotting

HERE = Path(__file__).resolve().parent      # <repo>/fig5/
INPUTS = HERE / 'build' / 'inputs'          # place the KiltHub feature tables / CLS caches here (or set env vars)

config = _S.make_config(
    HERE,
    # A (compounds): build inputs for the biotin barcode regenerator
    CMPD_WIDE=Path(os.environ.get('FIG5_CMPD_WIDE', INPUTS / 'cmpd_260522_collapsedWide.parquet')),
    CLEANDEL_WIDES=[INPUTS / 'cleanDel_260521_collapsedWide.parquet',
                    INPUTS / 'cleanDel_260522_collapsedWide.parquet'],
    REIM_WIDE=Path(os.environ.get('FIG5_REIM_WIDE', INPUTS / 'reimaging_collapsedWide.parquet')),
    # B (kleb): CLS embeddings + master frame (growth filter)
    KLEB_CLS=Path(os.environ.get('FIG5_KLEB_CLS', INPUTS / 'kleb_cls.npy')),
    KLEB_EMBIDX=Path(os.environ.get('FIG5_KLEB_EMBIDX', INPUTS / 'kleb_embIndex.csv')),
    KLEB_FRAME=Path(os.environ.get('FIG5_KLEB_FRAME', INPUTS / 'kleb_master_frame_features.csv')),
    # C (multispecies): CLS embeddings (10X) + index
    MULTI_CLS=Path(os.environ.get('FIG5_MULTI_CLS', INPUTS / 'multispecies_10X_cls.npy')),
    MULTI_EMBIDX=Path(os.environ.get('FIG5_MULTI_EMBIDX', INPUTS / 'multispecies_10X_embIndex.csv')),
)

# ── B: K. pneumoniae kleb strains (waaL dropped) ──
KLEB_FRAMES = list(range(9, 25))
KLEB_NO_GROWTH_FLOOR = 0.005
KLEB_PLATE_NV = {'250311_124651': 'NV_058', '250311_125358': 'NV_059', '250311_130104': 'NV_064',
                 '250311_130813': 'NV_065', '250311_131518': 'NV_066', '250311_132226': 'NV_070'}
KLEB_MUTANT_MAP = {'NV_058': 'WT', 'NV_059': 'WzcQ395K', 'NV_064': 'wcaJ', 'NV_065': 'galU', 'NV_066': 'waaL', 'NV_070': 'mrkA'}
KLEB_ORDER = ['WT', 'WzcQ395K', 'wcaJ', 'galU', 'mrkA']
KLEB_DISPLAY = {'WT': 'WT', 'WzcQ395K': r'$\mathit{wzc}^{Q395K}$', 'wcaJ': r'$\Delta \mathit{wcaJ}$',
                'galU': r'$\Delta \mathit{galU}$', 'mrkA': r'$\Delta \mathit{mrkA}$'}
KLEB_COLORS = {'WT': '#808588', 'WzcQ395K': '#ff6f00', 'wcaJ': '#180ae0', 'galU': '#00a2ff', 'mrkA': '#CF0000'}

# ── C: multispecies ──
SPECIES = ['A. baumannii', 'A. baylyi', 'E. coli', 'E. faecalis', 'K. pneumoniae', 'P. aeruginosa', 'S. aureus', 'V. cholerae']
SPECIES_COLORS = dict(zip(SPECIES, ['#e6194B', '#3cb44b', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#9A6324']))
SPECIES_DISPLAY = {s: '$\\mathit{' + s.replace(' ', '\\ ') + '}$' for s in SPECIES}
MULTI_FRAMES = list(range(9, 24))

# ── A: compound conditions (biotin view) ──
CMPD_MARKERS = {'WT_DMSO': 'o', 'WT_antiBio': 'D', 'bioD_biotin': 'v'}
CMPD_PRETTY = {'WT_DMSO': 'WT + DMSO', 'WT_antiBio': 'WT + anti-biotin', 'bioD_biotin': r'$\Delta bioD$ + biotin'}
