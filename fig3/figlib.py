"""fig3 figlib shim - per-figure paths + re-export of the shared library.

Panel scripts import: `from figlib import config, features, plotting`.
The shared, figure-agnostic code (colors, feature filter, plot style) lives in ../figlib_shared.py.
"""
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # -> repo root for figlib_shared
import figlib_shared as _S

features = _S.features
plotting = _S.plotting

HERE = Path(__file__).resolve().parent      # <repo>/fig3/
INPUTS = HERE / 'build' / 'inputs'          # place the KiltHub feature tables here (or set the env vars)

config = _S.make_config(
    HERE,
    # build-layer inputs (from KiltHub). MASTER_FRAME = the original full per-(well,frame) feature set
    # (step 0 builds WIDE from it). Override via env, e.g. FIG3_WIDE_TABLE=/path/to/reimaging_collapsedWide.parquet
    MASTER_FRAME=Path(os.environ.get('FIG3_MASTER_FRAME', INPUTS / 'master_frame_features.csv')),
    WIDE=Path(os.environ.get('FIG3_WIDE_TABLE', INPUTS / 'reimaging_collapsedWide.parquet')),
    INDEX=Path(os.environ.get('FIG3_REIMAGING_INDEX', INPUTS / 'reimagingIndex.csv')),
    REIMAGING_INDEX=Path(os.environ.get('FIG3_REIMAGING_INDEX', INPUTS / 'reimagingIndex.csv')),
    EMB=HERE / 'data' / 'reimaging_umapEmbeddings.parquet',   # bundled (small)
)
