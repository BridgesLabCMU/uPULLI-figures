"""figS7 figlib shim — per-figure paths + re-export of the shared library.

Panel scripts import: `from figlib import config, plotting`.

Supplemental Figure S7 = DINOv2-embedding views of the reimaging atlas (three panels):
  A  full dendrogram + heatmap on the top-15 embedding principal components (all 158 mutants),
  B  full reimaging UMAP from the DINOv2 embeddings, colored by functional annotation,
  C  the same embedding-PC dendrogram + heatmap for the functional-annotation subset (50 mutants).
The quantitative-feature dendrogram/heatmap now lives in Fig S6.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # -> paper-figures/ for figlib_shared
import figlib_shared as _S

features = _S.features
plotting = _S.plotting

HERE = Path(__file__).resolve().parent      # paper-figures/figS7/

config = _S.make_config(
    HERE,
    inputs=[
        'reimaging/embeddings/dendroPcaCentroids.csv',
        'reimaging/embeddings/dendroClusterOrder.csv',
    ],
)
