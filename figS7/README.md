# Supplemental Figure S7: DINOv2-Embedding Views of the Reimaging Atlas

The deep-embedding counterparts of the quantitative-feature clustering (Fig S6). Three panels.

## Panels

| Panel | View | `render/` | `data/` | `build/` |
|---|---|---|---|---|
| **A** | Full dendrogram + heatmap on the top-15 embedding PCs (all 158 mutants) | `S7A_embPCDendrogramHeatmap.py` (+ `_vertical`) | `embPC15_linkage.npy`, `embPC15_cluster_order.csv`, `embPC15_heatmap_matrix.csv` | `S7A_embPC15.py` |
| **B** | Full reimaging UMAP from DINOv2 embeddings, colored by functional annotation | `S7B_embeddingUmap.py` | `embUmap_pca50_nn10_md0.1_coords.csv` | (v2 `buildEmbeddingUmaps_pca50.py`) |
| **C** | Top-15 embedding-PC dendrogram + heatmap, functional-annotation subset (50 mutants) | `S7C_embPCDendrogramHeatmap_functional.py` | `embPC15func_linkage.npy`, `embPC15func_cluster_order.csv`, `embPC15func_heatmap_matrix.csv` | `S7C_embPC15_functional.py` |

- **A** clusters all 158 mutants (Ward) in the top-15 DINOv2-embedding principal-component space and
  shows those 15 PCs as a z-scored heatmap; leaves colored by functional annotation. Built by
  `S7A_embPC15.py` from the saved 50-PC embedding centroids. Horizontal + vertical renders.
- **B** is the embedding UMAP (CLS PCA-50; nn=10, md=0.1) colored by functional annotation; the pathways
  recover as distinct regions of the manifold.
- **C** is **A** restricted to the functionally-annotated subset (six pathways + WT), with a functional
  legend.

Renders A and C are large (all 158 mutants along one axis); the heatmap/strip are rasterized and the SVG
saved at dpi 200 so the embedded raster stays under the SVG 32767-px image limit. The quantitative-feature
dendrogram/heatmap is **Fig S6**; the low-biofilm chromosome heatmaps are **Fig S5**.
