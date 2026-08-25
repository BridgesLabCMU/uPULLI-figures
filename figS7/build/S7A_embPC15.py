"""figS7 build A: dendrogram + heatmap on the top-15 embedding PCs (source data for Panel A).

Embedding analogue of the numerical PCA-linkage dendrogram, restricted to the TOP 15 principal
components of the DINOv2-embedding per-mutant centroids. Reads the saved 50-PC embedding centroids
(one per reimaging mutant), takes PC1..PC15, Ward-clusters the mutants in that 15-PC space, and writes
the leaf order, linkage, and a z-scored PC x mutant heatmap matrix.

Reads:  reimaging/embeddings/dendroPcaCentroids.csv (build input; 158 mutants x 50 PCs),
        results/v2/reimaging/embeddings/dendogram/embDendro_cluster_order.csv (mutant -> annotation)
Writes: data/embPC15_linkage.npy, data/embPC15_cluster_order.csv, data/embPC15_heatmap_matrix.csv
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS7/ for figlib
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
from figlib import config, plotting

HERE = Path(__file__).resolve().parents[1]
CENTROIDS = config.input('reimaging/embeddings/dendroPcaCentroids.csv')
ANNOT = config.input('reimaging/embeddings/dendroClusterOrder.csv')
N_PC = 15

FC, BG = plotting.FUNCTION_COLORS, plotting.BACKGROUND_COLOR


def annotColor(a):
    if a == 'WT':
        return '#000000'
    return FC.get(a, BG)


cent = pd.read_csv(CENTROIDS).set_index('mutant')
pcCols = [f'PC{i}' for i in range(1, N_PC + 1)]
X = cent[pcCols]
print(f'{X.shape[0]} mutants x top {N_PC} PCs')

# Ward clustering in the top-15 PC space
Z = linkage(pdist(X.values, metric='euclidean'), method='ward')
leaves = dendrogram(Z, no_plot=True)['leaves']
ordered = X.index[leaves].tolist()

# leaf annotation + color
annot = pd.read_csv(ANNOT).set_index('mutant')['annotation'].to_dict()
orderDf = pd.DataFrame({'mutant': ordered})
orderDf['annotation'] = orderDf['mutant'].map(lambda m: annot.get(m, 'Other'))
orderDf['color'] = orderDf['annotation'].map(annotColor)

# heatmap matrix: top-15 PCs (rows) x mutants(ordered), z-scored per PC across mutants
z = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
heat = z.loc[ordered].T                       # PC x mutant
heat.index = pcCols

config.ensure(config.TABLES)
np.save(config.TABLES / 'embPC15_linkage.npy', Z)
orderDf.to_csv(config.TABLES / 'embPC15_cluster_order.csv', index=False)
heat.to_csv(config.TABLES / 'embPC15_heatmap_matrix.csv')
print('by annotation:', orderDf['annotation'].value_counts().to_dict())
print('Saved: embPC15_linkage.npy, embPC15_cluster_order.csv, embPC15_heatmap_matrix.csv')
