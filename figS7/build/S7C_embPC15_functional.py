"""figS7 build C: top-15 embedding-PC dendrogram + heatmap, FUNCTIONAL-ANNOTATION subset.

Same as build/S7C_embPC15.py but restricted to the functionally-annotated reimaging mutants (the six
highlight pathways + WT; annotation != 'Other'), mirroring the Fig-3D functional subset. Re-clusters
that subset in the top-15 embedding-PC space and writes its leaf order, linkage, and z-scored PC heatmap.

Reads:  reimaging/embeddings/dendroPcaCentroids.csv (build input; see inputs.json),
        results/v2/reimaging/embeddings/dendogram/embDendro_cluster_order.csv (mutant -> annotation)
Writes: data/embPC15func_linkage.npy, data/embPC15func_cluster_order.csv, data/embPC15func_heatmap_matrix.csv
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
annot = pd.read_csv(ANNOT).set_index('mutant')['annotation'].to_dict()
# functional-annotation subset: the six highlight pathways + WT (exclude 'Other')
func = [m for m in cent.index if annot.get(m, 'Other') != 'Other']
pcCols = [f'PC{i}' for i in range(1, N_PC + 1)]
X = cent.loc[func, pcCols]
print(f'functional subset: {X.shape[0]} mutants x top {N_PC} PCs')

Z = linkage(pdist(X.values, metric='euclidean'), method='ward')
ordered = X.index[dendrogram(Z, no_plot=True)['leaves']].tolist()

orderDf = pd.DataFrame({'mutant': ordered})
orderDf['annotation'] = orderDf['mutant'].map(lambda m: annot.get(m, 'Other'))
orderDf['color'] = orderDf['annotation'].map(annotColor)

z = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
heat = z.loc[ordered].T
heat.index = pcCols

config.ensure(config.TABLES)
np.save(config.TABLES / 'embPC15func_linkage.npy', Z)
orderDf.to_csv(config.TABLES / 'embPC15func_cluster_order.csv', index=False)
heat.to_csv(config.TABLES / 'embPC15func_heatmap_matrix.csv')
print('by annotation:', orderDf['annotation'].value_counts().to_dict())
print('Saved: embPC15func_linkage.npy, embPC15func_cluster_order.csv, embPC15func_heatmap_matrix.csv')
