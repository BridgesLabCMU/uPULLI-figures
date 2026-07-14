#!/usr/bin/env python3

import os
import re
import multiprocessing as mp
import numpy as np
import pandas as pd
import matplotlib as mpl

mpl.use('Agg')

import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist

import warnings
warnings.filterwarnings(
    'ignore',
    message='Mean of empty slice',
    category=RuntimeWarning
)

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig3/ for figlib
from figlib import config


mpl.rcParams.update({
    'font.family': 'Gillius ADF',
    'font.size': 44,
    'axes.titlesize': 44,
    'axes.labelsize': 44,
    'xtick.labelsize': 44,
    'ytick.labelsize': 44,
    'axes.linewidth': 2,
    'mathtext.fontset': 'stixsans'
})


# v2 port — functional-annotation subset (clusters only highlighted-loci mutants + WT).
# Faithful copy repointed at re-processed data + v2 output. RENDER_FRAMES gates per-frame render.
RENDER_FRAMES = False   # per-frame heatmaps not needed for the figure (peak-biomass-frame heatmap only)
dataPath = str(config.WIDE)                 # raw feature matrix (KiltHub; see build/README.md)
indexPath = str(config.INDEX)               # gene annotations (KiltHub)
outRoot = str(config.TABLES)                # regenerated tables land in ../data/
framesOutDir = str(config.FIGURES / 'dendrogram_frames')

os.makedirs(outRoot, exist_ok=True)
os.makedirs(framesOutDir, exist_ok=True)

clusterCsv = f'{outRoot}/functional_pcaLinkage_cluster_order.csv'
linkagePath = f'{outRoot}/functional_pcaLinkage_linkage.npy'
pcaCentroidPath = f'{outRoot}/functional_pcaLinkage_centroids.csv'
pcaExplainedPath = f'{outRoot}/functional_pcaLinkage_explainedVariance.csv'
peakMatrixPath = f'{outRoot}/functional_peakBiomass_featureMatrix.csv'
peakFrameCsv = f'{outRoot}/functional_peakBiomass_frames.csv'
peakHeatmapPng = str(config.FIGURES / 'functional_peakBiomass_dendoHeatmap.png')
legendPng = str(config.FIGURES / 'functional_legend.png')

plotMinFrame = 0
plotMaxFrame = 30

pcaMinFrame = 0
pcaMaxFrame = 30
pcaComponents = 50

peakMinFrame = 0
peakMaxFrame = 30

cellSize = 0.7
titleFontsize = 52

wtLabel = 'WT'
biomassFeatureBase = 'biomass'
wtPeakFraction = 0.5
lateThreshold = 0.005
lateTimepoint = 28
minReplicates = 5
excludedLoci = {'VC_0185', 'VC_1797', 'VC_2111', 'VC_A1031'}

labelOnlyHighlighted = False
maxWorkers = min(12, mp.cpu_count())

compoundKeywords = ['compound', 'drug', 'treatment', 'uM', 'um', 'µM']

highlightSets = {
    'Motility': [
        'VC_2059','VC_2066','VC_2067','VC_2069','VC_2120','VC_2121','VC_2122','VC_2123',
        'VC_2129','VC_2130','VC_2134','VC_2136','VC_2137','VC_2138','VC_2140','VC_2188',
        'VC_2191','VC_2196','VC_2197','VC_2198','VC_2200','VC_2203','VC_2204','VC_2206',
        'VC_2207','VC_2208'
    ],
    'O-Antigen Biosynthesis': [
        'VC_0212','VC_0223','VC_0239','VC_0241','VC_0242','VC_0245',
        'VC_0247','VC_0249','VC_0250','VC_0251','VC_0259','VC_0269'
    ],
    'Polyamine Import': ['VC_1424','VC_1426','VC_1427','VC_1428'],
    'Biotin Biosynthesis': ['VC_1111','VC_1113','VC_1114','VC_1115'],
    'Pyruvate Flux': ['VC_2413','VC_0943'],
    'Vibriobactin Biosynthesis': ['VC_0771','VC_0772']
}

functionColors = {
    'Motility': '#ff0004',
    'O-Antigen Biosynthesis': '#0096ff',
    'Polyamine Import': '#14f7f0',
    'Biotin Biosynthesis': '#ff9f1c',
    'Pyruvate Flux': '#39ff14',
    'Vibriobactin Biosynthesis': '#a200ff',
    'WT': '#000000',
    'Compound': '#d4af37',
    'Other': '#bdbdbd'
}

canonicalHaralickOrder = [
    'energy','contrast','correlation','variance','inverse_difference_moment',
    'sum_average','sum_variance','sum_entropy','entropy',
    'difference_variance','difference_entropy','imc1','imc2'
]

haralickPretty = {
    'energy': 'Energy',
    'contrast': 'Contrast',
    'correlation': 'Correlation',
    'variance': 'Variance',
    'inverse_difference_moment': 'Inv Diff Moment',
    'sum_average': 'Sum Average',
    'sum_variance': 'Sum Variance',
    'sum_entropy': 'Sum Entropy',
    'entropy': 'Entropy',
    'difference_variance': 'Diff Variance',
    'difference_entropy': 'Diff Entropy',
    'imc1': 'IMC 1',
    'imc2': 'IMC 2'
}

allowedColonyBases = [
    'nColonies',
    'colony_area_um2_mean',
    'colony_area_um2_std',
    'colony_bgCV',
    'colony_centroidOffset_um_mean',
    'colony_majorAxisLength_um_mean',
    'colony_meanIntensity_mean',
    'colony_meanIntensity_kurtosis',
    'colony_mstEdgeMax_um_mean',
    'colony_nnDistance1_um_mean',
    'colony_nnDistance1_um_std',
    'colony_eccentricity_mean'
]

workerDf = None
workerFeatureBases = None
workerOrdered = None
workerDend = None
workerDisplayLabels = None
workerMutantToFunc = None

def initWorker(df, featureBases, ordered, dend, displayLabels, mutantToFunc):
    global workerDf, workerFeatureBases, workerOrdered, workerDend, workerDisplayLabels, workerMutantToFunc
    workerDf = df
    workerFeatureBases = featureBases
    workerOrdered = ordered
    workerDend = dend
    workerDisplayLabels = displayLabels
    workerMutantToFunc = mutantToFunc

def inFrameRange(col, minFrame, maxFrame):
    m = re.search(r'_t(\d+)$', col)
    return m is not None and minFrame <= int(m.group(1)) <= maxFrame

def splitFeatureFrame(col):
    m = re.match(r'(.+)_t(\d+)$', col)
    return (m.group(1), int(m.group(2))) if m else (col, None)

def selectBiomass(cols, minFrame, maxFrame):
    return [c for c in cols if c.startswith('biomass_') and inFrameRange(c, minFrame, maxFrame)]

def selectWhole(cols, minFrame, maxFrame):
    return [
        c for c in cols
        if (c.startswith('whole_haralick_') or c.startswith('whole_entropy_'))
        and inFrameRange(c, minFrame, maxFrame)
    ]

def selectColony(cols, minFrame, maxFrame):
    out = []
    for c in cols:
        if not (c.startswith('colony') or c.startswith('nColonies')):
            continue
        if not inFrameRange(c, minFrame, maxFrame):
            continue
        if c.rsplit('_t', 1)[0] in allowedColonyBases:
            out.append(c)
    return out

def isCompoundLabel(label):
    s = str(label).lower()
    return any(k.lower() in s for k in compoundKeywords)

def applyGrowthFilter(df):
    biomassPattern = re.compile(rf'^{re.escape(biomassFeatureBase)}_t(\d+)$')
    biomassCols = {}

    for col in df.columns:
        m = biomassPattern.match(col)
        if m:
            biomassCols[int(m.group(1))] = col

    if not biomassCols:
        print('[WARN] No biomass columns found; skipping growth filter')
        return df

    allBio = [biomassCols[t] for t in sorted(biomassCols)]
    lateBio = [biomassCols[t] for t in sorted(biomassCols) if t >= lateTimepoint]

    repMax = df[allBio].max(axis=1)
    repLateMax = df[lateBio].max(axis=1) if lateBio else pd.Series(np.inf, index=df.index)

    wtMask = df['mutant'].astype(str) == wtLabel
    wtPeak = repMax[wtMask].median()
    threshold = wtPeakFraction * wtPeak

    growthMask = (repMax >= threshold) | (repLateMax > lateThreshold)

    print(f'[GROWTH] WT peak median: {wtPeak:.4f}')
    print(f'[GROWTH] Threshold: {threshold:.4f}')
    print(f'[GROWTH] Removing {(~growthMask).sum()}/{len(growthMask)} wells')

    return df[growthMask].reset_index(drop=True)

def applyMinReplicateFilter(df):
    counts = df['mutant'].value_counts()
    lowRep = counts[counts < minReplicates].index.tolist()
    lowRepNonCompounds = [x for x in lowRep if x != wtLabel and not isCompoundLabel(x)]

    if lowRepNonCompounds:
        keep = ~df['mutant'].isin(lowRepNonCompounds)
        print(f'[MIN-REP] Removing {len(lowRepNonCompounds)} mutants ({(~keep).sum()} rows)')
        return df[keep].reset_index(drop=True)

    print('[MIN-REP] No mutants removed')
    return df

def zscoreRows(mat):
    return mat.sub(mat.mean(axis=1), axis=0).div(
        mat.std(axis=1).replace(0, np.nan), axis=0
    ).fillna(0)

def featureGroup(base):
    if base == 'biomass':
        return 'biomass'
    if 'haralick' in base:
        return 'haralick'
    if 'entropy' in base:
        return 'entropy'
    return 'colony'

def getFeatureSortKey(base):
    if base == 'biomass':
        return (0, 0)
    if featureGroup(base) == 'colony':
        return (1, base)
    if 'entropy' in base:
        return (2, 0)
    if 'haralick' in base:
        key = re.sub(r'_(mean|std|var)$', '', base.replace('whole_haralick_', ''))
        idx = int(key) if key.isdigit() else (
            canonicalHaralickOrder.index(key) if key in canonicalHaralickOrder else 99
        )
        return (3, idx)
    return (4, base)

def prettyFeatureName(base):
    mapping = {
        'biomass': r'Biofilm Biomass (a.u.)',
        'nColonies': r'Colonies (count)',
        'colony_area_um2_mean': r'Area ($\mu$m$^2$)',
        'colony_area_um2_std': r'Area Variability ($\mu$m$^2$)',
        'colony_bgCV': r'Intensity Variability (CV)',
        'colony_centroidOffset_um_mean': r'Radial Offset ($\mu$m)',
        'colony_eccentricity_mean': r'Eccentricity',
        'colony_majorAxisLength_um_mean': r'Major Axis Length ($\mu$m)',
        'colony_meanIntensity_mean': r'Intensity (a.u.)',
        'colony_meanIntensity_kurtosis': r'Intensity Kurtosis',
        'colony_mstEdgeMax_um_mean': r'Max Distance ($\mu$m)',
        'colony_nnDistance1_um_mean': r'Nearest Neighbor ($\mu$m)',
        'colony_nnDistance1_um_std': r'NN Variability ($\mu$m)'
    }

    if base in mapping:
        return mapping[base]
    if base.startswith('whole_entropy'):
        return r'Global Entropy'
    if base.startswith('whole_haralick_'):
        key = re.sub(r'_(mean|std|var)$', '', base.replace('whole_haralick_', ''))
        if key.isdigit():
            idx = int(key)
            if idx < len(canonicalHaralickOrder):
                key = canonicalHaralickOrder[idx]
        return haralickPretty.get(key, key)
    return base

def cleanName(s):
    """Drop underscores (so mathmode doesn't subscript) + shorten fused loci 'VC_2718; VC_2719'->'VC2718-19'."""
    s = str(s)
    if ';' in s:
        parts = [p.strip().replace('_', '') for p in s.split(';')]
        if len(parts) >= 2 and parts[0].startswith('VC') and parts[1].startswith('VC'):
            return f'{parts[0]}-{parts[1][-2:]}'
    return s.replace('_', '')


def formatMutantLabel(name):
    if name == wtLabel:
        return 'WT'
    if isCompoundLabel(name):
        return cleanName(name)
    return f'${cleanName(name)}$'

def loadMetadata(df):
    if not os.path.exists(indexPath):
        print('[WARN] Index file not found; using mutant names only')
        mutantGene = {m: '' for m in df['mutant'].unique()}
        mutantDisplay = {m: m for m in df['mutant'].unique()}
        return mutantGene, mutantDisplay

    indexDf = pd.read_csv(indexPath)
    keepCols = [c for c in ['plateId', 'wellId', 'geneLocus', 'geneName'] if c in indexDf.columns]
    indexDf = indexDf[keepCols].drop_duplicates()

    meta = df[['plateId', 'wellId', 'mutant']].merge(indexDf, on=['plateId', 'wellId'], how='left')

    mutantGene = (
        meta.groupby('mutant')['geneLocus']
        .apply(lambda x: x.dropna().iloc[0] if x.dropna().shape[0] else '')
        .to_dict()
        if 'geneLocus' in meta.columns else {}
    )

    if 'geneName' in meta.columns:
        mutantDisplay = meta.groupby('mutant').apply(
            lambda g: g['geneName'].dropna().iloc[0]
            if g['geneName'].dropna().shape[0]
            else mutantGene.get(g.name, g.name)
        ).to_dict()
    else:
        mutantDisplay = {m: mutantGene.get(m, m) for m in df['mutant'].unique()}

    mutantDisplay[wtLabel] = wtLabel

    for m in df['mutant'].unique():
        if isCompoundLabel(m):
            mutantGene[m] = ''
            mutantDisplay[m] = m

    return mutantGene, mutantDisplay

def buildMutantFunctionMap(mutantGene):
    mutantToFunc = {}
    for mutant, locus in mutantGene.items():
        for func, genes in highlightSets.items():
            if locus in genes:
                mutantToFunc[mutant] = func
                break
    return mutantToFunc

def mutantColor(mutant, mutantToFunc):
    if mutant == wtLabel:
        return functionColors['WT']
    if isCompoundLabel(mutant):
        return functionColors['Compound']
    if mutant in mutantToFunc:
        return functionColors[mutantToFunc[mutant]]
    return functionColors['Other']

def mutantLabel(mutant, mutantDisplay, mutantToFunc):
    if not labelOnlyHighlighted:
        return formatMutantLabel(mutantDisplay.get(mutant, mutant))

    if mutant == wtLabel or isCompoundLabel(mutant) or mutant in mutantToFunc:
        return formatMutantLabel(mutantDisplay.get(mutant, mutant))

    return ''

def buildPeakFrameMap(df):
    biomassCols = selectBiomass(df.columns, peakMinFrame, peakMaxFrame)

    if not biomassCols:
        raise ValueError('No biomass columns available for peak-frame detection')

    frameByCol = {c: splitFeatureFrame(c)[1] for c in biomassCols}
    peakFrames = {}

    for mutant, g in df.groupby('mutant'):
        biomassTrajectory = {}
        for col in biomassCols:
            frame = frameByCol[col]
            biomassTrajectory[frame] = g[col].median()

        trajectory = pd.Series(biomassTrajectory).sort_index()
        peakFrames[mutant] = np.nan if trajectory.dropna().empty else int(trajectory.idxmax())

    peakDf = pd.DataFrame({'mutant': list(peakFrames.keys()), 'peakFrame': list(peakFrames.values())})
    peakDf.to_csv(peakFrameCsv, index=False)
    print('Saved:', peakFrameCsv)

    return peakFrames

def buildPeakMatrix(df, featureBases, peakFrames):
    grouped = list(df.groupby('mutant'))
    mutants = [m for m, _ in grouped]
    rows = []

    for base in featureBases:
        vals = []
        for mutant, g in grouped:
            peakFrame = peakFrames.get(mutant, np.nan)
            if pd.isna(peakFrame):
                vals.append(np.nan)
                continue

            col = f'{base}_t{int(peakFrame)}'
            vals.append(g[col].median() if col in g.columns else np.nan)

        rows.append(vals)

    return pd.DataFrame(rows, index=featureBases, columns=mutants)

def buildFrameMatrix(df, featureBases, frame):
    grouped = list(df.groupby('mutant'))
    mutants = [m for m, _ in grouped]
    rows = []

    for base in featureBases:
        col = f'{base}_t{frame}'
        vals = [g[col].median() if col in g.columns else np.nan for _, g in grouped]
        rows.append(vals)

    return pd.DataFrame(rows, index=featureBases, columns=mutants)

def buildPcaCentroids(df, featureCols, selected):
    meta = df[['plateId', 'wellId', 'mutant']].copy()
    X = df[featureCols].copy().fillna(0)

    zeroVarCols = X.columns[X.nunique(dropna=False) <= 1].tolist()
    if zeroVarCols:
        X = X.drop(columns=zeroVarCols)
        print(f'[PCA] Removed {len(zeroVarCols)} zero-variance columns')

    nComponents = min(pcaComponents, X.shape[0] - 1, X.shape[1])
    if nComponents < 1:
        raise ValueError('Not enough data/features for PCA')

    Xscaled = StandardScaler().fit_transform(X)

    pca = PCA(n_components=nComponents, random_state=42)
    Xpca = pca.fit_transform(Xscaled)

    explainedDf = pd.DataFrame({
        'pc': [f'PC{i + 1}' for i in range(nComponents)],
        'explainedVarianceRatio': pca.explained_variance_ratio_,
        'cumulativeExplainedVariance': np.cumsum(pca.explained_variance_ratio_)
    })
    explainedDf.to_csv(pcaExplainedPath, index=False)
    print('Saved:', pcaExplainedPath)
    print(f'[PCA] nComponents: {nComponents}')
    print(f'[PCA] cumulative explained variance: {explainedDf["cumulativeExplainedVariance"].iloc[-1]:.4f}')

    pcaDf = pd.DataFrame(
        Xpca,
        columns=[f'PC{i + 1}' for i in range(nComponents)]
    )
    pcaDf['mutant'] = meta['mutant'].values

    centroids = pcaDf.groupby('mutant').median()

    missing = [m for m in selected if m not in centroids.index]
    if missing:
        print(f'[PCA] Dropping {len(missing)} selected mutants absent from PCA centroids')
        selected = [m for m in selected if m in centroids.index]

    centroids = centroids.loc[selected]
    centroids.to_csv(pcaCentroidPath)
    print('Saved:', pcaCentroidPath)

    return centroids, selected

def drawBracket(ax, y0, y1, label):
    # bracket at the RIGHT edge of axBracket (adjacent to y-labels), label just left of it
    ax.plot([0.80, 0.80], [y0, y1], lw=3, color='black', clip_on=False)
    ax.plot([0.80, 0.97], [y0, y0], lw=3, color='black', clip_on=False)
    ax.plot([0.80, 0.97], [y1, y1], lw=3, color='black', clip_on=False)
    ax.text(
        0.42,
        (y0 + y1) / 2,
        label,
        rotation=90,
        va='center',
        ha='center',
        fontsize=44,
        transform=ax.transData,
        clip_on=False
    )
    
    
def drawFigure(heat, ordered, dend, labels, mutantToFunc, title, outPath):
    nF, nM = heat.shape
    # box is nM*cellSize x nF*cellSize inches -> imshow cells are square (cellSize x cellSize).
    # clamp scale so width*savefig_dpi stays under matplotlib's 65536-px Agg cap.
    scale = min(2.0, 60000.0 / (300.0 * max(nM * cellSize + 7, 1)))
    fw = (nM * cellSize + 7) * scale
    fh = (nF * cellSize + 6) * scale
    fig = plt.figure(figsize=(fw, fh), dpi=150)

    # ── inch-based margins (see fullReimaging variant) ── top→bottom: dendrogram, heatmap, mutant
    # x-labels (directly under heatmap), functional strip. left→right: bracket, y-labels, heatmap, cbar.
    _vis = lambda l: re.sub(r'\$|\\mathit|\{|\}|\\Delta|\^|_', '', str(l))
    maxFeatChars = max((len(prettyFeatureName(b)) for b in heat.index), default=10)
    maxMutChars = min(max((len(_vis(l)) for l in labels), default=6), 14)
    yLabelIn = maxFeatChars * (44 / 72.0) * 0.40
    bracketIn = 1.5
    mutLabelIn = maxMutChars * (44 / 72.0) * 0.55 * 0.71
    stripIn = 0.5
    padIn = 0.3

    leftIn = bracketIn + yLabelIn + padIn
    bottomIn = mutLabelIn + stripIn + 1.5 * padIn
    left, bottom = leftIn / fw, bottomIn / fh
    heatW = (nM * cellSize) / fw
    heatH = (nF * cellSize) / fh

    axH = fig.add_axes([left, bottom, heatW, heatH])
    axD = fig.add_axes([left, bottom + heatH + 0.003, heatW, 0.07])
    axD.set_title(title, pad=50, fontsize=titleFontsize, fontweight='bold')
    axC = fig.add_axes([left + heatW + 0.015, bottom, 0.012, heatH])
    stripBottom = max((bottomIn - mutLabelIn - padIn - stripIn) / fh, 0.004)
    axFunc = fig.add_axes([left, stripBottom, heatW, stripIn / fh])
    axBracket = fig.add_axes([0.4 / fw, bottom, bracketIn / fw, heatH])

    icoord, dcoord = dend['icoord'], dend['dcoord']
    xMin = min(min(i) for i in icoord)
    scaleX = (nM - 1) / (max(max(i) for i in icoord) - xMin)
    maxD = max(max(d) for d in dcoord)

    for xs, ys in zip(icoord, dcoord):
        axD.plot([(x - xMin) * scaleX for x in xs], ys, color='black', linewidth=3)

    for i, m in enumerate(ordered):
        color = mutantColor(m, mutantToFunc)
        markerSize = 24 if m == wtLabel else 20
        axD.plot(i, -0.05 * maxD, 'o', color=color, markersize=markerSize, clip_on=False)

    axD.set_xlim(-0.5, nM - 0.5)
    axD.set_ylim(-0.1 * maxD, maxD)
    axD.axis('off')

    im = axH.imshow(heat.values, cmap='RdBu_r', vmin=-3, vmax=3, interpolation='nearest', aspect='auto')

    axH.set_yticks(np.arange(nF))
    axH.set_yticklabels([prettyFeatureName(b) for b in heat.index])
    axH.set_xticks(np.arange(nM))
    axH.set_xticklabels(labels, rotation=45, ha='right', fontsize=44)
    # push labels below the axFunc strip (axFunc bottom is ~0.045*fh below the heatmap)
    axH.tick_params(axis='x', pad=4)  # labels directly under heatmap (strip moved below labels)


    cbar = fig.colorbar(im, cax=axC)
    cbar.set_label(r'Z-score')

    funcColors = [mutantColor(m, mutantToFunc) for m in ordered]
    axFunc.imshow([range(nM)], aspect='auto', cmap=mpl.colors.ListedColormap(funcColors))
    axFunc.axis('off')

    axBracket.set_xlim(0, 1)
    axBracket.set_ylim(nF - 0.5, -0.5)
    axBracket.axis('off')

    groups = [featureGroup(b) for b in heat.index]

    def getRange(name):
        idx = [i for i, g in enumerate(groups) if g == name]
        return (min(idx), max(idx)) if idx else None

    if (r := getRange('colony')):
        drawBracket(axBracket, r[0], r[1], 'Colony Segmentation-\nDerived Features')
    if (r := getRange('haralick')):
        drawBracket(axBracket, r[0], r[1], 'Whole Image\nHaralick Features')

    fig.savefig(outPath, dpi=300, bbox_inches='tight')
    fig.savefig(str(outPath).replace('.png', '.svg'), bbox_inches='tight')
    plt.close(fig)

def renderFrame(frame):
    tag = f't{frame:02d}h'
    outPath = f'{framesOutDir}/dendoHeatmap_{tag}.png'

    frameMatrix = zscoreRows(buildFrameMatrix(workerDf, workerFeatureBases, frame))
    frameMatrix = frameMatrix[workerOrdered]

    drawFigure(
        frameMatrix,
        workerOrdered,
        workerDend,
        workerDisplayLabels,
        workerMutantToFunc,
        f't = {frame}h',
        outPath
    )

    return outPath

def drawLegend():
    fig, ax = plt.subplots(figsize=(9, len(functionColors) * 0.7 + 0.8), dpi=180)
    ax.axis('off')

    for i, (label, color) in enumerate(functionColors.items()):
        y = 1 - (i + 0.5) / len(functionColors)
        ax.scatter(
            0.08,
            y,
            s=250,
            color=color,
            edgecolors='black',
            linewidths=1.2,
            transform=ax.transAxes,
            clip_on=False
        )
        ax.text(0.16, y, label, fontsize=36, va='center', transform=ax.transAxes)

    fig.savefig(legendPng, dpi=300, bbox_inches='tight')
    fig.savefig(str(legendPng).replace('.png', '.svg'), bbox_inches='tight')
    plt.close(fig)
    print('Saved:', legendPng)

def main():
    print('Loading data...')
    df = pd.read_parquet(dataPath)
    print(f'Initial shape: {df.shape}')

    df = applyGrowthFilter(df)
    df = applyMinReplicateFilter(df)
    print(f'Filtered shape: {df.shape}')

    mutantGene, mutantDisplay = loadMetadata(df)

    excludedMutants = {m for m, g in mutantGene.items() if g in excludedLoci}
    if excludedMutants:
        df = df[~df['mutant'].isin(excludedMutants)].reset_index(drop=True)
        for m in excludedMutants:
            mutantGene.pop(m, None)
            mutantDisplay.pop(m, None)
        print(f'[EXCLUDE] Removed {len(excludedMutants)} transposons: {sorted(excludedMutants)}')

    mutantToFunc = buildMutantFunctionMap(mutantGene)

    print('Selecting features...')
    featureCols = (
        selectBiomass(df.columns, plotMinFrame, plotMaxFrame)
        + selectWhole(df.columns, plotMinFrame, plotMaxFrame)
        + selectColony(df.columns, plotMinFrame, plotMaxFrame)
    )

    featureCols = [c for c in featureCols if c in df.columns]
    featureCols = [c for c in featureCols if df[c].notna().any()]
    featureCols = [c for c in featureCols if not (df[c].fillna(0) == 0).all()]
    featureCols = [c for c in featureCols if df[c].nunique(dropna=True) > 1]

    featureBases = sorted(set(splitFeatureFrame(c)[0] for c in featureCols), key=getFeatureSortKey)

    print(f'Selected feature-timepoint columns: {len(featureCols)}')
    print(f'Selected feature bases: {len(featureBases)}')

    print('Computing peak biomass frame per mutant...')
    peakFrames = buildPeakFrameMap(df)

    print('Building peak-biomass display matrix...')
    peakMatrix = zscoreRows(buildPeakMatrix(df, featureBases, peakFrames))
    peakMatrix.to_csv(peakMatrixPath)
    print('Saved:', peakMatrixPath)

    selected = list(peakMatrix.columns)

    allHighlightedLoci = set(g for genes in highlightSets.values() for g in genes)
    selected = [
        m for m in selected
        if m == wtLabel
        or isCompoundLabel(m)
        or mutantGene.get(m, '') in allHighlightedLoci
    ]

    if wtLabel in selected:
        selected = [m for m in selected if m != wtLabel] + [wtLabel]

    if len(selected) < 2:
        raise ValueError('Need at least two mutants/treatments for clustering')

    print(f'Building PCA centroids for functional/highlighted subset: {len(selected)} mutants/treatments')

    pcaFeatureCols = (
        selectBiomass(df.columns, pcaMinFrame, pcaMaxFrame)
        + selectWhole(df.columns, pcaMinFrame, pcaMaxFrame)
        + selectColony(df.columns, pcaMinFrame, pcaMaxFrame)
    )

    pcaFeatureCols = [c for c in pcaFeatureCols if c in df.columns]
    pcaFeatureCols = [c for c in pcaFeatureCols if df[c].notna().any()]
    pcaFeatureCols = [c for c in pcaFeatureCols if not (df[c].fillna(0) == 0).all()]
    pcaFeatureCols = [c for c in pcaFeatureCols if df[c].nunique(dropna=True) > 1]

    print(f'[PCA] selected feature-timepoint columns: {len(pcaFeatureCols)}')

    centroidsPca, selected = buildPcaCentroids(df, pcaFeatureCols, selected)

    print('Clustering PCA centroids with Ward linkage')
    dist = pdist(centroidsPca.values, metric='euclidean')
    Z = linkage(dist, method='ward')
    np.save(linkagePath, Z)
    print('Saved:', linkagePath)

    labelForSelected = [mutantLabel(m, mutantDisplay, mutantToFunc) for m in selected]
    dend = dendrogram(Z, labels=labelForSelected, no_plot=True)
    ordered = [selected[i] for i in dend['leaves']]

    clusterDf = pd.DataFrame({
        'mutant': ordered,
        'display': [mutantDisplay.get(m, m) for m in ordered],
        'gene': [mutantGene.get(m, '') for m in ordered],
        'peakFrame': [peakFrames.get(m, np.nan) for m in ordered],
        'annotation': [
            'WT' if m == wtLabel else
            'Compound' if isCompoundLabel(m) else
            mutantToFunc.get(m, 'Other')
            for m in ordered
        ],
        'color': [mutantColor(m, mutantToFunc) for m in ordered]
    })

    clusterDf.to_csv(clusterCsv, index=False)
    print('Saved:', clusterCsv)

    displayLabels = [mutantLabel(m, mutantDisplay, mutantToFunc) for m in ordered]

    print('Generating peak-biofilm-biomass dendrogram heatmap...')
    peakHeat = peakMatrix[ordered]

    drawFigure(
        peakHeat,
        ordered,
        dend,
        displayLabels,
        mutantToFunc,
        'Peak Biofilm Biomass Frame',
        peakHeatmapPng
    )
    print('Saved:', peakHeatmapPng)

    if RENDER_FRAMES:
        print(f'Generating dendrogram heatmaps for frames t{plotMinFrame:02d}h–t{plotMaxFrame:02d}h using {maxWorkers} workers...')

        frames = list(range(plotMinFrame, plotMaxFrame + 1))

        with mp.Pool(
            processes=maxWorkers,
            initializer=initWorker,
            initargs=(df, featureBases, ordered, dend, displayLabels, mutantToFunc)
        ) as pool:
            for outPath in pool.imap_unordered(renderFrame, frames):
                print('Saved:', outPath)
    else:
        print('[FRAMES] RENDER_FRAMES=False — skipping per-frame heatmaps (peak heatmap only)')

    print('Generating legend...')
    drawLegend()

    print('Summary:')
    print(f'  Total clustered: {len(selected)}')
    print(f'  Highlighted mutants: {sum(1 for m in selected if m in mutantToFunc)}')
    print(f'  Compounds: {sum(1 for m in selected if isCompoundLabel(m))}')
    print(f'  WT: {sum(1 for m in selected if m == wtLabel)}')
    print('Done. Saved all outputs to:', outRoot)

if __name__ == '__main__':
    main()