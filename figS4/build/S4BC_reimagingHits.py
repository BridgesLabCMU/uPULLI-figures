"""figS4 build S4B/S4C: per-gene reimaging biomass trajectories for the reimaging-selected mutants.

Builds the source tables for the Chromosome I / II heatmaps (panels B/C) from the AUTHORITATIVE inputs:
  * membership  = the 158-mutant Fig-3 reimaging set (paper-figures/fig3 landscape coords), non-WT (157).
  * phenotype   = ReimagingResults_10x.csv (the original reimaging analysis). Binary: 'Dispersal Defect'
                  as labeled there, everything else -> 'High Biofilm'.
  * trajectories= the reimaging feature data (data/v2/reimaging/reimaging_collapsedWide),
                  per-gene mean biomass over its ~24 replicate wells, normalized to the reimaging WT peak
                  mean (matching the Fig-3 "normalized to WT" convention).

Reads:  reimaging/collapsedWide.parquet (build input; see inputs.json),
        paper-figures/fig3/data/reimagingLandscape_nn10_md0.10_coords.csv,
        reimaging/reimagingResults_10x.csv (build input; see inputs.json)
Writes: data/reimHits_biomass_matrix.csv  (geneLocus x frame, mean normalized biomass, genome-ordered)
        data/reimHits_locus_meta.csv       (geneLocus, chromosome, phenotype)
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> figS4/ for figlib
import numpy as np
import pandas as pd
from figlib import config

HERE = Path(__file__).resolve().parents[1]
WIDE = config.input('reimaging/collapsedWide.parquet')
REIMRESULTS = config.input('reimaging/reimagingResults_10x.csv')
FIG3 = HERE.parent / 'fig3' / 'data' / 'reimagingLandscape_nn10_md0.10_coords.csv'   # bundled, sibling package


def locusSortKey(locus):
    first = str(locus).split('/')[0]
    m = re.match(r'VC_A?(\d+)', first)
    return (1 if 'VC_A' in first else 0, int(m.group(1)) if m else 0)


def chromosome(locus):
    return 'II' if 'VC_A' in str(locus).split('/')[0] else 'I'


# ── most-recent reimaging biomass, normalized to reimaging WT peak mean ──
wide = pd.read_parquet(WIDE)
bcols = sorted([c for c in wide.columns if re.match(r'^biomass_t\d+$', c)], key=lambda c: int(c.split('_t')[1]))
frames = [int(c.split('_t')[1]) for c in bcols]
wtPeak = wide.loc[wide['mutant'] == 'WT', bcols].max(axis=1).mean()
print(f'reimaging WT peak mean = {wtPeak:.5f} ({(wide["mutant"] == "WT").sum()} WT wells); frames {frames[0]}..{frames[-1]}')

# ── membership: the Fig-3 reimaging set (non-WT) ──
fig3 = pd.read_csv(FIG3, usecols=['geneLocus'])
loci158 = sorted(set(fig3['geneLocus'].unique()) - {'WT'}, key=locusSortKey)
print(f'Fig-3 reimaging-selected mutants (non-WT): {len(loci158)}')

# ── authoritative phenotype (binary): Dispersal Defect vs High Biofilm ──
reim = pd.read_csv(REIMRESULTS, usecols=['File', 'Well', 'Gene.Locus', 'Phenotype']).drop_duplicates(['File', 'Well'])
mode = reim.groupby('Gene.Locus')['Phenotype'].agg(lambda s: s.mode().iloc[0] if len(s.mode()) else s.iloc[0])
phenotype = {l: ('Dispersal Defect' if mode.get(l) == 'Dispersal Defect' else 'High Biofilm') for l in loci158}

# ── per-gene mean trajectory (normalized), genome-ordered ──
sub = wide[wide['geneLocus'].isin(loci158)]
mat = (sub.groupby('geneLocus')[bcols].mean() / wtPeak).reindex(loci158)
mat.columns = frames
meta = pd.DataFrame({'geneLocus': loci158})
meta['chromosome'] = meta['geneLocus'].map(chromosome)
meta['phenotype'] = meta['geneLocus'].map(phenotype)

nDisp = (meta['phenotype'] == 'Dispersal Defect').sum()
print(f'phenotype (binary): Dispersal Defect={nDisp}, High Biofilm={len(meta) - nDisp}')
print('by chromosome:', meta.groupby('chromosome')['phenotype'].value_counts().to_dict())
print(f'biomass matrix: {mat.shape} (max normalized = {np.nanmax(mat.values):.2f})')

config.ensure(config.TABLES)
mat.to_csv(config.TABLES / 'reimHits_biomass_matrix.csv')
meta.to_csv(config.TABLES / 'reimHits_locus_meta.csv', index=False)
print('Saved: reimHits_biomass_matrix.csv, reimHits_locus_meta.csv')
