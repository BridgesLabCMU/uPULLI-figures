#!/usr/bin/env python3
"""fig3 build 3A: the transposon screen's phenotype calls + class thresholds (source data for 3A).

Panel 3A ranks every screened mutant by peak biofilm biomass and colors it by its phenotype class, with
the B_max / B_min boundaries drawn. Both come from `TransposonResults_10x.csv`, the screen analysis of
record: those calls are what selected the reimaging library, so they are a fixed historical result and
are **read as given, never recomputed** here.

B_min / B_max / maxFinal are *recovered* from the class boundaries in that table — the WT reference set
that originally produced them is not archived alongside it — and each is reported as the midpoint of
the gap that brackets it. So the dashed lines in 3A are the criterion as it was actually applied.

This duplicates the same parse in figS4/build/S4_screen.py, which still needs it for its own
trajectory panels. That is deliberate: each figure package regenerates its own source data (as fig4 and
fig5 each carry their own landscape coords), so neither reaches into the other's data/.

Reads:  config.input('transposons/results_10x.csv')
Writes: data/tn_phenotype_summary.csv   (per well: plate, well, geneLocus, Peak, Final, Early, Phenotype)
        data/tn_thresholds.csv          (Bmax, Bmin, maxFinal, class counts)

Usage:  python fig3/build/3A_rankedTransposonBiomass.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # -> fig3/ for figlib
import pandas as pd
from figlib import config


def normLocus(s):
    """Multi-locus insertions are written `VC_0255; VC_A0373` here and `VC_0255/VC_A0373` in the
    screen index; normalize to the slash form so the two agree."""
    return (s.astype(str).str.replace(r'\s*;\s*', '/', regex=True).str.strip())


def midpoint(below, above, label):
    """Threshold implied by a class boundary: it lies in (max below, min above]. Report the midpoint."""
    lo, hi = float(below.max()), float(above.min())
    print(f'  {label}: in ({lo:.4f}, {hi:.4f}] -> {(lo + hi) / 2:.4f}')
    return (lo + hi) / 2


res = pd.read_csv(config.input('transposons/results_10x.csv'))
summary = (res[['Plate', 'Well', 'Gene.Locus', 'Peak', 'Final', 'Early', 'Phenotype']]
           .drop_duplicates()
           .rename(columns={'Well': 'well', 'Gene.Locus': 'geneLocus'}))
summary['geneLocus'] = normLocus(summary['geneLocus'])
summary['plate'] = summary['Plate'].apply(lambda n: f'TN-Plate{int(n):02d}')
summary = summary.drop(columns='Plate')
counts = summary['Phenotype'].value_counts()
print(f'Calls of record: {len(summary)} wells, {summary["geneLocus"].nunique()} loci')
print(counts.to_string())

print('Thresholds recovered from the class boundaries:')
Bmax = midpoint(summary.loc[summary['Phenotype'] != 'High Biofilm', 'Peak'],
                summary.loc[summary['Phenotype'] == 'High Biofilm', 'Peak'], 'B_max')
Bmin = midpoint(summary.loc[summary['Phenotype'] == 'Low Biofilm', 'Peak'],
                summary.loc[summary['Phenotype'] != 'Low Biofilm', 'Peak'], 'B_min')
inBand = summary[(summary['Peak'] >= Bmin) & (summary['Peak'] <= Bmax)]
maxFinal = midpoint(inBand.loc[inBand['Phenotype'] == 'Normal', 'Final'],
                    inBand.loc[inBand['Phenotype'] == 'Dispersal Defect', 'Final'], 'maxFinal')

config.ensure(config.TABLES)
summary[['plate', 'well', 'geneLocus', 'Peak', 'Final', 'Early', 'Phenotype']].to_csv(
    config.TABLES / 'tn_phenotype_summary.csv', index=False)
pd.DataFrame([{'Bmax': Bmax, 'Bmin': Bmin, 'maxFinal': maxFinal,
               'nWells': int(len(summary)), 'nLoci': int(summary['geneLocus'].nunique()),
               'nHigh': int(counts.get('High Biofilm', 0)), 'nLow': int(counts.get('Low Biofilm', 0)),
               'nDispersal': int(counts.get('Dispersal Defect', 0)),
               'nNormal': int(counts.get('Normal', 0))}]).to_csv(
    config.TABLES / 'tn_thresholds.csv', index=False)
print(f'Saved tn_phenotype_summary.csv + tn_thresholds.csv to {config.TABLES}')
