"""figS6 figlib shim — per-figure paths + re-export of the shared library.

Panel scripts import: `from figlib import config, plotting`.

Supplemental Figure S6 = the full reimaging atlas from the hand-engineered QUANTITATIVE features:
  A          all 158 mutants clustered: PCA-50 -> per-mutant centroid -> Ward dendrogram + z-scored
             peak-biomass feature heatmap (vertical render; a horizontal variant is also provided),
  B left     within- vs between-mutant pairwise-distance distributions,
  B right    paired within- vs between-mutant distance, one line per mutant,
  C left     observed vs permutation-null distribution of per-mutant within-replicate distance,
  C right    every mutant's within-replicate distance against its permutation null.
One claim per row, each shown pooled (left) then resolved per mutant (right): row B = replicates are
closer to each other than to other mutants; row C = replicates are closer than chance.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # -> paper-figures/ for figlib_shared
import figlib_shared as _S

features = _S.features
plotting = _S.plotting

HERE = Path(__file__).resolve().parent      # paper-figures/figS6/

config = _S.make_config(
    HERE,
    inputs=['reimaging/collapsedWide.parquet'],   # build layer for panels B and C only
)
