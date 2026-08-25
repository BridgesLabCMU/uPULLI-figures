"""movies figlib shim — per-figure paths + re-export of the shared library.

Build scripts import: `from figlib import config, plotting`.

Supplemental movies of the reimaging atlas. Each is rendered frame-by-frame with the same matplotlib
layout as its static counterpart, then encoded with ffmpeg to both MP4 (H.264) and AVI (MJPEG).

  Movie S3  the Figure 3D functional-annotation dendrogram + heatmap animated across the imaging
            timecourse (hourly frames), 8 fps, titled only "Time = X h".

`config.VIDEOS` is the output directory (videos/), alongside the usual TABLES (data/).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # -> paper-figures/ for figlib_shared
import figlib_shared as _S

features = _S.features
plotting = _S.plotting

HERE = Path(__file__).resolve().parent      # paper-figures/movies/

config = _S.make_config(HERE, VIDEOS=HERE / 'videos')
