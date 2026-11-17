#!/usr/bin/env python3
"""Turn deposited raw image stacks back into the per-frame files µPULLI-I reads.

The BioImage Archive deposit (S-BIAD3830) stores each well's raw brightfield timelapse as ONE
uint16 OME-TIFF stack, `raw-images/<organism>/…/<plate>/<well>_<obj>x[_<label>]_raw.ome.tiff`.
µPULLI-I instead reads a plate folder of single-frame Cytation 5 TIFFs,
`<well><mag>_1_1_Bright Field_<frame>.tif`, and takes the objective and pixel size from the
Cytation XML header in each file. This script writes exactly that: every frame as its own TIFF
with the original pixel values, and a minimal header carrying the three fields µPULLI-I reads
(ObjectiveSize, PixelWidth, ImageWidthMicrons). Nothing else in the header is used.

Output is grouped by plate (the stack's parent folder name), so a plate split across folders in
the deposit (the multispecies set is filed by species) is reassembled into one plate folder.

  python reproduce/deposit_to_cytation.py <deposit folders or .ome.tiff files…> --out raw_frames --mag _02

Then run µPULLI-I on each `raw_frames/<plate>` with the same `--mag`. Needs numpy + tifffile
(both installed with µPULLI-I). Verified: features from converted reimaging wells are identical to
the published tables.
"""
import argparse
import re
import sys
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np
import tifffile

NAME_RE = re.compile(r'^(?P<well>[A-P]\d{1,2})_(?P<obj>\d+)x(?:_.*)?_raw\.ome\.tiff?$')


def header(well, objective, width, height, umPerPx):
    """Minimal Cytation-style XML. Width in µm is chosen so that width_um / width_px == umPerPx exactly."""
    widthUm = next((f'{c:.{d}f}' for d in range(0, 10) for c in [round(umPerPx * width, d)]
                    if float(f'{c:.{d}f}') / width == umPerPx), repr(umPerPx * width))
    heightUm = repr(float(widthUm) / width * height)
    return ('<?xml version="1.0" encoding="UTF-8"?>\n<BTIImageMetaData>'
            f'<ImageReference><Well>{escape(well)}</Well></ImageReference>'
            '<ImageAcquisition>'
            f'<ObjectiveSize>{objective}</ObjectiveSize>'
            f'<PixelWidth>{width}</PixelWidth><PixelHeight>{height}</PixelHeight>'
            f'<ImageWidthMicrons>{widthUm}</ImageWidthMicrons><ImageHeightMicrons>{heightUm}</ImageHeightMicrons>'
            '</ImageAcquisition></BTIImageMetaData>')


def stacks(paths):
    for p in map(Path, paths):
        if p.is_dir():
            yield from sorted(p.rglob('*_raw.ome.tif*'))
        elif p.is_file():
            yield p
        else:
            sys.exit(f'not found: {p}')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('inputs', nargs='+', help='deposit folders (searched recursively) or *_raw.ome.tiff files')
    ap.add_argument('--out', required=True, help='output root; one sub-folder per plate')
    ap.add_argument('--mag', default='_02',
                    help='magnification token written into the file names; pass the same --mag to µPULLI-I (default _02)')
    a = ap.parse_args()
    if not re.fullmatch(r'_\d+', a.mag):
        sys.exit('--mag must look like _02')

    n = 0
    for f in stacks(a.inputs):
        m = NAME_RE.match(f.name)
        if not m:
            print(f'skip (unrecognised name): {f}')
            continue
        well, objective = m['well'], int(m['obj'])
        with tifffile.TiffFile(f) as t:
            om = t.ome_metadata or ''
            px = re.search(r'PhysicalSizeX="([^"]+)"', om)
            if not px:
                sys.exit(f'no PhysicalSizeX in OME metadata: {f}')
            data = t.asarray()
        if data.ndim == 2:
            data = data[None]
        if data.ndim != 3 or data.dtype != np.uint16:
            sys.exit(f'expected a uint16 T×Y×X stack, got {data.dtype} {data.shape}: {f}')
        T, H, W = data.shape
        xml = header(well, objective, W, H, float(px.group(1)))
        dst = Path(a.out) / f.parent.name
        dst.mkdir(parents=True, exist_ok=True)
        for t in range(T):
            out = dst / f'{well}{a.mag}_1_1_Bright Field_{t + 1:03d}.tif'
            tifffile.imwrite(out, data[t], description=xml, metadata=None, compression='zlib')
        n += 1
        print(f'{f.parent.name}/{well}: {T} frames, {objective}x, {float(px.group(1)):.6g} µm/px')
    print(f'converted {n} wells into {a.out}')


if __name__ == '__main__':
    main()
