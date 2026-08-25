#!/usr/bin/env python3
"""Stage build-layer inputs from local storage (maintainer helper).

Reads `inputs.local.json` — a machine-specific map of logical input name -> absolute path, which is
gitignored and never published — and symlinks each figure's declared inputs into
`fig*/build/inputs/<logical-name>`. Symlinks, so staging ~1 GB costs nothing.

Public users don't need this: `fetch_data.py` downloads the same logical names from the data deposit.
Both populate the identical layout, so the build scripts don't know or care which one ran.

Usage:
  python stage_inputs.py                 # stage every figure
  python stage_inputs.py figS4 fig3      # stage only these
  python stage_inputs.py --check         # report what is staged / missing, change nothing
"""
import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / 'inputs.json'
LOCAL = HERE / 'inputs.local.json'

ap = argparse.ArgumentParser()
ap.add_argument('figures', nargs='*', help='figure dirs to stage (default: all)')
ap.add_argument('--check', action='store_true', help='report status only')
args = ap.parse_args()

manifest = json.loads(MANIFEST.read_text())['inputs']
if not LOCAL.exists():
    raise SystemExit(f'{LOCAL.name} not found — this helper is for the maintainer machine only.\n'
                     f'Public users: run `python fetch_data.py` instead.')
local = json.loads(LOCAL.read_text())

wanted = set(args.figures) if args.figures else None
staged = missing = unresolved = 0
for name, entry in sorted(manifest.items()):
    for fig in entry['figures']:
        if wanted and fig not in wanted:
            continue
        dest = HERE / fig / 'build' / 'inputs' / name
        src = local.get(name)
        if src is None or not Path(src).exists():
            print(f'  [no source] {fig}: {name}')
            unresolved += 1
            continue
        if args.check:
            print(f'  {"[ok]      " if dest.exists() else "[missing] "}{fig}: {name}')
            staged += dest.exists()
            missing += not dest.exists()
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.is_symlink() or dest.exists():
            dest.unlink()
        dest.symlink_to(Path(src))
        staged += 1

if args.check:
    print(f'\n{staged} staged, {missing} missing, {unresolved} without a local source')
    sys.exit(1 if (missing or unresolved) else 0)
print(f'\nstaged {staged} input symlinks' + (f' ({unresolved} had no local source)' if unresolved else ''))
