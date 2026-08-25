#!/usr/bin/env python3
"""Download the build-layer inputs from the data deposit.

Only the `build/` layer needs these. The `render/` scripts that redraw each panel read the source-data
tables bundled in `fig*/data/`, so a fresh clone can reproduce every figure with no download at all.

Reads `inputs.json` (logical input name -> size, sha256, which figures use it) and populates
`fig*/build/inputs/<logical-name>`, verifying each checksum. Stdlib only — nothing to install beyond
the pinned requirements.

Usage:
  python fetch_data.py                   # everything (~1 GB)
  python fetch_data.py figS4             # only what figS4's build layer needs
  python fetch_data.py --list            # show the manifest and exit
  python fetch_data.py --verify          # re-check checksums of what is already present

Already have the deposit unpacked somewhere? Skip this and point the build layer at it:
  export UPULLI_DATA_ROOT=/path/to/deposit
"""
import argparse
import hashlib
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / 'inputs.json'
CHUNK = 1 << 22


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(CHUNK), b''):
            h.update(block)
    return h.hexdigest()


def human(n):
    return '—' if not n else (f'{n / 1e9:.2f} GB' if n >= 1e9 else f'{n / 1e6:.0f} MB')


ap = argparse.ArgumentParser()
ap.add_argument('figures', nargs='*', help='figure dirs to fetch for (default: all)')
ap.add_argument('--list', action='store_true', help='print the manifest and exit')
ap.add_argument('--verify', action='store_true', help='verify checksums of already-present inputs')
args = ap.parse_args()

spec = json.loads(MANIFEST.read_text())
deposit, inputs = spec['deposit'], spec['inputs']
wanted = set(args.figures) if args.figures else None

selected = {n: e for n, e in sorted(inputs.items())
            if not wanted or (wanted & set(e['figures']))}
if wanted:
    unknown = wanted - {f for e in inputs.values() for f in e['figures']}
    if unknown:
        raise SystemExit(f'no manifest entries for: {", ".join(sorted(unknown))}')

if args.list:
    print(f'{deposit["name"]}  (DOI: {deposit["doi"] or "pending"})\n')
    for name, e in selected.items():
        print(f'  {name}\n      {human(e["bytes"])}  used by {", ".join(e["figures"])}\n'
              f'      {e["description"]}')
    print(f'\n{len(selected)} inputs, {human(sum(e["bytes"] or 0 for e in selected.values()))} total')
    sys.exit(0)

if args.verify:
    bad = ok = absent = 0
    for name, e in selected.items():
        for fig in e['figures']:
            p = HERE / fig / 'build' / 'inputs' / name
            if not p.exists():
                absent += 1
                continue
            if p.is_dir():
                ok += 1          # directory inputs are checksummed as the deposited .zip, not in place
                continue
            if e['sha256'] and sha256(p) != e['sha256']:
                print(f'  [CHECKSUM MISMATCH] {fig}: {name}')
                bad += 1
            else:
                ok += 1
    print(f'\n{ok} verified, {bad} mismatched, {absent} not present')
    sys.exit(1 if bad else 0)

def depositUrl(name, entry):
    """Where a logical input lives in the deposit.

    Repositories like KiltHub/Figshare give every file its own opaque download URL and cannot serve
    nested paths, so each manifest entry may carry an explicit "url". Falling back to base_url + name
    only works on a plain file server that mirrors the logical layout."""
    if entry.get('url'):
        return entry['url']
    if not deposit.get('base_url'):
        return None
    return f"{deposit['base_url'].rstrip('/')}/{name}"


def unpackDirectory(zipPath, dest):
    """Directory inputs are deposited as a .zip (no repository serves a directory)."""
    import zipfile
    dest.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zipPath) as zf:
        zf.extractall(dest)
    zipPath.unlink()


if not any(e.get('url') for e in selected.values()) and not deposit.get('base_url'):
    raise SystemExit(
        'The data deposit URL has not been assigned yet.\n'
        f'  Once the deposit is published, "base_url" and "doi" in {MANIFEST.name} will point at it\n'
        '  and this script will fetch every input listed by `python fetch_data.py --list`.\n\n'
        '  In the meantime, if you already have the inputs, set:\n'
        '      export UPULLI_DATA_ROOT=/path/to/deposit\n'
        '  laid out with the logical names shown by `--list`.')

for name, e in selected.items():
    figs = [f for f in e['figures'] if not wanted or f in wanted]
    first = HERE / figs[0] / 'build' / 'inputs' / name
    isDir = bool(e.get('directory')) or name.endswith('/')
    if not first.exists():
        url = depositUrl(name, e)
        if not url:
            raise SystemExit(f'no deposit URL for {name} — add "url" to its manifest entry')
        print(f'  fetching {name}  ({human(e["bytes"])})')
        if isDir:
            tmp = first.parent / (first.name.rstrip('/') + '.zip')
            tmp.parent.mkdir(parents=True, exist_ok=True)
            try:
                urllib.request.urlretrieve(url, tmp)
            except urllib.error.URLError as exc:
                raise SystemExit(f'download failed for {name}\n  {url}\n  {exc}')
            if e['sha256'] and sha256(tmp) != e['sha256']:
                tmp.unlink()
                raise SystemExit(f'checksum mismatch for {name} — download discarded, please retry')
            unpackDirectory(tmp, first)
        else:
            first.parent.mkdir(parents=True, exist_ok=True)
            try:
                urllib.request.urlretrieve(url, first)
            except urllib.error.URLError as exc:
                raise SystemExit(f'download failed for {name}\n  {url}\n  {exc}')
            if e['sha256'] and sha256(first) != e['sha256']:
                first.unlink()
                raise SystemExit(f'checksum mismatch for {name} — download discarded, please retry')
    for fig in figs[1:]:                       # other figures share one copy via a symlink
        dest = HERE / fig / 'build' / 'inputs' / name
        if not dest.exists():
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.symlink_to(first)
print(f'\n{len(selected)} inputs ready under fig*/build/inputs/')
