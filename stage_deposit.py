#!/usr/bin/env python3
"""MAINTAINER: assemble everything that goes to the data deposit (KiltHub), ready to upload.

Repositories like KiltHub/Figshare store a flat list of files: filenames cannot contain `/`, and a
directory cannot be served as such. So this writes, into `deposit/`:

  files/<flat name>        every build input, one file each, `/` replaced by `__`
                           (directory inputs are zipped first — that is what the deposit stores)
  upulli-figure-inputs.zip the same inputs as one archive, in their LOGICAL layout, for anyone who
                           would rather download once and `export UPULLI_DATA_ROOT=<unzipped>`
  source-data-tables.zip   the bundled fig*/data tables, so the deposit is self-describing without
                           requiring a clone
  deposit_manifest.csv     deposited filename <-> logical name, bytes, sha256, which figures use it
  inputs.json              a copy of the manifest, as deposited

After uploading, paste each file's download URL into the matching entry of `inputs.json` as "url"
(and fill deposit.doi). `fetch_data.py` prefers that per-entry URL, which is the only scheme that
works with per-file repositories.

Reads inputs.local.json (the maintainer's absolute paths). Nothing here is published to GitHub —
`deposit/` should stay untracked.

Usage:
  python stage_deposit.py            # symlink the per-file copies (fast, no extra disk)
  python stage_deposit.py --copy     # real copies (needed if your upload tool will not follow links)
  python stage_deposit.py --no-zip   # skip the archives, just the flat files + manifest
"""
import argparse
import csv
import hashlib
import json
import os
import shutil
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEPOSIT = HERE / 'deposit'
CHUNK = 1 << 22

ap = argparse.ArgumentParser()
ap.add_argument('--copy', action='store_true', help='copy inputs instead of symlinking them')
ap.add_argument('--no-zip', action='store_true', help='skip building the archives')
args = ap.parse_args()


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(CHUNK), b''):
            h.update(block)
    return h.hexdigest()


def flatName(logical):
    return logical.rstrip('/').replace('/', '__')


# Never deposit these, wherever they turn up in a source directory: shell history, OS cruft, and
# superseded artifacts that would confuse a reader.
SKIP_NAMES = {'.Rhistory', '.DS_Store', 'Thumbs.db', 'Master_Volcano_Explorer.html'}


def keep(path):
    return path.name not in SKIP_NAMES and not path.name.endswith(('.pyc', '.Rhistory~'))


spec = json.loads((HERE / 'inputs.json').read_text())
local = json.loads((HERE / 'inputs.local.json').read_text())
extras = local.pop('extras', {})     # non-input material deposited alongside (e.g. robot protocols)
missing = [n for n in spec['inputs'] if n not in local]
if missing:
    raise SystemExit('inputs.local.json has no path for: ' + ', '.join(missing))

filesDir = DEPOSIT / 'files'
filesDir.mkdir(parents=True, exist_ok=True)
rows, staged = [], {}

for name, entry in sorted(spec['inputs'].items()):
    src = Path(local[name])
    if not src.exists():
        raise SystemExit(f'missing source for {name}: {src}')
    isDir = bool(entry.get('directory')) or name.endswith('/') or src.is_dir()
    out = filesDir / (flatName(name) + ('.zip' if isDir else ''))

    if isDir:
        if not out.exists():                       # a directory is deposited as one archive
            print(f'  zipping {name} -> {out.name}')
            with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as zf:
                for f in sorted(src.rglob('*')):
                    if f.is_file() and keep(f):
                        zf.write(f, f.relative_to(src))
    else:
        if out.exists() or out.is_symlink():
            out.unlink()
        if args.copy:
            shutil.copy2(src, out)
        else:
            out.symlink_to(os.path.realpath(src))

    real = Path(os.path.realpath(out))
    digest, size = sha256(real), real.stat().st_size
    staged[name] = out
    rows.append({'depositFile': out.name, 'logicalName': name, 'bytes': size, 'sha256': digest,
                 'isDirectoryArchive': isDir, 'figures': ' '.join(entry['figures']),
                 'description': entry['description']})
    flag = '' if (entry['sha256'] in (None, digest) or isDir) else '   [manifest sha differs!]'
    print(f'  {out.name:52s} {size/1e6:8.1f} MB{flag}')

with open(DEPOSIT / 'deposit_manifest.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
shutil.copy2(HERE / 'inputs.json', DEPOSIT / 'inputs.json')

if not args.no_zip:
    bundle = DEPOSIT / 'upulli-figure-inputs.zip'
    print(f'\n  building {bundle.name} (logical layout, for UPULLI_DATA_ROOT users)...')
    with zipfile.ZipFile(bundle, 'w', zipfile.ZIP_DEFLATED) as zf:
        for name in sorted(spec['inputs']):
            src = Path(os.path.realpath(local[name]))
            if src.is_dir():
                for f in sorted(src.rglob('*')):
                    if f.is_file() and keep(f):
                        zf.write(f, f'{name.rstrip("/")}/{f.relative_to(src)}')
            else:
                zf.write(src, name)
    print(f'  {bundle.name}: {bundle.stat().st_size/1e6:.0f} MB')

    tables = DEPOSIT / 'source-data-tables.zip'
    print(f'  building {tables.name} (the bundled fig*/data tables + their column dictionaries)...')
    with zipfile.ZipFile(tables, 'w', zipfile.ZIP_DEFLATED) as zf:
        for d in sorted(HERE.glob('fig*/data')) + [HERE / 'interactive' / 'data']:
            for f in sorted(d.rglob('*')):
                if f.is_file():
                    zf.write(f, str(f.relative_to(HERE)))
            # the figure's README IS the column dictionary for those tables — a download that
            # omitted it would leave the columns undefined
            readme = d.parent / 'README.md'
            if readme.exists():
                zf.write(readme, str(readme.relative_to(HERE)))
        for top in ('README.md', 'INPUTS.md'):
            if (HERE / top).exists():
                zf.write(HERE / top, top)
    print(f'  {tables.name}: {tables.stat().st_size/1e6:.0f} MB')

# ── extras: material that belongs in the deposit but is not a build input, so it is NOT in
# inputs.json (fetch_data.py must not try to download it). Zipped whole; a NOTES.txt, if present,
# is also copied out so it is readable on the landing page without downloading the archive.
for name, src in sorted(extras.items()):
    src = Path(src)
    if not src.exists():
        print(f'  [warn] extra {name}: missing source {src}')
        continue
    # A curated NOTES kept in this repo (protocols_NOTES.txt) supersedes whatever loose notes sit in
    # the source folder: it is the edited, figure-mapped version, and the source share is read-only.
    curated = HERE / 'protocols_NOTES.txt'
    out = DEPOSIT / f'{name}.zip'
    with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as zf:
        for f in sorted(src.rglob('*')):
            if f.is_file() and keep(f) and not (curated.exists() and f.name.upper().startswith('NOTES')):
                zf.write(f, f.relative_to(src))
        if curated.exists():
            zf.write(curated, 'NOTES.txt')
    print(f'  {out.name:52s} {out.stat().st_size/1e6:8.1f} MB   (extra: {src.name})')
    notes = curated if curated.exists() else next((f for f in src.glob('NOTES*')), None)
    if notes:
        dest = DEPOSIT / f'{name}_NOTES.txt'
        if dest.exists():
            dest.chmod(0o644)       # sources on the read-only share copy in as 0444
            dest.unlink()
        shutil.copy2(notes, dest)
        dest.chmod(0o644)
        print(f'  {dest.name}')

print(f'\n{len(rows)} inputs staged in {DEPOSIT}')
print('Next: upload deposit/, then paste each download URL into inputs.json as "url" (+ deposit.doi).')
