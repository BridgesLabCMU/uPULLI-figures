#!/usr/bin/env python3
"""MAINTAINER: upload the staged deposit to KiltHub (Figshare) over the API.

Stdlib only, like fetch_data.py. Reads the token from ~/.config/kilthub/token (mode 600) or
$KILTHUB_TOKEN — never from the repo, never on the command line.

The item is created as a PRIVATE draft and stays that way: there is deliberately no publish call in
this script. Publishing is a decision to make in the web UI, once.

  python kilthub_upload.py --create "uPULLI figure build inputs"   # -> prints the new article id
  python kilthub_upload.py --article 12345678                      # upload everything in deposit/
  python kilthub_upload.py --article 12345678 --verify             # remote md5/size vs local
  python kilthub_upload.py --article 12345678 --dry-run            # list what would upload
  python kilthub_upload.py --article 12345678 --dedupe             # drop duplicate uploads

Item metadata (title, authors, categories, license, funding, ...) is filled in the KiltHub submission
form, not here — this script only moves files.

Uploads are idempotent: a file already present remotely with a matching md5 is skipped, so an
interrupted run can simply be repeated.
"""
import argparse
import hashlib
import json
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path

API = 'https://api.figshare.com/v2'
HERE = Path(__file__).resolve().parent
DEPOSIT = HERE / 'deposit'
CHUNK = 1 << 22


def token():
    tok = os.environ.get('KILTHUB_TOKEN')
    if tok:
        return tok.strip()
    p = Path.home() / '.config' / 'kilthub' / 'token'
    if not p.exists():
        raise SystemExit(f'no token: set $KILTHUB_TOKEN or write it to {p} (chmod 600)')
    return p.read_text().strip()


TOKEN = token()


def call(method, url, body=None, raw=None, ctype='application/json'):
    data = raw if raw is not None else (json.dumps(body).encode() if body is not None else None)
    req = urllib.request.Request(url if url.startswith('http') else API + url, data=data, method=method)
    req.add_header('Authorization', f'token {TOKEN}')
    if data is not None:
        req.add_header('Content-Type', ctype)
    try:
        with urllib.request.urlopen(req) as r:
            payload = r.read()
            return json.loads(payload) if payload and r.headers.get_content_type() == 'application/json' else payload
    except urllib.error.HTTPError as e:
        raise SystemExit(f'{method} {url} -> HTTP {e.code}\n  {e.read()[:400].decode(errors="replace")}')


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(CHUNK), b''):
            h.update(block)
    return h.hexdigest()


def depositFiles():
    """Everything staged for upload: the flat per-input files plus the archives and manifests."""
    if not DEPOSIT.exists():
        raise SystemExit(f'{DEPOSIT} does not exist — run: python stage_deposit.py --copy')
    out = []
    for p in sorted(DEPOSIT.rglob('*')):
        if p.is_file() and not p.is_symlink():
            out.append(p)
    links = [p for p in DEPOSIT.rglob('*') if p.is_symlink()]
    if links:
        raise SystemExit(f'{len(links)} staged files are symlinks (e.g. {links[0].name}) — '
                         're-run: python stage_deposit.py --copy')
    return out


def remoteFileList(article):
    """Every file on the item. The endpoint paginates (10 per page by default) -- reading only the
    first page is what caused duplicate uploads, since unseen files looked absent."""
    out, page = [], 1
    while True:
        batch = call('GET', f'/account/articles/{article}/files?page={page}&page_size=100')
        if not batch:
            return out
        out.extend(batch)
        page += 1


def remoteFiles(article):
    """name -> file, keeping the newest entry when a name appears more than once."""
    return {f['name']: f for f in remoteFileList(article)}


def dedupe(article, apply=False):
    """Remove duplicate uploads of the same filename, keeping one copy (preferring a matching md5)."""
    byName = {}
    for f in remoteFileList(article):
        byName.setdefault(f['name'], []).append(f)
    local = {p.name: p for p in depositFiles()}
    removed = 0
    for name, fs in sorted(byName.items()):
        if len(fs) < 2:
            continue
        want = md5(local[name]) if name in local else None
        keep = next((f for f in fs if want and f.get('supplied_md5') == want), fs[0])
        for f in fs:
            if f['id'] == keep['id']:
                continue
            print(f"  {'deleting' if apply else 'would delete'} duplicate {name} (id {f['id']})")
            if apply:
                call('DELETE', f"/account/articles/{article}/files/{f['id']}")
            removed += 1
    print(f"{removed} duplicate(s) {'removed' if apply else 'found'}")
    return removed


ap = argparse.ArgumentParser()
ap.add_argument('--create', metavar='TITLE', help='create a new PRIVATE article and print its id')
ap.add_argument('--description', default='Build inputs for the uPULLI figure-reproduction code '
                                         '(github: BridgesLabCMU). See inputs.json for the manifest: '
                                         'logical name, size, sha256 and which figures use each file.')
ap.add_argument('--article', type=int, help='article id to upload into')
ap.add_argument('--verify', action='store_true', help='compare remote md5/size against local')
ap.add_argument('--dry-run', action='store_true')
ap.add_argument('--dedupe', action='store_true', help='delete duplicate uploads of the same filename')
args = ap.parse_args()

if args.create:
    res = call('POST', '/account/articles',
               {'title': args.create, 'description': args.description, 'defined_type': 'dataset'})
    articleId = int(str(res['location']).rstrip('/').split('/')[-1])
    print(f'created PRIVATE article {articleId}')
    print(f'  {API}/account/articles/{articleId}')
    print(f'  next: python kilthub_upload.py --article {articleId}')
    sys.exit(0)

if not args.article:
    raise SystemExit('need --create TITLE or --article ID')

files = depositFiles()
if args.dedupe:
    dedupe(args.article, apply=not args.dry_run)
    sys.exit(0)
remote = remoteFiles(args.article)

if args.verify:
    bad = ok = missing = 0
    for p in files:
        r = remote.get(p.name)
        if not r:
            print(f'  [MISSING]  {p.name}')
            missing += 1
        elif r.get('supplied_md5') and r['supplied_md5'] != md5(p):
            print(f'  [MD5 DIFF] {p.name}')
            bad += 1
        elif r['size'] != p.stat().st_size:
            print(f'  [SIZE DIFF] {p.name}  remote {r["size"]} vs local {p.stat().st_size}')
            bad += 1
        else:
            ok += 1
    print(f'\n{ok} verified, {bad} mismatched, {missing} missing')
    sys.exit(1 if (bad or missing) else 0)

total = sum(p.stat().st_size for p in files)
print(f'{len(files)} files, {total / 1e9:.2f} GB -> article {args.article}')
for p in files:
    size, digest = p.stat().st_size, md5(p)
    r = remote.get(p.name)
    if r and r.get('supplied_md5') == digest and r['size'] == size:
        print(f'  skip (already uploaded)  {p.name}')
        continue
    if args.dry_run:
        print(f'  would upload  {p.name}  ({size / 1e6:.1f} MB)')
        continue

    loc = call('POST', f'/account/articles/{args.article}/files',
               {'name': p.name, 'md5': digest, 'size': size})['location']
    info = call('GET', loc)
    parts = call('GET', info['upload_url'])['parts']
    print(f'  uploading {p.name}  ({size / 1e6:.1f} MB, {len(parts)} part(s))', flush=True)
    with open(p, 'rb') as fh:
        for part in parts:
            fh.seek(part['startOffset'])
            blob = fh.read(part['endOffset'] - part['startOffset'] + 1)
            call('PUT', f"{info['upload_url']}/{part['partNo']}", raw=blob, ctype='application/octet-stream')
    call('POST', loc)                      # complete
print('\ndone — the article is still PRIVATE (no publish call is made by this script)')
