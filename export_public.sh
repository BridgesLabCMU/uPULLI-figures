#!/usr/bin/env bash
# Export the public-facing copy of this package.
#
# The working copy (this directory) is where you edit and re-render. This script produces the
# publishable copy at $DEST containing EXACTLY the files git would track — nothing else. Anything
# .gitignore excludes (the manuscript PDF, SUPPLEMENTAL_LEGENDS.md, inputs.local.json, build inputs,
# __pycache__, rendered figures) is never copied, so the public tree cannot leak it even by accident.
#
#   bash export_public.sh                      # -> ~/multiphenotype-figs
#   bash export_public.sh /some/other/path
#
# Re-run any time. Existing files are overwritten and files no longer tracked are removed, so the
# export always matches the working copy. Edit only the working copy; treat $DEST as generated.
set -euo pipefail

SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="${1:-$HOME/multiphenotype-figs}"
STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT

echo "source : $SRC"
echo "dest   : $DEST"

# Resolve the tracked-file set by letting git apply .gitignore in a throwaway index.
rsync -a --exclude='.git/' --exclude='__pycache__/' "$SRC/" "$STAGE/"
git -C "$STAGE" init -q
git -C "$STAGE" add -A
git -C "$STAGE" ls-files -z > "$STAGE/.filelist"

mkdir -p "$DEST"
# --delete keeps the export in sync when files are renamed or removed upstream; .git is preserved.
rsync -a --delete --exclude='.git/' --from0 --files-from="$STAGE/.filelist" "$STAGE/" "$DEST/"

n=$(git -C "$STAGE" ls-files | wc -l)
csv=$(git -C "$STAGE" ls-files 'fig*/data/*.csv' | wc -l)
echo "exported $n files ($csv bundled source-data CSVs)"

# ── gates: things that must never appear in the public copy ──
fail=0
leak=$(grep -rIl -e '/mnt/' -e '/home/' -e '_MONO' --include='*.py' "$DEST" 2>/dev/null \
       | grep -v -e 'fetch_data.py' -e 'stage_inputs.py' || true)
[ -n "$leak" ] && { echo "FAIL absolute paths in: $leak"; fail=1; }
for f in "Multiphenotype Manuscript-2.pdf" SUPPLEMENTAL_LEGENDS.md inputs.local.json; do
  [ -e "$DEST/$f" ] && { echo "FAIL excluded file present: $f"; fail=1; }
done
find "$DEST" -name '__pycache__' -o -path '*/build/inputs/*' | grep -q . && {
  echo "FAIL build inputs or bytecode copied"; fail=1; }
[ "$csv" -lt 80 ] && { echo "FAIL only $csv data CSVs exported"; fail=1; }
[ "$fail" -eq 0 ] && echo "gates: clean" || { echo "gates: FAILED"; exit 1; }
