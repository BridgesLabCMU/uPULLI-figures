#!/usr/bin/env bash
# Export the public-facing copy of this package.
#
# The working copy (this directory) is where you edit and re-render. This script produces the
# publishable copy at $DEST containing EXACTLY the files git would track — nothing else. Anything
# .gitignore excludes (the manuscript PDF, SUPPLEMENTAL_LEGENDS.md, inputs.local.json, build inputs,
# __pycache__, rendered figures) is never copied, so the public tree cannot leak it even by accident.
#
#   bash export_public.sh                      # -> ~/uPULLI-figures (the published repo clone)
#   bash export_public.sh /some/other/path
#
# Re-run any time. Existing files are overwritten and files no longer tracked are removed, so the
# export always matches the working copy. Edit only the working copy; treat $DEST as generated.
set -euo pipefail

SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# The published repo is BridgesLabCMU/uPULLI-figures; ~/uPULLI-figures is its clone. Exporting
# straight into it means `git status` there shows exactly what this export changed, and the
# rsync below preserves .git. (Before 2026-08 this defaulted to ~/multiphenotype-figs, which had
# no remote -- exports went nowhere while the published repo drifted ahead of this working copy.)
DEST="${1:-$HOME/uPULLI-figures}"
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
#
# Every gate below checks the EXPORTED FILE SET ($STAGE/.filelist, resolved under $DEST), never
# $DEST as a whole. $DEST is a live clone of the published repo: people render figures in it and
# Python leaves __pycache__ there. Those are local, gitignored, and none of this script's
# business -- scanning the whole directory reports them as export failures.
fail=0
exported() { sed -z "s|^|$DEST/|" "$STAGE/.filelist"; }

# Absolute paths, in EVERY text file rather than only *.py: a mount point is just as leaky in a
# README, a JSON manifest, a CSV header or an interactive HTML page. -I skips binaries.
# export_public.sh is exempt because it necessarily contains these patterns itself (this line).
leak=$(exported | xargs -0 grep -Il -e '/mnt/' -e '/home/' -e '_MONO' 2>/dev/null \
       | grep -v -e 'fetch_data.py' -e 'stage_inputs.py' -e 'export_public.sh' || true)
[ -n "$leak" ] && { echo "FAIL absolute paths in: $leak"; fail=1; }

for f in "Multiphenotype Manuscript-2.pdf" SUPPLEMENTAL_LEGENDS.md inputs.local.json; do
  grep -qzxF "$f" "$STAGE/.filelist" && { echo "FAIL excluded file exported: $f"; fail=1; }
done

tr '\0' '\n' < "$STAGE/.filelist" | grep -qE '(^|/)__pycache__/|/build/inputs/' && {
  echo "FAIL build inputs or bytecode exported"; fail=1; }
[ "$csv" -lt 80 ] && { echo "FAIL only $csv data CSVs exported"; fail=1; }
[ "$fail" -eq 0 ] && echo "gates: clean" || { echo "gates: FAILED"; exit 1; }
