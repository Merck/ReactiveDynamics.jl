#!/usr/bin/env bash
# Assemble the standalone handoff folder (and a zip) from this demo: the talk page, model code, snapshot, the
# baseline results, and the agent guide. Run from anywhere after run_twin.jl has refreshed out/.
#   demo/nocobase_pharma_twin/handoff/build_handoff.sh [dest_dir]     (default: handoff/dist)
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"; demo="$(dirname "$here")"; repo="$(cd "$demo/../.." && pwd)"
dest="${1:-$here/dist}"; pkg="$dest/pharma_twin"
# The page must be pure ASCII (special characters as entities), so it renders right whatever encoding a viewer assumes.
if LC_ALL=C grep -q $'[\x80-\xff]' "$demo/talk.html"; then echo "talk.html contains non-ASCII bytes; escape them first" >&2; exit 1; fi
rm -rf "$pkg" "$dest/pharma_twin.zip"; mkdir -p "$pkg/data" "$pkg/vendor" "$pkg/out"
wrap() { { printf '<!doctype html>\n<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover"></head><body>\n'; cat "$1"; printf '\n</body></html>\n'; } > "$2"; }
wrap "$demo/talk.html" "$pkg/talk.html"
"$here/make_pdf.sh" "$pkg/talk.pdf" >/dev/null
cp "$demo/twin.jl" "$demo/run_twin.jl" "$pkg/"
cp "$demo/data/nocobase_snapshot.json" "$pkg/data/"
cp "$repo/agent/scripts/chartkit.jl" "$pkg/vendor/"; cp "$repo/LICENSE" "$pkg/LICENSE"   # the repo's MIT license covers this package, chartkit.jl included
cp "$here/Project.toml" "$here/AGENTS.md" "$here/CLAUDE.md" "$here/HANDOFF_PROMPT.md" "$pkg/"
[ -f "$here/Manifest.toml" ] && cp "$here/Manifest.toml" "$pkg/"
for f in "$demo"/out/*.json "$demo"/out/*.svg; do cp "$f" "$pkg/out/"; done
(cd "$dest" && zip -qr pharma_twin.zip pharma_twin)
echo "$pkg"; echo "$dest/pharma_twin.zip"
