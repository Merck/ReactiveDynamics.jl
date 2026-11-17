#!/usr/bin/env bash
# Print talk.html to a PDF with headless Chrome: A3 portrait, so the page keeps its desktop layout
# (A4 is narrower than its 880px breakpoint) and prints cleanly scaled onto A4. Collapsed sections
# are opened, since a PDF cannot expand them.
#   demo/nocobase_pharma_twin/handoff/make_pdf.sh [out.pdf]     (default: handoff/dist/talk.pdf)
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"; demo="$(dirname "$here")"
out="${1:-$here/dist/talk.pdf}"; mkdir -p "$(dirname "$out")"
chrome="${CHROME:-/Applications/Google Chrome.app/Contents/MacOS/Google Chrome}"
tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT
python3 -I - "$demo/talk.html" "$tmp/print.html" <<'PY'
import sys
src, dst = sys.argv[1], sys.argv[2]
s = open(src, encoding="ascii").read()
s = s.replace("<details>", "<details open>")
# Chrome evaluates the width media queries against a default sheet, not the A3 page, so turn the
# narrow-screen rules off and keep the desktop layout.
for q in ("880px", "760px", "480px"):
    assert s.count(f"@media (max-width: {q})") == 1, q
    s = s.replace(f"@media (max-width: {q})", "@media (max-width: 1px)")
s = s.replace("<summary>Show the month-by-month", "<summary>The month-by-month")
s = s.replace("<summary>Show the numbers behind the charts", "<summary>The numbers behind the charts")
css = """<style>
@page { size: A3 portrait; margin: 14mm 12mm; }
html, body { -webkit-print-color-adjust: exact; print-color-adjust: exact; }
#tip { display: none !important; }
nav { display: none; }
details { padding-top: 14px; }
details > summary { list-style: none; font-size: var(--t-md); }
details > summary::-webkit-details-marker { display: none; }
.decide li, .split > div, .step, .next li, .gaps li, .dl li, tr, figure, .cm, .vs, .headline, .flow { break-inside: avoid; }
h2, h3 { break-after: avoid; }
.feed { display: block; columns: 2; column-gap: 32px; }
.feed li { grid-template-columns: 50px 70px 72px minmax(0, 1fr); gap: 8px; padding: 6px 0; break-inside: avoid; }
.feed li .d { display: none; }
.feed .feed-break { column-span: all; }
#levers details { break-inside: avoid; }
section { break-inside: auto; }
section:not(#answer) { break-before: page; }
.sec-head, details > summary { break-after: avoid; }
</style>
"""
open(dst, "w", encoding="ascii").write("<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"></head><body>\n" + s + css + "\n</body></html>\n")
PY
"$chrome" --headless=new --disable-gpu --no-pdf-header-footer --force-prefers-reduced-motion \
    --virtual-time-budget=15000 --run-all-compositor-stages-before-draw \
    --print-to-pdf="$out" "file://$tmp/print.html" 2>/dev/null
echo "$out"
