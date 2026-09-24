#!/usr/bin/env julia
#
# render_deck.jl — render a slide deck document as one self-contained HTML file.
#
#   julia --project=. agent/scripts/render_deck.jl DECK.json [--out=DECK.html]
#
# The page renderer (`render_session.jl`) makes documents to be READ: long, scrolling, every claim next
# to its evidence. A deck is for a room — one idea per sheet, the number big, the evidence one click
# away on the companion page. So this is a second layout over the SAME grammar: it includes the page
# renderer for its inline/block markup, figure inliner and stat row, which means a table, a fold or a
# money literal renders identically in both, and a figure a study regenerates changes in both at once.
#
# On screen each slide fills the viewport (scaled from a fixed 1280×720 sheet, so the layout never
# reflows), with ←/→, PageUp/PageDown and space to move. On paper each slide is exactly one page, so
# printing to PDF from any browser — or headless Chrome — yields the deck as slides.
#
# Input schema — one JSON object:
#
#   { "title": …, "brand": …, "footer": "<short text on every slide>",
#     "slides": [ { "layout": "title" | "section" | …,   # optional; anything else is a content slide
#                   "kicker": …, "title": …, "subtitle": …,
#                   "stats": [ { "value": …, "label": …, "note": … } ],
#                   "text": "<the page grammar>",
#                   "columns": [ "<grammar>", … ],        # side-by-side bodies instead of `text`
#                   "figure": "<path to an .svg>", "caption": …,   # with `text`: text left, figure right
#                   "turns": [ { "role": "user"|"assistant", "text": …, "call": "<one tool line>" } ],
#                   "notes": "<speaker notes — kept in the HTML as a comment, never shown>" } ] }

include(joinpath(@__DIR__, "render_session.jl"))

const DECK_CSS = """
:root {
  --ink: #0A0A0A; --paper: #FFFFFF; --inset: #F4F4F2; --hero: #FAF8F3; --accent: #0A8A84;
  --accent-d: #076A65; --accent-lt: #3FB5AF; --amber: #B5740A; --rose: #E5243B;
  --grey: #6A6A68; --grey-lt: #9A9A98; --hair: #C7C7C4;
  --display: "Inter", "Helvetica Neue", Helvetica, Arial, system-ui, sans-serif;
  --mono: "JetBrains Mono", "JuliaMono", "SF Mono", ui-monospace, monospace;
}
* { box-sizing: border-box; }
html { scroll-snap-type: y mandatory; }
body { margin: 0; background: #1A1A1A; font-family: var(--display); color: var(--ink); -webkit-font-smoothing: antialiased; }
.frame { height: 100vh; display: flex; align-items: center; justify-content: center; scroll-snap-align: start; overflow: hidden; }
.slide {
  width: 1280px; height: 720px; flex: none; background: var(--paper); position: relative;
  padding: 58px 72px 64px; display: flex; flex-direction: column; transform-origin: center center;
}
.slide .kicker { font-family: var(--mono); font-size: 13px; letter-spacing: .12em; text-transform: uppercase; color: var(--accent-d); margin-bottom: 10px; }
.slide h2 { font-size: 38px; line-height: 1.12; letter-spacing: -.02em; font-weight: 700; margin: 0 0 22px; max-width: 1080px; }
.slide .subtitle { font-size: 20px; color: var(--grey); margin: -10px 0 22px; max-width: 1000px; line-height: 1.4; }
.slide .content { flex: 1; min-height: 0; display: flex; gap: 44px; font-size: 20px; line-height: 1.5; }
.slide .content > .col { flex: 1; min-width: 0; }
.slide .content > .col.fig { flex: 1.25; display: flex; flex-direction: column; justify-content: center; }
.slide .content.solo-fig { justify-content: center; }
.slide .content.solo-fig > .col.fig { flex: none; width: 100%; max-width: 1060px; }
.slide figure.fig { margin: 0; }
.slide figure.fig svg { width: 100%; height: auto; max-height: 480px; display: block; }
.slide .content.solo-fig figure.fig svg { max-height: 530px; }
.slide figcaption, .slide .caption { font-size: 14px; color: var(--grey); margin-top: 8px; line-height: 1.4; }
.slide .foot { position: absolute; left: 72px; right: 72px; bottom: 24px; display: flex; justify-content: space-between; font-size: 12px; color: var(--grey-lt); font-family: var(--mono); letter-spacing: .04em; }
.slide .foot b { color: var(--ink); font-family: var(--display); letter-spacing: 0; }
.slide::before { content: ""; position: absolute; left: 0; top: 0; width: 100%; height: 5px; background: var(--accent); }

/* the grammar's blocks, resized for a room */
.prose p { margin: 0 0 14px; } .prose ul, .prose ol { margin: 0 0 14px; padding-left: 24px; } .prose li { margin-bottom: 8px; }
.prose strong { font-weight: 650; }
.prose code { font-family: var(--mono); font-size: .86em; background: var(--inset); padding: .08em .3em; border-radius: 2px; }
.prose table { width: 100%; border-collapse: collapse; font-size: 16px; margin: 4px 0 12px; }
.prose th { text-align: left; font-weight: 650; border-bottom: 2px solid var(--ink); padding: 6px 10px 6px 0; }
.prose td { border-bottom: 1px solid var(--hair); padding: 7px 10px 7px 0; vertical-align: top; }
.prose .head2, .prose .head3 { font-size: 17px; font-weight: 700; color: var(--accent-d); margin: 0 0 8px; letter-spacing: .01em; }
.prose blockquote { margin: 0 0 14px; padding: 12px 18px; background: var(--inset); border-left: 4px solid var(--accent); }
pre.block { font-family: var(--mono); font-size: 13px; line-height: 1.5; background: var(--inset); border: 1px solid var(--hair); padding: 12px 14px; margin: 0 0 12px; white-space: pre-wrap; }
details.fold { display: none; }
.stats { display: flex; gap: 22px; margin: 0 0 24px; }
.stat { flex: 1; border-top: 3px solid var(--ink); padding-top: 10px; }
.stat .label { font-size: 13px; text-transform: uppercase; letter-spacing: .08em; color: var(--grey); }
.stat .value { font-size: 42px; font-weight: 700; letter-spacing: -.02em; line-height: 1.1; margin: 4px 0; }
.stat .note { font-size: 14px; color: var(--grey); line-height: 1.35; }

/* title and section sheets */
.slide.title { background: var(--ink); color: #F4F4F2; justify-content: center; }
.slide.title .kicker { color: var(--accent-lt); }
.slide.title h2 { font-size: 60px; max-width: 1100px; }
.slide.title .subtitle { color: #C7C7C4; font-size: 24px; margin-top: 0; }
.slide.title .content { color: #C7C7C4; font-size: 17px; flex: none; }
.slide.title .foot, .slide.section .foot { color: #6A6A68; } .slide.title .foot b, .slide.section .foot b { color: #C7C7C4; }
.slide.section { background: var(--hero); justify-content: center; }
.slide.section h2 { font-size: 50px; }

/* a dialogue excerpt */
.chat { display: flex; flex-direction: column; gap: 12px; font-size: 16.5px; line-height: 1.45; }
.chat .b { max-width: 88%; padding: 12px 16px; border-radius: 3px; }
.chat .b.user { align-self: flex-end; background: var(--hero); border: 1px solid #E6DFCF; }
.chat .b.assistant { align-self: flex-start; background: var(--paper); border: 1px solid var(--hair); border-left: 4px solid var(--accent); }
.chat .who { font-family: var(--mono); font-size: 11px; letter-spacing: .1em; text-transform: uppercase; color: var(--grey); margin-bottom: 4px; }
.chat .b p { margin: 0 0 6px; } .chat .b p:last-child { margin: 0; }
.chat .call { margin-top: 8px; font-family: var(--mono); font-size: 12.5px; color: var(--accent-d); background: #F7FBFB; border: 1px solid #BFDEDB; padding: 5px 9px; border-radius: 2px; white-space: pre-wrap; }

@media print {
  @page { size: 1280px 720px; margin: 0; }
  html { scroll-snap-type: none; } body { background: #fff; }
  .frame { height: 720px; width: 1280px; display: block; break-after: page; overflow: hidden; }
  .slide { transform: none !important; }
  * { -webkit-print-color-adjust: exact; print-color-adjust: exact; }
}
"""

# Fit each 1280×720 sheet to the viewport without reflowing it, and move by whole slides.
const DECK_JS = """
(function () {
  function fit() {
    var s = Math.min(window.innerWidth / 1280, window.innerHeight / 720) * 0.97;
    document.querySelectorAll('.slide').forEach(function (el) { el.style.transform = 'scale(' + s + ')'; });
  }
  window.addEventListener('resize', fit); fit();
  var frames = Array.prototype.slice.call(document.querySelectorAll('.frame'));
  function current() {
    var y = window.scrollY + window.innerHeight / 2;
    for (var i = 0; i < frames.length; i++) if (frames[i].offsetTop + frames[i].offsetHeight > y) return i;
    return frames.length - 1;
  }
  document.addEventListener('keydown', function (e) {
    var d = ({ArrowRight: 1, ArrowDown: 1, PageDown: 1, ' ': 1, ArrowLeft: -1, ArrowUp: -1, PageUp: -1})[e.key];
    if (!d) return;
    e.preventDefault();
    var i = Math.max(0, Math.min(frames.length - 1, current() + d));
    frames[i].scrollIntoView({behavior: 'instant', block: 'start'});
  });
})();
"""

function chat_html(turns)
    io = IOBuffer()
    println(io, "<div class=\"chat\">")
    for t in turns
        role = string(get(t, "role", "assistant"))
        println(io, "<div class=\"b $role\"><div class=\"who\">", esc(get(ROLE_LABEL, role, "")), "</div>")
        println(io, body_html(string(get(t, "text", ""))))
        call = string(get(t, "call", ""))
        isempty(call) || println(io, "<div class=\"call\">", esc(call), "</div>")
        println(io, "</div>")
    end
    println(io, "</div>")
    return String(take!(io))
end

function slide_html(s, n, total, footer)
    io = IOBuffer()
    layout = string(get(s, "layout", ""))
    cls = layout in ("title", "section") ? " $layout" : ""
    notes = string(get(s, "notes", ""))
    isempty(notes) || println(io, "<!-- notes: ", replace(notes, "--" => "—"), " -->")
    println(io, "<div class=\"frame\"><section class=\"slide$cls\">")
    k = string(get(s, "kicker", ""))
    isempty(k) || println(io, "<div class=\"kicker\">", inline(k), "</div>")
    println(io, "<h2>", inline(string(get(s, "title", ""))), "</h2>")
    sub = string(get(s, "subtitle", ""))
    isempty(sub) || println(io, "<div class=\"subtitle\">", inline(sub), "</div>")
    println(io, "<div class=\"prose\">", stats_html(get(s, "stats", [])), "</div>")

    text = string(get(s, "text", ""))
    cols = get(s, "columns", [])
    fig = string(get(s, "figure", ""))
    turns = get(s, "turns", [])
    cap = string(get(s, "caption", ""))
    figcol = isempty(fig) ? "" :
        "<div class=\"col fig\">" * figure_html(fig) *
        (isempty(cap) ? "" : "<div class=\"caption\">" * inline(cap) * "</div>") * "</div>"

    solo = !isempty(fig) && isempty(text) && isempty(cols) && isempty(turns)
    println(io, "<div class=\"content", solo ? " solo-fig" : "", "\">")
    if !isempty(cols)
        for c in cols
            println(io, "<div class=\"col prose\">", body_html(string(c)), "</div>")
        end
    else
        isempty(text) || println(io, "<div class=\"col prose\">", body_html(text), "</div>")
        isempty(turns) || println(io, "<div class=\"col\">", chat_html(turns), "</div>")
    end
    println(io, figcol)
    println(io, "</div>")

    println(io, "<div class=\"foot\"><span><b>ReactiveDynamics.jl</b> · ", inline(footer), "</span><span>", n, " / ", total, "</span></div>")
    println(io, "</section></div>")
    return String(take!(io))
end

function render_deck(doc, jsonpath)
    slides = get(doc, "slides", [])
    footer = string(get(doc, "footer", get(doc, "brand", "")))
    io = IOBuffer()
    println(io, "<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"utf-8\">")
    println(io, "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">")
    println(io, "<title>", esc(get(doc, "title", "deck")), "</title>")
    println(io, "<style>", DECK_CSS, "</style></head><body>")
    for (i, s) in enumerate(slides)
        println(io, slide_html(s, i, length(slides), footer))
    end
    println(io, "<script>", DECK_JS, "</script>")
    println(io, "<!-- generated by agent/scripts/render_deck.jl from ", esc(jsonpath), " -->")
    println(io, "</body></html>")
    return String(take!(io))
end

function deck_main(argv)
    isempty(argv) && (println(stderr, "usage: render_deck.jl DECK.json [--out=DECK.html]"); return 2)
    src = argv[1]
    isfile(src) || (println(stderr, "render_deck: no such file: $src"); return 2)
    out = ""
    for a in argv[2:end]
        startswith(a, "--out=") && (out = a[(length("--out=") + 1):end])
    end
    isempty(out) && (out = replace(src, r"\.json$" => "") * ".html")
    doc = JSON.parsefile(src)
    html = render_deck(doc, src)
    write(out, html)
    println("render_deck: wrote $out — $(length(get(doc, "slides", []))) slides, $(round(sizeof(html) / 1024, digits = 1)) KiB self-contained")
    return 0
end

abspath(PROGRAM_FILE) == (@__FILE__) && exit(deck_main(ARGS))
