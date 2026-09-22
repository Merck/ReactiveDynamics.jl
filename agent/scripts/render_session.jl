#!/usr/bin/env julia
#
# render_session.jl — render a captured rd-study session as a self-contained HTML transcript.
#
#   julia --project=. agent/scripts/render_session.jl SESSION.json [--out=SESSION.html]
#
# The point of the artifact: the four phases produce documents, but what SELLS the workflow is the
# shape of the conversation that produced them — the questions asked, the defaults proposed, the
# checks that ran, the number that came out. A transcript shows that; a directory of Markdown files
# does not. So this renders one session into a single HTML file with no external assets (fonts
# degrade to system stacks, figures are inlined), which means it can be mailed, committed, or opened
# from a USB stick five years from now and still look like itself.
#
# Input schema — one JSON object:
#
#   { "title": …, "subtitle": …, "headline": …, "headline_note": …,
#     "meta":  { "<label>": "<value>", … },          # rendered as the provenance strip
#     "turns": [ { "role": "user"|"assistant"|"note",
#                  "text": "<markdown-ish>",          # see `inline` for what is supported
#                  "phase": "FRAME"|…,                # optional; starts a new phase band
#                  "tools": [ { "cmd": …, "out": … } ],
#                  "files": [ "path (written)", … ],
#                  "figure": "<path to an .svg, inlined>" } ] }
#
# Deliberately NOT a Markdown renderer: the input is authored for this file, so it supports only the
# handful of inline constructs a transcript needs (`code`, **bold**, *italic*, bullet lists, tables,
# fenced blocks). A full Markdown dependency would buy formatting nobody needs and a dep nobody wants.

using JSON

const ROLE_LABEL = Dict("user" => "Operator", "assistant" => "Claude", "note" => "")

# ── escaping and the tiny inline grammar ──────────────────────────────────────────────────────────

esc(s) = replace(
    string(s),
    "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;"
)

# Inline markup, applied AFTER escaping so authored text can never inject markup. Order matters:
# code spans first, so `**` inside backticks is left alone.
function inline(s)
    out = esc(s)
    out = replace(out, r"`([^`]+)`" => s"<code>\1</code>")
    out = replace(out, r"\*\*([^*]+)\*\*" => s"<strong>\1</strong>")
    out = replace(out, r"(?<![\*\w])\*([^*\n]+)\*(?!\w)" => s"<em>\1</em>")
    out = replace(out, r"\$(\d[\d,.]*)" => s"$\1")   # keep money literal, not a math delimiter
    return out
end

# ── block layout ──────────────────────────────────────────────────────────────────────────────────

# Render one turn's body: paragraphs, bullet lists, pipe tables and fenced code, in source order.
# A transcript is mostly prose with occasional structure, so a single pass over lines is enough.
function body_html(text)
    io = IOBuffer()
    lines = split(rstrip(text), '\n')
    i = 1
    para = String[]

    flush_para!() = begin
        isempty(para) || println(io, "<p>", inline(join(para, " ")), "</p>")
        empty!(para)
    end

    while i <= length(lines)
        ln = lines[i]
        stripped = strip(ln)

        if startswith(stripped, "```")                      # fenced block
            flush_para!()
            fence = String[]
            i += 1
            while i <= length(lines) && !startswith(strip(lines[i]), "```")
                push!(fence, lines[i]); i += 1
            end
            println(io, "<pre class=\"block\">", esc(join(fence, "\n")), "</pre>")
            i += 1

        elseif startswith(stripped, "- ")                   # bullet list
            flush_para!()
            println(io, "<ul>")
            while i <= length(lines) && startswith(strip(lines[i]), "- ")
                println(io, "<li>", inline(strip(lines[i])[3:end]), "</li>")
                i += 1
            end
            println(io, "</ul>")

        elseif startswith(stripped, "> ")                   # blockquote — the pull-quote form
            flush_para!()
            quote_lines = String[]
            while i <= length(lines) && startswith(strip(lines[i]), "> ")
                push!(quote_lines, strip(lines[i])[3:end]); i += 1
            end
            println(io, "<blockquote>", inline(join(quote_lines, " ")), "</blockquote>")

        elseif startswith(stripped, "|") && endswith(stripped, "|")   # pipe table
            flush_para!()
            rows = String[]
            while i <= length(lines) && startswith(strip(lines[i]), "|")
                push!(rows, strip(lines[i])); i += 1
            end
            cells(r) = strip.(split(strip(r, '|'), '|'))
            println(io, "<table>")
            for (k, r) in enumerate(rows)
                all(c -> occursin(r"^:?-{2,}:?$", c), cells(r)) && continue   # separator row
                tag = k == 1 ? "th" : "td"
                println(io, "<tr>", join("<$tag>" * inline(c) * "</$tag>" for c in cells(r)), "</tr>")
            end
            println(io, "</table>")

        elseif isempty(stripped)
            flush_para!()
            i += 1
        else
            push!(para, stripped)
            i += 1
        end
    end
    flush_para!()
    return String(take!(io))
end

# A tool block: collapsed by default. The commands ARE the evidence — a reader who doubts a number
# opens this and re-runs it — but expanded by default they would drown the conversation.
function tools_html(tools, files)
    (isempty(tools) && isempty(files)) && return ""
    io = IOBuffer()
    n = length(tools)
    label = n == 0 ? "files written" : "$n command$(n == 1 ? "" : "s")"
    println(io, "<details class=\"tools\"><summary><span class=\"tk\">run</span> $label</summary>")
    for t in tools
        println(io, "<div class=\"cmd\"><span class=\"p\">\$</span> ", esc(get(t, "cmd", "")), "</div>")
        out = get(t, "out", "")
        isempty(out) || println(io, "<div class=\"res\">→ ", esc(out), "</div>")
    end
    isempty(files) || println(
        io, "<div class=\"files\">",
        join(("<span>" * esc(f) * "</span>" for f in files), ""), "</div>"
    )
    println(io, "</details>")
    return String(take!(io))
end

# Inline an SVG so the page stays self-contained. A missing figure is reported in place rather than
# silently dropped — a demo with a hole in it should look like one.
function figure_html(path)
    isempty(path) && return ""
    full = isabspath(path) ? path : joinpath(pwd(), path)
    isfile(full) || return "<div class=\"missing\">figure not found: $(esc(path))</div>"
    svg = read(full, String)
    svg = replace(svg, r"<\?xml[^>]*\?>" => "", r"<!DOCTYPE[^>]*>" => "")
    return "<figure class=\"fig\">$svg</figure>"
end

const CSS = """
:root {
  --ink: #0A0A0A; --paper: #FFFFFF; --inset: #F4F4F2; --accent: #0A8A84;
  --accent-d: #076A65; --accent-lt: #3FB5AF; --amber: #B5740A; --rose: #E5243B;
  --grey: #6A6A68; --grey-lt: #9A9A98; --hair: #C7C7C4;
  --display: "Inter", "Helvetica Neue", Helvetica, Arial, system-ui, sans-serif;
  --mono: "JetBrains Mono", "JuliaMono", "SF Mono", ui-monospace, monospace;
}
* { box-sizing: border-box; }
body {
  margin: 0; background: var(--inset); color: var(--ink);
  font: 16px/1.62 var(--display); -webkit-font-smoothing: antialiased;
}
.wrap { max-width: 50rem; margin: 0 auto; padding: 0 1.25rem 5rem; }

/* ── masthead ── */
header.mast { background: var(--ink); color: #F4F4F2; padding: 3rem 0 2.5rem; margin-bottom: 2.5rem; }
header.mast .wrap { padding-bottom: 0; }
.brand { font-size: .8rem; letter-spacing: .14em; text-transform: uppercase; color: var(--accent-lt); margin: 0 0 1.25rem; }
.brand b { font-weight: 600; color: #F4F4F2; }
h1 { font-size: 2.1rem; line-height: 1.18; letter-spacing: -.02em; margin: 0 0 .6rem; font-weight: 650; }
.sub { color: #9A9A98; margin: 0; font-size: 1.02rem; }
.headline {
  margin: 2rem 0 0; padding: 1.15rem 1.35rem; background: #161616;
  border-left: 3px solid var(--accent-lt); border-radius: 2px;
}
.headline .big { font-size: 1.7rem; font-weight: 650; letter-spacing: -.01em; color: #fff; }
.headline .note { display: block; margin-top: .4rem; color: #9A9A98; font-size: .93rem; }
.meta { display: flex; flex-wrap: wrap; gap: .4rem 1.75rem; margin: 1.5rem 0 0; font-size: .82rem; color: #9A9A98; }
.meta b { color: #C7C7C4; font-weight: 500; }
.meta code { font-family: var(--mono); font-size: .93em; color: var(--accent-lt); }

/* ── phase bands ── */
.phase { display: flex; align-items: center; gap: .85rem; margin: 2.75rem 0 1.5rem; }
.phase::after { content: ""; flex: 1; height: 1px; background: var(--hair); }
.phase span {
  font-size: .72rem; letter-spacing: .13em; text-transform: uppercase; font-weight: 650;
  color: var(--accent-d); background: var(--paper); border: 1px solid var(--hair);
  padding: .3rem .6rem; border-radius: 2px; white-space: nowrap;
}

/* ── turns ── */
.turn { margin: 0 0 1.5rem; }
.who { font-size: .72rem; letter-spacing: .1em; text-transform: uppercase; color: var(--grey); margin-bottom: .45rem; }
.bubble { padding: 1.1rem 1.3rem; border-radius: 3px; background: var(--paper); border: 1px solid var(--hair); }
.turn.user { max-width: 85%; margin-left: auto; }
.turn.user .who { text-align: right; color: var(--accent-d); }
.turn.user .bubble { background: #E8F3F2; border-color: #BFDEDB; }
.turn.assistant .bubble { border-left: 3px solid var(--accent); }
.turn.note .bubble { background: transparent; border: 1px dashed var(--hair); color: var(--grey); font-size: .91rem; }
.bubble p { margin: 0 0 .85rem; } .bubble p:last-child { margin-bottom: 0; }
.bubble ul { margin: 0 0 .85rem; padding-left: 1.2rem; } .bubble li { margin-bottom: .3rem; }
.bubble code { font-family: var(--mono); font-size: .88em; background: var(--inset); padding: .1em .32em; border-radius: 2px; }
.bubble table { width: 100%; border-collapse: collapse; margin: .35rem 0 1rem; font-size: .9rem; }
.bubble th { text-align: left; font-weight: 600; border-bottom: 1.5px solid var(--ink); padding: .4rem .6rem .35rem 0; }
.bubble td { border-bottom: 1px solid var(--hair); padding: .4rem .6rem .4rem 0; vertical-align: top; }
pre.block {
  font-family: var(--mono); font-size: .82rem; line-height: 1.55; background: var(--inset);
  border: 1px solid var(--hair); border-radius: 2px; padding: .8rem .9rem; overflow-x: auto; margin: 0 0 .9rem;
}
.bubble blockquote {
  margin: 1rem 0; padding: .7rem 1rem .7rem 1.1rem; background: var(--inset);
  border-left: 3px solid var(--accent); font-size: 1.02rem; line-height: 1.55;
}

/* ── tool blocks ── */
details.tools { margin: .9rem 0 0; border-top: 1px dashed var(--hair); padding-top: .75rem; }
details.tools summary {
  cursor: pointer; font-size: .78rem; color: var(--grey); letter-spacing: .02em; list-style: none;
}
details.tools summary::-webkit-details-marker { display: none; }
details.tools summary::before { content: "▸ "; color: var(--grey-lt); }
details.tools[open] summary::before { content: "▾ "; }
.tk {
  font-family: var(--mono); font-size: .72rem; text-transform: uppercase; letter-spacing: .08em;
  color: var(--accent-d); border: 1px solid #BFDEDB; background: #E8F3F2;
  padding: .08rem .35rem; border-radius: 2px; margin-right: .4rem;
}
.cmd { font-family: var(--mono); font-size: .79rem; line-height: 1.6; margin-top: .6rem; word-break: break-all; }
.cmd .p { color: var(--accent); font-weight: 600; margin-right: .35rem; }
.res { font-family: var(--mono); font-size: .79rem; color: var(--grey); padding-left: 1rem; word-break: break-word; }
.files { margin-top: .7rem; display: flex; flex-wrap: wrap; gap: .35rem; }
.files span {
  font-family: var(--mono); font-size: .72rem; color: var(--amber);
  border: 1px solid #E4D2B0; background: #FBF5EA; padding: .12rem .4rem; border-radius: 2px;
}

/* ── figures ── */
figure.fig { margin: 1rem 0 .3rem; padding: .9rem; background: var(--paper); border: 1px solid var(--hair); border-radius: 2px; overflow-x: auto; }
figure.fig svg { max-width: 100%; height: auto; display: block; margin: 0 auto; }
.missing { font-family: var(--mono); font-size: .78rem; color: var(--rose); padding: .6rem; border: 1px dashed var(--rose); border-radius: 2px; }

footer.foot { margin-top: 3.5rem; padding-top: 1.25rem; border-top: 1px solid var(--hair); font-size: .82rem; color: var(--grey); }
footer.foot code { font-family: var(--mono); font-size: .93em; }

@media print {
  body { background: #fff; } header.mast { background: #fff; color: #000; border-bottom: 2px solid #000; }
  h1, .headline .big { color: #000; } .headline { background: var(--inset); }
  .sub, .meta, .meta b { color: #444; } .brand { color: var(--accent-d); } .brand b { color: #000; }
  details.tools { display: block; } details.tools > *:not(summary) { display: block !important; }
  .turn { break-inside: avoid; }
}
"""

function render(doc, jsonpath)
    io = IOBuffer()
    println(io, "<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"utf-8\">")
    println(io, "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">")
    println(io, "<title>", esc(get(doc, "title", "rd-study session")), "</title>")
    println(io, "<style>", CSS, "</style></head><body>")

    # masthead
    println(io, "<header class=\"mast\"><div class=\"wrap\">")
    println(io, "<p class=\"brand\"><b>ReactiveDynamics.jl</b> &nbsp;·&nbsp; rd-study session</p>")
    println(io, "<h1>", inline(get(doc, "title", "")), "</h1>")
    sub = get(doc, "subtitle", "")
    isempty(sub) || println(io, "<p class=\"sub\">", inline(sub), "</p>")
    hl = get(doc, "headline", "")
    if !isempty(hl)
        println(io, "<div class=\"headline\"><span class=\"big\">", inline(hl), "</span>")
        hn = get(doc, "headline_note", "")
        isempty(hn) || println(io, "<span class=\"note\">", inline(hn), "</span>")
        println(io, "</div>")
    end
    meta = get(doc, "meta", Dict())
    if !isempty(meta)
        println(io, "<div class=\"meta\">")
        for (k, v) in meta
            println(io, "<div><b>", esc(k), "</b> ", inline(v), "</div>")
        end
        println(io, "</div>")
    end
    println(io, "</div></header>")

    # transcript
    println(io, "<main class=\"wrap\">")
    for turn in get(doc, "turns", [])
        ph = get(turn, "phase", "")
        isempty(ph) || println(io, "<div class=\"phase\"><span>", esc(ph), "</span></div>")

        role = get(turn, "role", "assistant")
        println(io, "<div class=\"turn $role\">")
        label = get(ROLE_LABEL, role, "")
        isempty(label) || println(io, "<div class=\"who\">", esc(label), "</div>")
        println(io, "<div class=\"bubble\">")
        println(io, body_html(get(turn, "text", "")))
        println(io, figure_html(get(turn, "figure", "")))
        println(io, tools_html(get(turn, "tools", []), get(turn, "files", [])))
        println(io, "</div></div>")
    end

    foot = get(doc, "footer", "")
    isempty(foot) || println(io, "<footer class=\"foot\">", body_html(foot), "</footer>")
    println(io, "</main>")
    println(io, "<!-- generated by agent/scripts/render_session.jl from ", esc(jsonpath), " -->")
    println(io, "</body></html>")
    return String(take!(io))
end

function main(argv)
    isempty(argv) && (println(stderr, "usage: render_session.jl SESSION.json [--out=SESSION.html]"); return 2)
    src = argv[1]
    isfile(src) || (println(stderr, "render_session: no such file: $src"); return 2)
    out = ""
    for a in argv[2:end]
        startswith(a, "--out=") && (out = a[(length("--out=") + 1):end])
    end
    isempty(out) && (out = replace(src, r"\.json$" => "") * ".html")

    doc = JSON.parsefile(src)
    html = render(doc, src)
    write(out, html)
    nturns = length(get(doc, "turns", []))
    ncmds = sum(length(get(t, "tools", [])) for t in get(doc, "turns", []); init = 0)
    println("render_session: wrote $out — $nturns turns, $ncmds commands, $(round(sizeof(html) / 1024, digits = 1)) KiB self-contained")
    return 0
end

exit(main(ARGS))
