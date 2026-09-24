#!/usr/bin/env julia
#
# render_session.jl — render a captured rd-study session as a self-contained HTML page.
#
#   julia --project=. agent/scripts/render_session.jl PAGE.json [--out=PAGE.html]
#
# The point of the artifact: the four phases produce documents, but what SELLS the workflow is the
# shape of the conversation that produced them — the questions asked, the defaults proposed, the
# checks that ran, the number that came out. A transcript shows that; a directory of Markdown files
# does not. So this renders one session into a single HTML file with no external assets (fonts
# degrade to system stacks, figures are inlined), which means it can be mailed, committed, or opened
# from a USB stick five years from now and still look like itself.
#
# A transcript is the spine, not the whole page: `sections` puts authored sheets BEFORE the dialogue,
# which is what turns a captured session into a capability overview that happens to end in a
# demonstration. `agent/overview.json` is that shape; `agent/sessions/*.json` are bare transcripts.
#
# Input schema — one JSON object:
#
#   { "title": …, "subtitle": …, "brand": …, "headline": …, "headline_note": …,
#     "compact": true,                               # on paper, let sheets flow instead of each
#                                                    #   starting a new page — for short documents
#     "meta":  { "<label>": "<value>", … },          # rendered as the provenance strip; write it as
#              [ { "<label>": "<value>" }, … ],       #   a list instead where the order matters
#     "hero":  { "stats":   [ { "value": …, "label": …, "note": … } ],   # every key optional
#                "text":    "<markdown-ish>",
#                "figures": [ { "svg": …, "caption": … } ] },
#     "sections": [ { "title": …, "text": … | "file": …, "stats": [ … ], "figures": [ … ] } ],
#     "turns": [ { "role": "user"|"assistant"|"note",
#                  "text": "<markdown-ish>",          # see the grammar note below
#                  "phase": "FRAME"|…,                # optional; starts a new phase band
#                  "stats": [ … ],                    # optional; a callout row above the text
#                  "tools": [ { "cmd": …, "out": …,     # terminal text: line breaks are preserved
#                               "exit": 0|1|2,         #   optional; rendered as an ok/findings/error badge
#                               "json": … } ],         #   optional; a result excerpt, shown as a code block
#                  "models": [ { "file": "<a .rdj.json>",  # the formalized model, as an expander: counts in
#                                "label": …,              #   the summary, the net drawing, then the document
#                                "svg": "<net .svg>",     #   itself — read at render time, never pasted — or,
#                                "base": "<a .rdj.json>" } ],  # given `base`, only its diff against that arm
#                  "open": true,                      # optional; leave the tool block expanded
#                  "files": [ "path (written)", … ],
#                  "figures": [ { "svg": …, "caption": … } ],
#                  "figure": "<path to an .svg>" } ],  # one-element alias for `figures`, no caption
#     "turns": { "file": "<another page>.json",         # …or import another document's transcript,
#                "pick": [1, 2, 6, 7],                   #   optionally an excerpt of it (1-based)
#                "elisions": { "6": "…", "end": "…" } },  #   what each cut run contained, keyed by the
#                                                        #   pick it precedes ("end" for the tail)
#     "appendix": [ { "title": …, "text": "<markdown-ish>" | "file": "<path to a .md>",
#                     "stats": [ … ], "figures": [ … ] } ],
#     "footer": "<markdown-ish>" }
#
# `hero` is the skim card between the masthead and the first turn, laid out stats → text → figures: a
# reader who never opens a turn must be able to state the recommendation from it alone. `sections` and
# `appendix` render identically — a band heavier than a phase over a full-width sheet, never collapsed —
# because both are authored divisions rather than dialogue; they differ only in where they sit.
#
# Any node that takes `text` takes `file` instead: a path to a Markdown deliverable, imported at render
# time (see `md_import`). Prefer it over pasting the document into the JSON, because a report
# transcribed by hand is a report that can quietly disagree with the one the study shipped. `turns`
# imports the same way, from another page document — so a page can present a captured session it does
# not own, and there is still exactly one copy of that session under version control.
#
# The block grammar in every `text`: paragraphs, `- ` bullets, `1. ` ordered lists, `## `/`### `
# sub-heads, `| pipe |` tables, ``` fences, `> ` blockquotes, `!figure path.svg | caption` (consecutive
# lines share a row), and an explicit collapsible — a line `+++ summary text` opens one, a bare `+++`
# closes it, and the body between goes back through the grammar, so collapsibles nest. Inline: `code`,
# **bold**, *italic*, money literals, and [links](https://…) — absolute http(s) only, because a
# relative link 404s the moment the page is mailed, so anything else stays plain text.
#
# Deliberately NOT a Markdown renderer: the input is authored for this file, so it supports only the
# handful of constructs a transcript and its appendix need. A full Markdown dependency would buy
# formatting nobody needs and a dep nobody wants. Raw HTML is never passed through — everything is
# escaped BEFORE markup is applied, which is why the collapsible is a fence and not a `<details>`.

using JSON

const ROLE_LABEL = Dict("user" => "Operator", "assistant" => "Claude", "note" => "")

# ── escaping and the tiny inline grammar ──────────────────────────────────────────────────────────

esc(s) = replace(
    string(s),
    "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;"
)

# Inline markup, applied AFTER escaping so authored text can never inject markup. Order matters:
# code spans first, so `**` inside backticks is left alone, then links, so a `*` in a URL survives the
# emphasis pass. The link href needs no separate attribute escape — `esc` already ran over it, and the
# pattern only admits absolute http(s), so a relative target degrades to the literal `[text](…)`
# rather than to a link that 404s the moment the page leaves this machine.
function inline(s)
    out = esc(s)
    out = replace(out, r"`([^`]+)`" => s"<code>\1</code>")
    out = replace(
        out, r"\[([^\]\n]+)\]\((https?://[^\s)]+)\)" =>
            s"<a href=\"\2\" target=\"_blank\" rel=\"noopener\">\1</a>"
    )
    out = replace(out, r"\*\*([^*]+)\*\*" => s"<strong>\1</strong>")
    out = replace(out, r"(?<![\*\w])\*([^*\n]+)\*(?!\w)" => s"<em>\1</em>")
    out = replace(out, r"\$(\d[\d,.]*)" => s"$\1")   # keep money literal, not a math delimiter
    return out
end

# ── block layout ──────────────────────────────────────────────────────────────────────────────────

# Render one body — a turn, a hero blurb, an appendix entry — in source order. A transcript is mostly
# prose with occasional structure, so a single pass over lines is enough; the one construct that needs
# more is the `+++` collapsible, which re-enters this function on its own contents and so nests.
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

        elseif startswith(stripped, "+++")                  # explicit collapsible
            flush_para!()
            summary = strip(stripped[4:end])
            inner = String[]
            depth = 1
            i += 1
            while i <= length(lines)                        # count depth so folds can nest
                s = strip(lines[i])
                if s == "+++"
                    depth -= 1
                    depth == 0 && break
                elseif startswith(s, "+++")
                    depth += 1
                end
                push!(inner, lines[i]); i += 1
            end
            i += 1
            println(io, "<details class=\"fold\"><summary>", inline(isempty(summary) ? "details" : summary), "</summary>")
            println(io, "<div class=\"fold-body\">", body_html(join(inner, "\n")), "</div></details>")

        elseif startswith(stripped, "### ") || startswith(stripped, "## ")   # sub-heads
            flush_para!()
            deep = startswith(stripped, "### ")
            tag, cls = deep ? ("h4", "head3") : ("h3", "head2")
            println(io, "<$tag class=\"$cls\">", inline(strip(stripped[(deep ? 5 : 4):end])), "</$tag>")
            i += 1

        elseif startswith(stripped, "!figure ")             # an SVG from inside the prose grammar
            flush_para!()
            figs = String[]
            while i <= length(lines) && startswith(strip(lines[i]), "!figure ")
                spec = strip(strip(lines[i])[9:end])
                path, cap = occursin('|', spec) ? strip.(split(spec, '|', limit = 2)) : (spec, "")
                push!(figs, figure_html(path, cap))
                i += 1
            end
            println(io, "<div class=\"figs\">", join(figs, "\n"), "</div>")

        elseif startswith(stripped, "- ")                   # bullet list
            flush_para!()
            println(io, "<ul>")
            while i <= length(lines) && startswith(strip(lines[i]), "- ")
                println(io, "<li>", inline(strip(lines[i])[3:end]), "</li>")
                i += 1
            end
            println(io, "</ul>")

        elseif occursin(r"^\d+[.)] ", stripped)             # ordered list — a decision rule, in order
            flush_para!()
            println(io, "<ol>")
            while i <= length(lines) && occursin(r"^\d+[.)] ", strip(lines[i]))
                item = replace(strip(lines[i]), r"^\d+[.)] " => "")
                println(io, "<li>", inline(item), "</li>")
                i += 1
            end
            println(io, "</ol>")

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

# ── importing a Markdown deliverable ──────────────────────────────────────────────────────────────

# `file` on a text-bearing node publishes a real deliverable — a study's FINDINGS.md — instead of a
# copy of it pasted into the JSON, for the same reason the tool blocks carry real commands: a report
# transcribed by hand is a report that can disagree with the one the study shipped. The conversion is
# deliberately narrow, because the block grammar above already covers what these documents use. Five
# things differ:
#
#   `<details>`/`<summary>`  → the `+++` fence, since raw HTML is escaped rather than passed through
#   `![caption](fig.svg)`    → `!figure`, resolved against the MARKDOWN FILE's own directory, which is
#                              what a relative path in a Markdown document means
#   `[text](relative)`       → `text`, because only absolute http(s) links survive being mailed
#   `#`/`####` heads         → the `##`/`###` the grammar has; a leading H1 is dropped, as the band
#                              above the entry already carries the title
#   `\$`                     → `$`, an escape a Markdown document needs so a pair of dollars is not
#                              read as math, and this grammar does not: it keeps money literal already
#
# Anything else — tables, fences, bullets, quotes, emphasis — is already common to both grammars and
# passes through untouched.
function md_import(path)
    full = isabspath(path) ? path : joinpath(pwd(), path)
    isfile(full) || return "*(deliverable not found: `$(path)`)*"
    dir = dirname(abspath(full))
    out = String[]
    pending = false     # a `<details>` whose `<summary>` is on the next line
    body = false        # anything non-blank seen yet, so the leading H1 can be told from a later one
    fenced = false

    open_fold(summary) = (push!(out, "", "+++ " * summary); pending = false)

    for raw in eachline(full)
        ln = rstrip(raw)
        s = strip(ln)

        if startswith(s, "```")                  # code passes through verbatim, markers and all
            fenced = !fenced
            push!(out, ln)
            continue
        elseif fenced
            push!(out, ln)
            continue
        end

        if occursin("\\\$", ln)              # outside fences only: in a fence the backslash is code
            ln = replace(ln, "\\\$" => "\$")
            s = strip(ln)
        end

        # a `<details>` still waiting for its summary, interrupted by content: open it unlabelled
        pending && !isempty(s) && !startswith(s, "<summary>") && open_fold("")

        if startswith(s, "<details")
            m = match(r"<summary>(.*?)</summary>", s)
            m === nothing ? (pending = true) : open_fold(strip(m[1]))
        elseif startswith(s, "<summary>")
            m = match(r"^<summary>(.*?)</summary>", s)
            open_fold(m === nothing ? "" : strip(m[1]))
        elseif startswith(s, "</details>")
            push!(out, "", "+++")
        elseif (m = match(r"^!\[(.*?)\]\((\S+?)\)$", s)) !== nothing
            src = String(m[2])
            local_src = startswith(src, "http") || isabspath(src) ? src : joinpath(dir, src)
            push!(out, "!figure " * local_src * (isempty(m[1]) ? "" : " | " * m[1]))
        elseif (m = match(r"^(#+)\s*(.*)$", s)) !== nothing
            level = length(m[1])
            level == 1 && !body && continue      # the title is already on the band above
            push!(out, (level >= 3 ? "### " : "## ") * strip(m[2]))
        else
            push!(out, plain_local_links(ln))
        end
        isempty(s) || (body = true)
    end
    pending && open_fold("")
    return join(out, "\n")
end

# A relative Markdown link degrades to its own text. The inline grammar would leave it as the literal
# `[text](path)`, which reads as noise in a published deliverable; dropping the target is honest about
# the fact that the link cannot be followed from a mailed page.
function plain_local_links(s)
    pat = r"(?<!!)\[([^\]]*)\]\((\S+?)\)"
    return replace(
        s, pat => function (whole)
            m = match(pat, whole)
            return startswith(m[2], "http://") || startswith(m[2], "https://") ? whole : m[1]
        end
    )
end

# The text of a node: authored inline, or imported from the Markdown file it publishes. Paths resolve
# against the working directory, the same rule the figures use, so a session renders from the repo root.
node_text(node) = let f = string(get(node, "file", ""))
    isempty(f) ? string(get(node, "text", "")) : md_import(f)
end

# A tool block: collapsed by default. The commands ARE the evidence — a reader who doubts a number
# opens this and re-runs it — but expanded by default they would drown the conversation. `open` on a
# turn overrides that, for the one or two turns where the commands are the claim and not a footnote.
function tools_html(tools, files, expanded = false)
    (isempty(tools) && isempty(files)) && return ""
    io = IOBuffer()
    n = length(tools)
    label = n == 0 ? "files written" : "$n command$(n == 1 ? "" : "s")"
    println(io, "<details class=\"tools\"", expanded ? " open" : "", "><summary><span class=\"tk\">run</span> $label</summary>")
    for t in tools
        code = get(t, "exit", nothing)
        badge = code === nothing ? "" :
            " <span class=\"exit x$(Int(code))\">exit $(Int(code)) · $(get(EXIT_LABEL, Int(code), "?"))</span>"
        println(io, "<div class=\"cmd\"><span class=\"p\">\$</span> ", esc(get(t, "cmd", "")), badge, "</div>")
        out = get(t, "out", "")
        isempty(out) || println(io, "<div class=\"res\">→ ", esc(out), "</div>")
        js = get(t, "json", "")
        isempty(js) || println(io, "<pre class=\"block json\">", esc(js isa AbstractString ? js : JSON.json(js, 2)), "</pre>")
    end
    isempty(files) || println(
        io, "<div class=\"files\">",
        join(("<span>" * esc(f) * "</span>" for f in files), ""), "</div>"
    )
    println(io, "</details>")
    return String(take!(io))
end

const EXIT_LABEL = Dict(0 => "ok", 1 => "findings", 2 => "error")

# The formalized model, as an expander. The point of showing it is that it is DATA: a reader who opens
# it sees a document with a counted shape, not a program. So the summary carries the counts, the body
# the net drawing and then the document itself, read from disk at render time for the same reason the
# appendix imports FINDINGS — a pasted copy is a copy that can disagree with the one that ran. An arm is
# shown as its diff against `base`, because the size of that diff IS the claim: one lever, a few lines.
function models_html(models)
    isempty(models) && return ""
    io = IOBuffer()
    for m in models
        f = string(get(m, "file", ""))
        isfile(f) || (println(io, "<div class=\"missing\">model not found: $(esc(f))</div>"); continue)
        d = JSON.parsefile(f)
        n(k) = length(get(d, k, []))
        counts = join(
            (
                "$(n(k)) $(n(k) == 1 ? chop(k) : k)" for k in ("places", "transitions", "arcs", "rules")
                    if n(k) > 0
            ), " · "
        )
        label = string(get(m, "label", basename(f)))
        base = string(get(m, "base", ""))
        println(
            io, "<details class=\"model\"><summary><span class=\"tk\">model</span> ", inline(label),
            " <span class=\"mc\">", esc(counts), "</span></summary><div class=\"model-body\">"
        )
        svg = string(get(m, "svg", ""))
        isempty(svg) || println(io, figure_html(svg))
        if !isempty(base) && isfile(base)
            diff = read(ignorestatus(`diff -u -L $(basename(base)) -L $(basename(f)) $base $f`), String)
            lines = split(diff, '\n')
            body = join(
                (
                    startswith(l, "+") && !startswith(l, "+++") ? "<ins>$(esc(l))</ins>" :
                        startswith(l, "-") && !startswith(l, "---") ? "<del>$(esc(l))</del>" : esc(l)
                    for l in lines
                ), "\n"
            )
            println(io, "<div class=\"mlabel\">the whole difference from <code>", esc(basename(base)), "</code></div>")
            println(io, "<pre class=\"block diff\">", body, "</pre>")
        else
            println(io, "<div class=\"mlabel\"><code>", esc(f), "</code> — the complete document, as it ran</div>")
            println(io, "<pre class=\"block json model-src\">", esc(rstrip(read(f, String))), "</pre>")
        end
        println(io, "</div></details>")
    end
    return String(take!(io))
end

# Inline an SVG so the page stays self-contained. A missing figure is reported in place rather than
# silently dropped — a demo with a hole in it should look like one.
function figure_html(path, caption = "")
    isempty(path) && return ""
    full = isabspath(path) ? path : joinpath(pwd(), path)
    isfile(full) || return "<div class=\"missing\">figure not found: $(esc(path))</div>"
    svg = read(full, String)
    svg = replace(svg, r"<\?xml[^>]*\?>" => "", r"<!DOCTYPE[^>]*>" => "")
    cap = isempty(caption) ? "" : "<figcaption>$(inline(caption))</figcaption>"
    return "<figure class=\"fig\">$svg$cap</figure>"
end

# A figure row. `figures` is the list form; `figure` is the one-element string alias that predates it
# and that the captured sessions still use. Two fit side by side when the measure allows.
function figures_html(figures, alias = "")
    items = [(string(get(f, "svg", "")), string(get(f, "caption", ""))) for f in figures]
    isempty(items) && !isempty(alias) && push!(items, (alias, ""))
    isempty(items) && return ""
    io = IOBuffer()
    println(io, "<div class=\"figs\">")
    for (path, cap) in items
        println(io, figure_html(path, cap))
    end
    println(io, "</div>")
    return String(take!(io))
end

# A row of big-number callouts. The value is the only thing a skimming reader is guaranteed to read,
# so it carries all the weight and the label and note stay deliberately quiet. Three fit across the
# measure; narrower than that they stack.
function stats_html(stats)
    isempty(stats) && return ""
    io = IOBuffer()
    println(io, "<div class=\"stats\">")
    for s in stats
        println(io, "<div class=\"stat\">")
        for (key, cls) in (("label", "label"), ("value", "value"), ("note", "note"))
            v = string(get(s, key, ""))
            isempty(v) || println(io, "<div class=\"$cls\">", inline(v), "</div>")
        end
        println(io, "</div>")
    end
    println(io, "</div>")
    return String(take!(io))
end

# The hero: the skim path. A reader who never opens a turn should still be able to state the
# recommendation, so the card gets the numbers, the one-paragraph why and the charts, on paper distinct
# from both the masthead above it and the bubbles below.
function hero_html(hero)
    isempty(hero) && return ""
    io = IOBuffer()
    println(io, "<section class=\"hero prose\">")
    println(io, stats_html(get(hero, "stats", [])))
    txt = node_text(hero)
    isempty(txt) || println(io, "<div class=\"hero-text\">", body_html(txt), "</div>")
    println(io, figures_html(get(hero, "figures", [])))
    println(io, "</section>")
    return String(take!(io))
end

# Sheets are the page's authored divisions — the overview sections before the dialogue, the
# deliverables after it. A band heavier than a phase, a full-width sheet, the same block grammar as a
# turn — and never collapsed, because a folded deliverable reads as an afterthought nobody was meant
# to open.
function sheets_html(entries, default_title = "Appendix")
    io = IOBuffer()
    for e in entries
        println(io, "<div class=\"phase band\"><span>", esc(string(get(e, "title", default_title))), "</span></div>")
        println(io, "<section class=\"appx prose\">")
        println(io, stats_html(get(e, "stats", [])))
        println(io, body_html(node_text(e)))
        println(io, figures_html(get(e, "figures", [])))
        println(io, "</section>")
    end
    return String(take!(io))
end

# `turns` is normally the transcript itself. It may instead name another page document, whose turns are
# spliced in — the `file` idiom the appendix already uses, for the same reason: a page that presents a
# captured session should not carry a second copy of it that can drift from the one being maintained.
#
# `pick` takes an excerpt instead of the whole thing, for a page that shows what the workflow is LIKE
# without reprinting an hour of it. An excerpt has to read as one: every cut run leaves a marker carrying
# its own count, so the page can never imply the session was four exchanges long, and `elisions` lets the
# author say what was in the cut — the count is computed, only the description is written by hand.
function resolve_turns(t)
    t isa AbstractVector && return t
    t isa AbstractDict || return Any[]
    f = string(get(t, "file", ""))
    isempty(f) && return Any[]
    isfile(f) || error("render_session: turns file not found: $f")
    turns = get(JSON.parsefile(f), "turns", Any[])
    pick = get(t, "pick", nothing)
    pick === nothing && return turns

    notes = get(t, "elisions", Dict())
    gap(n, key) = Dict(
        "role" => "gap", "text" => string(n, n == 1 ? " exchange omitted" : " exchanges omitted"),
        "note" => string(get(notes, key, "")),
    )
    out, prev = Any[], 0
    for i in pick
        idx = Int(i)
        1 <= idx <= length(turns) ||
            error("render_session: turns pick $idx is outside 1:$(length(turns)) in $f")
        idx > prev && idx > prev + 1 && push!(out, gap(idx - prev - 1, string(idx)))
        push!(out, turns[idx])
        prev = idx
    end
    prev < length(turns) && push!(out, gap(length(turns) - prev, "end"))
    return out
end

const CSS = """
:root {
  --ink: #0A0A0A; --paper: #FFFFFF; --inset: #F4F4F2; --hero: #FAF8F3; --accent: #0A8A84;
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

/* ── hero: the skim card ── */
.hero {
  background: var(--hero); border: 1px solid var(--hair); border-top: 3px solid var(--ink);
  border-radius: 2px; padding: 2.1rem 1.9rem 1.9rem; margin: 0 0 .5rem;
}
.hero-text { margin-top: 1.5rem; font-size: 1.04rem; }
.hero-text p, .hero-text ul, .hero-text ol { max-width: 38rem; }
.hero .figs { margin-top: 1.6rem; }

/* ── stat cards ── */
.stats { display: flex; flex-wrap: wrap; gap: .75rem; }
.stat {
  flex: 1 1 8.5rem; min-width: 8rem; background: var(--paper); border: 1px solid var(--hair);
  border-top: 2px solid var(--accent); border-radius: 2px; padding: .9rem 1rem 1rem;
}
.stat .label { font-size: .67rem; letter-spacing: .13em; text-transform: uppercase; color: var(--grey); }
.stat .value { font-size: 2.45rem; font-weight: 650; letter-spacing: -.035em; line-height: 1.08; margin-top: .15rem; }
.stat .note { font-size: .78rem; line-height: 1.45; color: var(--grey); margin-top: .45rem; }
.bubble .stats, .appx .stats { margin-bottom: 1.15rem; }

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

/* An elision, drawn as a seam across the column rather than a bubble: it is not part of the dialogue,
   it is the place where the dialogue was cut. */
.gap { display: flex; align-items: baseline; gap: .8rem; margin: 1.4rem 0 1.9rem; color: var(--grey-lt); }
.gap::before, .gap::after { content: ""; flex: 1; height: 1px; background: var(--hair); }
.gap b { font-family: var(--mono); font-size: .7rem; letter-spacing: .08em; text-transform: uppercase; white-space: nowrap; }
.gap i { font-style: normal; font-size: .8rem; color: var(--grey); max-width: 60%; }

/* ── prose: shared by turn bubbles, the hero blurb and appendix entries ── */
.prose p { margin: 0 0 .85rem; } .prose p:last-child { margin-bottom: 0; }
.prose ul, .prose ol { margin: 0 0 .85rem; padding-left: 1.2rem; } .prose li { margin-bottom: .3rem; }
.prose code { font-family: var(--mono); font-size: .88em; background: var(--inset); padding: .1em .32em; border-radius: 2px; }
.prose a { color: var(--accent-d); text-decoration: underline; text-underline-offset: .15em; }
.prose table { width: 100%; border-collapse: collapse; margin: .35rem 0 1rem; font-size: .9rem; }
.prose th { text-align: left; font-weight: 600; border-bottom: 1.5px solid var(--ink); padding: .4rem .6rem .35rem 0; }
.prose td { border-bottom: 1px solid var(--hair); padding: .4rem .6rem .4rem 0; vertical-align: top; }
pre.block {
  font-family: var(--mono); font-size: .82rem; line-height: 1.55; background: var(--inset);
  border: 1px solid var(--hair); border-radius: 2px; padding: .8rem .9rem; overflow-x: auto; margin: 0 0 .9rem;
}
.prose blockquote {
  margin: 1rem 0; padding: .7rem 1rem .7rem 1.1rem; background: var(--inset);
  border-left: 3px solid var(--accent); font-size: 1.02rem; line-height: 1.55;
}
.prose .head2 {
  font-size: 1.12rem; font-weight: 650; letter-spacing: -.01em; margin: 1.7rem 0 .6rem;
  padding-bottom: .3rem; border-bottom: 1px solid var(--hair);
}
.prose .head3 { font-size: .95rem; font-weight: 650; letter-spacing: .015em; color: var(--accent-d); margin: 1.35rem 0 .45rem; }
.prose > .head2:first-child, .prose > .head3:first-child { margin-top: 0; }

/* ── collapsibles authored in the text (`+++`), as opposed to the tool blocks ── */
details.fold { margin: 1.1rem 0; border: 1px solid var(--hair); border-radius: 2px; background: var(--paper); }
details.fold > summary {
  cursor: pointer; list-style: none; padding: .55rem .85rem; font-size: .86rem; font-weight: 600;
  color: var(--accent-d); background: var(--inset); border-bottom: 1px solid transparent;
}
details.fold[open] > summary { border-bottom-color: var(--hair); }
details.fold > summary::-webkit-details-marker { display: none; }
details.fold > summary::before { content: "▸ "; color: var(--grey-lt); font-weight: 400; }
details.fold[open] > summary::before { content: "▾ "; }
.fold-body { padding: .95rem 1rem 1rem; }

/* ── appendix: a deliverable, so a heavier band and a full sheet under it ── */
.phase.band { margin: 3.5rem 0 1.4rem; }
.phase.band::after { height: 2px; background: var(--ink); }
.phase.band span {
  font-size: .78rem; letter-spacing: .11em; padding: .4rem .8rem;
  color: #F4F4F2; background: var(--ink); border-color: var(--ink);
}
.appx {
  background: var(--paper); border: 1px solid var(--hair); border-radius: 2px;
  padding: 1.7rem 1.9rem; line-height: 1.72;
}
.appx p, .appx ul, .appx ol, .appx blockquote, .appx .head2, .appx .head3 { max-width: 40rem; }
.appx > p:first-of-type { font-size: 1.05rem; }

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
/* terminal text is line-oriented: a captured command or its output keeps its own line breaks, or a
   column-aligned table of results reflows into one paragraph and stops being readable as evidence */
.cmd { font-family: var(--mono); font-size: .79rem; line-height: 1.6; margin-top: .6rem; white-space: pre-wrap; overflow-wrap: break-word; }
.cmd .p { color: var(--accent); font-weight: 600; margin-right: .35rem; }
.res { font-family: var(--mono); font-size: .79rem; color: var(--grey); padding-left: 1rem; white-space: pre-wrap; overflow-wrap: break-word; }
.exit { font-size: .68rem; letter-spacing: .04em; padding: .04rem .32rem; border-radius: 2px; margin-left: .45rem; white-space: nowrap; }
.exit.x0 { color: var(--accent-d); background: #E8F3F2; border: 1px solid #BFDEDB; }
.exit.x1 { color: var(--amber); background: #FBF5EA; border: 1px solid #E4D2B0; }
.exit.x2 { color: var(--rose); background: #FDF3F4; border: 1px solid #F3C2C8; }
pre.block.json { font-size: .74rem; max-height: 22rem; margin: .45rem 0 .2rem 1rem; }
details.model { margin: 1rem 0 .2rem; border: 1px solid var(--accent); border-radius: 2px; background: #F7FBFB; }
details.model > summary {
  cursor: pointer; list-style: none; padding: .55rem .85rem; font-size: .86rem; font-weight: 600; color: var(--ink);
}
details.model > summary::-webkit-details-marker { display: none; }
details.model > summary::before { content: "▸ "; color: var(--accent); font-weight: 400; }
details.model[open] > summary::before { content: "▾ "; }
details.model .mc { font-family: var(--mono); font-size: .74rem; font-weight: 400; color: var(--grey); margin-left: .5rem; }
.model-body { padding: .2rem .9rem .9rem; background: var(--paper); border-top: 1px solid #BFDEDB; }
.mlabel { font-size: .78rem; color: var(--grey); margin: .7rem 0 .35rem; }
pre.block.model-src { max-height: 28rem; font-size: .72rem; }
pre.block.diff { font-size: .74rem; max-height: 28rem; }
pre.block.diff ins { text-decoration: none; color: var(--accent-d); background: #E8F3F2; }
pre.block.diff del { text-decoration: none; color: var(--rose); background: #FDF3F4; }
.files { margin-top: .7rem; display: flex; flex-wrap: wrap; gap: .35rem; }
.files span {
  font-family: var(--mono); font-size: .72rem; color: var(--amber);
  border: 1px solid #E4D2B0; background: #FBF5EA; padding: .12rem .4rem; border-radius: 2px;
}

/* ── figures ── */
.figs { display: flex; flex-wrap: wrap; gap: .8rem; margin: 1rem 0 .3rem; }
.figs > * { flex: 1 1 17rem; min-width: 0; }
figure.fig { margin: 0; padding: .9rem; background: var(--paper); border: 1px solid var(--hair); border-radius: 2px; overflow-x: auto; }
figure.fig svg { max-width: 100%; height: auto; display: block; margin: 0 auto; }
figure.fig figcaption { margin-top: .6rem; font-size: .78rem; line-height: 1.5; color: var(--grey); }
.missing { font-family: var(--mono); font-size: .78rem; color: var(--rose); padding: .6rem; border: 1px dashed var(--rose); border-radius: 2px; }

@media (max-width: 34rem) {
  .stat, .figs > * { flex-basis: 100%; }
  .hero, .appx { padding-left: 1.15rem; padding-right: 1.15rem; }
}

footer.foot { margin-top: 3.5rem; padding-top: 1.25rem; border-top: 1px solid var(--hair); font-size: .82rem; color: var(--grey); }
footer.foot code { font-family: var(--mono); font-size: .93em; }

@media print {
  /* On paper the measure is the page, not 50rem: a 50rem wrap is wider than Letter's printable width,
     so any full-width box (the hero, a stat row, a sheet) used to run off the right edge. Hand the
     margin to `@page` and let the wrap fill what is left — which also gets four stat cards back onto
     one row, since they no longer have to wrap inside a box that is itself overflowing. */
  @page { margin: 14mm 12mm; }
  .wrap { max-width: none; padding: 0; }
  /* A fence scrolls sideways on screen; on paper there is nowhere to scroll, so an overlong line would
     be silently cut — the one failure mode a page arguing for auditability cannot have. Wrap instead.
     `.cmd`/`.res` already wrap, which is why captured output only ever soft-wraps. */
  pre.block { white-space: pre-wrap; overflow-wrap: break-word; }
  body { background: #fff; } header.mast { background: #fff; color: #000; border-bottom: 2px solid #000; }
  h1, .headline .big { color: #000; } .headline { background: var(--inset); }
  .sub, .meta, .meta b { color: #444; } .brand { color: var(--accent-d); } .brand b { color: #000; }
  /* force every collapsible open on paper: the `display` pair covers engines that hide the children
     themselves, `::details-content` the ones (Chrome 131+) that hide the slot instead. */
  details.tools, details.fold, details.model { display: block; }
  details.tools > *:not(summary), details.fold > *:not(summary), details.model > *:not(summary) { display: block !important; }
  /* …except a model's full source: on paper it would be pages of JSON nobody reads there. The drawing
     and an arm's diff stay, and the summary line says where the document lives. */
  pre.block.model-src, pre.block.json { max-height: none; } pre.block.model-src { display: none; }
  details.tools::details-content, details.fold::details-content, details.model::details-content {
    content-visibility: visible !important; block-size: auto !important;
  }
  /* Keep the short things whole. An assistant turn with its commands expanded is routinely taller than
     a sheet, and `break-inside: avoid` on that just pushes it to the next page and leaves this one
     half blank — so let it flow across the break instead. */
  .turn.user, .turn.note, .stats, figure.fig { break-inside: avoid; }
  .hero { background: #fff; border-top-color: #000; }
  .stat { border-top-color: #000; }
  .figs { display: block; } .figs > * + * { margin-top: .8rem; }
  .phase.band { break-before: page; } .phase.band span { color: #fff; background: #000; }
  .appx { border: none; padding: 0; }
}
"""

function render(doc, jsonpath)
    io = IOBuffer()
    println(io, "<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"utf-8\">")
    println(io, "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">")
    println(io, "<title>", esc(get(doc, "title", "rd-study session")), "</title>")
    println(io, "<style>", CSS, "</style>")
    # A sheet normally starts a fresh sheet of paper, which is right for a long document and wrong for a
    # short one: four brief sections become four mostly-empty pages. `"compact": true` lets them flow.
    get(doc, "compact", false) === true &&
        println(io, "<style>@media print { .phase.band { break-before: auto; } }</style>")
    println(io, "</head><body>")

    # masthead
    println(io, "<header class=\"mast\"><div class=\"wrap\">")
    kicker = get(doc, "brand", "rd-study session")
    println(io, "<p class=\"brand\"><b>ReactiveDynamics.jl</b> &nbsp;·&nbsp; ", esc(kicker), "</p>")
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
        # A JSON object parses to an unordered Dict, so a strip written `what it is` first can render
        # last. Where the order is part of the argument, author `meta` as a list of one-pair objects.
        for (k, v) in (meta isa AbstractVector ? [only(pairs(e)) for e in meta] : collect(pairs(meta)))
            println(io, "<div><b>", esc(k), "</b> ", inline(v), "</div>")
        end
        println(io, "</div>")
    end
    println(io, "</div></header>")

    # hero, then the authored sections, then the transcript, then the deliverables
    println(io, "<main class=\"wrap\">")
    println(io, hero_html(get(doc, "hero", Dict())))
    secs = sheets_html(get(doc, "sections", []), "Overview")
    isempty(secs) || println(io, secs)
    for turn in resolve_turns(get(doc, "turns", []))
        role = get(turn, "role", "assistant")
        if role == "gap"
            note = string(get(turn, "note", ""))
            println(io, "<div class=\"gap\"><b>", esc(node_text(turn)), "</b>")
            isempty(note) || println(io, "<i>", inline(note), "</i>")
            println(io, "</div>")
            continue
        end

        ph = get(turn, "phase", "")
        isempty(ph) || println(io, "<div class=\"phase\"><span>", esc(ph), "</span></div>")

        println(io, "<div class=\"turn $role\">")
        label = get(ROLE_LABEL, role, "")
        isempty(label) || println(io, "<div class=\"who\">", esc(label), "</div>")
        println(io, "<div class=\"bubble prose\">")
        println(io, stats_html(get(turn, "stats", [])))
        println(io, body_html(node_text(turn)))
        println(io, figures_html(get(turn, "figures", []), string(get(turn, "figure", ""))))
        println(io, models_html(get(turn, "models", [])))
        println(
            io, tools_html(
                get(turn, "tools", []), get(turn, "files", []), get(turn, "open", false) == true
            )
        )
        println(io, "</div></div>")
    end
    println(io, sheets_html(get(doc, "appendix", [])))

    foot = get(doc, "footer", "")
    isempty(foot) || println(io, "<footer class=\"foot\">", body_html(foot), "</footer>")
    println(io, "</main>")
    println(io, "<!-- generated by agent/scripts/render_session.jl from ", esc(jsonpath), " -->")
    println(io, "</body></html>")
    return String(take!(io))
end

function main(argv)
    isempty(argv) && (println(stderr, "usage: render_session.jl PAGE.json [--out=PAGE.html]"); return 2)
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
    turns = resolve_turns(get(doc, "turns", []))
    nsec = length(get(doc, "sections", []))
    ncmds = sum(length(get(t, "tools", [])) for t in turns; init = 0)
    ngap = count(t -> get(t, "role", "") == "gap", turns)
    nturn = length(turns) - ngap
    println("render_session: wrote $out — $nsec sections, $nturn turns$(ngap > 0 ? " (excerpt, $ngap elisions)" : ""), $ncmds commands, $(round(sizeof(html) / 1024, digits = 1)) KiB self-contained")
    return 0
end

# Run as a script, render; included (as `render_deck.jl` does, for the grammar and the figure inliner),
# define and stay quiet.
abspath(PROGRAM_FILE) == (@__FILE__) && exit(main(ARGS))
