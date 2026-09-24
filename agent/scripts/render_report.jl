#!/usr/bin/env julia
#
# render_report.jl — assemble a finished study into ONE shareable report, laid out to the rd-report
# signature, as a self-contained HTML file (print it to PDF from any browser or headless Chrome).
#
#   julia --project=. agent/scripts/render_report.jl STUDY_DIR --arm=ARM.rdj.json
#         [--base=model.rdj.json] [--exec-map=figures/MAP.svg [--exec-map-from=DOC]] [--hero=figures/CHART.svg] [--out=PATH]
#
# A study ends in a directory the modeler understands and a reader outside the team does not: a
# FINDINGS.md that cites `BRIEF.md` and `runs/`, a model that is 8 000 lines of JSON, figures in a
# subfolder. This script is what the modeler sends instead. It is GENERATED, never edited, so it
# cannot disagree with the study it was made from: the prose comes from FINDINGS.md, the operator-
# language element table from BRIEF.md's `## Slots`, and everything else is computed from the files.
#
# The six slots of the rd-report signature (agent/skills/rd-report/SKILL.md), and where each comes from:
#
#   1. headline — number and verdict      FINDINGS.md's first `##` heading and its body
#   2. decision rule, caveats as scope    FINDINGS.md's remaining sections, in their own order
#   3. why a spreadsheet cannot           FINDINGS.md's section whose title says "spreadsheet"
#   4. the model                          GENERATED: BRIEF.md's slot table, an inventory of the
#                                          document (pools, timed and state-dependent transitions,
#                                          decision rules), the arm diff, the full element table
#   5. hero visual                        the execution map, `rdx draw ARM --exec-map --seed=1`
#   6. how we got it, provenance          FINDINGS.md's how-it-was-made section, plus GENERATED:
#                                          commit, content hashes, the runs on disk, the replay rule
#
# A report opens with a strip that says, slot by slot, whether it is present and where it came from —
# a missing slot is shown as missing rather than papered over, the same rule the figure inliner follows.

include(joinpath(@__DIR__, "render_session.jl"))

using SHA: sha256
using Dates: today

# ── FINDINGS.md, split into its `##` sections ─────────────────────────────────────────────────────

# md_import has already converted the deliverable to the page grammar, so a section is everything
# from one `## ` line to the next. Fences are tracked so a `## ` inside a code block is not a heading.
function split_sections(text)
    secs = Pair{String, Vector{String}}[]
    pre = String[]
    fenced = false
    for ln in split(text, '\n')
        startswith(strip(ln), "```") && (fenced = !fenced)
        if !fenced && startswith(ln, "## ")
            push!(secs, strip(ln[4:end]) => String[])
        elseif isempty(secs)
            push!(pre, ln)
        else
            push!(last(secs[end]), ln)
        end
    end
    return join(pre, "\n"), [k => join(v, "\n") for (k, v) in secs]
end

strip_bold(s) = strip(replace(s, "**" => ""))

# ── the model document, summarized from its own structure ─────────────────────────────────────────

# Every name an expression reads, by kind. The IR is closed, so this walk is total: a `ref` node is
# the only way an expression can see the state.
function refs!(acc, x)
    if x isa AbstractDict
        get(x, "node", "") == "ref" && push!(get!(acc, string(x["kind"]), Set{String}()), string(x["name"]))
        get(x, "node", "") == "timeref" && push!(get!(acc, "time", Set{String}()), "t")
        foreach(v -> refs!(acc, v), values(x))
    elseif x isa AbstractVector
        foreach(v -> refs!(acc, v), x)
    end
    return acc
end
refs(x) = refs!(Dict{String, Set{String}}(), x)

# What a rule's action writes: markings and params, from the closed verb family.
function writes!(acc, x)
    if x isa AbstractDict
        v = get(x, "verb", "")
        v == "set_marking" && push!(get!(acc, "place", Set{String}()), string(x["name"]))
        v == "set_params" && foreach(a -> push!(get!(acc, "param", Set{String}()), string(a["name"])), x["assigns"])
        foreach(y -> writes!(acc, y), values(x))
    elseif x isa AbstractVector
        foreach(y -> writes!(acc, y), x)
    end
    return acc
end
writes(x) = writes!(Dict{String, Set{String}}(), x)

names_str(s; max = 6) = let v = sort!(collect(s))
    isempty(v) ? "—" : join(("`$n`" for n in first(v, max)), ", ") * (length(v) > max ? " +$(length(v) - max)" : "")
end

num(x) = x isa Real ? (isinteger(x) && abs(x) < 1.0e12 ? string(Int(x)) : string(round(x; sigdigits = 4))) : "expr"

# A rate is read by kind, because the kind is the modeling claim: a constant is an assumption, a
# time-dependent rate is a schedule, a rate that reads a place is FEEDBACK — the state steers the flow.
rate_kind(r) = r isa Real ? "constant $(num(r))" : let a = refs(r)
    haskey(a, "place") ? "state-dependent" : haskey(a, "time") ? "time-dependent" : "parametric"
end

# An arc names a place, creates a structured token, or selects live tokens by predicate (ADR 0008).
arc_target(a) = haskey(a, "place") ? "`$(a["place"])`" :
    haskey(a, "structured") ? "creates a `$(a["structured"]["kind"])` token" :
    haskey(a, "predicate") ? "selects `$(a["predicate"]["kind"])` tokens" : "?"

function modality_str(m)
    m isa AbstractDict || return "consumed"
    return get(m, "return", "") == "conserved" ? "held and returned" : string(get(m, "allocation", ""), "/", get(m, "return", ""))
end

function inventory(d)
    places, trans, arcs, rules, params = (get(d, k, Any[]) for k in ("places", "transitions", "arcs", "rules", "params"))
    held = Set(string(a["place"]) for a in arcs if haskey(a, "place") && get(something(get(a, "modality", nothing), Dict()), "return", "") == "conserved")
    ledger = [p for p in places if haskey(p, "cost") || haskey(p, "reward")]
    feedback = [t for t in trans if rate_kind(t["rate"]) == "state-dependent"]
    timed = [t for t in trans if get(t, "cycletime", 0) isa Real && get(t, "cycletime", 0) > 0]
    risky = [t for t in trans if get(t, "prob_of_success", 1.0) isa Real && get(t, "prob_of_success", 1.0) < 1]
    reads_state = [r for r in rules if haskey(refs(r["guard"]), "place")]
    stochastic = count(r -> occursin("\"sample\"", JSON.json(r)), rules)
    kinds = Set(string(a[k]["kind"]) for a in arcs for k in ("structured", "predicate") if haskey(a, k))
    m = get(d, "meta", Dict())
    rows = [
        "| What | In this model |",
        "| --- | --- |",
        "| Horizon and tick | `tspan = $(num(get(m, "tspan", 0)))`, `dt = $(num(get(m, "dt", 1)))` |",
        "| Places | $(length(places)): $(length(held)) are pools a transition holds and returns ($(names_str(held))); $(length(ledger)) carry a cost or reward into the ledger ($(names_str(Set(p["name"] for p in ledger)))) |",
        "| Transitions | $(length(trans)): $(length(timed)) take time (cycle time > 0), $(length(risky)) can fail (success probability < 1), **$(length(feedback)) have a rate that reads the state** ($(names_str(Set(t["id"] for t in feedback)))) |",
        "| Decision rules | **$(length(rules))**, $(count(r -> r["fire_mode"] == "once", rules)) firing once and $(count(r -> r["fire_mode"] == "every_tick", rules)) every tick; $(length(reads_state)) decide on the live state, $(stochastic) draw a random input |",
        "| Parameters | $(length(params)), every one a named value an arm can change |",
    ]
    isempty(kinds) || insert!(rows, 6, "| Structured tokens | $(names_str(kinds)): entities with attributes, created and selected by predicate as they move through their lifecycle |")
    return join(rows, "\n")
end

# One row per element, in the formal vocabulary — the reviewer's table, folded by default.
function element_table(d)
    arcs = get(d, "arcs", Any[])
    io = IOBuffer()
    println(io, "| Element | Kind | Formal role | Reads | Writes / emits |")
    println(io, "| --- | --- | --- | --- | --- |")
    for p in get(d, "places", Any[])
        n = string(p["name"])
        uses = [a for a in arcs if get(a, "place", nothing) == n]
        role = isempty(uses) ? "state or tally, set by rules" :
            join(unique(string(a["side"] == "lhs" ? "input to" : "output of", " `", a["transition"], "` (", a["side"] == "lhs" ? modality_str(get(a, "modality", nothing)) : "emitted", ")") for a in uses), "; ")
        led = join(filter(!isempty, [haskey(p, "cost") ? "cost $(num(p["cost"]))/unit" : "", haskey(p, "reward") ? "reward $(num(p["reward"]))/unit" : ""]), ", ")
        println(io, "| `", n, "` | place | ", haskey(p, "init") ? "init $(num(p["init"])); " : "", role, " | — | ", isempty(led) ? "—" : led, " |")
    end
    for t in get(d, "transitions", Any[])
        a = refs(t["rate"])
        out = [arc_target(x) for x in arcs if x["transition"] == t["id"] && x["side"] == "rhs"]
        ins = [arc_target(x) for x in arcs if x["transition"] == t["id"] && x["side"] == "lhs" && !haskey(x, "place")]
        println(
            io, "| `", t["id"], "` | transition | rate ", rate_kind(t["rate"]), ", cycle ", num(get(t, "cycletime", 0)),
            ", success ", num(get(t, "prob_of_success", 1.0)), " | ", join(filter(!isempty, [names_str(get(a, "place", Set{String}())) == "—" ? "" : names_str(get(a, "place", Set{String}())); ins]), "; ") |> x -> isempty(x) ? "—" : x, " | ",
            isempty(out) ? "—" : join(out, ", "), " |"
        )
    end
    for r in get(d, "rules", Any[])
        g, w = refs(r["guard"]), writes(r["action"])
        println(
            io, "| `", r["id"], "` | rule | fires ", replace(string(r["fire_mode"]), "_" => " "), " when its guard holds | ",
            names_str(union(get(g, "place", Set{String}()), get(g, "param", Set{String}()))), " | ",
            names_str(union(get(w, "place", Set{String}()), get(w, "param", Set{String}()))), " |"
        )
    end
    return String(take!(io))
end

# BRIEF.md's `## Slots` section: the operator-language table the study was framed in. Taken verbatim,
# because it is the one place the business words and the formal targets sit side by side and were
# signed off by the person who owns the decision.
function brief_slots(study)
    f = joinpath(study, "BRIEF.md")
    isfile(f) || return ""
    _, secs = split_sections(md_import(f))
    i = findfirst(s -> lowercase(first(s)) == "slots", secs)
    return i === nothing ? "" : last(secs[i])
end

# ── provenance, from the repository and the files on disk ─────────────────────────────────────────

git(args...) = strip(read(ignorestatus(`git $args`), String))

sha12(f) = bytes2hex(open(sha256, f))[1:12]

# `runs/` names are `<arm>__<world>__s<seed>.json` by the study convention; anything else is counted
# but not grouped, so an unconventional study still gets an honest total.
function runs_table(study)
    dirs = filter(isdir, [joinpath(study, "runs"), joinpath(study, "runs_ens")])
    isempty(dirs) && return "No `runs/` directory: the report's numbers cannot be traced to stored runs."
    groups = Dict{Tuple{String, String}, Vector{Int}}()
    other = 0
    for dir in dirs, f in readdir(dir)
        m = match(r"^(.+?)__(.+?)__s(\d+)\.json$", f)
        m === nothing ? (other += 1) : push!(get!(groups, (m[1], m[2]), Int[]), parse(Int, m[3]))
    end
    io = IOBuffer()
    seedstr(v) = (sort!(v); v == collect(first(v):last(v)) ? (length(v) == 1 ? "$(first(v))" : "$(first(v))–$(last(v))") : join(v, ","))
    println(io, "| Document | Paths | Scenario: seeds |")
    println(io, "| --- | --- | --- |")
    for arm in sort!(unique(first.(keys(groups))))
        ws = sort!([(w, v) for ((a, w), v) in groups if a == arm]; by = first)
        println(io, "| `", arm, "` | ", sum(length(v) for (_, v) in ws), " | ", join(("$w: $(seedstr(v))" for (w, v) in ws), " · "), " |")
    end
    total = sum(length, values(groups); init = 0)
    other > 0 && println(io, "\nPlus $other stored run file(s) outside the naming convention.")
    return "$total stored single-seed runs, each the verbatim JSON output of one `rdx run`:\n\n" * String(take!(io))
end

function provenance(study, docs, argv)
    head = git("rev-parse", "--short", "HEAD")
    last = git("log", "-1", "--format=%h %cs", "--", study)
    dirty = !isempty(git("status", "--porcelain", "--", study, ":!$study/REPORT.*"))     # the report itself does not count
    rows = [
        "| | |",
        "| --- | --- |",
        "| Repository commit | `$head`$(dirty ? " — **with uncommitted changes in the study**" : "") |",
        "| Study last changed | `$last` |",
        "| Report generated | $(today()), by `agent/scripts/render_report.jl` |",
        "| Engine | ReactiveDynamics.jl $(something(match(r"version = \"(.*?)\"", read("Project.toml", String)), ["?"])[1]) |",
    ]
    for f in docs
        push!(rows, "| `$(relpath(f))` | sha256 `$(sha12(f))…` |")
    end
    cmd = "julia --project=. agent/scripts/render_report.jl " * join(argv, " ")
    return join(rows, "\n") * "\n\n**The replay rule.** A run is determined by its document and its seed: the same pair gives the same output, bit for bit, on any machine. To check a number, run the document named above with the seed in the table below and compare against the stored run.\n\n" *
        runs_table(study) * "\n\n**Regenerate this report** from the study as it stands:\n\n```bash\n$cmd\n```"
end

# ── assembling the page ───────────────────────────────────────────────────────────────────────────

const REPORT_CSS = """
.execmap { max-height: 78vh; overflow: auto; border: 1px solid #C7C7C4; background: #fff; padding: 8px; }
.execmap figure.fig { border: 0; padding: 0; }
@media print {
  .screen-only { display: none !important; }
  .execmap { max-height: none; overflow: visible; border: 0; padding: 0; break-inside: avoid; }
  .execmap figure.fig svg { width: auto !important; max-width: 100%; height: 8.6in !important; }
}
"""

function build(study, opts, argv)
    findings = joinpath(study, "FINDINGS.md")
    isfile(findings) || error("render_report: no FINDINGS.md in $study — the report publishes it, it does not write it")
    base = joinpath(study, get(opts, "base", "model.rdj.json"))
    arm = joinpath(study, opts["arm"])
    isfile(arm) || error("render_report: no such arm: $arm")
    d = JSON.parsefile(arm)

    title_line = something(findfirst(l -> startswith(l, "# "), readlines(findings)), 0)
    title = title_line == 0 ? basename(study) : strip(readlines(findings)[title_line][3:end])
    pre, secs = split_sections(md_import(findings))
    isempty(secs) && error("render_report: FINDINGS.md has no `##` sections to lay out")

    headline = strip_bold(first(secs[1]))
    is_sheet(pat) = s -> occursin(pat, lowercase(first(s)))
    sheet = findfirst(is_sheet("spreadsheet"), secs)
    how = findfirst(s -> occursin(r"how (the numbers|we got)|provenance"i, first(s)), secs)
    rest = [s for (i, s) in enumerate(secs) if i != 1 && i != sheet && i != how]

    execmap = get(opts, "exec-map", "")
    execmap = isempty(execmap) ? "" : joinpath(study, execmap)
    hero = get(opts, "hero", "")
    slots = brief_slots(study)
    docs = filter(isfile, [base, arm])

    have(x) = x === nothing || x == "" || (x isa AbstractString && !isfile(x)) ? "**missing**" : "present"
    strip_rows = [
        "| Slot | Status | Source |",
        "| --- | --- | --- |",
        "| 1 · Headline: number and verdict | present | FINDINGS.md, first section |",
        "| 2 · Decision rule, caveats as scope | $(isempty(rest) ? "**missing**" : "present") | FINDINGS.md, $(length(rest)) section(s) |",
        "| 3 · Why a spreadsheet cannot produce this | $(sheet === nothing ? "**missing**" : "present") | FINDINGS.md |",
        "| 4 · The model: elements, units, the arm diff | present | generated from `$(basename(arm))`$(isempty(slots) ? "" : " and BRIEF.md's slot table") |",
        "| 5 · Hero visual: the execution map | $(have(execmap)) | $(isempty(execmap) ? "run `rdx draw ARM --exec-map --seed=1` and pass `--exec-map`" : "`$(relpath(execmap))`") |",
        "| 6 · How we got it, and provenance | present | $(how === nothing ? "" : "FINDINGS.md, plus ")generated from the repository and `runs/` |",
    ]

    sections = Any[]
    push!(
        sections, Dict(
            "title" => "About this report",
            "text" => "This report is generated from the study directory `$(relpath(study))` and is laid out to the six-part report template every study follows. It is safe to forward: every figure is embedded, nothing links back into the repository, and nothing in it was written by hand after the study closed. The prose is the study's own briefing, the model section is computed from the model document that produced the numbers, and the provenance section names the commit and the runs.\n\n" * join(strip_rows, "\n") * (isempty(strip(pre)) ? "" : "\n\n" * strip(pre)),
        )
    )
    push!(sections, Dict("title" => "1 · The answer and the decision rule", "text" => last(secs[1])))
    for (t, body) in rest
        push!(sections, Dict("title" => "2 · " * strip_bold(t), "text" => body))
    end
    sheet === nothing || push!(sections, Dict("title" => "3 · " * strip_bold(first(secs[sheet])), "text" => last(secs[sheet])))

    model_text = "The model is data: a document in the engine's closed vocabulary, checked by a validator before it runs, with no code in it. What follows is computed from `$(basename(arm))`, the recommended option, so it describes the model that produced the numbers above.\n\n" * inventory(d)
    isempty(slots) || (model_text *= "\n\n### The elements in the language of the business\n\nThe table the study was framed in (BRIEF.md), signed off before the model was built: each slot in the operator's words, and the formal element it became.\n\n" * slots)
    nrows = sum(length(get(d, k, [])) for k in ("places", "transitions", "rules"))
    model_text *= "\n\n### Every element\n\nThe reviewer's table: each place, transition and rule in the formal vocabulary, with what it reads and what it writes. It is folded in the HTML version and left out of print, where its $nrows rows would outweigh the report.\n\n+++ The full element table ($nrows rows)\n" * element_table(d) * "\n+++"
    model_text *= "\n\n### What the option changes\n\nThe recommended option is the base model with the lines below changed — that diff *is* the lever. Open it to see it whole."
    push!(
        sections, Dict(
            "title" => "4 · The model",
            "text" => model_text,
            "models" => isfile(base) && base != arm ? [Dict("file" => arm, "base" => base, "label" => "`$(basename(arm))` against `$(basename(base))`")] : Any[],
        )
    )
    push!(
        sections, Dict(
            "title" => "5 · Hero visual: the execution map",
            "text" => isempty(execmap) ? "*(no execution map was passed — run `rdx draw $(basename(arm)) --exec-map --seed=1` and pass `--exec-map`)*" :
                "The network of the recommended option after one run (seed 1): places are pools and states, boxes are transitions, and a highlighted place is one that ran dry at some point — where the contention was. Drawn by `rdx draw $(get(opts, "exec-map-from", basename(arm))) --exec-map --seed=1`.\n\n@@EXECMAP@@",
        )
    )
    how_text = (how === nothing ? "" : last(secs[how]) * "\n\n") * "### Provenance\n\n" * provenance(study, docs, argv)
    push!(sections, Dict("title" => "6 · " * (how === nothing ? "Provenance" : strip_bold(first(secs[how]))), "text" => how_text))

    doc = Dict(
        "brand" => "study report · generated",
        "title" => title,
        "subtitle" => "A ReactiveDynamics.jl study, generated from `$(relpath(study))` and laid out to the six-part report template.",
        "headline" => headline,
        "meta" => [
            Dict("study" => "`$(relpath(study))`"),
            Dict("recommended option" => "`$(basename(arm))`"),
            Dict("generated" => "$(today()) at commit `$(git("rev-parse", "--short", "HEAD"))`"),
        ],
        "hero" => isempty(hero) ? Dict() : Dict("figures" => [Dict("svg" => joinpath(study, hero), "caption" => "The study's headline chart.")]),
        "sections" => sections,
        "compact" => true,
        "footer" => "Generated by `agent/scripts/render_report.jl` from the files of the study; edit the study, not this page, and regenerate. The template is the one in `agent/skills/rd-report/SKILL.md`.",
    )
    html = render(doc, relpath(findings))
    html = replace(html, "</head>" => "<style>$REPORT_CSS</style></head>"; count = 1)
    html = replace(html, "<details class=\"fold\"><summary>The full element table" => "<details class=\"fold screen-only\"><summary>The full element table")
    isempty(execmap) || (html = replace(html, "<p>@@EXECMAP@@</p>" => "<div class=\"execmap\">" * figure_html(execmap) * "</div>"))
    return html, length(sections)
end

function report_main(argv)
    pos = filter(a -> !startswith(a, "--"), argv)
    opts = Dict(String(m[1]) => String(m[2]) for m in (match(r"^--([^=]+)=(.*)$", a) for a in argv) if m !== nothing)
    (length(pos) == 1 && haskey(opts, "arm")) ||
        (println(stderr, "usage: render_report.jl STUDY_DIR --arm=ARM.rdj.json [--base=model.rdj.json] [--exec-map=SVG [--exec-map-from=DOC]] [--hero=SVG] [--out=PATH]"); return 2)
    study = String(rstrip(pos[1], '/'))
    isdir(study) || (println(stderr, "render_report: no such study: $study"); return 2)
    out = get(opts, "out", joinpath(study, "REPORT.html"))
    html, nsec = try
        build(study, opts, argv)
    catch e
        e isa ErrorException || rethrow()
        println(stderr, e.msg)
        return 2
    end
    write(out, html)
    println("render_report: wrote $out — $nsec sections, $(round(sizeof(html) / 1024, digits = 1)) KiB self-contained")
    return 0
end

abspath(PROGRAM_FILE) == (@__FILE__) && exit(report_main(ARGS))
