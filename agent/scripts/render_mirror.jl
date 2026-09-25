#!/usr/bin/env julia
#
# render_mirror.jl — show how a captured conversation became the model, from the files themselves.
#
#   julia --project=. agent/scripts/render_mirror.jl STUDY/mirror.json
#
# A reader who is told "the interview became a model" has to take it on faith. This makes it
# checkable, in two artifacts written next to the spec:
#
#   MIRROR.md             for each thing the operator said, the verbatim quote, the model element it
#                         became — pretty-printed from the model document, never pasted — and why
#   figures/model_map.svg the model at the level a reader can hold: every place, transition and rule
#                         falls in exactly one named group, and the arrows are the model's own arcs,
#                         aggregated, plus the rules' writes and the rates that read the live state
#   figures/model_net.svg the same model with every element drawn, inside the same groups: the
#                         engine's `network_graph` extraction, styled, plus the rules from the document
#
# Two refusals make it honest rather than illustrative. A quote that is not a substring of that
# operator turn in the session transcript is an error, so a mirror cannot put words in the operator's
# mouth. An element the groups do not cover (or cover twice) is an error, so the map cannot quietly
# leave out the part of the model that would complicate it.
#
# Spec — one JSON object:
#
#   { "model": "<study-relative .rdj.json>", "base": "<optional, the arm the diff is against>",
#     "session": "<repo-relative session .json>", "intro": "<markdown-ish>", "caption": "…",
#     "groups": [ { "id": …, "label": …, "kind": "places"|"transitions"|"rules",
#                   "match": ["line", "k_base_*", …] } ],          # `*` is the only wildcard
#     "rows": [ { "turn": 3, "quote": "…", "elements": ["place:line", "rule:t_file", "param:exit_at_k"],
#                 "gloss": "…" } ],
#     "open": 3 }                                                     # rows shown before the fold
#
# The grammar the output uses is render_session.jl's (`> ` quotes, fences, `+++` folds, `!figure`),
# so MIRROR.md is imported into a page section with `"file"` and renders there unchanged.

using JSON
using ReactiveDynamics
const RD = ReactiveDynamics

include(joinpath(@__DIR__, "chartkit.jl"))
const MK = ChartKit

# ── the closed expression and action vocabulary, printed as a reader would write it ───────────────

const INFIX = Set(["+", "-", "*", "/", ">", "<", ">=", "<=", "==", "&&", "||"])

fmt(v::Real) = isinteger(v) ? string(Int(v)) : string(round(v; sigdigits = 4))
fmt(v) = string(v)

function ex(n)
    n isa Dict || return fmt(n)
    k = get(n, "node", "")
    k == "const" && return fmt(n["value"])
    k == "ref" && return n["name"]
    k == "timeref" && return "t"
    k == "sample" && return string(n["dist"], "(", join(ex.(get(n, "args", [])), ", "), ")")
    if k == "call"
        a = ex.(n["args"])
        length(a) == 2 && n["op"] in INFIX && return "($(a[1]) $(n["op"]) $(a[2]))"
        length(a) == 1 && n["op"] == "!" && return "!" * a[1]
        return string(n["op"], "(", join(a, ", "), ")")
    end
    return JSON.json(n)
end

# the outermost parentheses carry no information in a guard or a rate
unwrap(s) = startswith(s, "(") && endswith(s, ")") && balanced(s[2:(end - 1)]) ? s[2:(end - 1)] : s
function balanced(s)
    d = 0
    for c in s
        d += c == '(' ? 1 : c == ')' ? -1 : 0
        d < 0 && return false
    end
    return d == 0
end

function act(a)
    v = a["verb"]
    v == "seq" && return reduce(vcat, act.(a["stmts"]); init = String[])
    v == "set_params" && return ["$(x["name"]) := $(unwrap(ex(x["value"])))" for x in a["assigns"]]
    v == "set_marking" && return ["$(a["name"]) $(get(a, "mode", "set") == "set" ? ":=" : get(a, "mode", "")) $(unwrap(ex(a["value"])))"]
    return [JSON.json(a)]
end

# ── the model's element index ─────────────────────────────────────────────────────────────────────

struct Model
    d::Dict
    places::Dict{String, Dict}
    trans::Dict{String, Dict}
    rules::Dict{String, Dict}
    params::Dict{String, Any}
end

function Model(d::Dict)
    return Model(
        d, Dict(p["name"] => p for p in d["places"]), Dict(t["id"] => t for t in d["transitions"]),
        Dict(r["id"] => r for r in get(d, "rules", [])), Dict(p["name"] => p["value"] for p in get(d, "params", []))
    )
end

arc_str(a) = string(a["place"], fmt(get(a, "multiplicity", 1)) == "1" ? "" : "×" * fmt(a["multiplicity"]))

function show_element(m::Model, base, spec)
    kind, id = split(spec, ':'; limit = 2)
    if kind == "place"
        p = get(m.places, id, nothing)
        p === nothing && error("render_mirror: no place `$id` in the model")
        bits = ["place $id"]
        haskey(p, "init") && push!(bits, "starts at " * fmt(p["init"]))
        get(p, "cost", nothing) === nothing || push!(bits, "cost " * fmt(p["cost"]) * " per unit")
        get(p, "reward", nothing) === nothing || push!(bits, "reward " * fmt(p["reward"]) * " per unit")
        return join(bits, " · ")
    elseif kind == "transition"
        t = get(m.trans, id, nothing)
        t === nothing && error("render_mirror: no transition `$id` in the model")
        arcs = [a for a in m.d["arcs"] if get(a, "transition", "") == id && haskey(a, "place")]
        lhs = [arc_str(a) * (haskey(a, "modality") ? " [$(a["modality"]["return"])]" : "") for a in arcs if a["side"] == "lhs"]
        rhs = [arc_str(a) for a in arcs if a["side"] == "rhs"]
        out = ["transition $id", "  takes    " * join(lhs, " + "), "  gives    " * (isempty(rhs) ? "nothing (a pure cost)" : join(rhs, " + "))]
        push!(out, "  rate     " * unwrap(ex(t["rate"])))
        get(t, "cycletime", 0) == 0 || push!(out, "  holds    " * fmt(t["cycletime"]) * " month(s)")
        ps = get(t, "prob_of_success", 1.0)
        ps == 1.0 || push!(out, "  succeeds " * unwrap(ex(ps)))
        haskey(t, "priority") && push!(out, "  priority " * fmt(t["priority"]))
        return join(out, "\n")
    elseif kind == "rule"
        r = get(m.rules, id, nothing)
        r === nothing && error("render_mirror: no rule `$id` in the model")
        a = act(r["action"])
        return join(
            vcat(
                ["rule $id  ($(replace(get(r, "fire_mode", "once"), "_" => " ")))", "  when  " * unwrap(ex(r["guard"]))],
                ["  " * (i == 1 ? "then  " : "      ") * s for (i, s) in enumerate(a)]
            ), "\n"
        )
    elseif kind == "param"
        haskey(m.params, id) || error("render_mirror: no parameter `$id` in the model")
        v = fmt(m.params[id])
        b = base === nothing ? nothing : get(base.params, id, nothing)
        return b === nothing || fmt(b) == v ? "param $id = $v" : "param $id = $v    (the base model has $(fmt(b)): this is the option)"
    end
    return error("render_mirror: element kind `$kind` is not place/transition/rule/param")
end

# ── quotes: verbatim or refused ───────────────────────────────────────────────────────────────────

plain(s) = replace(s, "**" => "", r"\s+" => " ")

function check_quote(turns, n, q)
    1 <= n <= length(turns) || error("render_mirror: turn $n is outside the transcript (1:$(length(turns)))")
    t = turns[n]
    t["role"] == "user" || error("render_mirror: turn $n is not an operator turn")
    occursin(plain(q), plain(t["text"])) || error("render_mirror: the quote is not in operator turn $n: \"$q\"")
    return nothing
end

# ── the group map ─────────────────────────────────────────────────────────────────────────────────

globre(p) = Regex("^" * replace(p, "*" => ".*") * "\$")

function assign(groups, kind, ids)
    owner = Dict{String, String}()
    for g in groups
        g["kind"] == kind || continue
        for p in g["match"], id in ids
            occursin(globre(p), id) || continue
            haskey(owner, id) && owner[id] != g["id"] && error("render_mirror: `$id` is in both `$(owner[id])` and `$(g["id"])`")
            owner[id] = g["id"]
        end
    end
    missing_ = [id for id in ids if !haskey(owner, id)]
    isempty(missing_) || error("render_mirror: $kind not in any group: $(join(missing_, ", "))")
    return owner
end

# the places (or params) an expression reads, and the places (or params) an action writes
function refs(n, kind, acc = Set{String}())
    if n isa Dict
        get(n, "node", "") == "ref" && get(n, "kind", "") == kind && push!(acc, n["name"])
        foreach(v -> refs(v, kind, acc), values(n))
    elseif n isa AbstractVector
        foreach(v -> refs(v, kind, acc), n)
    end
    return acc
end

function writes(a, kind, acc = Set{String}())
    a["verb"] == "seq" && foreach(s -> writes(s, kind, acc), a["stmts"])
    kind == "place" && a["verb"] == "set_marking" && push!(acc, a["name"])
    kind == "param" && a["verb"] == "set_params" && foreach(x -> push!(acc, x["name"]), a["assigns"])
    return acc
end

dotesc(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")

function members_label(g, ids)
    shown = length(ids) <= 5 ? ids : vcat(ids[1:4], ["+ $(length(ids) - 4) more"])
    rows = join(["<TR><TD ALIGN=\"LEFT\"><FONT POINT-SIZE=\"10\" FACE=\"Menlo\" COLOR=\"$(MK.GREY)\">$(dotesc(i))</FONT></TD></TR>" for i in shown])
    noun = g["kind"] == "places" ? "place" : g["kind"] == "transitions" ? "transition" : "rule"
    n = length(ids)
    return "<<TABLE BORDER=\"0\" CELLSPACING=\"0\" CELLPADDING=\"1\">" *
        "<TR><TD ALIGN=\"LEFT\"><B><FONT POINT-SIZE=\"14\">$(dotesc(g["label"]))</FONT></B></TD></TR>" *
        "<TR><TD ALIGN=\"LEFT\"><FONT POINT-SIZE=\"10\" COLOR=\"$(MK.ACCENT_D)\">$n $noun$(n == 1 ? "" : "s")</FONT></TD></TR>" *
        rows * "</TABLE>>"
end

function model_map(m::Model, groups)
    po = assign(groups, "places", collect(keys(m.places)))
    to = assign(groups, "transitions", collect(keys(m.trans)))
    ro = assign(groups, "rules", collect(keys(m.rules)))
    byg = Dict(g["id"] => sort([id for (id, o) in merge(po, to, ro) if o == g["id"]]) for g in groups)

    flow = Dict{Tuple{String, String}, Int}()
    for a in m.d["arcs"]
        haskey(a, "place") || continue
        pg, tg = po[a["place"]], to[a["transition"]]
        key = a["side"] == "lhs" ? (pg, tg) : (tg, pg)
        flow[key] = get(flow, key, 0) + 1
    end
    reads = Set{Tuple{String, String}}()
    for (id, t) in m.trans, p in refs(t["rate"], "place")
        haskey(po, p) && push!(reads, (po[p], to[id]))
    end
    # a rule reaches the flow two ways: it sets a marking, or it sets a parameter a rate reads
    sets = Set{Tuple{String, String}}()
    for (id, r) in m.rules
        for p in writes(r["action"], "place")
            haskey(po, p) && push!(sets, (ro[id], po[p]))
        end
        pw = writes(r["action"], "param")
        for (tid, t) in m.trans
            isempty(intersect(pw, refs(t["rate"], "param"))) || push!(sets, (ro[id], to[tid]))
        end
    end

    io = IOBuffer()
    println(io, "digraph model {")
    println(io, "  graph [rankdir=TB, nodesep=0.3, ranksep=0.45, fontname=\"Helvetica\", bgcolor=\"white\", pad=0.2];")
    println(io, "  node [fontname=\"Helvetica\", margin=\"0.12,0.06\", penwidth=1.2];")
    println(io, "  edge [fontname=\"Helvetica\", fontsize=11, arrowsize=0.6, color=\"$(MK.GREY_LT)\", fontcolor=\"$(MK.GREY)\"];")
    for g in groups
        style = g["kind"] == "places" ? "shape=box, style=\"rounded\", color=\"$(MK.ACCENT)\"" :
            g["kind"] == "transitions" ? "shape=box, style=\"filled\", fillcolor=\"$(MK.INSET)\", color=\"$(MK.HAIR)\"" :
            "shape=note, style=\"filled,dashed\", fillcolor=\"#FBF6EC\", color=\"$(MK.AMBER)\""
        println(io, "  \"$(g["id"])\" [$style, label=$(members_label(g, byg[g["id"]]))];")
    end
    for ((a, b), n) in sort(collect(flow))
        println(io, "  \"$a\" -> \"$b\" [label=\"$(n == 1 ? "" : string(n))\", penwidth=$(round(1 + log(n); digits = 2)), color=\"$(MK.INK)\"];")
    end
    for (a, b) in sort(collect(reads))
        println(io, "  \"$a\" -> \"$b\" [style=dotted, color=\"$(MK.ACCENT)\", constraint=false];")
    end
    for (a, b) in sort(collect(sets))
        println(io, "  \"$a\" -> \"$b\" [style=dashed, color=\"$(MK.AMBER)\"];")
    end
    println(io, "}")
    return String(take!(io)), (length(m.places), length(m.trans), length(m.rules))
end

# ── the full net: every place, transition and rule, laid out inside the same groups ───────────────
#
# The structure is the engine's own Layer A extraction (`network_graph`, src/visualize.jl): the
# places, the transitions and every arc with its multiplicity and modality, read off the built model.
# Only the styling is ours. `draw_network` renders the same extraction unstyled and without the rules,
# which carry the model's calendar and decisions, so they are added here from the document.

function net_map(ng, m::Model, groups)
    po = assign(groups, "places", [string(p.name) for p in ng.places])
    to = assign(groups, "transitions", [string(t.name) for t in ng.transitions])
    ro = assign(groups, "rules", collect(keys(m.rules)))
    owner = merge(po, to, ro)
    node(id) = "\"n_$(id)\""

    function tlabel(id)
        t = m.trans[id]
        bits = String[]
        ps = get(t, "prob_of_success", 1.0)
        ps isa Real && ps != 1.0 && push!(bits, "p " * fmt(ps))
        get(t, "cycletime", 0) in (0, 1) || push!(bits, fmt(t["cycletime"]) * " mo")
        haskey(t, "priority") && push!(bits, "prio " * fmt(t["priority"]))
        sub = isempty(bits) ? "" : "<BR/><FONT POINT-SIZE=\"8\" COLOR=\"$(MK.GREY)\">$(join(bits, " · "))</FONT>"
        return "<<FONT FACE=\"Menlo\" POINT-SIZE=\"10\">$(dotesc(id))</FONT>$sub>"
    end
    function plabel(id)
        p = m.places[id]
        init = get(p, "init", 0)
        sub = init isa Real && init != 0 && abs(init) < 1.0e9 ? "<BR/><FONT POINT-SIZE=\"8\" COLOR=\"$(MK.GREY)\">starts $(fmt(init))</FONT>" : ""
        return "<<FONT FACE=\"Menlo\" POINT-SIZE=\"10\">$(dotesc(id))</FONT>$sub>"
    end

    io = IOBuffer()
    println(io, "digraph net {")
    println(io, "  graph [rankdir=LR, nodesep=0.12, ranksep=0.5, fontname=\"Helvetica\", bgcolor=\"white\", pad=0.2, newrank=true];")
    println(io, "  node [fontname=\"Helvetica\", margin=\"0.08,0.03\", penwidth=1.1, height=0.3];")
    println(io, "  edge [arrowsize=0.5, penwidth=0.9, fontname=\"Helvetica\", fontsize=9, fontcolor=\"$(MK.GREY)\"];")
    for g in groups
        ids = sort([id for (id, o) in owner if o == g["id"]])
        isempty(ids) && continue
        frame = g["kind"] == "places" ? MK.ACCENT : g["kind"] == "transitions" ? MK.GREY_LT : MK.AMBER
        println(io, "  subgraph \"cluster_$(g["id"])\" {")
        println(io, "    label=<<B><FONT POINT-SIZE=\"12\">$(dotesc(g["label"]))</FONT></B>>; labeljust=l; style=\"rounded$(g["kind"] == "rules" ? ",dashed" : "")\"; color=\"$frame\"; fontcolor=\"$(MK.INK)\";")
        for id in ids
            if g["kind"] == "places"
                println(io, "    $(node(id)) [shape=box, style=\"rounded\", color=\"$(MK.ACCENT)\", label=$(plabel(id))];")
            elseif g["kind"] == "transitions"
                println(io, "    $(node(id)) [shape=box, style=\"filled\", fillcolor=\"$(MK.INSET)\", color=\"$(MK.HAIR)\", label=$(tlabel(id))];")
            else
                println(io, "    $(node(id)) [shape=note, style=\"filled\", fillcolor=\"#FBF6EC\", color=\"$(MK.AMBER)\", label=<<FONT FACE=\"Menlo\" POINT-SIZE=\"10\">$(dotesc(id))</FONT>>];")
            end
        end
        println(io, "  }")
    end
    seen = Set{Tuple{String, String}}()
    for a in ng.arcs
        f, t = string(a.from), string(a.to)
        (f, t) in seen && continue
        push!(seen, (f, t))
        conserved = :conserved in a.modality
        lbl = a.multiplicity == 1.0 ? "" : ", label=\"$(fmt(a.multiplicity))\""
        println(io, "  $(node(f)) -> $(node(t)) [color=\"$(conserved ? MK.ACCENT : MK.INK)\"$lbl];")
    end
    for (id, r) in sort(collect(m.rules); by = first)
        for p in sort(collect(writes(r["action"], "place")))
            haskey(po, p) && println(io, "  $(node(id)) -> $(node(p)) [style=dashed, color=\"$(MK.AMBER)\"];")
        end
    end
    println(io, "}")
    return String(take!(io))
end

function run_dot(src, out)
    svg = read(pipeline(IOBuffer(src), `dot -Tsvg`), String)
    svg = replace(svg, r"<\?xml.*?\?>\s*"s => "", r"<!DOCTYPE.*?>\s*"s => "", r"<!--.*?-->\s*"s => "")
    write(out, svg)
    return out
end

# ── assemble ──────────────────────────────────────────────────────────────────────────────────────

function mirror_main(args)
    length(args) == 1 || (println(stderr, "usage: render_mirror.jl STUDY/mirror.json"); return 2)
    specf = args[1]
    study = dirname(specf)
    spec = JSON.parsefile(specf)
    m = Model(JSON.parsefile(joinpath(study, spec["model"])))
    base = haskey(spec, "base") ? Model(JSON.parsefile(joinpath(study, spec["base"]))) : nothing
    turns = JSON.parsefile(spec["session"])["turns"]

    mkpath(joinpath(study, "figures"))
    svgf = joinpath(study, "figures", "model_map.svg")
    src, (np, nt, nr) = model_map(m, spec["groups"])
    run_dot(src, svgf)
    netf = joinpath(study, "figures", "model_net.svg")
    modelf = joinpath(study, spec["model"])
    ng = RD.network_graph(RD.from_json_model(read(modelf, String); seed = 1))
    run_dot(net_map(ng, m, spec["groups"]), netf)

    io = IOBuffer()
    println(io, spec["intro"], "\n")
    println(
        io, "!figure ", svgf, " | ", spec["caption"], " Generated from `", spec["model"], "`: all ", np, " places, ", nt,
        " transitions and ", nr, " rules, each in exactly one box. Solid arrows are the model's arcs (the number is how many), dotted arrows are rates that read the live state, dashed arrows are rules setting a marking, or a parameter that a rate reads.\n"
    )
    println(io, "+++~ The same model, every element: ", np, " places, ", nt, " transitions and ", nr, " rules\n")
    println(
        io, "!figure ", netf, " | Every place (teal outline), transition (grey) and rule (amber note) of `", spec["model"],
        "`, inside the boxes of the map above. The places, transitions and arcs are the engine's own extraction (`network_graph`); a teal arc is a resource the transition returns when it completes (`conserved`), a black arc is consumed or produced, and the number is the multiplicity. A transition's second line is its success probability, duration when longer than a month, and priority. Dashed arrows are rules setting a marking.\n"
    )
    println(io, "+++\n")
    nopen = get(spec, "open", length(spec["rows"]))
    for (i, r) in enumerate(spec["rows"])
        check_quote(turns, r["turn"], r["quote"])
        i == nopen + 1 && println(io, "+++~ ", length(spec["rows"]) - nopen, " more: what else the operator said, and where it went\n")
        println(io, "> *Operator, turn ", r["turn"], ":* “", r["quote"], "”\n")
        println(io, "```")
        println(io, join([show_element(m, base, e) for e in r["elements"]], "\n\n"))
        println(io, "```\n")
        println(io, r["gloss"], "\n")
    end
    length(spec["rows"]) > nopen && println(io, "+++\n")
    println(io, "*Generated by `julia --project=. agent/scripts/render_mirror.jl ", specf, "`. Every quote is checked verbatim against the transcript, and every element is printed from the model document.*")
    mdf = joinpath(study, "MIRROR.md")
    write(mdf, String(take!(io)))
    println("render_mirror: wrote $mdf ($(length(spec["rows"])) quotes, all verbatim), $svgf and $netf")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(mirror_main(ARGS))
end
