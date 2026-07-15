# Visualization — the network "exec map" (ADR 0014 §B / CONTRACT §15.2).
#
# Three SEPARABLE layers, so a user can extract/inspect the graph with no plotting deps, render it
# with Graphviz, and decorate it with run results — each step optional:
#
#   Layer A  network_graph(prob) -> NetworkGraph   — dependency-free, pure data, no plotting dep.
#   Layer B  to_graphviz(g) / draw_network(prob)   — Petri-net DOT, rendered via AA's run_graphviz,
#                                                     composable with AA's wiring_diagram.
#   Layer C  exec_map(prob_or_ens; highlight)      — decorate Layer A with §14 run/ensemble stats
#                                                     (bottlenecks hot, starvation, token paths).
#
# The result-plot RECIPES (§15.1) are NOT here — they live in ext/RDPlotsExt.jl (they need Plots).
# This file is dependency-free for Layer A; Layer B emits a DOT string and only touches Graphviz
# through AA's `run_graphviz` when actually rendering; Layer C reads finished-run statistics.

export NetworkGraph, network_graph, to_graphviz, draw_network, exec_map

# ════════════════════════════════════════════════════════════════════════════════════════
# Layer A — structure extraction (dep-free, pure data)
# ════════════════════════════════════════════════════════════════════════════════════════

# A species (place) node: its name and whether it is a structured/agentic species.
struct SpeciesNode
    name::Symbol
    structured::Bool
end

# A transition node: its name (or a synthesized `transition_<i>` when unnamed) and the attributes
# worth annotating on the diagram (rate/priority/cycletime as display strings — kept as strings so
# Layer A stays a pure value with no live attribute eval).
struct TransitionNode
    name::Symbol
    index::Int
    label::String
end

# An arc: LHS species → transition (`:in`) or transition → RHS species (`:out`), with stoichiometry
# and the modality set (the §1 truth-table tags consumed/conserved/nonblock/rate) that styles it.
struct Arc
    from::Symbol
    to::Symbol
    dir::Symbol                 # :in (species→transition) or :out (transition→species)
    stoich::Float64
    modality::Set{Symbol}
end

"""
    NetworkGraph

A plain, inspectable Petri-net view of a model (ADR 0014 Layer A / §15.2 Invariant 1): species
(place) nodes, transition nodes, and the arcs between them with stoichiometry + modality. Built by
`network_graph` with NO plotting/Graphviz dependency and NO simulation — a pure function of the
model — so the diagram is authoring-time documentation as well as a run artifact.
"""
struct NetworkGraph
    species::Vector{SpeciesNode}
    transitions::Vector{TransitionNode}
    arcs::Vector{Arc}
end

# Normalize a reactant's `species` to a Symbol node id. A plain species is already a Symbol; a
# structured/parameterized product can be an `Expr` (FoldedReactant.species is `Union{Expr,Symbol}`)
# — render it as a Symbol of its source text so the place node is still well-defined and stable.
_species_sym(s::Symbol) = s
_species_sym(s) = Symbol(string(s))

# The DOT node id for transition index `i` — the single source of truth shared by `network_graph`
# (which builds the nodes) and `exec_map` (which must reference the SAME id when highlighting a
# token's path). A named transition uses its `transName`; an unnamed one uses `transition_<i>`. NB
# this is NOT the per-instance name `past_bonds` carries (`"<transName>_@<t>"`, solvers.jl) — the
# instance suffix must be dropped, which is why highlighting maps through the transition INDEX
# (`Transition.i`) rather than the bond's instance name.
function _transition_node_name(net, i)
    tname = net[i, :transName]
    return (tname === nothing || tname === missing) ? Symbol("transition_$i") : Symbol(tname)
end

# Display label for a transition node — its rate/priority/cycletime, read off the net recipe columns
# as strings (no eval). Falls back gracefully when a column is unset.
function _transition_label(net, i)
    name = net[i, :transName]
    base = (name === nothing || name === missing) ? "transition_$i" : string(name)
    parts = String[base]
    return join(parts, "\n")
end

"""
    network_graph(prob::ReactionNetworkProblem) -> NetworkGraph

Extract the Petri-net structure of a constructed model (Layer A). Walks the species table for the
place nodes (flagging structured/agentic species) and the transition incidence for the arcs — today
via `transLHS`/`transRHS` (the parsed reactant lists + the RHS expression). Because that incidence is
realized by `sample_transitions!` (which draws stoichiometries through the RNG), this runs on a
`deepcopy` of `prob` so the caller's `state.rng` is NOT perturbed — `network_graph` is observationally
pure (no simulation, Invariant 1). The extraction simplifies (typed FKs, no reactant re-parse) when
the ADR 0003 `ReactantSpec` table lands; this is the `transLHS`/`transRHS` form noted in §15.2.
"""
function network_graph(prob::ReactionNetworkProblem)
    net = prob.network
    species = SpeciesNode[
        SpeciesNode(net[i, :specName], net[i, :specStructured] === true) for i in row_ids(net, :S)
    ]

    # Realize the incidence on a copy so the original RNG is untouched.
    work = deepcopy(prob)
    sample_transitions!(work)

    transitions = TransitionNode[]
    arcs = Arc[]
    lhs = work.transitions[:transLHS]
    rhs = work.transitions[:transRHS]
    known_species = Set(net[i, :specName] for i in row_ids(net, :S))
    for i in eachindex(lhs)
        tnode_name = _transition_node_name(net, i)
        push!(transitions, TransitionNode(tnode_name, i, _transition_label(net, i)))

        # LHS reactants → transition (consumed/blocking/nonblock arcs). The structured LHS species of
        # this transition (if any) is the @advance/@move target's true place — the produced token IS
        # the bound program (identity preserved, ADR 0008 §D), so an @advance RHS resolves back to it.
        struct_lhs = nothing
        for r in lhs[i]
            sp = _species_sym(r.species)
            (sp in known_species && net[find_index(sp, work), :specStructured] === true) && (struct_lhs = sp)
            push!(arcs, Arc(sp, tnode_name, :in,
                r.stoich isa Real ? Float64(r.stoich) : 1.0, r.modality))
        end
        # transition → RHS products. The RHS expr is parsed by extract_reactants on the copy.
        rprods = try
            extract_reactants(rhs[i], work)
        catch
            []
        end
        for r in rprods
            modality = r.modality isa Set ? r.modality : Set{Symbol}()
            stoich = hasproperty(r, :stoich) && r.stoich isa Real ? Float64(r.stoich) : 1.0
            sp = _species_sym(r.species)
            # An @advance/@move RHS is a macro Expr, not a plain species; its destination place is the
            # transition's structured LHS species (phase is an attribute, the kind is unchanged). Map
            # such a non-species RHS node back to that place so the arc connects to a real place rather
            # than a synthetic node named after the raw macro text.
            sp in known_species || (struct_lhs === nothing || (sp = struct_lhs))
            push!(arcs, Arc(tnode_name, sp, :out, stoich, modality))
        end
    end
    return NetworkGraph(species, transitions, arcs)
end

# ════════════════════════════════════════════════════════════════════════════════════════
# Layer B — static rendering (Petri-net DOT; render via AA's run_graphviz)
# ════════════════════════════════════════════════════════════════════════════════════════

# Arc color by modality (the §1 truth-table legend, ADR 0014 §15.2 / §15.4 Invariant — documented as
# the overlay legend): conserved arcs (returned) distinct from consumed; nonblock distinct again.
function _arc_color(modality::Set{Symbol})
    :conserved in modality && return "darkgreen"
    :nonblock in modality && return "orange"
    :rate in modality && return "gray40"
    return "black"            # plain consumed
end

# Quote + escape a DOT node id/label.
_dotstr(s) = "\"" * replace(string(s), "\"" => "\\\"") * "\""

"""
    to_graphviz(g::NetworkGraph; highlight_species = Symbol[], highlight_arcs = Tuple{Symbol,Symbol}[]) -> String

Emit Graphviz DOT for the Petri net (Layer B): species as circles, transitions as boxes, arcs with
stoichiometry labels and color by §1 modality. `highlight_species`/`highlight_arcs` paint a subset
(used by Layer C's overlay). Returns a DOT digraph STRING — rendering is deferred to `draw_network`
(via AA's `run_graphviz`), so emitting the structure needs no Graphviz backend (Invariant 2). Valid
DOT for any model; the smoke tests check `dot` accepts it for SIR/toy-pharma.
"""
function to_graphviz(g::NetworkGraph;
                     highlight_species::AbstractVector = Symbol[],
                     highlight_arcs::AbstractVector = Tuple{Symbol,Symbol}[])
    hs = Set(Symbol.(highlight_species))
    ha = Set(highlight_arcs)
    io = IOBuffer()
    println(io, "digraph \"reactive_network\" {")
    println(io, "  rankdir=LR;")
    println(io, "  node [fontsize=9];")
    # species = circles (double circle / filled if highlighted)
    for s in g.species
        shape = s.structured ? "doublecircle" : "circle"
        fill = s.name in hs ? ", style=filled, fillcolor=gold" : ""
        println(io, "  $(_dotstr(s.name)) [shape=$shape$fill];")
    end
    # transitions = boxes
    for t in g.transitions
        println(io, "  $(_dotstr(t.name)) [shape=box, label=$(_dotstr(t.label))];")
    end
    # arcs
    for a in g.arcs
        color = _arc_color(a.modality)
        key = (a.from, a.to)
        pen = key in ha ? ", penwidth=3.0" : ""
        lbl = a.stoich == 1.0 ? "" : ", label=$(_dotstr(string(a.stoich)))"
        println(io, "  $(_dotstr(a.from)) -> $(_dotstr(a.to)) [color=$color$lbl$pen];")
    end
    println(io, "}")
    return String(take!(io))
end

"""
    draw_network(prob; format = "svg", prog = :dot, path = nothing, kwargs...)

Render the Petri net of `prob` (Layer B): builds the `network_graph`, emits DOT (`to_graphviz`, with
any `highlight_*` kwargs forwarded), and renders it through AlgebraicAgents' `run_graphviz` (which
uses `Graphviz_jll` if present, else a system `dot`). With `path` given, writes the rendered output
there and returns the path; otherwise returns the rendered bytes as a `String`. Rendering is reuse,
not reinvention — no new graph library (Invariant 2). If no Graphviz backend is available the DOT
string is still obtainable via `to_graphviz(network_graph(prob))`.
"""
function draw_network(prob::ReactionNetworkProblem; format::AbstractString = "svg",
                      prog::Symbol = :dot, path = nothing, kwargs...)
    dot = to_graphviz(network_graph(prob); kwargs...)
    if path === nothing
        io = IOBuffer()
        AlgebraicAgents.run_graphviz(io, dot; prog = prog, format = format)
        return String(take!(io))
    else
        AlgebraicAgents.run_graphviz(String(path), dot; prog = prog, format = format)
        return path
    end
end

# ════════════════════════════════════════════════════════════════════════════════════════
# Layer C — the result overlay ("exec map" / inefficiency view)
# ════════════════════════════════════════════════════════════════════════════════════════

# Pool trough (lowest level reached) per species over a run — the starvation signal.
function _pool_troughs(prob::ReactionNetworkProblem)
    troughs = Dict{Symbol,Float64}()
    for s in prob.network[:, :specName]
        col = string(s)
        if col in names(prob.sol)
            troughs[s] = minimum(prob.sol[!, col])
        end
    end
    return troughs
end

"""
    exec_map(prob; highlight = nothing, format = "svg", path = nothing, prog = :dot) -> String or path

The result-decorated "exec map" (Layer C / §15.2): the Petri net of `prob` with run statistics
painted on — species nodes filled where their pool ran to a trough (starvation), and, when a
`highlight::TokenPredicate` (a `@select` set, ADR 0008 / §9.5 — Invariant 4) is given, the matching
cohort's `past_bonds` path THROUGH the net drawn as thickened arcs ("where did these programs go").
Decorates Layer A with finished-run statistics ONLY — it never mutates state or re-runs dynamics
(Invariant 3). Returns the rendered output (or the `path` written to); the underlying DOT is always
available via `to_graphviz`. The overlay is a styling pass over §14 data, not a new computation.
"""
function exec_map(prob::ReactionNetworkProblem; highlight = nothing,
                  format::AbstractString = "svg", path = nothing, prog::Symbol = :dot)
    g = network_graph(prob)

    # Starvation: species whose pool hit (near) zero at its trough.
    troughs = _pool_troughs(prob)
    starved = Symbol[s for (s, v) in troughs if v <= 0.0]

    # Token-path highlighting from past_bonds, scoped by the @select predicate. A bond is a
    # `(species, t, transition)` triple; the transition node id must be the SAME `_transition_node_name`
    # the graph uses — derived from the transition's INDEX (`Transition.i`), NOT the bond's per-instance
    # name `"<transName>_@<t>"` (which would never match a graph node). Each bond highlights the
    # species→transition arc the token traversed.
    hi_arcs = Tuple{Symbol,Symbol}[]
    if highlight isa TokenPredicate
        for tok in select_tokens(prob, highlight)
            for (species, _t, transition) in tok.past_bonds
                push!(hi_arcs, (species, _transition_node_name(prob.network, transition.i)))
            end
        end
    end

    dot = to_graphviz(g; highlight_species = starved, highlight_arcs = hi_arcs)
    if path === nothing
        io = IOBuffer()
        AlgebraicAgents.run_graphviz(io, dot; prog = prog, format = format)
        return String(take!(io))
    else
        AlgebraicAgents.run_graphviz(String(path), dot; prog = prog, format = format)
        return path
    end
end
