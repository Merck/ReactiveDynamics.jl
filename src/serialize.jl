# Single-JSON model serialization (ADR 0005), eval-free. A model is one JSON object with `meta`
# + top-level arrays `params[]`, `places[]`, `transitions[]`, `arcs[]`, `observables[]`,
# `events[]`. Every expression is a node-tagged `ExprNode` dict (never a Julia source string);
# `from_json_model` parses → validates → lowers each node via `to_expr` into the same `Expr`
# columns the `@reaction_network` DSL fills, then constructs a `ReactionNetworkProblem`.
# Uses the existing JSON.jl dependency with hand-rolled node-tagged (de)serialization.

import JSON

export node_to_dict, node_from_dict, model_to_dict, build_network_from_dict
export from_json_model, to_json_model

# ── ADR 0017 Tier 3: the retired wire spellings, accepted for ONE release ───────────────
# The rename reached the format last and on its own commit: `species[]` → `places[]`,
# the arc array `reactants` → `arcs`, an arc's `species` → `place`, a population entry's `species` →
# `place`, the action verb `set_species` → `set_marking`, a ref kind `species` → `place`.
#
# The WRITER emits only the new spellings. The READER takes either for one release and warns
# once per site on the old one, retiring together with the Tier-1 name shims. There is no
# permanent alias: the maintainer resolved (2026-08-22) that no saved document exists outside
# this repository, and a second accepted spelling kept forever would just rebuild the
# two-vocabulary problem ADR 0017 exists to remove. Re-export a document to update it.
_legacy_key(d::AbstractDict, new::String, old::String, default) = if haskey(d, new)
    d[new]
elseif haskey(d, old)
    Base.depwarn(
        "the serialized key `\"$old\"` is deprecated (ADR 0017); use `\"$new\"` — " *
            "re-export the document to update it.", :build_network_from_dict,
    )
    d[old]
else
    default
end

# Which spelling a document used for one of the renamed arrays, for a diagnostic PATH string —
# `validate` reports the path the reader actually walked, so a message stays greppable in the
# document it describes.
_legacy_path(d::AbstractDict, new::String, old::String) = haskey(d, new) || !haskey(d, old) ? new : old

# A `NodeRef` kind off the wire: the retired `species` spelling normalizes to `place`.
_ref_kind(k::Symbol) = if k === :species
    Base.depwarn(
        "the serialized ref kind `\"species\"` is deprecated (ADR 0017); use `\"place\"`.",
        :node_from_dict,
    )
    :place
else
    k
end

# ── ExprNode ⟷ JSON dict (the recursive node-tagged union) ──────────────────────────────
"""
    node_to_dict(n::ExprNode) -> Dict{String, Any}

Serialize one [`ExprNode`](@ref) to its node-tagged JSON dict — the recursive half of the eval-free (de)serialization. Each dict carries a `"node"` tag (`"const"`/`"ref"`/`"call"`/`"sample"`/`"timeref"`/`"choose"`/`"field"`/`"externalref"`) plus that node's fields, with child nodes serialized recursively; a `Symbol`-valued `Const` is flagged so [`node_from_dict`](@ref) can recover it as a Symbol rather than a string. Never emits a Julia source string. The inverse is [`node_from_dict`](@ref).
"""
node_to_dict(n::Const) = Dict{String, Any}(
    "node" => "const",
    "value" => n.value isa Symbol ? string(n.value) : n.value,
    # tag a Symbol-valued const so node_from_dict can recover it (vs a string param name)
    "symbol" => n.value isa Symbol,
)
node_to_dict(n::NodeRef) =
    Dict{String, Any}("node" => "ref", "kind" => string(n.kind), "name" => string(n.name))
node_to_dict(n::Call) =
    Dict{String, Any}("node" => "call", "op" => string(n.op), "args" => map(node_to_dict, n.args))
node_to_dict(n::Sample) =
    Dict{String, Any}("node" => "sample", "dist" => string(n.dist), "args" => map(node_to_dict, n.args))
node_to_dict(::TimeRef) = Dict{String, Any}("node" => "timeref")
node_to_dict(n::Choose) = Dict{String, Any}(
    "node" => "choose",
    "alts" => [Dict{String, Any}("weight" => w, "value" => node_to_dict(v)) for (w, v) in n.alts],
)
node_to_dict(n::Field) = Dict{String, Any}("node" => "field", "name" => string(n.name))
# ExternalRef (ADR 0012 §B2): a declared inputs[] port read. The JSON carries only the port NAME;
# the foreign-agent topology that fills it lives host-side in add_wire! (Invariant 4, eval-free).
node_to_dict(n::ExternalRef) = Dict{String, Any}("node" => "externalref", "port" => string(n.port))

"""
    node_from_dict(d::AbstractDict) -> ExprNode

Parse a node-tagged JSON dict back to the typed [`ExprNode`](@ref) it denotes — the inverse of [`node_to_dict`](@ref) and the recursive half of the eval-free deserialization. Dispatches on the `"node"` tag; child nodes are parsed recursively. NEVER `Meta.parse`/`eval`s — a `"const"` string is recovered as a Symbol only when its `"symbol"` flag is set (else Int/Float64/Bool per the JSON scalar), and an unknown tag is a hard `error`. An out-of-whitelist op/dist is NOT caught here (it round-trips as a Symbol) — `validate` is the whitelist gate.
"""
function node_from_dict(d::AbstractDict)
    tag = d["node"]
    if tag == "const"
        v = d["value"]
        get(d, "symbol", false) === true && return Const(Symbol(v))
        # JSON numbers: keep Int vs Float64 (an integer-valued Float64 from a `1.0` literal stays
        # Float64 here; validate() coerces where integrality is required, ADR 0005:166).
        v isa Bool && return Const(v)
        v isa Integer && return Const(Int(v))
        return Const(Float64(v))
    elseif tag == "ref"
        return NodeRef(_ref_kind(Symbol(d["kind"])), Symbol(d["name"]))
    elseif tag == "call"
        return Call(Symbol(d["op"]), ExprNode[node_from_dict(a) for a in d["args"]])
    elseif tag == "sample"
        return Sample(Symbol(d["dist"]), ExprNode[node_from_dict(a) for a in d["args"]])
    elseif tag == "timeref"
        return TimeRef()
    elseif tag == "choose"
        return Choose([(Float64(a["weight"]), node_from_dict(a["value"])) for a in d["alts"]])
    elseif tag == "field"
        return Field(Symbol(d["name"]))
    elseif tag == "externalref"
        return ExternalRef(Symbol(d["port"]))
    else
        error("node_from_dict: unknown node tag $(tag)")
    end
end

# A scalar attribute may be authored as a bare JSON number/bool/string OR a node dict; normalize
# to an ExprNode either way (Const for a bare literal, the tagged node otherwise).
_attr_node(x::AbstractDict) = node_from_dict(x)
_attr_node(x::Bool) = Const(x)
_attr_node(x::Integer) = Const(Int(x))
_attr_node(x::Real) = Const(Float64(x))
_attr_node(x::AbstractString) = Const(Symbol(x))   # a string scalar is a literal symbol (e.g. a phase)

# ── Build a ReactionNetwork from a parsed model dict (E2: scalar attrs) ─────
# Lowers each ExprNode to the Expr column the constructor consumes. Transitions are assembled
# from arcs[] into the :trans reaction line (E4); for E2 a transition may carry an explicit
# `reaction` string-free node-list, but the minimal path supports params + place + a transition
# whose arcs[] are plain (no modality/predicate) — assembled by assemble_reaction_line (E4).
"""
    build_network_from_dict(d::AbstractDict; registry = Dict{Symbol, Any}()) -> ReactionNetwork

Build a static [`ReactionNetwork`](@ref) store from a parsed model dict `d`, eval-free. Reads the top-level `params[]`/`places[]`/`transitions[]`/`arcs[]`/`observables[]` arrays, lowering each attribute [`ExprNode`](@ref) (or bare literal) to the `Expr` column the constructor consumes via [`to_expr`](@ref) — param values are JSON numbers taken verbatim (never `eval`'d), and a transition's `arcs[]` assemble into its `:trans` reaction-line Expr. The retired ADR-0017 key spellings (`species[]`, `reactants[]`, an arc's `species`) are still read for one release, with a deprecation warning. Assumes `d` has already passed `validate`; it is the structural core of [`from_json_model`](@ref) and the inverse of [`model_to_dict`](@ref). Host functions referenced by an action/genesis kind are resolved by NAME through `registry`, never carried in `d`.
"""
function build_network_from_dict(d::AbstractDict; registry = Dict{Symbol, Any}())
    net = ReactionNetwork()

    # params[] → :P (values are JSON numbers, never eval'd — replaces loadsave.jl:65)
    for pr in get(d, "params", [])
        add_row!(net, :P; prmName = Symbol(pr["name"]), prmVal = pr["value"])
    end

    # places[] → :S (placeInitVal/placeCost/… are scalar Const ExprNodes → literals)
    for pl in _legacy_key(d, "places", "species", [])
        i = add_row!(net, :S; placeName = Symbol(pl["name"]))
        haskey(pl, "init") && (net[i, :placeInitVal] = to_expr(_attr_node(pl["init"])))
        haskey(pl, "cost") && (net[i, :placeCost] = to_expr(_attr_node(pl["cost"])))
        haskey(pl, "reward") && (net[i, :placeReward] = to_expr(_attr_node(pl["reward"])))
        haskey(pl, "valuation") && (net[i, :placeValuation] = to_expr(_attr_node(pl["valuation"])))
        get(pl, "structured", false) === true && (net[i, :placeStructured] = true)
        # modality 3-axis → Set{Symbol} (E6); default empty set = row 1
        haskey(pl, "modality") && (net[i, :placeModality] = modality_from_dict(pl["modality"]))
    end

    # transitions[] → :T. arcs[] for this transition assemble into the :trans reaction line.
    arcs_by_tr = Dict{String, Vector{Any}}()
    for r in _legacy_key(d, "arcs", "reactants", [])
        push!(get!(arcs_by_tr, string(r["transition"]), []), r)
    end
    for tr in get(d, "transitions", [])
        id = string(tr["id"])
        line = assemble_reaction_line(get(arcs_by_tr, id, []))
        rate_node = _attr_node(tr["rate"])
        rate_mode = Symbol(get(tr, "rate_mode", "poisson"))
        i = add_row!(
            net,
            :T;
            trans = line,
            transRate = lower_rate(rate_node, rate_mode),
            transName = haskey(tr, "name") ? String(tr["name"]) : missing,
        )
        for (jsonkey, col) in (
                "cycletime" => :transCycleTime,
                "prob_of_success" => :transProbOfSuccess,
                "capacity" => :transCapacity,
                "priority" => :transPriority,
                "max_lifetime" => :transMaxLifeTime,
                "multiplier" => :transMultiplier,
            )
            haskey(tr, jsonkey) && (net[i, col] = to_expr(_attr_node(tr[jsonkey])))
        end
    end

    # observables[] (E6) and events[] (E5) are added by their step's loaders.
    haskey(d, "observables") && _load_observables!(net, d["observables"])

    assign_defaults!(net)
    return net
end

# Deprecated alias (ADR 0015 Tier 2): `build_acs_from_dict` → `build_network_from_dict`.
@deprecate build_acs_from_dict(d::AbstractDict; registry = Dict{Symbol, Any}()) build_network_from_dict(d; registry = registry)

# meta[] → keywords. A few string-valued meta keys are symbolized for backward compatibility.
# `alloc_strategy`/`strategy` are now accepted-and-ignored (ADR 0002 makes priority-weighted
# progressive filling the single allocation policy — there is no longer a :weighted/:greedy
# switch), so symbolizing them is harmless; `schedule` is likewise a legacy no-op key.
const _SYMBOL_META = (:alloc_strategy, :strategy, :schedule)
function _meta_kwargs(d::AbstractDict)
    m = get(d, "meta", Dict{String, Any}())
    kw = Dict{Symbol, Any}()
    for (k, v) in m
        key = Symbol(k)
        kw[key] = (key in _SYMBOL_META && v isa AbstractString) ? Symbol(v) : v
    end
    return kw
end

# Load a model FRAGMENT from a JSON file as a static `ReactionNetwork` (NOT a constructed
# ReactionNetworkProblem). This is the target of `@join`'s file-include branch (joins.jl): when a
# `@join` argument is a path/macrocall, the macro expands to `include_model(path)`, whose result is
# fed to `merge_networks!`. Unlike `from_json_model`, it does NOT require `meta.tspan` — a joined fragment
# carries no simulation horizon of its own; the composed whole supplies it. The loader is eval-free
# (build_network_from_dict / validate — the ADR-0005 typed IR), so the file-include path carries no RCE.
function include_model(path::AbstractString; registry = Dict{Symbol, Any}())
    d = JSON.parse(read(path, String))
    diags = validate(d; registry = registry)
    isempty(diags) ||
        error("include_model: $(repr(path)) failed validation:\n" * join(string.(diags), "\n"))
    return build_network_from_dict(d; registry = registry)
end

# ── from_json / to_json (the model envelope) ────────────────────────────────────────────
"""
    from_json_model(json::AbstractString; seed = nothing, registry = Dict{Symbol, Any}(), population = []) -> ReactionNetworkProblem

Parse a single-JSON model document into a runnable [`ReactionNetworkProblem`](@ref) — the public IMPORT half of the round-trip (its inverse is [`to_json_model`](@ref)). Pipeline: `JSON.parse` → [`validate`](@ref) (gates on `isempty(diags)`, else a hard `error` listing the diagnostics) → [`build_network_from_dict`](@ref) → construct. Eval-free throughout: a model IS data. `meta.tspan` (the simulation horizon) is REQUIRED; `meta.seed` supplies `seed` when the kwarg is `nothing`. `rules[]` (ADR 0010) become typed `Rule`s and `inputs[]` (ADR 0012 §B1) become the declared external read ports + pre-wire defaults. Host functions named by an action/genesis `kind` are resolved through `registry`; the AA wiring that fills `inputs[]` ports is host-side, never in the document.
"""
function from_json_model(json::AbstractString; seed = nothing, registry = Dict{Symbol, Any}(), population = [])
    d = JSON.parse(json)
    diags = validate(d; registry = registry)
    isempty(diags) || error("from_json_model: model failed validation:\n" * join(string.(diags), "\n"))
    net = build_network_from_dict(d; registry = registry)
    kw = _meta_kwargs(d)
    haskey(kw, :tspan) || error(
        "from_json_model: meta.tspan is required (the simulation horizon) — add e.g. " *
            "\"meta\": { \"tspan\": 100.0, \"dt\": 1.0 } to the model document.",
    )
    seed === nothing && haskey(kw, :seed) && (seed = kw[:seed])
    # rules[] (ADR 0010) → typed Rule structs passed to the constructor (the endogenous channel).
    rules = Any[rule_from_dict(r) for r in get(d, "rules", [])]
    # inputs[] (ADR 0012 §B1) → the declared external read ports + their pre-wire defaults. The
    # default seeds `state.external_inputs[port]` at construction so a port read before any wire
    # delivers a value (or in a standalone run with no wires) is still well-defined (§B3).
    external_inputs = inputs_from_dict(get(d, "inputs", []))
    return ReactionNetworkProblem(
        net;
        seed = seed,
        registry = registry,
        population = population,
        rules = rules,
        external_inputs = external_inputs,
        filter(p -> p.first ∉ (:seed,), kw)...,
    )
end

# The declared place + param NAME sets, typed Set{Symbol} (an empty comprehension would infer
# Set{Any}, which from_expr/rate_from_expr reject). Threaded into every from_expr call so a stored
# attribute Expr's bare symbols classify back to the right NodeRef kind.
function _name_sets(net::ReactionNetwork)
    places = Set{Symbol}(net[i, :placeName] for i in row_ids(net, :S))
    params = Set{Symbol}(net[i, :prmName] for i in row_ids(net, :P) if !isnothing(net[i, :prmName]))
    return places, params
end

# model_to_dict is the EXPORT envelope — the structural inverse of build_network_from_dict
# (~line 73). It emits every top-level array build_network_from_dict reads back: params[], places[],
# transitions[], arcs[], observables[], plus rules[] when the caller passes a constructed
# model's typed Rule vector (a DSL/loaded model → JSON → from_json_model → equivalent model).
# The place/param NAME SETS are threaded into every from_expr call below so a stored attribute
# Expr's bare symbols are classified back to the right NodeRef kind (place vs param), matching
# how the authoring DSL named them — exactly the inverse of to_expr's name→state.u[i]/state.p[:k]
# substitution (ADR 0005 §66).
"""
    model_to_dict(net::ReactionNetwork; meta = Dict{String, Any}(), rules = [], inputs = Dict{Symbol, Any}()) -> Dict{String, Any}

The EXPORT envelope: emit the JSON dict of a static [`ReactionNetwork`](@ref) — every top-level array [`build_network_from_dict`](@ref) reads back (`params[]`/`places[]`/`transitions[]`/`arcs[]`, plus `observables[]`/`rules[]`/`inputs[]` when present), and hence its structural inverse. Only the current ADR-0017 key spellings are written; the retired ones are read, never emitted. Each stored attribute `Expr` is lowered back to an [`ExprNode`](@ref) dict via [`from_expr`](@ref), with the net's place/param NAME sets threaded through so a bare symbol classifies to the right `NodeRef` kind (matching how the DSL named it). `rules`/`inputs` are passed through by `to_json_model(::ReactionNetworkProblem)` since they live on the problem, not the net; legacy `:E` event rows and `RawExpr` actions are intentionally not emitted (not JSON-serializable). Called by [`to_json_model`](@ref).
"""
function model_to_dict(
        net::ReactionNetwork; meta = Dict{String, Any}(), rules = [],
        inputs = Dict{Symbol, Any}()
    )
    places, params = _name_sets(net)
    d = Dict{String, Any}(
        "rd_format" => "reactive-dynamics-model",
        "version" => "1.0",
        "meta" => meta,
        "params" => [
            Dict{String, Any}("name" => string(net[i, :prmName]), "value" => net[i, :prmVal])
                for i in row_ids(net, :P) if !isnothing(net[i, :prmName])
        ],
        "places" => [_place_to_dict(net, i) for i in row_ids(net, :S)],
        "transitions" => [_transition_to_dict(net, i; places, params) for i in row_ids(net, :T)],
        "arcs" => _arcs_to_dict(net),
    )
    # observables[] (the inverse of _load_observables!, ~line 413) — emit only if any :obs row.
    obs = [obs_to_dict(net[i, :obsName], net[i, :obsOpts]) for i in row_ids(net, :obs)]
    isempty(obs) || (d["observables"] = obs)
    # rules[] (the endogenous channel, ADR 0010) — a constructed model carries its typed Rules on
    # the ReactionNetworkProblem (prob.rules), so to_json_model(::ReactionNetworkProblem) passes
    # them through here. A bare schema net has none. Legacy :E event rows are NOT emitted: they are
    # lifted to RawExpr-action Rules at construction (solvers.jl ~line 671), and RawExpr is the
    # non-typed bridge that is intentionally not JSON-serializable (stmt_to_dict(::RawExpr) errors).
    rule_dicts = [rule_to_dict(r) for r in rules if r.action isa ActionStmt && !(r.action isa RawExpr)]
    isempty(rule_dicts) || (d["rules"] = rule_dicts)
    # inputs[] (ADR 0012 §B1) — the inverse of inputs_from_dict (~line 430). Declared external read
    # ports + their pre-wire literal defaults live on the ReactionNetworkProblem (external_input_-
    # defaults), so to_json_model(::ReactionNetworkProblem) passes them through. A default is a
    # literal value (inputs_from_dict requires a Const), so it round-trips as a `const` node. Sorted
    # by port name for deterministic output (the buffer is an unordered Dict).
    input_dicts = [
        Dict{String, Any}("port" => string(p), "default" => node_to_dict(Const(inputs[p])))
            for p in sort!(collect(keys(inputs)))
    ]
    isempty(input_dicts) || (d["inputs"] = input_dicts)
    return d
end

"""
    to_json_model(net::ReactionNetwork; meta = Dict{String, Any}()) -> String
    to_json_model(prob::ReactionNetworkProblem; meta = Dict{String, Any}()) -> String

Serialize a model to a single JSON string — the public EXPORT half of the round-trip (its inverse is [`from_json_model`](@ref); `to_json_model` ∘ `from_json_model` reconstructs an equivalent model). Delegates to [`model_to_dict`](@ref) then `JSON.json`. Given a constructed [`ReactionNetworkProblem`](@ref) it also round-trips the problem-level state the net does not carry: the typed `Rule`s (via `rules[]`), the declared external-input ports + defaults (via `inputs[]`), and — reconstructed into `meta` from the solver fields when the caller supplies none — `tspan`/`dt`/`seed`, so the exported document is COMPLETE and re-importable on its own (`from_json_model` requires `meta.tspan`). An explicit `meta` always wins. Host functions are referenced by NAME through the registry, never embedded.
"""
to_json_model(net::ReactionNetwork; meta = Dict{String, Any}()) =
    JSON.json(model_to_dict(net; meta = meta))
# A constructed model also carries its typed Rules — round-trip them through rules[] — and its
# solver settings (tspan/dt/seed), which live on the ReactionNetworkProblem (not the net) and merge
# into the meta bag at construction (solvers.jl ~line 544). When the caller supplies no `meta`, we
# reconstruct it from those fields so the exported document is COMPLETE and re-importable on its own
# (from_json_model requires meta.tspan). An explicit `meta` always takes precedence.
function to_json_model(prob::ReactionNetworkProblem; meta = Dict{String, Any}())
    full = _meta_from_prob(prob)
    merge!(full, meta)   # caller-supplied keys win
    return JSON.json(
        model_to_dict(
            prob.network; meta = full, rules = prob.rules,
            inputs = prob.external_input_defaults
        )
    )
end

# Reconstruct the meta bag from a constructed model's solver fields. `tspan` is stored as a
# (t0, tend) tuple; the JSON `tspan` is the horizon `tend` (the scalar from_json_model passes on).
function _meta_from_prob(prob::ReactionNetworkProblem)
    m = Dict{String, Any}("tspan" => float(prob.tspan[2]), "dt" => float(prob.dt))
    prob.seed === nothing || (m["seed"] = prob.seed)
    return m
end

# ── Rate lowering + unwrapping (E3) ─────────────────────────────────────────────────────
# A JSON rate carries the BARE intensity; the engine's expand_rate wraps it (Poisson per tick,
# or used as-is under @deterministic). lower_rate reproduces that so the :transRate column matches.
function lower_rate(rate_node::ExprNode, rate_mode::Symbol)
    bare = to_expr(rate_node)
    if rate_mode === :deterministic
        return bare                                   # used as-is (create.jl @deterministic path)
    else
        return :(rand(state.rng, Poisson(max(state.dt * $bare, 0))))   # matches expand_rate
    end
end

# The inverse: recover (bare-intensity ExprNode, rate_mode) from a stored :transRate Expr, so a
# DSL-authored model can be serialized to JSON. A wrapped `rand(state.rng, Poisson(max(state.dt *
# <bare>, 0)))` ⇒ (<bare>, :poisson); anything else ⇒ (rate, :deterministic).
function rate_from_expr(rate; places = Set{Symbol}(), params = Set{Symbol}())
    if rate isa Expr && rate.head == :call && rate.args[1] == :rand
        distcall = rate.args[end]
        if distcall isa Expr && distcall.head == :call && distcall.args[1] == :Poisson
            maxcall = distcall.args[2]                       # max(state.dt * <bare>, 0)
            if maxcall isa Expr && maxcall.head == :call && maxcall.args[1] == :max
                prod = maxcall.args[2]                       # state.dt * <bare>
                if prod isa Expr && prod.head == :call && prod.args[1] == :*
                    bare = prod.args[3]                      # the bare intensity (state.dt is args[2])
                    return from_expr(bare; places, params), :poisson
                end
            end
        end
    end
    return from_expr(rate; places, params), :deterministic
end

# ── Reaction-line assembly (E4) ─────────────────────────────────────────────────────────
# Assemble a transition's arcs[] (grouped lhs/rhs) into the single :trans reaction-line Expr
# `LHS --> RHS` that merge_network!/the runtime parser consume. An arc row may carry: `place`
# (or a `predicate` for @select on the LHS), `stoich`, `modality` (3-axis, LHS), or `advance`
# (field+value, an RHS @advance) — assembled into exactly the macrocall Expr shapes the parser
# (reaction_parser.jl / create.jl) expects.
const _LN = LineNumberNode(0, :none)

# The place "atom" of an arc: a bare place symbol, or a @select(Kind, clauses) macrocall
# (LHS predicate), or a @advance(field, value)/@structured/@move macrocall (RHS).
function _arc_atom(r::AbstractDict)
    if haskey(r, "predicate")                       # @select(Kind, clause && clause …)
        pd = r["predicate"]
        kind = Symbol(pd["kind"])
        clauses = get(pd, "clauses", [])
        if isempty(clauses)
            return Expr(:macrocall, Symbol("@select"), _LN, kind)
        end
        clause_exprs = [_clause_expr(c) for c in clauses]
        conj = foldl((a, b) -> Expr(:(&&), a, b), clause_exprs)
        return Expr(:macrocall, Symbol("@select"), _LN, kind, conj)
    elseif haskey(r, "advance")                     # @advance(field, value) — RHS field write
        adv = r["advance"]
        return Expr(
            :macrocall,
            Symbol("@advance"),
            _LN,
            Symbol(adv["field"]),
            to_expr(_attr_node(adv["value"])),
        )
    elseif haskey(r, "structured")                  # @structured(:Kind, field = value, …) — RHS genesis
        st = r["structured"]
        # The named, eval-free genesis product (ADR 0005 §39): the JSON carries the registry KIND
        # plus field-value nodes; the host constructor is resolved by name at firing time. Assemble
        # the exact macrocall the runtime's named branch consumes (structured_rhs, solvers.jl) — a
        # quoted kind symbol followed by `field = <lowered node>` kwargs, one per fields[] entry.
        kwargs = [Expr(:(=), Symbol(f["name"]), to_expr(_attr_node(f["value"]))) for f in st["fields"]]
        return Expr(:macrocall, Symbol("@structured"), _LN, QuoteNode(Symbol(st["kind"])), kwargs...)
    else
        return Symbol(_legacy_key(r, "place", "species", nothing))
    end
end

# A predicate clause `field op value` (op ∈ PRED_OP_WHITELIST), value lowered via to_expr.
function _clause_expr(c)
    field = Symbol(c[1])
    op = Symbol(c[2])
    val = to_expr(_attr_node(c[3]))
    return Expr(:call, op, field, val)
end

# Wrap an LHS atom in its modality macros (@conserved/@rate/@nonblock) per the 3-axis modality.
# The parser unions these macro names into the arc's modality Set (reaction_parser.jl).
function _apply_modality(atom, m)
    m === nothing && return atom
    s = m isa AbstractDict ? modality_from_dict(m) : m   # Set{Symbol}
    out = atom
    # nest so the innermost wraps the place; order is irrelevant (the parser unions a Set)
    :conserved in s && (out = Expr(:macrocall, Symbol("@conserved"), _LN, out))
    :rate in s && (out = Expr(:macrocall, Symbol("@rate"), _LN, out))
    :nonblock in s && (out = Expr(:macrocall, Symbol("@nonblock"), _LN, out))
    return out
end

# A full arc term: optional integer stoich coefficient × the (modality-wrapped) atom.
function _arc_term(r::AbstractDict; lhs::Bool)
    atom = _arc_atom(r)
    lhs && haskey(r, "modality") && (atom = _apply_modality(atom, r["modality"]))
    stv = to_expr(_attr_node(get(r, "stoich", 1)))
    return (stv == 1 || stv === 1.0) ? atom : Expr(:call, :*, stv, atom)
end

function _sum_terms(terms)
    isempty(terms) && return :∅
    length(terms) == 1 && return terms[1]
    return foldl((a, b) -> Expr(:call, :+, a, b), terms)
end

function assemble_reaction_line(arcs)
    lhs = [_arc_term(r; lhs = true) for r in arcs if String(r["side"]) == "lhs"]
    rhs = [_arc_term(r; lhs = false) for r in arcs if String(r["side"]) == "rhs"]
    return Expr(:call, :→, _sum_terms(lhs), _sum_terms(rhs))
end

# ── E5: action statement + predicate (de)serialization ──────────────────────────────────
# Action VALUES are ExprNodes in JSON, lowered via to_expr to the Expr the apply_action!/
# _eval_value path consumes (actions.jl). RawExpr is intentionally NOT serializable (the legacy
# bridge for non-typed Exprs — a JSON model uses typed verbs only; ADR 0005 open question).
#
# NOTE (serialize-direction fidelity): the to_dict path lowers a stored action-value Expr back to
# a node via from_expr WITHOUT a place/params context, so a bare symbol classifies to its
# default NodeRef(:place, …). This is runtime-harmless — both NodeRef kinds lower to the same
# bare symbol via to_expr, so a from_json-loaded model is unaffected — and only mislabels the JSON
# `ref.kind` tag when re-serializing a model whose actions were hand-built from raw Exprs. JSON-
# authored actions carry typed nodes and never round-trip through from_expr, so they are exact.
stmt_to_dict(s::SetMarking) = Dict{String, Any}(
    "verb" => "set_marking", "name" => string(s.name),
    "value" => node_to_dict(from_expr(s.value)), "mode" => string(s.mode)
)
stmt_to_dict(s::SetParams) = Dict{String, Any}(
    "verb" => "set_params",
    "assigns" => [Dict("name" => string(n), "value" => node_to_dict(from_expr(v))) for (n, v) in s.assigns]
)
stmt_to_dict(s::SetField) = Dict{String, Any}(
    "verb" => "set_field", "field" => string(s.field), "value" => node_to_dict(from_expr(s.value))
)
stmt_to_dict(s::SetTokens) = Dict{String, Any}(
    "verb" => "set_tokens", "predicate" => pred_to_dict(s.predicate),
    "assigns" => [Dict("name" => string(n), "value" => node_to_dict(from_expr(v))) for (n, v) in s.assigns]
)
stmt_to_dict(s::AddToken) = Dict{String, Any}(
    "verb" => "add_token", "kind" => string(s.kind),
    "fields" => [Dict("name" => string(n), "value" => node_to_dict(from_expr(v))) for (n, v) in s.fields]
)
stmt_to_dict(s::Activate) = Dict{String, Any}("verb" => "activate", "transition" => string(s.transition))
stmt_to_dict(s::Deactivate) = Dict{String, Any}("verb" => "deactivate", "transition" => string(s.transition))
stmt_to_dict(s::Invoke) =
    Dict{String, Any}("verb" => "invoke", "fn" => string(s.fn), "args" => [node_to_dict(from_expr(a)) for a in s.args])
stmt_to_dict(s::Log) = Dict{String, Any}("verb" => "log", "msg" => s.msg isa Union{Expr, Symbol} ? node_to_dict(from_expr(s.msg)) : s.msg)
stmt_to_dict(s::Seq) = Dict{String, Any}("verb" => "seq", "stmts" => [stmt_to_dict(x) for x in s.stmts])
stmt_to_dict(::RawExpr) =
    error("RawExpr is not JSON-serializable (the legacy non-typed bridge) — re-express with typed action verbs")

# A value-expr in an action field arrives as a node dict (typed) or a bare literal.
_stmt_value(x) = to_expr(_attr_node(x))

function stmt_from_dict(d::AbstractDict)
    verb = d["verb"]
    if verb == "set_marking" || verb == "set_species"
        verb == "set_species" && Base.depwarn(
            "the serialized action verb `\"set_species\"` is deprecated (ADR 0017); use " *
                "`\"set_marking\"` — re-export the document to update it.", :stmt_from_dict,
        )
        return SetMarking(Symbol(d["name"]), _stmt_value(d["value"]), Symbol(get(d, "mode", "set")))
    elseif verb == "set_params"
        return SetParams([Symbol(a["name"]) => _stmt_value(a["value"]) for a in d["assigns"]])
    elseif verb == "set_field"
        return SetField(Symbol(d["field"]), _stmt_value(d["value"]))
    elseif verb == "set_tokens"
        return SetTokens(
            pred_from_dict(d["predicate"]),
            [Symbol(a["name"]) => _stmt_value(a["value"]) for a in d["assigns"]]
        )
    elseif verb == "add_token"
        return AddToken(
            Symbol(d["kind"]),
            [Symbol(a["name"]) => _stmt_value(a["value"]) for a in d["fields"]]
        )
    elseif verb == "activate"
        return Activate(Symbol(d["transition"]))
    elseif verb == "deactivate"
        return Deactivate(Symbol(d["transition"]))
    elseif verb == "invoke"
        return Invoke(Symbol(d["fn"]), Any[_stmt_value(a) for a in get(d, "args", [])])
    elseif verb == "log"
        return Log(d["msg"] isa AbstractDict ? _stmt_value(d["msg"]) : d["msg"])
    elseif verb == "seq"
        return Seq(ActionStmt[stmt_from_dict(x) for x in d["stmts"]])
    else
        error("stmt_from_dict: unknown action verb $(verb)")
    end
end

# TokenPredicate ⟷ JSON (the @select predicate). Clause value is an ExprNode (or literal).
pred_to_dict(p::TokenPredicate) = Dict{String, Any}(
    "kind" => string(p.kind),
    "clauses" => [[string(c.field), string(c.op), node_to_dict(from_expr(c.value))] for c in p.clauses]
)
function pred_from_dict(d::AbstractDict)
    clauses = Clause[]
    for c in get(d, "clauses", [])
        op = Symbol(c[2])
        op in PRED_OP_WHITELIST || error("pred_from_dict: op $op ∉ PRED_OP_WHITELIST")
        push!(clauses, Clause(Symbol(c[1]), op, _stmt_value(c[3])))
    end
    return TokenPredicate(Symbol(d["kind"]), clauses)
end

# A Rule ⟷ JSON (ADR 0010): id, guard ExprNode, action stmt, fire_mode.
rule_to_dict(r::Rule) = Dict{String, Any}(
    "id" => string(r.id), "guard" => node_to_dict(from_expr(r.guard)),
    "action" => stmt_to_dict(r.action), "fire_mode" => string(r.fire_mode)
)
rule_from_dict(d::AbstractDict) = Rule(
    Symbol(d["id"]), _stmt_value(d["guard"]), stmt_from_dict(d["action"]);
    fire_mode = Symbol(get(d, "fire_mode", "every_tick"))
)

# ── E6: modality 3-axis ⟷ Set{Symbol} — the bijective 5-row translation (CONTRACT §1.1/§1.3) ──
# The 5 legal rows (CONTRACT §1.3). An illegal combination is rejected with a citing message.
function to_set(allocation::Symbol, ret::Symbol, blocking::Symbol)
    # row 5: nonblock ⇒ consumed (CONTRACT §1.4 illegal: nonblock+conserved)
    blocking === :nonblock && ret === :conserved &&
        error("illegal modality (CONTRACT §1.4): blocking=nonblock requires return=consumed")
    s = Set{Symbol}()
    allocation === :perstep && push!(s, :rate)
    ret === :conserved && push!(s, :conserved)
    blocking === :nonblock && push!(s, :nonblock)
    return s
end
function from_set(s::Set{Symbol})
    return (
        allocation = (:rate in s ? :perstep : :upfront),
        return_ = (:conserved in s ? :conserved : :consumed),
        blocking = (:nonblock in s ? :nonblock : :block),
    )
end

function modality_from_dict(m::AbstractDict)
    return to_set(
        Symbol(get(m, "allocation", "upfront")),
        Symbol(get(m, "return", "consumed")),
        Symbol(get(m, "blocking", "block")),
    )
end

# ── inputs[] loader (ADR 0012 §B1): declared external read ports + pre-wire defaults ────
# A model-local `inputs[]` array names the OBSERVABLE-level ports the network may read; each port
# may carry an optional `default` — the value used before any AA wire has delivered (and in a
# standalone, wire-less run). The default is a LITERAL value (a bare JSON scalar or a Const node),
# not a live expression: it must be a concrete `state.external_inputs[port]` seed, evaluable with
# no `state`/RNG. Returns a `Dict{Symbol,Any}(port => default_value)`; a port without a default is
# omitted (a read before its wire delivers then KeyErrors, surfacing the missing-default — exactly
# as a missing param would). The wiring itself is host-side (`add_wire!`), never in this document.
function inputs_from_dict(inputs)
    seed = Dict{Symbol, Any}()
    for inp in inputs
        port = Symbol(inp["port"])
        if haskey(inp, "default")
            node = _attr_node(inp["default"])
            node isa Const ||
                error(
                "inputs_from_dict: port `$port` default must be a literal value (a bare " *
                    "scalar or a Const node), not a live expression — got $(typeof(node))"
            )
            seed[port] = node.value
        end
    end
    return seed
end

# ── observables[] loader (E6): structured FoldedObservable, eval-free ───────────────────
# {name, every, on:[ExprNode], range:[{weight, value:ExprNode}]} → an :obs row (no eval).
function _load_observables!(net, obs)
    for o in obs
        every = Float64(get(o, "every", Inf))
        on = SampleableValues[to_expr(_attr_node(e)) for e in get(o, "on", [])]
        range = SampleableRange[]
        for r in get(o, "range", [])
            w = Float64(get(r, "weight", 1.0))
            v = to_expr(_attr_node(r["value"]))
            push!(range, (w, v))
        end
        fo = FoldedObservable(range, every, on)
        add_row!(net, :obs; obsName = Symbol(o["name"]), obsOpts = fo)
    end
    return net
end

function obs_to_dict(name, o::FoldedObservable)
    return Dict{String, Any}(
        "name" => string(name),
        "every" => o.every,
        "on" => [node_to_dict(from_expr(e)) for e in o.on],
        "range" => [
            Dict(
                    "weight" => (r isa Tuple ? r[1] : 1.0),
                    "value" => node_to_dict(from_expr(r isa Tuple ? r[2] : r))
                ) for r in o.range
        ],
    )
end

# ── E7: validate — the eval-free pre-load self-check (ADR 0005 §70 + ADR 0007/0008 rules) ──
# A PURE walk over the parsed Dict: NO node_from_dict, NO to_expr, NO eval. An LLM runs this to
# self-check a model before from_json. Returns diagnostics; from_json_model gates on isempty.
struct Diagnostic
    severity::Symbol   # :error | :warn
    path::String
    msg::String
end
Base.string(d::Diagnostic) = "[$(d.severity)] $(d.path): $(d.msg)"

# Walk an ExprNode dict, collecting op/dist/ref/arity diagnostics (rule 1). `allow_sample` governs
# whether a Sample node is legal here (false in a predicate clause value — no RNG, §9.5); `allow_field`
# governs whether a Field (@field) node is legal here (true only in a SetField/@advance value,
# ADR 0008 §D — a Field in a predicate clause would crash at runtime since @field is a macro).
function _validate_node!(
        diags, d, path; places, params, obs, ports = Set{Symbol}(),
        allow_sample = true, allow_field = true
    )
    d isa AbstractDict || return diags        # a bare literal scalar — fine
    tag = get(d, "node", nothing)
    if tag == "const"
        # ok (integrality/range checked by the consuming attribute, not here)
    elseif tag == "ref"
        # the retired `species` kind still validates for one release (ADR 0017 Tier 3)
        kind = Symbol(get(d, "kind", "")) === :species ? :place : Symbol(get(d, "kind", ""))
        nm = Symbol(get(d, "name", ""))
        kind in REF_KINDS || push!(diags, Diagnostic(:error, path, "ref kind $kind ∉ $REF_KINDS"))
        pool = kind === :place ? places : kind === :param ? params : obs
        nm in pool || push!(diags, Diagnostic(:error, path, "ref to undeclared $kind `$nm`"))
    elseif tag == "call"
        op = Symbol(get(d, "op", ""))
        op in OP_WHITELIST || push!(diags, Diagnostic(:error, path, "op $op ∉ OP_WHITELIST"))
        for (i, a) in enumerate(get(d, "args", []))
            _validate_node!(diags, a, "$path.args[$i]"; places, params, obs, ports, allow_sample, allow_field)
        end
    elseif tag == "sample"
        allow_sample || push!(diags, Diagnostic(:error, path, "Sample (RNG) is not 𝓕ₜ-measurable here (no draws in a predicate)"))
        Symbol(get(d, "dist", "")) in DIST_WHITELIST ||
            push!(diags, Diagnostic(:error, path, "dist $(get(d, "dist", "")) ∉ DIST_WHITELIST"))
        for (i, a) in enumerate(get(d, "args", []))
            _validate_node!(diags, a, "$path.args[$i]"; places, params, obs, ports, allow_sample, allow_field)
        end
    elseif tag == "field"
        allow_field || push!(
            diags, Diagnostic(
                :error, path,
                "Field (@field) is legal only in a SetField/@advance value, not here (ADR 0008 §D)"
            )
        )
    elseif tag == "externalref"
        # rule 8 (ADR 0012 §B2): an ExternalRef's port must be a declared inputs[] port. The node
        # is eval-free and 𝓕ₜ-measurable everywhere (it reads the latched buffer, never the RNG),
        # so it is legal in any value context — only an UNDECLARED port is flagged.
        Symbol(get(d, "port", "")) in ports ||
            push!(diags, Diagnostic(:error, path, "ExternalRef port `$(get(d, "port", ""))` is not a declared inputs[] port"))
    elseif tag == "timeref"
        # ok
    elseif tag == "choose"
        for (i, alt) in enumerate(get(d, "alts", []))
            _validate_node!(diags, get(alt, "value", nothing), "$path.alts[$i]"; places, params, obs, ports, allow_sample, allow_field)
        end
    else
        push!(diags, Diagnostic(:error, path, "unknown node tag `$tag`"))
    end
    return diags
end

# Is a JSON attribute value a literal (a number/bool/string or a Const node)? Used for rule 5
# (a TVE=no attribute must be a literal, not a non-trivial tree).
_is_literal(x) = !(x isa AbstractDict) || get(x, "node", "") == "const"

"""
    validate(d::AbstractDict; registry = Dict{Symbol, Any}()) -> Vector{Diagnostic}

Statically check a parsed model dict `d` against the closed, eval-free IR and return a `Vector{Diagnostic}` — it never runs a model field, so it is the trust boundary [`from_json_model`](@ref) gates on (construction proceeds only when the returned vector is empty). Checks include: every [`NodeRef`](@ref) names a declared place/param/observable, every [`Call`](@ref) op is in [`OP_WHITELIST`](@ref) and every [`Sample`](@ref) distribution in [`DIST_WHITELIST`](@ref), reference kinds are in [`REF_KINDS`](@ref), and a non-time-varying attribute is a bare literal rather than a non-trivial tree. Host functions named by an action/genesis `kind` are resolved by name through `registry`. Unexported; call it as `ReactiveDynamics.validate(dict)` to inspect a document directly.
"""
function validate(d::AbstractDict; registry = Dict{Symbol, Any}())
    diags = Diagnostic[]
    # ADR 0017 Tier 3: read either spelling, and report the path the document actually uses.
    place_rows = _legacy_key(d, "places", "species", [])
    arc_rows = _legacy_key(d, "arcs", "reactants", [])
    place_path = _legacy_path(d, "places", "species")
    arc_path = _legacy_path(d, "arcs", "reactants")
    places = Set(Symbol(s["name"]) for s in place_rows)
    structured = Set(Symbol(s["name"]) for s in place_rows if get(s, "structured", false) === true)
    params = Set(Symbol(p["name"]) for p in get(d, "params", []))
    obs = Set(Symbol(o["name"]) for o in get(d, "observables", []))
    # rule 8 (ADR 0012 §B2): the declared external input ports an ExternalRef may reference.
    ports = Set(Symbol(inp["port"]) for inp in get(d, "inputs", []))
    regnames = Set(keys(registry))

    # rule 5 (TVE policy): init/cost-type place attrs must be literals; structured/modality too
    for (i, s) in enumerate(place_rows)
        for k in ("init", "cost", "reward", "valuation")
            haskey(s, k) && !_is_literal(s[k]) &&
                push!(diags, Diagnostic(:error, "$place_path[$i].$k", "must be a literal (TVE=no, §5 A3)"))
        end
        # rule 4: modality must be a legal 5-row combo (§1.4)
        if haskey(s, "modality")
            m = s["modality"]
            try
                to_set(
                    Symbol(get(m, "allocation", "upfront")), Symbol(get(m, "return", "consumed")),
                    Symbol(get(m, "blocking", "block"))
                )
            catch e
                push!(diags, Diagnostic(:error, "$place_path[$i].modality", sprint(showerror, e)))
            end
        end
    end

    # rule 1 + rule 3: transition attr nodes + ranges/integrality
    ids = Set{String}()
    for (i, tr) in enumerate(get(d, "transitions", []))
        push!(ids, string(tr["id"]))
        haskey(tr, "rate") && _validate_node!(diags, tr["rate"], "transitions[$i].rate"; places, params, obs, ports)
        for (k, lo, hi) in (("prob_of_success", 0.0, 1.0),)
            if haskey(tr, k) && _is_literal(tr[k])
                v = tr[k] isa AbstractDict ? get(tr[k], "value", nothing) : tr[k]
                v isa Real && !(lo <= v <= hi) &&
                    push!(diags, Diagnostic(:error, "transitions[$i].$k", "must be in [$lo,$hi], got $v"))
            end
        end
        for k in ("cycletime", "capacity", "max_lifetime")
            if haskey(tr, k) && _is_literal(tr[k])
                v = tr[k] isa AbstractDict ? get(tr[k], "value", nothing) : tr[k]
                v isa Real && v < 0 &&
                    push!(diags, Diagnostic(:error, "transitions[$i].$k", "must be ≥ 0, got $v"))
            end
        end
    end

    # rule 2: dangling arc FKs (+ §9.5 predicate well-formedness, rule 7)
    for (i, r) in enumerate(arc_rows)
        string(get(r, "transition", "")) in ids ||
            push!(diags, Diagnostic(:error, "$arc_path[$i].transition", "dangling FK `$(get(r, "transition", ""))`"))
        pl_key = haskey(r, "place") ? "place" : haskey(r, "species") ? "species" : nothing
        if pl_key !== nothing
            Symbol(r[pl_key]) in places ||
                push!(diags, Diagnostic(:error, "$arc_path[$i].$pl_key", "undeclared place `$(r[pl_key])`"))
        end
        if haskey(r, "predicate")   # rule 7: predicate over a structured kind; clause values 𝓕ₜ-measurable
            pd = r["predicate"]
            Symbol(get(pd, "kind", "")) in structured ||
                push!(diags, Diagnostic(:error, ".predicate.kind", "kind `$(get(pd, "kind", ""))` is not a structured place"))
            for (j, c) in enumerate(get(pd, "clauses", []))
                Symbol(c[2]) in PRED_OP_WHITELIST ||
                    push!(diags, Diagnostic(:error, ".predicate.clauses[$j]", "op `$(c[2])` ∉ PRED_OP_WHITELIST"))
                # the clause VALUE must be 𝓕ₜ-measurable (no Sample)
                c[3] isa AbstractDict &&
                    _validate_node!(diags, c[3], ".predicate.clauses[$j].value"; places, params, obs, ports, allow_sample = false, allow_field = false)
            end
        end
        if haskey(r, "structured")   # named @structured(:Kind, field = value, …) genesis product
            st = r["structured"]
            k = Symbol(get(st, "kind", ""))
            k in regnames ||
                push!(diags, Diagnostic(:error, ".structured.kind", "kind `$k` not in registry"))
            k in structured ||
                push!(diags, Diagnostic(:error, ".structured.kind", "kind `$k` is not a structured place"))
            # field VALUES are genesis attributes: a Sample draw is legal (like an AddToken field),
            # but a bound-token @field read is not (there is no bound token at a genesis product).
            for (j, f) in enumerate(get(st, "fields", []))
                haskey(f, "value") && f["value"] isa AbstractDict &&
                    _validate_node!(diags, f["value"], ".structured.fields[$j].value"; places, params, obs, ports, allow_field = false)
            end
        end
    end

    # rules[] / events[]: guard nodes + action verb + AddToken.kind/Invoke.fn registry resolution
    for (i, r) in enumerate(get(d, "rules", []))
        haskey(r, "guard") && _validate_node!(diags, r["guard"], "rules[$i].guard"; places, params, obs, ports)
        haskey(r, "action") && _validate_action!(diags, r["action"], "rules[$i].action"; places, params, obs, ports, structured, regnames, in_rule = true)
    end

    # rule 6: population[] well-formedness (ADR 0007)
    for (i, pe) in enumerate(get(d, "population", []))
        pl_key = haskey(pe, "place") ? "place" : "species"      # ADR 0017 Tier 3 legacy spelling
        haskey(pe, pl_key) && Symbol(pe[pl_key]) in structured ||
            push!(diags, Diagnostic(:error, "population[$i].$pl_key", "must be a declared structured place"))
        haskey(pe, "kind") && !(Symbol(pe["kind"]) in regnames) &&
            push!(diags, Diagnostic(:error, "population[$i].kind", "kind `$(pe["kind"])` not in registry"))
    end

    return diags
end

# Validate an action statement dict (rule 1 over its values + verb/registry checks). `ports`
# (ADR 0012 §B2) lets rule 8 flag an undeclared ExternalRef inside an action VALUE (e.g. a
# `set_params` value driven by an external signal) — actions are among the "value" contexts §B2
# enumerates.
function _validate_action!(diags, a, path; places, params, obs, ports = Set{Symbol}(), structured, regnames, in_rule)
    a isa AbstractDict || return diags
    # Walk every node-valued field of the statement (rule 1 + rule 8 over its values).
    _walk_action_values!(diags, a, path; places, params, obs, ports)
    verb = get(a, "verb", nothing)
    if verb == "set_field" && in_rule
        push!(diags, Diagnostic(:error, path, "SetField is illegal in a Rule (no bound token, ADR 0010 §C)"))
    elseif verb == "add_token"
        Symbol(get(a, "kind", "")) in regnames ||
            push!(diags, Diagnostic(:error, "$path.kind", "AddToken kind `$(get(a, "kind", ""))` not in registry"))
    elseif verb == "invoke"
        Symbol(get(a, "fn", "")) in regnames ||
            push!(diags, Diagnostic(:error, "$path.fn", "Invoke fn `$(get(a, "fn", ""))` not in registry"))
    elseif verb == "seq"
        for (i, s) in enumerate(get(a, "stmts", []))
            _validate_action!(diags, s, "$path.stmts[$i]"; places, params, obs, ports, structured, regnames, in_rule)
        end
    elseif verb === nothing
        push!(diags, Diagnostic(:error, path, "action missing `verb`"))
    elseif !(Symbol(verb) in ACTION_VERBS || Symbol(verb) in _LEGACY_ACTION_VERBS)
        push!(diags, Diagnostic(:error, path, "unknown action verb `$verb`"))
    end
    return diags
end

# Validate the node-valued fields of an action statement (the `value`/`assigns[].value`/`fields[].
# value`/`args[]` slots). Used by rule 1 (op/dist/ref) and rule 8 (ExternalRef port declared).
# A `seq`'s nested stmts are walked by `_validate_action!` itself; we skip them here.
function _walk_action_values!(diags, a, path; places, params, obs, ports)
    haskey(a, "value") &&
        _validate_node!(diags, a["value"], "$path.value"; places, params, obs, ports)
    for key in ("assigns", "fields")
        for (i, asg) in enumerate(get(a, key, []))
            haskey(asg, "value") &&
                _validate_node!(diags, asg["value"], "$path.$key[$i].value"; places, params, obs, ports)
        end
    end
    for (i, arg) in enumerate(get(a, "args", []))
        _validate_node!(diags, arg, "$path.args[$i]"; places, params, obs, ports)
    end
    return diags
end

# ── to_json helpers ─────────────────────────────────────────────────────────────────────
# Emit only NON-DEFAULT scalar attrs (mirroring defargs in ReactiveDynamics.jl), so the JSON stays
# clean and re-import reconstructs the same value via assign_defaults!. A literal Const lowered by
# the loader is a bare Number here, so we emit the bare number (the loader's _attr_node wraps it).
function _place_to_dict(net, i)
    pl = Dict{String, Any}("name" => string(net[i, :placeName]))
    iv = net[i, :placeInitVal]
    iv isa Number && iv != 0 && (pl["init"] = iv)
    # cost/reward/valuation default to 0.0 (defargs[:S]); emit only when set (TVE=no literals).
    for (col, key) in (:placeCost => "cost", :placeReward => "reward", :placeValuation => "valuation")
        v = net[i, col]
        v isa Number && v != 0 && (pl[key] = v)
    end
    net[i, :placeStructured] && (pl["structured"] = true)
    isempty(net[i, :placeModality]) || (pl["modality"] = modality_to_dict(net[i, :placeModality]))
    return pl
end

# Emit a transition's `id`/`name`, its `rate`(+`rate_mode`), and every non-default attr column
# (cycletime/prob_of_success/capacity/priority/max_lifetime/multiplier) as an ExprNode dict, plus
# pre/post actions when they are typed ActionStmts. The structural inverse of the transitions[]
# loop in build_network_from_dict (~line 98): that loop reads `id` as the arc FK, `rate`+`rate_mode`
# (lowered by lower_rate), and the jsonkey→col attrs (lowered by to_expr); we recover each.
#
# `id` is the transition's `transName` (the FK the BD model.rdj.json uses, e.g. "adv_discovery") —
# falling back to a positional "t<i>" for an unnamed transition so the arc FKs still resolve.
function _transition_to_dict(net, i; places = Set{Symbol}(), params = Set{Symbol}())
    name = net[i, :transName]
    tr = Dict{String, Any}("id" => string(ismissing(name) || isnothing(name) ? Symbol("t", i) : name))
    (ismissing(name) || isnothing(name)) || (tr["name"] = string(name))

    # rate: rate_from_expr Poisson-unwraps the stored :transRate Expr back to (bare-intensity node,
    # rate_mode) — the inverse of lower_rate (~line 190). A bare value ⇒ :deterministic.
    rate_node, rate_mode = rate_from_expr(net[i, :transRate]; places, params)
    tr["rate"] = node_to_dict(rate_node)
    tr["rate_mode"] = string(rate_mode)

    # the ExprNode-valued attrs: emit only when present and DIFFERENT from the construction default
    # (defargs[:T]) so the JSON matches what import expects as the default. Each is lowered back to
    # a node by from_expr (the inverse of to_expr on import).
    for (col, key, default) in (
            (:transCycleTime, "cycletime", 0.0),
            (:transProbOfSuccess, "prob_of_success", 1),
            (:transCapacity, "capacity", Inf),
            (:transPriority, "priority", 1),
            (:transMaxLifeTime, "max_lifetime", Inf),
            (:transMultiplier, "multiplier", 1),
        )
        v = net[i, col]
        (isnothing(v) || _is_default_attr(v, default)) && continue
        tr[key] = node_to_dict(from_expr(v; places, params))
    end

    # pre/post actions: a JSON-authored model holds these as typed ActionStmts (the action family,
    # ADR 0010/0011), which serialize losslessly via stmt_to_dict. A DSL-built model instead stores
    # a raw Expr here (the default empty `:(())` for no action) — that has no typed verb, so it is
    # NOT serializable; we emit it only when it is a typed ActionStmt, and otherwise drop the empty
    # default. A non-empty raw-Expr action would be lost (documented limitation, mirrors RawExpr).
    for (col, key) in (:transPreAction => "pre_action", :transPostAction => "post_action")
        a = net[i, col]
        a isa ActionStmt && (tr[key] = stmt_to_dict(a))
    end
    return tr
end

# Is a stored attribute value the construction default (so it can be omitted)? Numeric compare with
# Inf/Int/Float tolerance; `:(())` (the empty pre/post action) is handled separately above.
_is_default_attr(v, default) = v isa Number && default isa Number && (v == default)

# ── _arcs_to_dict — the EXPORT inverse of assemble_reaction_line (~line 289) ─────────
# assemble_reaction_line turns arcs[] dicts INTO a transition's :trans reaction-line Expr;
# this is the exact inverse: it decomposes each transition's stored :trans Expr back into the flat
# arcs[] list (one dict per arc term, tagged with its transition FK + side). It is the
# highest-risk export piece — it must reproduce exactly the arc dicts the loader consumes.
#
# Method: split the stored `LHS → RHS` line (the same `prune_r_line` split the runtime does, but
# eval-free — we only need the LHS/RHS arms, not the @choose resolution), then walk each arm with
# the runtime `recursive_find_arcs!` (reaction_parser.jl:74). That walker is PURE on a raw
# Expr (no state) and yields the same FoldedArc structs the runtime extracts — place/kind,
# integer-or-Expr stoich, the 3-axis modality Set (from @conserved/@rate/@nonblock wrappers), and
# the @select TokenPredicate. From each FoldedArc we emit the inverse of _arc_atom /
# _apply_modality / the stoich coefficient.
function _arcs_to_dict(net)
    places, params = _name_sets(net)
    out = Dict{String, Any}[]
    for i in row_ids(net, :T)
        name = net[i, :transName]
        id = string(ismissing(name) || isnothing(name) ? Symbol("t", i) : name)
        line = net[i, :trans]
        lhs, rhs = _split_reaction_line(line)
        for r in _static_arcs(lhs)
            push!(out, _arc_to_dict(r, id, "lhs"; places, params))
        end
        for r in _static_arcs(rhs)
            push!(out, _arc_to_dict(r, id, "rhs"; places, params))
        end
    end
    return out
end

# Split a stored reaction line into (LHS, RHS) Exprs. The line is `Expr(:call, arrow, LHS, RHS)`;
# forward arrows keep (args[2], args[3]) and backward arrows flip — the eval-free counterpart of
# prune_r_line (state.jl:192). `∅` (empty_set) arms yield no arcs via _static_arcs.
function _split_reaction_line(line)
    if line isa Expr && line.head == :call && line.args[1] in fwd_arrows
        return line.args[2], line.args[3]
    elseif line isa Expr && line.head == :call && line.args[1] in bwd_arrows
        return line.args[3], line.args[2]
    else
        error(
            "_arcs_to_dict: unexpected reaction line shape $(repr(line)) " *
                "(a @choose/bidirectional line is not yet supported by the export path)"
        )
    end
end

# Walk one arm of the reaction line into FoldedArcs, reusing the runtime parser unchanged.
# An `∅`/`0` arm contributes nothing (recursive_find_arcs! drops it).
_static_arcs(arm) =
    recursive_find_arcs!(arm, 1.0, Set{Symbol}(), Vector{FoldedArc}())

# A single FoldedArc → its arcs[] dict. Three shapes, mirroring _arc_atom's three
# branches in reverse:
#   • a @select LHS  → {side, predicate:{kind, clauses}}  (FoldedArc.predicate ≠ nothing)
#   • an @advance/@structured/@move RHS → {side, advance:{field, value}} / {side, structured/move:…}
#     (FoldedArc.place is a macrocall Expr)
#   • a plain place → {side, place, stoich, modality}
function _arc_to_dict(r::FoldedArc, id, side; places, params)
    d = Dict{String, Any}("transition" => id, "side" => side)
    if r.predicate !== nothing
        # @select(Kind, clauses) — the inverse of _arc_atom's predicate branch. pred_to_dict
        # emits {kind, clauses:[[field, op, value-node], …]} exactly as the loader's predicate{} reads.
        d["predicate"] = pred_to_dict(r.predicate)
    elseif r.place isa Expr && isexpr(r.place, :macrocall)
        _emit_macro_arc!(d, r.place; places, params)
    else
        # a plain place term: name, integer-or-Expr stoich (omit the default 1), 3-axis modality.
        d["place"] = string(r.place)
        _emit_stoich!(d, r.stoich)
        isempty(r.modality) || (d["modality"] = modality_to_dict(r.modality))
    end
    return d
end

# Decompose a structured RHS macrocall (@advance / @structured / @move) back to its JSON form —
# the inverse of _arc_atom's @advance branch and the structured_rhs (solvers.jl:352) shapes.
function _emit_macro_arc!(d, mc::Expr; places, params)
    name = macroname(mc)
    if name === :advance
        # @advance(field, value) — an RHS lifecycle field-write (ADR 0008 §D).
        field = mc.args[3]
        valex = mc.args[4]
        d["advance"] = Dict{String, Any}(
            "field" => string(field isa QuoteNode ? field.value : field),
            "value" => node_to_dict(from_expr(valex; places, params)),
        )
    elseif name === :move
        # @move(from, to) — place relabel (ADR 0006). Emit the two place symbols.
        d["move"] = Dict{String, Any}(
            "from" => string(_macro_sym(mc.args[3])),
            "to" => string(_macro_sym(mc.args[4])),
        )
    elseif name === :structured
        # NAMED form `@structured(:Kind, field = node, …)` — the eval-free-serializable genesis
        # product (ADR 0005 §39), and the ONLY @structured form the engine accepts (the raw
        # constructor form was removed and is rejected at construction, create.jl). Emit
        # {kind, fields:[{name, value-node}, …]}, the inverse of _arc_atom's structured branch;
        # each field value lowers through from_expr exactly like an @advance value or an AddToken
        # field. args[3] is always a QuoteNode here (construction guarantees it).
        d["structured"] = Dict{String, Any}(
            "kind" => string(mc.args[3].value),
            "fields" => [
                Dict{String, Any}(
                        "name" => string(kw.args[1]),
                        "value" => node_to_dict(from_expr(kw.args[2]; places, params)),
                    ) for kw in @view mc.args[4:end]
            ],
        )
    else
        error("_arcs_to_dict: unsupported arc macrocall @$(name)")
    end
    return d
end

_macro_sym(x) = x isa QuoteNode ? x.value : x

# Emit a stoich coefficient, omitting the default 1. The runtime parser carries stoich as a Float
# multiplier (multiplex), so an integer authored as `2` comes back as `2.0`; coerce an integral
# Float back to Int so the re-imported reaction line is the SAME Expr (`2 * X`, not `2.0 * X`).
function _emit_stoich!(d, stoich)
    if stoich isa Number
        (stoich == 1) && return d                      # default coefficient — omit
        s = (stoich isa AbstractFloat && isinteger(stoich)) ? Int(stoich) : stoich
        d["stoich"] = s
    else
        # an expression-valued stoich (rare) — emit as a node; the loader lowers it via to_expr.
        d["stoich"] = node_to_dict(from_expr(stoich))
    end
    return d
end

modality_to_dict(s::Set{Symbol}) = Dict{String, Any}(
    "allocation" => (:rate in s ? "perstep" : "upfront"),
    "return" => (:conserved in s ? "conserved" : "consumed"),
    "blocking" => (:nonblock in s ? "nonblock" : "block"),
)

# ── ADR 0003 Phase 2: populate the promoted ArcSpec incidence table from `:trans` ────────
# Derive `net.arcs` from the authoritative `:trans` column, reusing the SAME eval-free static
# decomposition the JSON exporter uses (`_split_reaction_line` + `_static_arcs`, ~line 844/857),
# so the table is exactly the arc set the runtime/exporter see. Each FoldedArc becomes one
# ArcSpec row: a plain place term gets an integer `place` FK (`find_index` into :S) and its
# static stoich/modality; a DYNAMIC term (a @select predicate, an @advance/@move/@structured/@choose
# macrocall, or a place not found in :S) gets `place = 0` and stashes its term Expr in `expr`
# (the ADR escape-hatch). Idempotent: clears and rebuilds. Lines the static splitter cannot handle
# (a raw @choose or bidirectional arrow at top level) are left un-promoted for that transition —
# their arcs stay Expr-only in `:trans`, which is the escape-hatch at the whole-line grain.
function populate_arcs!(net::ReactionNetwork)
    empty!(net.arcs)
    for t in row_ids(net, :T)
        line = net[t, :trans]
        lhs, rhs = try
            _split_reaction_line(line)
        catch
            continue    # @choose / bidirectional / non-standard line: leave this transition Expr-only
        end
        for (arm, side) in ((lhs, :lhs), (rhs, :rhs))
            for r in _static_arcs(arm)
                if r.predicate !== nothing || (r.place isa Expr)
                    # dynamic: a @select predicate or an @advance/@move/@structured macrocall term.
                    push!(
                        net.arcs, ArcSpec(
                            t, 0, r.stoich, side, r.modality,
                            r.place isa Union{Expr, Symbol} ? r.place : nothing
                        )
                    )
                else
                    pl = r.place isa Symbol ? r.place : Symbol(r.place)
                    j = find_index(pl, net)
                    if j === nothing
                        push!(net.arcs, ArcSpec(t, 0, r.stoich, side, r.modality, pl))
                    else
                        push!(net.arcs, ArcSpec(t, j, r.stoich, side, r.modality, nothing))
                    end
                end
            end
        end
    end
    return net.arcs
end
