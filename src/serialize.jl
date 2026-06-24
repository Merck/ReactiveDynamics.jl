# Single-JSON model serialization (ADR 0005), eval-free. A model is one JSON object with `meta`
# + top-level arrays `params[]`, `species[]`, `transitions[]`, `reactants[]`, `observables[]`,
# `events[]`. Every expression is a node-tagged `ExprNode` dict (never a Julia source string);
# `from_json_model` parses → validates → lowers each node via `to_expr` into the same `Expr`
# columns the `@ReactionNetworkSchema` DSL fills, then constructs a `ReactionNetworkProblem`.
# Uses the existing JSON.jl dependency with hand-rolled node-tagged (de)serialization.

import JSON

export node_to_dict, node_from_dict, model_to_dict, build_acs_from_dict
export from_json_model, to_json_model

# ── ExprNode ⟷ JSON dict (the recursive node-tagged union) ──────────────────────────────
node_to_dict(n::Const) = Dict{String,Any}(
    "node" => "const",
    "value" => n.value isa Symbol ? string(n.value) : n.value,
    # tag a Symbol-valued const so node_from_dict can recover it (vs a string param name)
    "symbol" => n.value isa Symbol,
)
node_to_dict(n::NodeRef) =
    Dict{String,Any}("node" => "ref", "kind" => string(n.kind), "name" => string(n.name))
node_to_dict(n::Call) =
    Dict{String,Any}("node" => "call", "op" => string(n.op), "args" => map(node_to_dict, n.args))
node_to_dict(n::Sample) =
    Dict{String,Any}("node" => "sample", "dist" => string(n.dist), "args" => map(node_to_dict, n.args))
node_to_dict(::TimeRef) = Dict{String,Any}("node" => "timeref")
node_to_dict(n::Choose) = Dict{String,Any}(
    "node" => "choose",
    "alts" => [Dict{String,Any}("weight" => w, "value" => node_to_dict(v)) for (w, v) in n.alts],
)
node_to_dict(n::Field) = Dict{String,Any}("node" => "field", "name" => string(n.name))

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
        return NodeRef(Symbol(d["kind"]), Symbol(d["name"]))
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

# ── Build a ReactionNetworkSchema acset from a parsed model dict (E2: scalar attrs) ─────
# Lowers each ExprNode to the Expr column the constructor consumes. Transitions are assembled
# from reactants[] into the :trans reaction line (E4); for E2 a transition may carry an explicit
# `reaction` string-free node-list, but the minimal path supports params + species + a transition
# whose reactants[] are plain (no modality/predicate) — assembled by assemble_reaction_line (E4).
function build_acs_from_dict(d::AbstractDict; registry = Dict{Symbol,Any}())
    acs = ReactionNetworkSchema()

    # params[] → :P (values are JSON numbers, never eval'd — replaces loadsave.jl:65)
    for pr in get(d, "params", [])
        add_part!(acs, :P; prmName = Symbol(pr["name"]), prmVal = pr["value"])
    end

    # species[] → :S (specInitVal/specCost/… are scalar Const ExprNodes → literals)
    for sp in get(d, "species", [])
        i = add_part!(acs, :S; specName = Symbol(sp["name"]))
        haskey(sp, "init") && (acs[i, :specInitVal] = to_expr(_attr_node(sp["init"])))
        haskey(sp, "cost") && (acs[i, :specCost] = to_expr(_attr_node(sp["cost"])))
        haskey(sp, "reward") && (acs[i, :specReward] = to_expr(_attr_node(sp["reward"])))
        haskey(sp, "valuation") && (acs[i, :specValuation] = to_expr(_attr_node(sp["valuation"])))
        get(sp, "structured", false) === true && (acs[i, :specStructured] = true)
        # modality 3-axis → Set{Symbol} (E6); default empty set = row 1
        haskey(sp, "modality") && (acs[i, :specModality] = modality_from_dict(sp["modality"]))
    end

    # transitions[] → :T. reactants[] for this transition assemble into the :trans reaction line.
    reactants_by_tr = Dict{String,Vector{Any}}()
    for r in get(d, "reactants", [])
        push!(get!(reactants_by_tr, string(r["transition"]), []), r)
    end
    for tr in get(d, "transitions", [])
        id = string(tr["id"])
        line = assemble_reaction_line(get(reactants_by_tr, id, []))
        rate_node = _attr_node(tr["rate"])
        rate_mode = Symbol(get(tr, "rate_mode", "poisson"))
        i = add_part!(
            acs,
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
            haskey(tr, jsonkey) && (acs[i, col] = to_expr(_attr_node(tr[jsonkey])))
        end
    end

    # observables[] (E6) and events[] (E5) are added by their step's loaders.
    haskey(d, "observables") && _load_observables!(acs, d["observables"])

    assign_defaults!(acs)
    return acs
end

# meta[] → keywords. Returns (acs, meta_kwargs, population) for from_json_model.
function _meta_kwargs(d::AbstractDict)
    m = get(d, "meta", Dict{String,Any}())
    kw = Dict{Symbol,Any}()
    for (k, v) in m
        kw[Symbol(k)] = v
    end
    return kw
end

# ── from_json / to_json (the model envelope) ────────────────────────────────────────────
function from_json_model(json::AbstractString; seed = nothing, registry = Dict{Symbol,Any}(), population = [])
    d = JSON.parse(json)
    diags = validate(d; registry = registry)
    isempty(diags) || error("from_json_model: model failed validation:\n" * join(string.(diags), "\n"))
    acs = build_acs_from_dict(d; registry = registry)
    kw = _meta_kwargs(d)
    haskey(kw, :tspan) || error(
        "from_json_model: meta.tspan is required (the simulation horizon) — add e.g. " *
        "\"meta\": { \"tspan\": 100.0, \"dt\": 1.0 } to the model document.",
    )
    seed === nothing && haskey(kw, :seed) && (seed = kw[:seed])
    return ReactionNetworkProblem(
        acs;
        seed = seed,
        registry = registry,
        population = population,
        filter(p -> p.first ∉ (:seed,), kw)...,
    )
end

function model_to_dict(acs::ReactionNetworkSchema; meta = Dict{String,Any}())
    return Dict{String,Any}(
        "rd_format" => "reactive-dynamics-model",
        "version" => "1.0",
        "meta" => meta,
        "params" => [
            Dict{String,Any}("name" => string(acs[i, :prmName]), "value" => acs[i, :prmVal])
            for i in parts(acs, :P) if !isnothing(acs[i, :prmName])
        ],
        "species" => [_species_to_dict(acs, i) for i in parts(acs, :S)],
        "transitions" => [_transition_to_dict(acs, i) for i in parts(acs, :T)],
        "reactants" => _reactants_to_dict(acs),
    )
end

to_json_model(acs::ReactionNetworkSchema; meta = Dict{String,Any}()) =
    JSON.json(model_to_dict(acs; meta = meta))
to_json_model(prob::ReactionNetworkProblem; meta = Dict{String,Any}()) =
    to_json_model(prob.acs; meta = meta)

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
function rate_from_expr(rate; species = Set{Symbol}(), params = Set{Symbol}())
    if rate isa Expr && rate.head == :call && rate.args[1] == :rand
        distcall = rate.args[end]
        if distcall isa Expr && distcall.head == :call && distcall.args[1] == :Poisson
            maxcall = distcall.args[2]                       # max(state.dt * <bare>, 0)
            if maxcall isa Expr && maxcall.head == :call && maxcall.args[1] == :max
                prod = maxcall.args[2]                       # state.dt * <bare>
                if prod isa Expr && prod.head == :call && prod.args[1] == :*
                    bare = prod.args[3]                      # the bare intensity (state.dt is args[2])
                    return from_expr(bare; species, params), :poisson
                end
            end
        end
    end
    return from_expr(rate; species, params), :deterministic
end

# ── Reaction-line assembly (E4) ─────────────────────────────────────────────────────────
# Assemble a transition's reactants[] (grouped lhs/rhs) into the single :trans reaction-line Expr
# `LHS --> RHS` that merge_acs!/the runtime parser consume. A reactant row may carry: `species`
# (or a `predicate` for @select on the LHS), `stoich`, `modality` (3-axis, LHS), or `advance`
# (field+value, an RHS @advance) — assembled into exactly the macrocall Expr shapes the parser
# (reaction_parser.jl / create.jl) expects.
const _LN = LineNumberNode(0, :none)

# The species "atom" of a reactant: a bare species symbol, or a @select(Kind, clauses) macrocall
# (LHS predicate), or a @advance(field, value)/@structured/@move macrocall (RHS).
function _reactant_atom(r::AbstractDict)
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
    else
        return Symbol(r["species"])
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
# The parser unions these macro names into the reactant's modality Set (reaction_parser.jl).
function _apply_modality(atom, m)
    m === nothing && return atom
    s = m isa AbstractDict ? modality_from_dict(m) : m   # Set{Symbol}
    out = atom
    # nest so the innermost wraps the species; order is irrelevant (the parser unions a Set)
    :conserved in s && (out = Expr(:macrocall, Symbol("@conserved"), _LN, out))
    :rate in s && (out = Expr(:macrocall, Symbol("@rate"), _LN, out))
    :nonblock in s && (out = Expr(:macrocall, Symbol("@nonblock"), _LN, out))
    return out
end

# A full reactant term: optional integer stoich coefficient × the (modality-wrapped) atom.
function _reactant_term(r::AbstractDict; lhs::Bool)
    atom = _reactant_atom(r)
    lhs && haskey(r, "modality") && (atom = _apply_modality(atom, r["modality"]))
    stv = to_expr(_attr_node(get(r, "stoich", 1)))
    return (stv == 1 || stv === 1.0) ? atom : Expr(:call, :*, stv, atom)
end

function _sum_terms(terms)
    isempty(terms) && return :∅
    length(terms) == 1 && return terms[1]
    return foldl((a, b) -> Expr(:call, :+, a, b), terms)
end

function assemble_reaction_line(reactants)
    lhs = [_reactant_term(r; lhs = true) for r in reactants if String(r["side"]) == "lhs"]
    rhs = [_reactant_term(r; lhs = false) for r in reactants if String(r["side"]) == "rhs"]
    return Expr(:call, :→, _sum_terms(lhs), _sum_terms(rhs))
end

# ── modality 3-axis ⟷ Set{Symbol} (full bijection + illegal-row rejection is E6) ─────────
function modality_from_dict(m::AbstractDict)
    s = Set{Symbol}()
    get(m, "allocation", "upfront") == "perstep" && push!(s, :rate)
    get(m, "return", "consumed") == "conserved" && push!(s, :conserved)
    get(m, "blocking", "block") == "nonblock" && push!(s, :nonblock)
    return s
end

# ── observables[] loader (full structured form is E6) ───────────────────────────────────
_load_observables!(acs, obs) = acs   # E6

# ── validate (full 7-rule pass is E7; E2 ships a permissive stub so from_json_model works) ──
struct Diagnostic
    severity::Symbol
    path::String
    msg::String
end
Base.string(d::Diagnostic) = "[$(d.severity)] $(d.path): $(d.msg)"
validate(d::AbstractDict; registry = Dict{Symbol,Any}()) = Diagnostic[]   # E7 fills the rules

# ── to_json helpers ─────────────────────────────────────────────────────────────────────
function _species_to_dict(acs, i)
    sp = Dict{String,Any}("name" => string(acs[i, :specName]))
    iv = acs[i, :specInitVal]
    iv isa Number && iv != 0 && (sp["init"] = iv)
    acs[i, :specStructured] && (sp["structured"] = true)
    isempty(acs[i, :specModality]) || (sp["modality"] = modality_to_dict(acs[i, :specModality]))
    return sp
end

function _transition_to_dict(acs, i)
    tr = Dict{String,Any}("id" => string(coalesce(acs[i, :transName], Symbol("t", i))))
    !ismissing(acs[i, :transName]) && (tr["name"] = string(acs[i, :transName]))
    return tr   # rate/attrs round-trip is exercised at the node level; full emit is E8
end

_reactants_to_dict(acs) = Dict{String,Any}[]   # full emit assembled in E8

modality_to_dict(s::Set{Symbol}) = Dict{String,Any}(
    "allocation" => (:rate in s ? "perstep" : "upfront"),
    "return" => (:conserved in s ? "conserved" : "consumed"),
    "blocking" => (:nonblock in s ? "nonblock" : "block"),
)
