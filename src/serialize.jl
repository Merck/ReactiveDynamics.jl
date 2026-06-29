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
# ExternalRef (ADR 0012 §B2): a declared inputs[] port read. The JSON carries only the port NAME;
# the foreign-agent topology that fills it lives host-side in add_wire! (Invariant 4, eval-free).
node_to_dict(n::ExternalRef) = Dict{String,Any}("node" => "externalref", "port" => string(n.port))

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

# meta[] → keywords. A few string-valued meta keys are symbolized for backward compatibility.
# `alloc_strategy`/`strategy` are now accepted-and-ignored (ADR 0002 makes priority-weighted
# progressive filling the single allocation policy — there is no longer a :weighted/:greedy
# switch), so symbolizing them is harmless; `schedule` is likewise a legacy no-op key.
const _SYMBOL_META = (:alloc_strategy, :strategy, :schedule)
function _meta_kwargs(d::AbstractDict)
    m = get(d, "meta", Dict{String,Any}())
    kw = Dict{Symbol,Any}()
    for (k, v) in m
        key = Symbol(k)
        kw[key] = (key in _SYMBOL_META && v isa AbstractString) ? Symbol(v) : v
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
    # rules[] (ADR 0010) → typed Rule structs passed to the constructor (the endogenous channel).
    rules = Any[rule_from_dict(r) for r in get(d, "rules", [])]
    # inputs[] (ADR 0012 §B1) → the declared external read ports + their pre-wire defaults. The
    # default seeds `state.external_inputs[port]` at construction so a port read before any wire
    # delivers a value (or in a standalone run with no wires) is still well-defined (§B3).
    external_inputs = inputs_from_dict(get(d, "inputs", []))
    return ReactionNetworkProblem(
        acs;
        seed = seed,
        registry = registry,
        population = population,
        rules = rules,
        external_inputs = external_inputs,
        filter(p -> p.first ∉ (:seed,), kw)...,
    )
end

# The declared species + param NAME sets, typed Set{Symbol} (an empty comprehension would infer
# Set{Any}, which from_expr/rate_from_expr reject). Threaded into every from_expr call so a stored
# attribute Expr's bare symbols classify back to the right NodeRef kind.
function _name_sets(acs::ReactionNetworkSchema)
    species = Set{Symbol}(acs[i, :specName] for i in parts(acs, :S))
    params = Set{Symbol}(acs[i, :prmName] for i in parts(acs, :P) if !isnothing(acs[i, :prmName]))
    return species, params
end

# model_to_dict is the EXPORT envelope — the structural inverse of build_acs_from_dict
# (~line 73). It emits every top-level array build_acs_from_dict reads back: params[], species[],
# transitions[], reactants[], observables[], plus rules[] when the caller passes a constructed
# model's typed Rule vector (a DSL/loaded model → JSON → from_json_model → equivalent model).
# The species/param NAME SETS are threaded into every from_expr call below so a stored attribute
# Expr's bare symbols are classified back to the right NodeRef kind (species vs param), matching
# how the authoring DSL named them — exactly the inverse of to_expr's name→state.u[i]/state.p[:k]
# substitution (ADR 0005 §66).
function model_to_dict(acs::ReactionNetworkSchema; meta = Dict{String,Any}(), rules = [],
                       inputs = Dict{Symbol,Any}())
    species, params = _name_sets(acs)
    d = Dict{String,Any}(
        "rd_format" => "reactive-dynamics-model",
        "version" => "1.0",
        "meta" => meta,
        "params" => [
            Dict{String,Any}("name" => string(acs[i, :prmName]), "value" => acs[i, :prmVal])
            for i in parts(acs, :P) if !isnothing(acs[i, :prmName])
        ],
        "species" => [_species_to_dict(acs, i) for i in parts(acs, :S)],
        "transitions" => [_transition_to_dict(acs, i; species, params) for i in parts(acs, :T)],
        "reactants" => _reactants_to_dict(acs),
    )
    # observables[] (the inverse of _load_observables!, ~line 413) — emit only if any :obs row.
    obs = [obs_to_dict(acs[i, :obsName], acs[i, :obsOpts]) for i in parts(acs, :obs)]
    isempty(obs) || (d["observables"] = obs)
    # rules[] (the endogenous channel, ADR 0010) — a constructed model carries its typed Rules on
    # the ReactionNetworkProblem (prob.rules), so to_json_model(::ReactionNetworkProblem) passes
    # them through here. A bare schema acs has none. Legacy :E event rows are NOT emitted: they are
    # lifted to RawExpr-action Rules at construction (solvers.jl ~line 671), and RawExpr is the
    # non-typed bridge that is intentionally not JSON-serializable (stmt_to_dict(::RawExpr) errors).
    rule_dicts = [rule_to_dict(r) for r in rules if r.action isa ActionStmt && !(r.action isa RawExpr)]
    isempty(rule_dicts) || (d["rules"] = rule_dicts)
    # inputs[] (ADR 0012 §B1) — the inverse of inputs_from_dict (~line 430). Declared external read
    # ports + their pre-wire literal defaults live on the ReactionNetworkProblem (external_input_-
    # defaults), so to_json_model(::ReactionNetworkProblem) passes them through. A default is a
    # literal value (inputs_from_dict requires a Const), so it round-trips as a `const` node. Sorted
    # by port name for deterministic output (the buffer is an unordered Dict).
    input_dicts = [Dict{String,Any}("port" => string(p), "default" => node_to_dict(Const(inputs[p])))
                   for p in sort!(collect(keys(inputs)))]
    isempty(input_dicts) || (d["inputs"] = input_dicts)
    return d
end

to_json_model(acs::ReactionNetworkSchema; meta = Dict{String,Any}()) =
    JSON.json(model_to_dict(acs; meta = meta))
# A constructed model also carries its typed Rules — round-trip them through rules[] — and its
# solver settings (tspan/dt/seed), which live on the ReactionNetworkProblem (not the acs) and merge
# into the meta bag at construction (solvers.jl ~line 544). When the caller supplies no `meta`, we
# reconstruct it from those fields so the exported document is COMPLETE and re-importable on its own
# (from_json_model requires meta.tspan). An explicit `meta` always takes precedence.
function to_json_model(prob::ReactionNetworkProblem; meta = Dict{String,Any}())
    full = _meta_from_prob(prob)
    merge!(full, meta)   # caller-supplied keys win
    return JSON.json(model_to_dict(prob.acs; meta = full, rules = prob.rules,
                                   inputs = prob.external_input_defaults))
end

# Reconstruct the meta bag from a constructed model's solver fields. `tspan` is stored as a
# (t0, tend) tuple; the JSON `tspan` is the horizon `tend` (the scalar from_json_model passes on).
function _meta_from_prob(prob::ReactionNetworkProblem)
    m = Dict{String,Any}("tspan" => float(prob.tspan[2]), "dt" => float(prob.dt))
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

# ── E5: action statement + predicate (de)serialization ──────────────────────────────────
# Action VALUES are ExprNodes in JSON, lowered via to_expr to the Expr the apply_action!/
# _eval_value path consumes (actions.jl). RawExpr is intentionally NOT serializable (the legacy
# bridge for non-typed Exprs — a JSON model uses typed verbs only; ADR 0005 open question).
#
# NOTE (serialize-direction fidelity): the to_dict path lowers a stored action-value Expr back to
# a node via from_expr WITHOUT a species/params context, so a bare symbol classifies to its
# default NodeRef(:species, …). This is runtime-harmless — both NodeRef kinds lower to the same
# bare symbol via to_expr, so a from_json-loaded model is unaffected — and only mislabels the JSON
# `ref.kind` tag when re-serializing a model whose actions were hand-built from raw Exprs. JSON-
# authored actions carry typed nodes and never round-trip through from_expr, so they are exact.
stmt_to_dict(s::SetSpecies) = Dict{String,Any}(
    "verb" => "set_species", "name" => string(s.name),
    "value" => node_to_dict(from_expr(s.value)), "mode" => string(s.mode))
stmt_to_dict(s::SetParams) = Dict{String,Any}(
    "verb" => "set_params",
    "assigns" => [Dict("name" => string(n), "value" => node_to_dict(from_expr(v))) for (n, v) in s.assigns])
stmt_to_dict(s::SetField) = Dict{String,Any}(
    "verb" => "set_field", "field" => string(s.field), "value" => node_to_dict(from_expr(s.value)))
stmt_to_dict(s::SetTokens) = Dict{String,Any}(
    "verb" => "set_tokens", "predicate" => pred_to_dict(s.predicate),
    "assigns" => [Dict("name" => string(n), "value" => node_to_dict(from_expr(v))) for (n, v) in s.assigns])
stmt_to_dict(s::AddToken) = Dict{String,Any}(
    "verb" => "add_token", "kind" => string(s.kind),
    "fields" => [Dict("name" => string(n), "value" => node_to_dict(from_expr(v))) for (n, v) in s.fields])
stmt_to_dict(s::Activate) = Dict{String,Any}("verb" => "activate", "transition" => string(s.transition))
stmt_to_dict(s::Deactivate) = Dict{String,Any}("verb" => "deactivate", "transition" => string(s.transition))
stmt_to_dict(s::Invoke) =
    Dict{String,Any}("verb" => "invoke", "fn" => string(s.fn), "args" => [node_to_dict(from_expr(a)) for a in s.args])
stmt_to_dict(s::Log) = Dict{String,Any}("verb" => "log", "msg" => s.msg isa Union{Expr,Symbol} ? node_to_dict(from_expr(s.msg)) : s.msg)
stmt_to_dict(s::Seq) = Dict{String,Any}("verb" => "seq", "stmts" => [stmt_to_dict(x) for x in s.stmts])
stmt_to_dict(::RawExpr) =
    error("RawExpr is not JSON-serializable (the legacy non-typed bridge) — re-express with typed action verbs")

# A value-expr in an action field arrives as a node dict (typed) or a bare literal.
_stmt_value(x) = to_expr(_attr_node(x))

function stmt_from_dict(d::AbstractDict)
    verb = d["verb"]
    if verb == "set_species"
        return SetSpecies(Symbol(d["name"]), _stmt_value(d["value"]), Symbol(get(d, "mode", "set")))
    elseif verb == "set_params"
        return SetParams([Symbol(a["name"]) => _stmt_value(a["value"]) for a in d["assigns"]])
    elseif verb == "set_field"
        return SetField(Symbol(d["field"]), _stmt_value(d["value"]))
    elseif verb == "set_tokens"
        return SetTokens(pred_from_dict(d["predicate"]),
            [Symbol(a["name"]) => _stmt_value(a["value"]) for a in d["assigns"]])
    elseif verb == "add_token"
        return AddToken(Symbol(d["kind"]),
            [Symbol(a["name"]) => _stmt_value(a["value"]) for a in d["fields"]])
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
pred_to_dict(p::TokenPredicate) = Dict{String,Any}(
    "kind" => string(p.kind),
    "clauses" => [[string(c.field), string(c.op), node_to_dict(from_expr(c.value))] for c in p.clauses])
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
rule_to_dict(r::Rule) = Dict{String,Any}(
    "id" => string(r.id), "guard" => node_to_dict(from_expr(r.guard)),
    "action" => stmt_to_dict(r.action), "fire_mode" => string(r.fire_mode))
rule_from_dict(d::AbstractDict) = Rule(
    Symbol(d["id"]), _stmt_value(d["guard"]), stmt_from_dict(d["action"]);
    fire_mode = Symbol(get(d, "fire_mode", "every_tick")))

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
    seed = Dict{Symbol,Any}()
    for inp in inputs
        port = Symbol(inp["port"])
        if haskey(inp, "default")
            node = _attr_node(inp["default"])
            node isa Const ||
                error("inputs_from_dict: port `$port` default must be a literal value (a bare " *
                      "scalar or a Const node), not a live expression — got $(typeof(node))")
            seed[port] = node.value
        end
    end
    return seed
end

# ── observables[] loader (E6): structured FoldedObservable, eval-free ───────────────────
# {name, every, on:[ExprNode], range:[{weight, value:ExprNode}]} → an :obs row (no eval).
function _load_observables!(acs, obs)
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
        add_part!(acs, :obs; obsName = Symbol(o["name"]), obsOpts = fo)
    end
    return acs
end

function obs_to_dict(name, o::FoldedObservable)
    return Dict{String,Any}(
        "name" => string(name),
        "every" => o.every,
        "on" => [node_to_dict(from_expr(e)) for e in o.on],
        "range" => [Dict("weight" => (r isa Tuple ? r[1] : 1.0),
                         "value" => node_to_dict(from_expr(r isa Tuple ? r[2] : r))) for r in o.range],
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
function _validate_node!(diags, d, path; species, params, obs, ports = Set{Symbol}(),
        allow_sample = true, allow_field = true)
    d isa AbstractDict || return diags        # a bare literal scalar — fine
    tag = get(d, "node", nothing)
    if tag == "const"
        # ok (integrality/range checked by the consuming attribute, not here)
    elseif tag == "ref"
        kind = Symbol(get(d, "kind", ""))
        nm = Symbol(get(d, "name", ""))
        kind in REF_KINDS || push!(diags, Diagnostic(:error, path, "ref kind $kind ∉ $REF_KINDS"))
        pool = kind === :species ? species : kind === :param ? params : obs
        nm in pool || push!(diags, Diagnostic(:error, path, "ref to undeclared $kind `$nm`"))
    elseif tag == "call"
        op = Symbol(get(d, "op", ""))
        op in OP_WHITELIST || push!(diags, Diagnostic(:error, path, "op $op ∉ OP_WHITELIST"))
        for (i, a) in enumerate(get(d, "args", []))
            _validate_node!(diags, a, "$path.args[$i]"; species, params, obs, ports, allow_sample, allow_field)
        end
    elseif tag == "sample"
        allow_sample || push!(diags, Diagnostic(:error, path, "Sample (RNG) is not 𝓕ₜ-measurable here (no draws in a predicate)"))
        Symbol(get(d, "dist", "")) in DIST_WHITELIST ||
            push!(diags, Diagnostic(:error, path, "dist $(get(d,"dist","")) ∉ DIST_WHITELIST"))
        for (i, a) in enumerate(get(d, "args", []))
            _validate_node!(diags, a, "$path.args[$i]"; species, params, obs, ports, allow_sample, allow_field)
        end
    elseif tag == "field"
        allow_field || push!(diags, Diagnostic(:error, path,
            "Field (@field) is legal only in a SetField/@advance value, not here (ADR 0008 §D)"))
    elseif tag == "externalref"
        # rule 8 (ADR 0012 §B2): an ExternalRef's port must be a declared inputs[] port. The node
        # is eval-free and 𝓕ₜ-measurable everywhere (it reads the latched buffer, never the RNG),
        # so it is legal in any value context — only an UNDECLARED port is flagged.
        Symbol(get(d, "port", "")) in ports ||
            push!(diags, Diagnostic(:error, path, "ExternalRef port `$(get(d,"port",""))` is not a declared inputs[] port"))
    elseif tag == "timeref"
        # ok
    elseif tag == "choose"
        for (i, alt) in enumerate(get(d, "alts", []))
            _validate_node!(diags, get(alt, "value", nothing), "$path.alts[$i]"; species, params, obs, ports, allow_sample, allow_field)
        end
    else
        push!(diags, Diagnostic(:error, path, "unknown node tag `$tag`"))
    end
    return diags
end

# Is a JSON attribute value a literal (a number/bool/string or a Const node)? Used for rule 5
# (a TVE=no attribute must be a literal, not a non-trivial tree).
_is_literal(x) = !(x isa AbstractDict) || get(x, "node", "") == "const"

function validate(d::AbstractDict; registry = Dict{Symbol,Any}())
    diags = Diagnostic[]
    species = Set(Symbol(s["name"]) for s in get(d, "species", []))
    structured = Set(Symbol(s["name"]) for s in get(d, "species", []) if get(s, "structured", false) === true)
    params = Set(Symbol(p["name"]) for p in get(d, "params", []))
    obs = Set(Symbol(o["name"]) for o in get(d, "observables", []))
    # rule 8 (ADR 0012 §B2): the declared external input ports an ExternalRef may reference.
    ports = Set(Symbol(inp["port"]) for inp in get(d, "inputs", []))
    regnames = Set(keys(registry))

    # rule 5 (TVE policy): init/cost-type species attrs must be literals; structured/modality too
    for (i, s) in enumerate(get(d, "species", []))
        for k in ("init", "cost", "reward", "valuation")
            haskey(s, k) && !_is_literal(s[k]) &&
                push!(diags, Diagnostic(:error, "species[$i].$k", "must be a literal (TVE=no, §5 A3)"))
        end
        # rule 4: modality must be a legal 5-row combo (§1.4)
        if haskey(s, "modality")
            m = s["modality"]
            try
                to_set(Symbol(get(m, "allocation", "upfront")), Symbol(get(m, "return", "consumed")),
                    Symbol(get(m, "blocking", "block")))
            catch e
                push!(diags, Diagnostic(:error, "species[$i].modality", sprint(showerror, e)))
            end
        end
    end

    # rule 1 + rule 3: transition attr nodes + ranges/integrality
    ids = Set{String}()
    for (i, tr) in enumerate(get(d, "transitions", []))
        push!(ids, string(tr["id"]))
        haskey(tr, "rate") && _validate_node!(diags, tr["rate"], "transitions[$i].rate"; species, params, obs, ports)
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

    # rule 2: dangling reactant FKs (+ §9.5 predicate well-formedness, rule 7)
    for (i, r) in enumerate(get(d, "reactants", []))
        string(get(r, "transition", "")) in ids ||
            push!(diags, Diagnostic(:error, "reactants[$i].transition", "dangling FK `$(get(r,"transition",""))`"))
        if haskey(r, "species")
            Symbol(r["species"]) in species ||
                push!(diags, Diagnostic(:error, "reactants[$i].species", "undeclared species `$(r["species"])`"))
        end
        if haskey(r, "predicate")   # rule 7: predicate over a structured kind; clause values 𝓕ₜ-measurable
            pd = r["predicate"]
            Symbol(get(pd, "kind", "")) in structured ||
                push!(diags, Diagnostic(:error, "reactants[$i].predicate.kind", "kind `$(get(pd,"kind",""))` is not a structured species"))
            for (j, c) in enumerate(get(pd, "clauses", []))
                Symbol(c[2]) in PRED_OP_WHITELIST ||
                    push!(diags, Diagnostic(:error, "reactants[$i].predicate.clauses[$j]", "op `$(c[2])` ∉ PRED_OP_WHITELIST"))
                # the clause VALUE must be 𝓕ₜ-measurable (no Sample)
                c[3] isa AbstractDict &&
                    _validate_node!(diags, c[3], "reactants[$i].predicate.clauses[$j].value"; species, params, obs, ports, allow_sample = false, allow_field = false)
            end
        end
    end

    # rules[] / events[]: guard nodes + action verb + AddToken.kind/Invoke.fn registry resolution
    for (i, r) in enumerate(get(d, "rules", []))
        haskey(r, "guard") && _validate_node!(diags, r["guard"], "rules[$i].guard"; species, params, obs, ports)
        haskey(r, "action") && _validate_action!(diags, r["action"], "rules[$i].action"; species, params, obs, ports, structured, regnames, in_rule = true)
    end

    # rule 6: population[] well-formedness (ADR 0007)
    for (i, pe) in enumerate(get(d, "population", []))
        haskey(pe, "species") && Symbol(pe["species"]) in structured ||
            push!(diags, Diagnostic(:error, "population[$i].species", "must be a declared structured species"))
        haskey(pe, "kind") && !(Symbol(pe["kind"]) in regnames) &&
            push!(diags, Diagnostic(:error, "population[$i].kind", "kind `$(pe["kind"])` not in registry"))
    end

    return diags
end

# Validate an action statement dict (rule 1 over its values + verb/registry checks). `ports`
# (ADR 0012 §B2) lets rule 8 flag an undeclared ExternalRef inside an action VALUE (e.g. a
# `set_params` value driven by an external signal) — actions are among the "value" contexts §B2
# enumerates.
function _validate_action!(diags, a, path; species, params, obs, ports = Set{Symbol}(), structured, regnames, in_rule)
    a isa AbstractDict || return diags
    # Walk every node-valued field of the statement (rule 1 + rule 8 over its values).
    _walk_action_values!(diags, a, path; species, params, obs, ports)
    verb = get(a, "verb", nothing)
    if verb == "set_field" && in_rule
        push!(diags, Diagnostic(:error, path, "SetField is illegal in a Rule (no bound token, ADR 0010 §C)"))
    elseif verb == "add_token"
        Symbol(get(a, "kind", "")) in regnames ||
            push!(diags, Diagnostic(:error, "$path.kind", "AddToken kind `$(get(a,"kind",""))` not in registry"))
    elseif verb == "invoke"
        Symbol(get(a, "fn", "")) in regnames ||
            push!(diags, Diagnostic(:error, "$path.fn", "Invoke fn `$(get(a,"fn",""))` not in registry"))
    elseif verb == "seq"
        for (i, s) in enumerate(get(a, "stmts", []))
            _validate_action!(diags, s, "$path.stmts[$i]"; species, params, obs, ports, structured, regnames, in_rule)
        end
    elseif verb === nothing
        push!(diags, Diagnostic(:error, path, "action missing `verb`"))
    elseif !(Symbol(verb) in ACTION_VERBS)
        push!(diags, Diagnostic(:error, path, "unknown action verb `$verb`"))
    end
    return diags
end

# Validate the node-valued fields of an action statement (the `value`/`assigns[].value`/`fields[].
# value`/`args[]` slots). Used by rule 1 (op/dist/ref) and rule 8 (ExternalRef port declared).
# A `seq`'s nested stmts are walked by `_validate_action!` itself; we skip them here.
function _walk_action_values!(diags, a, path; species, params, obs, ports)
    haskey(a, "value") &&
        _validate_node!(diags, a["value"], "$path.value"; species, params, obs, ports)
    for key in ("assigns", "fields")
        for (i, asg) in enumerate(get(a, key, []))
            haskey(asg, "value") &&
                _validate_node!(diags, asg["value"], "$path.$key[$i].value"; species, params, obs, ports)
        end
    end
    for (i, arg) in enumerate(get(a, "args", []))
        _validate_node!(diags, arg, "$path.args[$i]"; species, params, obs, ports)
    end
    return diags
end

# ── to_json helpers ─────────────────────────────────────────────────────────────────────
# Emit only NON-DEFAULT scalar attrs (mirroring defargs in ReactiveDynamics.jl), so the JSON stays
# clean and re-import reconstructs the same value via assign_defaults!. A literal Const lowered by
# the loader is a bare Number here, so we emit the bare number (the loader's _attr_node wraps it).
function _species_to_dict(acs, i)
    sp = Dict{String,Any}("name" => string(acs[i, :specName]))
    iv = acs[i, :specInitVal]
    iv isa Number && iv != 0 && (sp["init"] = iv)
    # cost/reward/valuation default to 0.0 (defargs[:S]); emit only when set (TVE=no literals).
    for (col, key) in (:specCost => "cost", :specReward => "reward", :specValuation => "valuation")
        v = acs[i, col]
        v isa Number && v != 0 && (sp[key] = v)
    end
    acs[i, :specStructured] && (sp["structured"] = true)
    isempty(acs[i, :specModality]) || (sp["modality"] = modality_to_dict(acs[i, :specModality]))
    return sp
end

# Emit a transition's `id`/`name`, its `rate`(+`rate_mode`), and every non-default attr column
# (cycletime/prob_of_success/capacity/priority/max_lifetime/multiplier) as an ExprNode dict, plus
# pre/post actions when they are typed ActionStmts. The structural inverse of the transitions[]
# loop in build_acs_from_dict (~line 98): that loop reads `id` as the reactant FK, `rate`+`rate_mode`
# (lowered by lower_rate), and the jsonkey→col attrs (lowered by to_expr); we recover each.
#
# `id` is the transition's `transName` (the FK the BD model.rdj.json uses, e.g. "adv_discovery") —
# falling back to a positional "t<i>" for an unnamed transition so the reactant FKs still resolve.
function _transition_to_dict(acs, i; species = Set{Symbol}(), params = Set{Symbol}())
    name = acs[i, :transName]
    tr = Dict{String,Any}("id" => string(ismissing(name) || isnothing(name) ? Symbol("t", i) : name))
    (ismissing(name) || isnothing(name)) || (tr["name"] = string(name))

    # rate: rate_from_expr Poisson-unwraps the stored :transRate Expr back to (bare-intensity node,
    # rate_mode) — the inverse of lower_rate (~line 190). A bare value ⇒ :deterministic.
    rate_node, rate_mode = rate_from_expr(acs[i, :transRate]; species, params)
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
        v = acs[i, col]
        (isnothing(v) || _is_default_attr(v, default)) && continue
        tr[key] = node_to_dict(from_expr(v; species, params))
    end

    # pre/post actions: a JSON-authored model holds these as typed ActionStmts (the action family,
    # ADR 0010/0011), which serialize losslessly via stmt_to_dict. A DSL-built model instead stores
    # a raw Expr here (the default empty `:(())` for no action) — that has no typed verb, so it is
    # NOT serializable; we emit it only when it is a typed ActionStmt, and otherwise drop the empty
    # default. A non-empty raw-Expr action would be lost (documented limitation, mirrors RawExpr).
    for (col, key) in (:transPreAction => "pre_action", :transPostAction => "post_action")
        a = acs[i, col]
        a isa ActionStmt && (tr[key] = stmt_to_dict(a))
    end
    return tr
end

# Is a stored attribute value the construction default (so it can be omitted)? Numeric compare with
# Inf/Int/Float tolerance; `:(())` (the empty pre/post action) is handled separately above.
_is_default_attr(v, default) = v isa Number && default isa Number && (v == default)

# ── _reactants_to_dict — the EXPORT inverse of assemble_reaction_line (~line 289) ─────────
# assemble_reaction_line turns reactants[] dicts INTO a transition's :trans reaction-line Expr;
# this is the exact inverse: it decomposes each transition's stored :trans Expr back into the flat
# reactants[] list (one dict per reactant term, tagged with its transition FK + side). It is the
# highest-risk export piece — it must reproduce exactly the reactant dicts the loader consumes.
#
# Method: split the stored `LHS → RHS` line (the same `prune_r_line` split the runtime does, but
# eval-free — we only need the LHS/RHS arms, not the @choose resolution), then walk each arm with
# the runtime `recursive_find_reactants!` (reaction_parser.jl:74). That walker is PURE on a raw
# Expr (no state) and yields the same FoldedReactant structs the runtime extracts — species/kind,
# integer-or-Expr stoich, the 3-axis modality Set (from @conserved/@rate/@nonblock wrappers), and
# the @select TokenPredicate. From each FoldedReactant we emit the inverse of _reactant_atom /
# _apply_modality / the stoich coefficient.
function _reactants_to_dict(acs)
    species, params = _name_sets(acs)
    out = Dict{String,Any}[]
    for i in parts(acs, :T)
        name = acs[i, :transName]
        id = string(ismissing(name) || isnothing(name) ? Symbol("t", i) : name)
        line = acs[i, :trans]
        lhs, rhs = _split_reaction_line(line)
        for r in _static_reactants(lhs)
            push!(out, _reactant_to_dict(r, id, "lhs"; species, params))
        end
        for r in _static_reactants(rhs)
            push!(out, _reactant_to_dict(r, id, "rhs"; species, params))
        end
    end
    return out
end

# Split a stored reaction line into (LHS, RHS) Exprs. The line is `Expr(:call, arrow, LHS, RHS)`;
# forward arrows keep (args[2], args[3]) and backward arrows flip — the eval-free counterpart of
# prune_r_line (state.jl:192). `∅` (empty_set) arms yield no reactants via _static_reactants.
function _split_reaction_line(line)
    if line isa Expr && line.head == :call && line.args[1] in fwd_arrows
        return line.args[2], line.args[3]
    elseif line isa Expr && line.head == :call && line.args[1] in bwd_arrows
        return line.args[3], line.args[2]
    else
        error("_reactants_to_dict: unexpected reaction line shape $(repr(line)) " *
              "(a @choose/bidirectional line is not yet supported by the export path)")
    end
end

# Walk one arm of the reaction line into FoldedReactants, reusing the runtime parser unchanged.
# An `∅`/`0` arm contributes nothing (recursive_find_reactants! drops it).
_static_reactants(arm) =
    recursive_find_reactants!(arm, 1.0, Set{Symbol}(), Vector{FoldedReactant}())

# A single FoldedReactant → its reactants[] dict. Three shapes, mirroring _reactant_atom's three
# branches in reverse:
#   • a @select LHS  → {side, predicate:{kind, clauses}}  (FoldedReactant.predicate ≠ nothing)
#   • an @advance/@structured/@move RHS → {side, advance:{field, value}} / {side, structured/move:…}
#     (FoldedReactant.species is a macrocall Expr)
#   • a plain species → {side, species, stoich, modality}
function _reactant_to_dict(r::FoldedReactant, id, side; species, params)
    d = Dict{String,Any}("transition" => id, "side" => side)
    if r.predicate !== nothing
        # @select(Kind, clauses) — the inverse of _reactant_atom's predicate branch. pred_to_dict
        # emits {kind, clauses:[[field, op, value-node], …]} exactly as the loader's predicate{} reads.
        d["predicate"] = pred_to_dict(r.predicate)
    elseif r.species isa Expr && isexpr(r.species, :macrocall)
        _emit_macro_reactant!(d, r.species; species, params)
    else
        # a plain species term: name, integer-or-Expr stoich (omit the default 1), 3-axis modality.
        d["species"] = string(r.species)
        _emit_stoich!(d, r.stoich)
        isempty(r.modality) || (d["modality"] = modality_to_dict(r.modality))
    end
    return d
end

# Decompose a structured RHS macrocall (@advance / @structured / @move) back to its JSON form —
# the inverse of _reactant_atom's @advance branch and the structured_rhs (solvers.jl:352) shapes.
function _emit_macro_reactant!(d, mc::Expr; species, params)
    name = macroname(mc)
    if name === :advance
        # @advance(field, value) — an RHS lifecycle field-write (ADR 0008 §D).
        field = mc.args[3]
        valex = mc.args[4]
        d["advance"] = Dict{String,Any}(
            "field" => string(field isa QuoteNode ? field.value : field),
            "value" => node_to_dict(from_expr(valex; species, params)),
        )
    elseif name === :move
        # @move(from, to) — species relabel (ADR 0006). Emit the two species symbols.
        d["move"] = Dict{String,Any}(
            "from" => string(_macro_sym(mc.args[3])),
            "to" => string(_macro_sym(mc.args[4])),
        )
    elseif name === :structured
        # @structured(expr) / @structured(token, species) — a host-built RHS token. The body is
        # host Julia (an Expr), not a typed node, so it cannot round-trip eval-free; surface that.
        error("_reactants_to_dict: @structured RHS carries a host Expr body and is not " *
              "JSON-serializable (use @advance / a typed AddToken rule instead)")
    else
        error("_reactants_to_dict: unsupported reactant macrocall @$(name)")
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

modality_to_dict(s::Set{Symbol}) = Dict{String,Any}(
    "allocation" => (:rate in s ? "perstep" : "upfront"),
    "return" => (:conserved in s ? "conserved" : "consumed"),
    "blocking" => (:nonblock in s ? "nonblock" : "block"),
)
