# Endogenous decision channel — the closed action family + the Rule record (ADR 0010/0011,
# CONTRACT §12). An action is a declarative, eval-free statement; a Rule is a (guard, action,
# fire_mode) triple evaluated once per tick at a fixed point in `_step!` (step 10). A transition
# also carries a stateless `guard` AND-ed with its latching `transActivated` gate.
#
# Action VALUES are held as raw `Expr`/literals and evaluated through the existing
# RNG-threaded `context_eval(state, transition, state.wrap_fun(expr))` closure path (§4 D5) —
# the same verified hot path the rest of the engine uses. The typed-ExprNode IR + JSON
# lowering (ADR 0005) is a later stage; until then the action structs ARE the IR.

# ── The closed action type family (ADR 0010 §C, extended by ADR 0011) ───────────────────
# {SetSpecies, SetParams, SetField, SetTokens, AddToken, Activate, Deactivate, Invoke, Log, Seq}

export Rule, ActionStmt
export SetSpecies, SetParams, SetField, SetTokens, AddToken, Activate, Deactivate, Invoke, Log, Seq
export apply_action!, fire_rules!, activate!, deactivate!, set_guard!

abstract type ActionStmt end

# Write a plain-species pool column (`state.u`). mode ∈ {:set, :inc}.
struct SetSpecies <: ActionStmt
    name::Symbol
    value::Any            # Expr / literal, evaluated via context_eval
    mode::Symbol          # :set | :inc
end
SetSpecies(name, value) = SetSpecies(name, value, :set)

# Write one or more model params.
struct SetParams <: ActionStmt
    assigns::Vector{Pair{Symbol,Any}}   # name => value-expr
end

# Write a field on the FIRING transition instance's bound token(s) (ADR 0008 §D).
# Transition post-action ONLY — a standalone Rule has no bound token (validated at build).
struct SetField <: ActionStmt
    field::Symbol
    value::Any
end

# Write a field over a SELECTED token population (ADR 0011 §A). Carries its own predicate,
# so it is legal in a Rule. `predicate` is a TokenPredicate (ADR 0008, Stage C); until Stage C
# lands it may also be a (kind, clauses) tuple — `apply_action!` dispatches on what is present.
struct SetTokens <: ActionStmt
    predicate::Any
    assigns::Vector{Pair{Symbol,Any}}
end

# Create a structured token of a registered kind (ADR 0006 §C); the acquisition lever.
struct AddToken <: ActionStmt
    kind::Symbol
    fields::Vector{Pair{Symbol,Any}}
end

# Toggle a transition line (ADR 0004 soft gate). `transition` is matched by transName/hash.
struct Activate <: ActionStmt
    transition::Symbol
end
struct Deactivate <: ActionStmt
    transition::Symbol
end

# General-code escape hatch (ADR 0011 §B): lowers to `registry[fn](state, transition, args…)`.
# The file carries only the NAME; the body is host Julia (trusted, obligations O1–O4).
struct Invoke <: ActionStmt
    fn::Symbol
    args::Vector{Any}
end
Invoke(fn) = Invoke(fn, Any[])

struct Log <: ActionStmt
    msg::Any
end

# Compose actions; run in order.
struct Seq <: ActionStmt
    stmts::Vector{ActionStmt}
end

# Legacy/bridge: run a raw action expression (e.g. `B += 100`) through the compiled closure.
# This is how legacy `:E` event actions (authored as `cond && action` blocks) are carried as
# Rules until they are re-expressed with the typed verbs. Side-effecting; result discarded.
struct RawExpr <: ActionStmt
    expr::Any
end

const ACTION_VERBS =
    (:set_species, :set_params, :set_field, :set_tokens, :add_token, :activate, :deactivate, :invoke, :log, :seq)

# ── A Rule is the repaired Event (ADR 0010 §A) ──────────────────────────────────────────
# guard resolves to Bool (fire once) or numeric v (fire rand(rng, Poisson(v)) times).
# fire_mode ∈ {:every_tick, :once}; `enabled` is run-state reset by _reinit! (§4 D7).
mutable struct Rule
    id::Symbol
    guard::Any           # Expr / literal, resolves to Bool or numeric
    action::ActionStmt
    fire_mode::Symbol    # :every_tick | :once
    enabled::Bool
end
Rule(id, guard, action; fire_mode = :every_tick, enabled = true) =
    Rule(id, guard, action, fire_mode, enabled)

# ── Transition line toggles (ADR 0004 soft gate; create the missing activate!/deactivate!) ─
# Find a live transition's recipe index by its transName/transHash key.
function _transition_index(state::ReactionNetworkProblem, t::Symbol)
    hashes = state.transition_recipes[:transHash]
    ix = findfirst(==(t), hashes)
    isnothing(ix) || return ix
    # fall back to matching the transName column on the schema
    return findfirst(i -> state.acs[i, :transName] == t, parts(state, :T))
end

function activate!(state::ReactionNetworkProblem, t::Symbol)
    ix = _transition_index(state, t)
    isnothing(ix) && error("activate!: no transition named $t")
    return state.transition_recipes[:transActivated][ix] = true
end

function deactivate!(state::ReactionNetworkProblem, t::Symbol)
    ix = _transition_index(state, t)
    isnothing(ix) && error("deactivate!: no transition named $t")
    return state.transition_recipes[:transActivated][ix] = false
end

# Attach a stateless per-tick guard to a transition (ADR 0010 §B). `guard` is an Expr/literal
# (e.g. `:(cash >= phase3_cost)`); it is compiled to the seeded closure and AND-ed with the
# latching transActivated gate in sample_transitions!. A transition with a false guard makes
# no genesis proposal that tick.
function set_guard!(state::ReactionNetworkProblem, t::Symbol, guard)
    ix = _transition_index(state, t)
    isnothing(ix) && error("set_guard!: no transition named $t")
    return state.transition_recipes[:transGuard][ix] = state.wrap_fun(guard)
end

# Live-phase guard (ADR 0007 §A): species identification (equalize!) reindexes the :S table via
# rem_parts! (operators/equalize.jl), which would invalidate the construction-frozen, position-
# indexed compiled closures (ADR 0004 INV-2). It is an AUTHORING-only op; calling it on a
# constructed/live ReactionNetworkProblem must refuse rather than corrupt the closures. (The
# authoring-phase equalize!(::ReactionNetworkSchema, …) stays unrestricted.)
function equalize!(state::ReactionNetworkProblem, args...)
    return error(
        "equalize! reindexes the species table (rem_parts!) and is illegal on a live, constructed " *
        "model (ADR 0004 INV-2 / ADR 0007 §A): the position-indexed compiled closures are frozen " *
        "at construction. Identify species at AUTHORING time, before ReactionNetworkProblem(...).",
    )
end

# ── Evaluate an action value in (state, transition) context via the seeded closure path ──
# Mirrors how rate/stoich/action exprs are evaluated elsewhere (context_eval + wrap_fun). A bare
# QuoteNode (a literal symbol like `:Phase2` in an action field) is the symbol it wraps —
# wrap_fun/context_eval pass QuoteNodes through unevaluated, so normalize here.
function _eval_value(state::ReactionNetworkProblem, transition, v)
    v isa QuoteNode && return v.value
    r = context_eval(state, transition, state.wrap_fun(v))
    return r isa QuoteNode ? r.value : r
end

# ── apply_action! — the lowering table (ADR 0010 §C / ADR 0011 §C), all eval-free ───────
function apply_action!(state::ReactionNetworkProblem, transition, a::SetSpecies)
    ix = find_index(a.name, state)
    isnothing(ix) && error("SetSpecies: unknown species $(a.name)")
    v = _eval_value(state, transition, a.value)
    if a.mode === :inc
        state.u[ix] += v
    else
        state.u[ix] = v
    end
    return state.u[ix]
end

function apply_action!(state::ReactionNetworkProblem, transition, a::SetParams)
    for (p, vex) in a.assigns
        state.p[p] = _eval_value(state, transition, vex)
    end
    return state.p
end

function apply_action!(state::ReactionNetworkProblem, transition, a::SetField)
    transition === nothing &&
        error("SetField is a transition post-action only — no bound token in a Rule (ADR 0010 §C)")
    for tok in transition.bound_structured_agents
        setproperty!(tok, a.field, _eval_value(state, transition, a.value))
    end
    return nothing
end

function apply_action!(state::ReactionNetworkProblem, transition, a::SetTokens)
    # Iterate matched tokens in the (species, creation_index) total order (§9.2) and write
    # each field. The predicate evaluator + total order arrive with Stage C; this dispatches
    # to `select_tokens` once that exists.
    for tok in select_tokens(state, a.predicate)
        for (f, vex) in a.assigns
            setproperty!(tok, f, _eval_value(state, tok.bound_transition, vex))
        end
    end
    return nothing
end

function apply_action!(state::ReactionNetworkProblem, transition, a::AddToken)
    haskey(state.registry, a.kind) ||
        error("AddToken: kind $(a.kind) not in the network registry (ADR 0006 §C)")
    ctor = state.registry[a.kind]
    fieldvals = Dict(f => _eval_value(state, transition, vex) for (f, vex) in a.fields)
    token = ctor(state, fieldvals)
    add_structured_token!(state, token)
    return token
end

apply_action!(state::ReactionNetworkProblem, transition, a::Activate) =
    activate!(state, a.transition)
apply_action!(state::ReactionNetworkProblem, transition, a::Deactivate) =
    deactivate!(state, a.transition)

function apply_action!(state::ReactionNetworkProblem, transition, a::Invoke)
    haskey(state.registry, a.fn) ||
        error("Invoke: fn $(a.fn) not in the network registry (ADR 0006 §C)")
    args = map(arg -> _eval_value(state, transition, arg), a.args)
    state.registry[a.fn](state, transition, args...)   # return discarded (statement position)
    return nothing
end

apply_action!(state::ReactionNetworkProblem, transition, a::Log) =
    log(state, a.msg isa Union{Expr,Symbol} ? _eval_value(state, transition, a.msg) : a.msg)

function apply_action!(state::ReactionNetworkProblem, transition, a::Seq)
    for s in a.stmts
        apply_action!(state, transition, s)
    end
    return nothing
end

# Evaluate the raw action expression for its side effects (the legacy `:E` event-action path).
apply_action!(state::ReactionNetworkProblem, transition, a::RawExpr) =
    (_eval_value(state, transition, a.expr); nothing)

# ── Fire all rules once (ADR 0010 §A/§D) — invoked at _step! step 10 ────────────────────
function fire_rules!(state::ReactionNetworkProblem)
    for r in state.rules
        r.enabled || continue
        v = _eval_value(state, nothing, r.guard)
        q = v isa Bool ? (v ? 1 : 0) : (v isa Number ? rand(state.rng, Poisson(v)) : 0)
        for _ = 1:q
            apply_action!(state, nothing, r.action)
        end
        # `once` rules latch off after they first fire (run-state, reset by _reinit!, §4 D7).
        (r.fire_mode === :once && q > 0) && (r.enabled = false)
    end
    return state
end
