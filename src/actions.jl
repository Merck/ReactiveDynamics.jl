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
# {SetMarking, SetParams, SetField, SetTokens, AddToken, Activate, Deactivate, Invoke, Log, Seq}

export Rule, ActionStmt
export SetMarking, SetParams, SetField, SetTokens, AddToken, Activate, Deactivate, Invoke, Log, Seq
export apply_action!, fire_rules!, activate!, deactivate!, set_guard!

"""
Abstract supertype of the closed, serializable action whitelist (ADR 0010 §C / ADR 0011, CONTRACT §12). The concrete family is `{SetMarking, SetParams, SetField, SetTokens, AddToken, Activate, Deactivate, Invoke, Log, Seq}` (plus the internal `RawExpr` legacy bridge). An `ActionStmt` is a declarative, eval-free record of a state mutation; `apply_action!` dispatches on the concrete type to perform it, and the closed set is the trust boundary for eval-free (de)serialization. Actions are the payload of a `Rule` (the endogenous decision channel) or a transition post-action; `Seq` composes them.
"""
abstract type ActionStmt end

"""
    SetMarking(name, value, mode = :set)

Action (`ActionStmt`) that writes a plain-place pool column of `state.u`: for place `name`, evaluate `value` (an `Expr`/literal, through the seeded closure path) and either set it (`mode = :set`) or increment it (`mode = :inc`).
"""
struct SetMarking <: ActionStmt
    name::Symbol
    value::Any            # Expr / literal, evaluated via context_eval
    mode::Symbol          # :set | :inc
end
SetMarking(name, value) = SetMarking(name, value, :set)

"""
    SetParams(assigns)

Action (`ActionStmt`) that writes one or more model parameters in `state.p`. `assigns` is a vector of `name => value-expr` pairs; each value is evaluated through the seeded closure path.
"""
struct SetParams <: ActionStmt
    assigns::Vector{Pair{Symbol, Any}}   # name => value-expr
end

"""
    SetField(field, value)

Action (`ActionStmt`) that writes `field` on the FIRING transition instance's bound token(s) (ADR 0008 §D), evaluating `value` in the firing context. This is a transition post-action ONLY — a standalone `Rule` has no bound token, so `apply_action!` errors on it; use `SetTokens` for a Rule (validated at build).
"""
struct SetField <: ActionStmt
    field::Symbol
    value::Any
end

"""
    SetTokens(predicate, assigns)

Action (`ActionStmt`) that writes `assigns` (`field => value-expr` pairs) over the token population selected by `predicate` — the population generalization of `SetField` (ADR 0011 §A). Because it carries its own predicate it is legal in a `Rule`. Matched tokens are iterated in the `(place, creation_index)` total order (§9.2), and each value is evaluated IN THE SELECTED TOKEN's context (so `@field(name)` reads that token's own current attribute). `predicate` is a `TokenPredicate` (ADR 0008); a `(kind, clauses)` tuple form is also accepted.
"""
struct SetTokens <: ActionStmt
    predicate::Any
    assigns::Vector{Pair{Symbol, Any}}
end

"""
    AddToken(kind, fields)

Action (`ActionStmt`) that creates a structured token of a registry-registered `kind` (ADR 0006 §C) — the acquisition lever. `fields` is a vector of `name => value-expr` pairs; the values are evaluated, passed to the kind's registered constructor, and the new token is entangled into the network's structured container.
"""
struct AddToken <: ActionStmt
    kind::Symbol
    fields::Vector{Pair{Symbol, Any}}
end

"""
    Activate(transition)

Action (`ActionStmt`) that soft-activates a transition line (ADR 0004 soft gate) by setting its latching `transActivated` flag true. `transition` is matched by `transName`/`transHash`. Lowers to `activate!`.
"""
struct Activate <: ActionStmt
    transition::Symbol
end
"""
    Deactivate(transition)

Action (`ActionStmt`) that soft-deactivates a transition line (ADR 0004 soft gate) by clearing its latching `transActivated` flag. `transition` is matched by `transName`/`transHash`. Lowers to `deactivate!`.
"""
struct Deactivate <: ActionStmt
    transition::Symbol
end

"""
    Invoke(fn, args = Any[])

Action (`ActionStmt`) — the general-code escape hatch (ADR 0011 §B). Lowers to `registry[fn](state, transition, args…)`: the serialized action carries only the NAME `fn`, resolved against the per-network host-function registry (never `eval`'d); the body is trusted host Julia bound by obligations O1–O4. `args` value-exprs are evaluated before the call; the result is discarded.
"""
struct Invoke <: ActionStmt
    fn::Symbol
    args::Vector{Any}
end
Invoke(fn) = Invoke(fn, Any[])

"""
    Log(msg)

Action (`ActionStmt`) that appends `msg` to the state log (via `log`). If `msg` is an `Expr`/`Symbol` it is evaluated in context first, otherwise it is logged as-is.
"""
struct Log <: ActionStmt
    msg::Any
end

"""
    Seq(stmts)

Action (`ActionStmt`) that composes a vector of `ActionStmt`s, applied in order. The composite that lets one `Rule` or transition post-action perform several mutations.
"""
struct Seq <: ActionStmt
    stmts::Vector{ActionStmt}
end

# Legacy/bridge: run a raw action expression (e.g. `B += 100`) through the compiled closure.
# This is how legacy `:E` event actions (authored as `cond && action` blocks) are carried as
# Rules until they are re-expressed with the typed verbs. Side-effecting; result discarded.
struct RawExpr <: ActionStmt
    expr::Any
end

# The SERIALIZED verb tags (ADR 0005 §E5 wire vocabulary), one per concrete `ActionStmt`. This is
# the CURRENT vocabulary; the retired `:set_species` spelling of `:set_marking` is accepted by the
# loader for one release (ADR 0017 Tier 3, `_LEGACY_ACTION_VERBS`) but is not a member here.
const ACTION_VERBS =
    (:set_marking, :set_params, :set_field, :set_tokens, :add_token, :activate, :deactivate, :invoke, :log, :seq)

# Retired ADR-0017 verb spellings a document may still carry for one release.
const _LEGACY_ACTION_VERBS = (:set_species,)

# ── A Rule is the repaired Event (ADR 0010 §A) ──────────────────────────────────────────
"""
    Rule(id, guard, action; fire_mode = :every_tick, enabled = true)

The endogenous decision channel (ADR 0010 §A, CONTRACT §12) — a `(guard, action, fire_mode)` triple evaluated once per tick at a fixed point in `_step!` (step 10, via `fire_rules!`). `guard` is an `Expr`/literal resolving to a `Bool` (fire once when true) or a numeric `v` (fire `rand(rng, Poisson(v))` times, so the RNG-threaded rate stays deterministic); `action` is any `ActionStmt`. `fire_mode ∈ {:every_tick, :once}` — an `:once` rule latches its `enabled` flag off after it first fires. `enabled` is run-state reset by `_reinit!` (§4 D7). Unlike a transition's stateless per-tick `guard`, a Rule has no bound token, so its action must be one of the population/global verbs (`SetTokens`, not `SetField`).
"""
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
    return findfirst(i -> state.network[i, :transName] == t, row_ids(state, :T))
end

"""
    activate!(state, t::Symbol)

Soft-activate the transition line named `t` (matched by `transName`/`transHash`) on a live `state`, by setting its latching `transActivated` gate true (ADR 0004 soft gate). Errors if no transition matches. The imperative twin of the `Activate` action.
"""
function activate!(state::ReactionNetworkProblem, t::Symbol)
    ix = _transition_index(state, t)
    isnothing(ix) && error("activate!: no transition named $t")
    return state.transition_recipes[:transActivated][ix] = true
end

"""
    deactivate!(state, t::Symbol)

Soft-deactivate the transition line named `t` (matched by `transName`/`transHash`) on a live `state`, by clearing its latching `transActivated` gate (ADR 0004 soft gate). Errors if no transition matches. The imperative twin of the `Deactivate` action.
"""
function deactivate!(state::ReactionNetworkProblem, t::Symbol)
    ix = _transition_index(state, t)
    isnothing(ix) && error("deactivate!: no transition named $t")
    return state.transition_recipes[:transActivated][ix] = false
end

# Attach a stateless per-tick guard to a transition (ADR 0010 §B).
"""
    set_guard!(state, t::Symbol, guard)

Attach a stateless per-tick `guard` to the transition named `t` (ADR 0010 §B). `guard` is an `Expr`/literal (e.g. `:(cash >= phase3_cost)`) compiled to the seeded closure and AND-ed with the latching `transActivated` gate in `sample_transitions!` — a transition whose guard evaluates false makes no genesis proposal that tick. Errors if no transition matches.
"""
function set_guard!(state::ReactionNetworkProblem, t::Symbol, guard)
    ix = _transition_index(state, t)
    isnothing(ix) && error("set_guard!: no transition named $t")
    return state.transition_recipes[:transGuard][ix] = state.wrap_fun(guard)
end

# Live-phase guard (ADR 0007 §A): place identification (equalize!) reindexes the :S table via
# rem_parts! (operators/equalize.jl), which would invalidate the construction-frozen, position-
# indexed compiled closures (ADR 0004 INV-2). It is an AUTHORING-only op; calling it on a
# constructed/live ReactionNetworkProblem must refuse rather than corrupt the closures. (The
# authoring-phase equalize!(::ReactionNetwork, …) stays unrestricted.)
function equalize!(state::ReactionNetworkProblem, args...)
    return error(
        "equalize! reindexes the place table (rem_parts!) and is illegal on a live, constructed " *
            "model (ADR 0004 INV-2 / ADR 0007 §A): the position-indexed compiled closures are frozen " *
            "at construction. Identify place at AUTHORING time, before ReactionNetworkProblem(...).",
    )
end

# ── Evaluate an action value in (state, transition) context via the seeded closure path ──
# Mirrors how rate/multiplicity/action exprs are evaluated elsewhere (context_eval + wrap_fun). A bare
# QuoteNode (a literal symbol like `:Phase2` in an action field) is the symbol it wraps —
# wrap_fun/context_eval pass QuoteNodes through unevaluated, so normalize here.
function _eval_value(state::ReactionNetworkProblem, transition, v)
    v isa QuoteNode && return v.value
    r = context_eval(state, transition, state.wrap_fun(v))
    return r isa QuoteNode ? r.value : r
end

# ── apply_action! — the lowering table (ADR 0010 §C / ADR 0011 §C), all eval-free ───────
"""
    apply_action!(state, transition, a::ActionStmt)

Perform the action `a` against `state`, dispatching on the concrete `ActionStmt` type — the eval-free lowering table for the closed action family (ADR 0010 §C / ADR 0011 §C). `transition` is the firing transition instance for a transition post-action, or `nothing` when the action comes from a `Rule`; it carries through to the seeded-closure value eval (params/observables/time/Sample). `SetField` requires a non-`nothing` `transition` (it writes bound tokens); `AddToken`/`Invoke` resolve their name against the per-network registry and never `eval`. `Seq` applies its statements in order.
"""
function apply_action!(state::ReactionNetworkProblem, transition, a::SetMarking)
    ix = find_index(a.name, state)
    isnothing(ix) && error("SetMarking: unknown place $(a.name)")
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
    for tok in transition.bound_tokens
        setproperty!(tok, a.field, _eval_value(state, transition, a.value))
    end
    return nothing
end

function apply_action!(state::ReactionNetworkProblem, transition, a::SetTokens)
    # Iterate matched tokens in the (place, creation_index) total order (§9.2) and write each
    # field. SetTokens is the population generalization of SetField (ADR 0011 §A), so a value is
    # evaluated IN THE SELECTED TOKEN's context via `eval_with_token` (NOT `_eval_value`): that
    # rewrites every `@field(name)` to a literal read of the token's own current attribute before
    # the seeded-closure eval — so `pos_remaining => @field(pos_remaining) * 0.9` writes each token
    # down by 10%. (`@field` is a syntactic marker, not a real macro; sending it through wrap_fun
    # directly would error at macro-expansion.) `transition` is nothing in a Rule, the firing
    # instance in a post-action — it carries through to context_eval for params/obs/time/Sample.
    for tok in select_tokens(state, a.predicate)
        for (f, vex) in a.assigns
            setproperty!(tok, f, eval_with_token(state, transition, tok, vex))
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
    log(state, a.msg isa Union{Expr, Symbol} ? _eval_value(state, transition, a.msg) : a.msg)

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
"""
    fire_rules!(state)

Evaluate every enabled `Rule` in `state.rules` once, in order — the endogenous decision channel, invoked at `_step!` step 10 (ADR 0010 §A/§D). For each rule, the guard is evaluated (`nothing` transition context): a `Bool` fires the action 0/1 times, a numeric `v` fires it `rand(state.rng, Poisson(v))` times (RNG-threaded for determinism). A `:once` rule that fired latches its `enabled` flag off (reset by `_reinit!`, §4 D7). Returns `state`.
"""
function fire_rules!(state::ReactionNetworkProblem)
    for r in state.rules
        r.enabled || continue
        v = _eval_value(state, nothing, r.guard)
        q = v isa Bool ? (v ? 1 : 0) : (v isa Number ? rand(state.rng, Poisson(v)) : 0)
        for _ in 1:q
            apply_action!(state, nothing, r.action)
        end
        # `once` rules latch off after they first fire (run-state, reset by _reinit!, §4 D7).
        (r.fire_mode === :once && q > 0) && (r.enabled = false)
    end
    return state
end
