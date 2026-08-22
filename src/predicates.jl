# Token filtration (ADR 0008, CONTRACT §9.5) — select agentic tokens by an 𝓕ₜ-measurable
# predicate, not kind alone. A `TokenPredicate{kind, clauses}` is the BIND analogue of the
# ADR-0006 `TokenAgg` READ node: it narrows the candidate set a structured LHS arc binds
# (in front of the unchanged priority/creation-index sort + integer take), and it powers the
# population-level `SetTokens` write (ADR 0011). Phase is an ATTRIBUTE (maintainer-canonical):
# one `Project` kind with a `phase` field, `@select(Project, phase==:Phase2)` selects, and
# `@advance(phase,:Phase3)` (a `SetField`) advances — no per-phase kinds.

# A conjunctive clause: read token field `field`, compare with `op` against `value` (an
# Expr/literal evaluated through the seeded closure path, 𝓕ₜ-measurable — no Sample).
struct Clause
    field::Symbol
    op::Symbol
    value::Any
end

struct TokenPredicate
    kind::Symbol
    clauses::Vector{Clause}
end
TokenPredicate(kind) = TokenPredicate(kind, Clause[])

const PRED_OP_WHITELIST = (:(==), :(!=), :(<), :(<=), :(>), :(>=), :in)

# Apply a whitelisted comparison op to (token-field-value, rhs-value).
function _apply_op(op::Symbol, lhs, rhs)
    op === :(==) && return lhs == rhs
    op === :(!=) && return lhs != rhs
    op === :(<) && return lhs < rhs
    op === :(<=) && return lhs <= rhs
    op === :(>) && return lhs > rhs
    op === :(>=) && return lhs >= rhs
    op === :in && return lhs in rhs
    error("predicate op $op not in PRED_OP_WHITELIST")
end

# Read a token field (protocol `place` or a host-struct extra like `phase`/`npv_estimate`),
# guarded so a typo is a clear error rather than AA's silent-false swallow (ADR 0008 §A).
function _token_field(token, field::Symbol)
    field === :species && return get_place(token)
    hasproperty(token, field) ||
        error("token of kind $(get_place(token)) has no field $field (predicate/SetField)")
    return getproperty(token, field)
end

# Evaluate a clause/SetField RHS value through the seeded closure path. A bare QuoteNode
# (a literal `:Phase2` written in a clause) is the symbol it wraps — `wrap_fun`/`context_eval`
# pass QuoteNodes through unevaluated, so normalize here to compare against a Symbol field.
function _eval_pred_value(state, transition, v)
    v isa QuoteNode && return v.value
    r = context_eval(state, transition, state.wrap_fun(v))
    return r isa QuoteNode ? r.value : r
end

# Does `token` satisfy the predicate, evaluated at the allocation-point snapshot (ADR 0008 §C)?
# An empty/nothing predicate matches any token of the right kind (degenerate = today's behavior).
matches(::Nothing, token, state, transition) = true
function matches(pred::TokenPredicate, token, state, transition)
    get_place(token) == pred.kind || return false
    for c in pred.clauses
        lhs = _token_field(token, c.field)
        rhs = _eval_pred_value(state, transition, c.value)
        _apply_op(c.op, lhs, rhs) || return false
    end
    return true
end

# All active (unblocked) tokens of the predicate's kind that match, in the deterministic
# (place, creation_index) total order (ADR 0006 inv 5 / §9.2). Used by population-level
# writes (SetTokens, ADR 0011) where there is no per-transition priority to order by.
function select_tokens(state::ReactionNetworkProblem, pred::TokenPredicate)
    toks = collect(values(inners(getagent(state, "structured"))))
    sel = filter(a -> !isblocked(a) && matches(pred, a, state, nothing), toks)
    return sort!(sel; by = a -> token_sortkey(state, a))
end

# Evaluate a `SetField`/`@advance` value expression with a bound token in scope, so it can read
# the token's own current fields via `@field(name)` (ADR 0008 §D Field leaf). `@field(name)`
# rewrites to a literal read of `getproperty(token, name)` before the usual seeded-closure eval;
# everything else (params, obs, time, arithmetic, Sample draws) goes through context_eval as normal.
function eval_with_token(state::ReactionNetworkProblem, transition, token, valex)
    valex isa QuoteNode && return valex.value
    rewritten = _substitute_fields(valex, token)
    r = context_eval(state, transition, state.wrap_fun(rewritten))
    return r isa QuoteNode ? r.value : r
end

# Replace every `@field(name)` macrocall in `ex` with the token's current field value (a literal),
# so the resulting expression no longer depends on the token at closure-eval time.
function _substitute_fields(ex, token)
    if ex isa Expr
        if isexpr(ex, :macrocall) && macroname(ex) == :field
            name = ex.args[3]
            name = name isa QuoteNode ? name.value : name
            return _token_field(token, Symbol(name))
        end
        return Expr(ex.head, map(a -> _substitute_fields(a, token), ex.args)...)
    end
    return ex
end

# Deterministic total order over tokens: (place, creation_index) — the engine-owned monotonic
# creation index (assigned in add_structured_token!, ADR 0006 §E). uuid is the stable final
# tie-break so the order is total even before any index is assigned.
function token_sortkey(state::ReactionNetworkProblem, a)
    ci = get(state.creation_index, AlgebraicAgents.getname(a), 0)
    return (string(get_place(a)), ci, string(AlgebraicAgents.getname(a)))
end
