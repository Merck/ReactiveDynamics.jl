# Typed, eval-free expression IR (ADR 0005). Every time-varying attribute (rate, stoich,
# cycletime, prob_of_success, priority, …) and action value is authored — by a human or an LLM —
# as a closed tagged-union `ExprNode` tree, NEVER a Julia source string. `to_expr` lowers a tree
# to EXACTLY the `Expr` the existing DSL produces (which `wrap_fun`/`compile_attrs` then compile
# to a (state, transition) closure ONCE at construction); `from_expr` is the structural inverse.
# The runtime hot path is untouched — this is purely the authoring/serialization boundary, so the
# closed whitelists below are the only Julia ever produced from an inert model document.

abstract type ExprNode end

# Scalar leaves and the arithmetic/comparison core.
struct Const <: ExprNode
    value::Union{Float64, Int, Bool, Symbol}   # Symbol ⇒ a literal like :Phase2 (lowers to a QuoteNode)
end
struct NodeRef <: ExprNode
    kind::Symbol   # ∈ REF_KINDS — what `name` refers to
    name::Symbol
end
struct Call <: ExprNode
    op::Symbol     # ∈ OP_WHITELIST
    args::Vector{ExprNode}
end

# Distribution draw, time reference, and weighted choice (the stochastic / dynamic core, E3).
struct Sample <: ExprNode
    dist::Symbol   # ∈ DIST_WHITELIST
    args::Vector{ExprNode}
end
struct TimeRef <: ExprNode end
struct Choose <: ExprNode
    alts::Vector{Tuple{Float64, ExprNode}}   # (weight, value) alternatives
end

# Bound-token field read (ADR 0008 §D) — resolved against the firing instance's bound token at
# apply time; legal only in a SetField/@advance value context (validate rule 7).
struct Field <: ExprNode
    name::Symbol
end

# External-agent read (ADR 0012 §B2) — a NEW closed leaf, sibling of NodeRef/TimeRef. Reads the
# LATCHED value of a declared `inputs[]` port from the per-tick external-input buffer. The port is
# a NAME only; the foreign-agent topology that fills it lives host-side in `add_wire!` calls, never
# in the RD document (Invariant 4, eval-free coupling). `ExternalRef` is itself 𝓕ₜ-measurable and
# RNG-free: it returns the value pinned at `_prestep!` (the source's previous-tick-boundary
# projection), so a rate and a guard reading the same port in one tick agree (Invariant 2). An
# undeclared port is a `validate` rule-8 diagnostic, never an eval.
struct ExternalRef <: ExprNode
    port::Symbol   # a declared inputs[] port name; reads state.external_inputs[port]
end

# Structural equality/hash over the node trees (the default is identity for these non-isbits
# structs), so round-trip tests and dedup compare by value.
Base.:(==)(a::ExprNode, b::ExprNode) =
    typeof(a) === typeof(b) && all(getfield(a, f) == getfield(b, f) for f in fieldnames(typeof(a)))
Base.hash(n::ExprNode, h::UInt) =
    foldr((f, acc) -> hash(getfield(n, f), acc), fieldnames(typeof(n)); init = hash(typeof(n), h))

const OP_WHITELIST = (
    :+, :-, :*, :/, :^, :>, :<, :>=, :<=, :(==), :!=, :&&, :||, :!,
    :min, :max, :exp, :log, :floor, :ceil, :abs,
)
const DIST_WHITELIST =
    (:Poisson, :Binomial, :Normal, :Uniform, :Exponential, :Bernoulli, :LogNormal, :Gamma, :Beta)
const REF_KINDS = (:species, :param, :obs)
# Short-circuit boolean ops use an Expr HEAD (not a :call); the rest are ordinary calls.
const _SHORTCIRCUIT_OPS = (:&&, :||)

export ExprNode, Const, NodeRef, Call, Sample, TimeRef, Choose, Field, ExternalRef
export OP_WHITELIST, DIST_WHITELIST, REF_KINDS, to_expr, from_expr

# ── to_expr: lower a typed node to the exact Expr the DSL/compilers consume ──────────────
to_expr(n::Const) = n.value isa Symbol ? QuoteNode(n.value) : n.value
to_expr(n::NodeRef) = n.name   # a bare Symbol; wrap_fun substitutes species→state.u[i], param→state.p[:name]

function to_expr(n::Call)
    n.op in OP_WHITELIST || error("Call op $(n.op) ∉ OP_WHITELIST")
    if n.op in _SHORTCIRCUIT_OPS
        length(n.args) == 2 || error("$(n.op) takes exactly 2 args, got $(length(n.args))")
        return Expr(n.op, to_expr(n.args[1]), to_expr(n.args[2]))
    end
    return Expr(:call, n.op, map(to_expr, n.args)...)
end

function to_expr(n::Sample)
    n.dist in DIST_WHITELIST || error("Sample dist $(n.dist) ∉ DIST_WHITELIST")
    # Matches context_eval's Sampleable path + the explicit-rng spawn shape (create.jl): the draw
    # routes through the state-owned RNG (§4 D5).
    return Expr(:call, :rand, :(state.rng), Expr(:call, n.dist, map(to_expr, n.args)...))
end

# @t() → rewritten to t(state) by wrap_fun (compilers.jl reserved-name rewrite).
to_expr(::TimeRef) = Expr(:macrocall, Symbol("@t"), LineNumberNode(0, :none))

# @choose((w1, v1), (w2, v2), …) — the macrocall shape recursively_choose/prune_r_line consume.
to_expr(n::Choose) = Expr(
    :macrocall,
    Symbol("@choose"),
    LineNumberNode(0, :none),
    (Expr(:tuple, w, to_expr(v)) for (w, v) in n.alts)...,
)

# @field(name) — the bound-token field read (predicates.jl _substitute_fields rewrites it).
to_expr(n::Field) = Expr(:macrocall, Symbol("@field"), LineNumberNode(0, :none), QuoteNode(n.name))

# ExternalRef(port) → a read of the per-tick external-input buffer (ADR 0012 §B2). Lowers to the
# literal index `state.external_inputs[:port]`. wrap_fun/compile_attrs leave this `state.<field>`
# access untouched (it is not a species/param name, so neither the varmap substitution nor the
# ref/dot-escaping passes rewrite it — verified against compilers.jl), so the compiled (state,
# transition) closure reads the buffer `_prestep!` filled this tick. The `state.dt`/`state.rng`
# accesses the engine already lowers elsewhere are the precedent for a bare `state.<field>` leaf.
to_expr(n::ExternalRef) = :(state.external_inputs[$(QuoteNode(n.port))])

# ── from_expr: structural inverse, classifying bare symbols via the known species/param sets ──
# Used to lower a DSL-authored attribute Expr back to a typed tree (for to_json of a DSL model).
function from_expr(ex; species::Set{Symbol} = Set{Symbol}(), params::Set{Symbol} = Set{Symbol}())
    if ex isa Bool
        return Const(ex)
    elseif ex isa Union{Int, Float64}
        return Const(ex)
    elseif ex isa QuoteNode
        return Const(ex.value)   # a literal symbol
    elseif ex isa Symbol
        ex in params && return NodeRef(:param, ex)
        ex in species && return NodeRef(:species, ex)
        return NodeRef(:species, ex)   # default: an unclassified bare symbol is a species ref
    elseif ex isa Expr
        return _from_expr_compound(ex; species = species, params = params)
    else
        error("from_expr: cannot classify $(ex) :: $(typeof(ex))")
    end
end

# Does `ex` match the ExternalRef lowering `state.external_inputs[:port]`?
# i.e. Expr(:ref, Expr(:., :state, QuoteNode(:external_inputs)), <port literal>).
function _is_external_input_ref(ex::Expr)
    ex.head === :ref && length(ex.args) == 2 || return false
    dot = ex.args[1]
    return dot isa Expr && dot.head === :. && dot.args[1] === :state &&
        dot.args[2] isa QuoteNode && dot.args[2].value === :external_inputs
end

function _from_expr_compound(ex::Expr; species, params)
    # short-circuit boolean head: Expr(:&&, a, b) / Expr(:||, a, b)
    if ex.head in _SHORTCIRCUIT_OPS
        return Call(ex.head, [from_expr(a; species, params) for a in ex.args])
    elseif ex.head == :ref && _is_external_input_ref(ex)
        # state.external_inputs[:port] ⇒ ExternalRef(:port) — the structural inverse of the
        # ADR 0012 §B2 lowering, so a stored attr Expr round-trips back to the typed leaf.
        return ExternalRef(ex.args[2] isa QuoteNode ? ex.args[2].value : Symbol(ex.args[2]))
    elseif ex.head == :call
        head = ex.args[1]
        # rand(state.rng, Dist(args...))  OR  rand(Dist(args...))  ⇒ Sample (E3 unwrapping)
        if head == :rand
            distcall = ex.args[end]
            if distcall isa Expr && distcall.head == :call && distcall.args[1] in DIST_WHITELIST
                return Sample(distcall.args[1], [from_expr(a; species, params) for a in distcall.args[2:end]])
            end
        end
        head in OP_WHITELIST || error("from_expr: call head $head ∉ OP_WHITELIST")
        return Call(head, [from_expr(a; species, params) for a in ex.args[2:end]])
    elseif ex.head == :macrocall
        m = ex.args[1]
        if m == Symbol("@t")
            return TimeRef()
        elseif m == Symbol("@field")
            nm = ex.args[end]
            return Field(nm isa QuoteNode ? nm.value : Symbol(nm))
        elseif m == Symbol("@choose")
            alts = Tuple{Float64, ExprNode}[]
            for a in ex.args[3:end]
                a isa Expr && a.head == :tuple ||
                    error("from_expr: @choose alt must be a (weight, value) tuple")
                push!(alts, (Float64(a.args[1]), from_expr(a.args[2]; species, params)))
            end
            return Choose(alts)
        end
        error("from_expr: unsupported macrocall $m")
    else
        error("from_expr: unsupported Expr head $(ex.head)")
    end
end
