using MacroTools: prewalk

# compilers.jl — the expr→closure back end (Catalyst-lineage; NOT the Phase-1 rework).
#
# The Phase-1 rework added the closed, serializable `ExprNode` IR (exprnode.jl); this file is the
# UNCHANGED final stage that consumes what the IR lowers to. `ExprNode.to_expr` emits EXACTLY the
# `Expr` the `@reaction_network` DSL always produced, and [`compile_attrs`](@ref) / [`wrap_expr`](@ref)
# here compile that `Expr` into a `(state, firing)` closure via `eval`. That `eval` is NOT a
# breach of the eval-free rule (ADR 0005/0006): the trust boundary — `validate` + the closed
# `ExprNode`/action whitelist — sits UPSTREAM. By the time an expr reaches this file it is proven
# closed-vocabulary; compiling a proven-closed expr to a callable is the intended terminal step. No
# model FIELD is ever `Meta.parse`/`eval`'d here — only the DSL/IR-produced attribute Exprs are.

"""
Recursively collect every bare `Symbol` (a candidate model variable) appearing in `exs` into `r`.
Descends into `Expr`s, skipping the head/callee of a `:call` (its args 2:end are the operands) so an
operator name is not mistaken for a variable. Used by [`get_contained_params`](@ref) to find which
declared params an attribute expression reads.
"""
function recursively_find_vars!(r, exs...)
    for ex in exs
        if ex isa Symbol
            push!(r, ex)
        else
            (
                ex isa Expr &&
                    recursively_find_vars!(r, ex.args[(!isexpr(ex, :call) ? 1 : 2):end]...)
            )
        end
    end
    return r
end

"""
The declared params (`prms`) actually referenced by expression `ex` — the intersection of the variables
[`recursively_find_vars!`](@ref) collects with the param-name set. [`wrap_expr`](@ref) uses it to emit a
minimal `let` prologue binding only the params a given expression reads.
"""
function get_contained_params(ex, prms)
    r = []
    return recursively_find_vars!(r, ex) ∩ prms
end

"""
Recursively rewrite `ex` in place, replacing every `Symbol` that is a key of `varmap` with its
mapped value. The core name-substitution primitive: [`wrap_expr`](@ref) uses it to point place/param
names at their state-space slots (`S → state.u[i]`, `β → state.p[:β]`), and the composition operators
([`equalize!`](@ref)/`prepend!`) use it to apply a place-rename map. Mutates and returns `ex`.
"""
function recursively_substitute_vars!(varmap, ex)
    if ex isa Symbol
        return haskey(varmap, ex) ? varmap[ex] : ex
    elseif ex isa Expr
        for i in 1:length(ex.args)
            if ex.args[i] isa Expr
                ex.args[i] = recursively_substitute_vars!(varmap, ex.args[i])
            else
                if ex.args[i] isa Symbol && haskey(varmap, ex.args[i])
                    ex.args[i] = varmap[ex.args[i]]
                end
            end
        end
    end

    return ex
end

"""
Recursively fold dotted notation `A.B` into the single flattened place name `Symbol("A.B")`, but
ONLY where that flattened name is a real place (present in `vars`); other dotted expressions (e.g. a
genuine field access) are left untouched. Mutates and returns `ex`. Complements [`escape_ref`](@ref),
which does the analogous folding for indexed `A[1]` place names.
"""
function recursively_expand_dots_in_ex!(ex, vars)
    if isexpr(ex, :.)
        expanded = recursively_expand_dots(ex)
        return if expanded isa Symbol && expanded in vars
            expanded
        else
            recursively_expand_dots_in_ex!.(ex.args, Ref(vars))
            ex
        end
    end
    ex isa Expr && for i in 1:length(ex.args)
        ex.args[i] isa Union{Expr, Symbol} &&
            (ex.args[i] = recursively_expand_dots_in_ex!(ex.args[i], vars))
    end

    return ex
end

# The query-metalanguage macro names `wrap_expr` rewrites into `(state, …)` calls: `@t()` → `t(state)`,
# `@obs(x)` → `obs(state, x)`, etc. Any `@name` in this set is lowered to a plain call on `state`;
# everything else is left for the IR/DSL to have already resolved.
reserved_names = [:t, :obs, :resample, :solverarg, :take, :log, :periodic, :set_params]

"""
Fold an indexed place reference `A[1]` into the single flattened name `Symbol("A[1]")` wherever that
flattened name is a real place (present in `place`), leaving other `:ref` expressions unchanged.
The indexed-name counterpart of [`recursively_expand_dots_in_ex!`](@ref); run first in [`wrap_expr`](@ref)
so array-style place names survive as atomic symbols before variable substitution.
"""
function escape_ref(ex, places)
    return if ex isa Symbol
        ex
    else
        prewalk(
            ex ->
            isexpr(ex, :ref) && Symbol(string(ex)) ∈ places ? Symbol(string(ex)) : ex,
            ex,
        )
    end
end

"""
Compile one attribute expression `fex` into a `(state, firing) -> value` closure — the terminal
expr→callable step every rate/multiplicity/cost/action/guard passes through. A non-expression `fex` (a bare
literal) is returned unchanged.

The rewrite pipeline, in order: fold indexed/dotted place names to atomic symbols
([`escape_ref`](@ref), [`recursively_expand_dots_in_ex!`](@ref)); lower the query metalanguage —
`@t()`/`@obs(x)`/… (see [`reserved_names`](@ref)) become `state`-threaded calls, `@firing`/`@state`
become the closure's own arguments (`@transition` is accepted as a legacy spelling of `@firing`, from
before ADR 0018 named the in-flight instance); then substitute names for state-space slots via `varmap`
(place → `state.u[i]`, param → `state.p[:name]`, read-only). Params actually read by `fex` are bound
in a `let` prologue.

The trailing `eval` builds the closure from an inert, ALREADY-VALIDATED expression: the eval-free trust
boundary (ADR 0005/0006) is `validate` + the closed `ExprNode`/action whitelist upstream, so no model
field is parsed or eval'd here — only a proven-closed DSL/IR Expr is turned into a callable.
"""
function wrap_expr(fex, place_names, prm_names, varmap)
    !isa(fex, Union{Expr, Symbol}) && return fex
    # escape refs in place names: A[1] -> Symbol("A[1]")
    fex = escape_ref(fex, place_names)
    # escape dots in place names: A.B -> Symbol("A.B")
    fex = deepcopy(fex)
    fex = recursively_expand_dots_in_ex!(fex, place_names)

    # prepare the function's body
    letex = :(
        let
        end
    )

    # expression walking (MacroTools): visit each expression, subsitute with the body's return value
    fex = prewalk(fex) do x
        # here we convert the query metalanguage: @t() -> time(state) etc.
        if isexpr(x, :macrocall) && (macroname(x) ∈ reserved_names)
            Expr(:call, macroname(x), :state, x.args[3:end]...)
        elseif isexpr(x, :macrocall) && (macroname(x) ∈ (:firing, :transition))
            :firing
        elseif isexpr(x, :macrocall) && (macroname(x) == :state)
            :state
        else
            x
        end
    end

    # substitute the place names with "pointers" into the state space: S -> state.u[1]
    fex = recursively_substitute_vars!(varmap, fex)
    # substitute the params names with "pointers" into the parameter space: β -> state.p[:β]
    # params can't be mutated!
    foreach(
        v -> push!(letex.args[1].args, :($v = state.p[$(QuoteNode(v))])),
        get_contained_params(fex, prm_names),
    )
    push!(letex.args[2].args, fex)

    # the function shall be a function of the dynamic ReactionNetwork structure: letex -> :(state -> $letex)
    # eval the expression to a Julia function, save that function into the "compiled" network

    return eval(
        quote
            function (state, firing)
                return $letex
            end
        end
    )
end

"""
`true` if attribute `attr` must be carried VERBATIM rather than compiled to a closure: the name/hash
columns (`*Name`, `obs`, `meta`), the raw `trans` reaction-line column (re-parsed per tick, not a
scalar-valued expr), and `placeRole` (a closed `Symbol` tag, ADR 0009 §A). [`compile_attrs`](@ref)
routes these around [`wrap_expr`](@ref).
"""
function skip_compile(attr)
    return any(contains.(Ref(string(attr)), ("Name", "obs", "meta"))) ||
        (string(attr) == "trans") ||
        (attr === :placeRole)          # ADR 0009 §A — a closed Symbol tag, never a compiled expr
end

"""
Compile a static [`ReactionNetwork`](@ref) into the runtime's closure tables — the bridge from the inert
IR store to the executable model consumed by [`ReactionNetworkProblem`](@ref). Every attribute column is
mapped through [`wrap_expr`](@ref) (except the [`skip_compile`](@ref) columns, carried verbatim),
splitting into `attrs` (place/obs/param/meta columns) and `transitions` (the `trans*` columns — the
static transition table, stored as `state.transitions`; the per-tick realized values live separately in
`state.sampled_transitions`), plus
the shared `wrap_fun = ex -> wrap_expr(ex, …)` closure the step loop reuses for per-tick exprs (stashed
as `state.wrap_fun`). Also seeds the runtime-only transition columns: `transActivated` (the latching
per-transition gate, all `true`), `transToSpawn` (pending-spawn counts, zero), `transHash` (per-row
identity), and `transGuard` (the stateless per-tick guard, default `true`; ADR 0010 §B). Returns
`(attrs, transitions, wrap_fun)`.
"""
function compile_attrs(net::ReactionNetwork, structured_token)
    place_names = collect(net[:, :placeName])

    prm_names = collect(net[:, :prmName])
    varmap = Dict([name => :(state.u[$i]) for (i, name) in enumerate(place_names)])
    for name in prm_names
        push!(varmap, name => :(state.p[$(QuoteNode(name))]))
    end
    wrap_fun = ex -> wrap_expr(ex, place_names, prm_names, varmap)
    attrs = Dict{Symbol, Vector}()
    transitions = Dict{Symbol, Vector}()
    for attr in propertynames(net.columns)
        attrs_ = column(net, attr)
        if !contains(string(attr), "trans")
            (
                attrs[attr] = map(
                    i -> skip_compile(attr) ? attrs_[i] : wrap_fun(attrs_[i]),
                    1:length(attrs_),
                )
            )
        else
            (
                transitions[attr] = map(
                    i -> skip_compile(attr) ? attrs_[i] : wrap_fun(attrs_[i]),
                    1:length(attrs_),
                )
            )
        end
    end
    transitions[:transActivated] = fill(true, nrows(net, :T))
    transitions[:transToSpawn] = zeros(nrows(net, :T))
    transitions[:transHash] =
        [coalesce(net[i, :transName], gensym()) for i in row_ids(net, :T)]
    # Stateless per-tick guard (ADR 0010 §B): default `true` (compiled away). AND-ed with the
    # latching transActivated gate in sample_transitions!. Authored via @conditional / the
    # `guard =>` transition attr; carried as a compiled closure like the other trans attrs.
    transitions[:transGuard] = Any[true for _ in row_ids(net, :T)]

    return attrs, transitions, wrap_fun
end

"""
Hoist inline `@register(expr)` computed-value declarations out of a network's attribute expressions into
first-class named observables. Walks every attribute column with [`normalize_pcs!`](@ref), which strips
each `@register` site (rewriting it to a `@take` of a generated observable name) and collects the lifted
declarations, then attaches them via `add_obs!`. Returns a deep copy — the input `net` is not mutated —
so it is called once at construction ([`ReactionNetworkProblem`](@ref)) before [`compile_attrs`](@ref).

(Named for historical reasons around the `@choose`/computed-value machinery; it no longer touches
`@choose`.)
"""
function register_observables(net::ReactionNetwork)
    net = deepcopy(net)
    pcs = []
    for attr in propertynames(net.columns)
        attrs_ = column(net, attr)
        foreach(
            i ->
            !isnothing(attrs_[i]) &&
                attrs_[i] isa Expr &&
                (attrs_[i] = normalize_pcs!(pcs, attrs_[i])),
            1:length(attrs_),
        )
    end

    add_obs!(net, pcs)
    return net
end
