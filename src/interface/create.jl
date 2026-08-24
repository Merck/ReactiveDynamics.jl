# reaction network DSL: CREATE part; reaction line and event parsing

export @reaction_network
# Deprecated compatibility alias (ADR 0015 Tier 1) — exported so existing caller code keeps
# working (with a depwarn) for one release; removed in a follow-up.
export @ReactionNetworkSchema
export @append_transitions

using MacroTools: prewalk, postwalk, striplines, isexpr
using Symbolics: build_function, get_variables

empty_set = Set{Symbol}([:∅])
fwd_arrows = Set{Symbol}([:>, :→, :↣, :↦, :⇾, :⟶, :⟼, :⥟, :⥟, :⇀, :⇁, :⇒, :⟾])
bwd_arrows =
    Set{Symbol}([:<, :←, :↢, :↤, :⇽, :⟵, :⟻, :⥚, :⥞, :↼, :↽, :⇐, :⟽, Symbol("<--")])
double_arrows = Set{Symbol}([:↔, :⟷, :⇄, :⇆, :⇌, :⇋, :⇔, :⟺, Symbol("<-->")])

arrows = fwd_arrows ∪ bwd_arrows ∪ double_arrows ∪ [:-->]
ifs = [:&&, :if]

reserved_sampling_macros = [:register, :sample, :take]

struct ArcTerm
    place::Symbol
    multiplicity::SampleableValues
    modality::Set{Symbol}
end

struct FoldedReactionStruct
    rate::SampleableValues
    reaction::SampleableValues
end

struct Event
    trigger::SampleableValues
    action::SampleableValues
end

# Declares symbols which may neither be used as parameters not varriables.
forbidden_symbols = [:t, :π, :pi, :ℯ, :im, :nothing, :∅]

"""
Macro that takes an expression corresponding to a reaction network and outputs a `ReactionNetwork`
(the static struct-of-columns model), ready to pass to `ReactionNetworkProblem` for simulation.

Most arrows accepted (both right, left, and bi-drectional arrows). Use 0 or ∅ for annihilation/creation to/from nothing.

Custom functions and sampleable objects can be used as numeric parameters. Note that these have to be accessible from ReactiveDynamics's source code.

# Examples

```julia
net = @reaction_network begin
    1.0, X ⟶ Y
    1.0, X ⟶ Y, priority => 6.0, prob => 0.7, capacity => 3.0
    1.0, ∅ --> (Poisson(0.3γ)X, Poisson(0.5)Y)
    (XY > 100) && (XY -= 1)
end
@push net 1.0 X ⟶ Y
@prob_init net X = 1 Y = 2 XY = α
@prob_params net γ = 1 α = 4
```
"""
macro reaction_network end

macro reaction_network()
    return make_ReactionNetwork(:())
end

macro reaction_network(ex)
    return make_ReactionNetwork(ex; eval_module = __module__)
end

macro reaction_network(ex, args...)
    # Multiple positional statements (a plain multi-statement authoring call): fold them into one
    # block expression and parse. The GeneratedExpressions brace-comprehension pass was retired
    # (ADR 0015 Tier 3) — the block is handed to the parser directly.
    return make_ReactionNetwork(Expr(:block, ex, args...); eval_module = __module__)
end

# Deprecated compatibility shim (ADR 0015 Tier 1): the macro was renamed `@ReactionNetworkSchema`
# → `@reaction_network`. Forward to the new macro after emitting a depwarn.
macro ReactionNetworkSchema(args...)
    Base.depwarn(
        "`@ReactionNetworkSchema` is deprecated, use `@reaction_network` instead.",
        Symbol("@ReactionNetworkSchema"),
    )
    return esc(
        Expr(
            :macrocall, GlobalRef(ReactiveDynamics, Symbol("@reaction_network")),
            __source__, args...
        )
    )
end

"""
Build the expression that constructs a [`ReactionNetwork`](@ref) from a `@reaction_network` body `ex`:
flatten the block ([`unblock_shallow!`](@ref)) and splice the parsed `(transitions, arcs, obs,
events)` from [`get_data`](@ref) into the constructor call. The core of the [`@reaction_network`](@ref)
macro.
"""
function make_ReactionNetwork(ex::Expr; eval_module = @__MODULE__)
    blockex = unblock_shallow!(ex)

    return :(ReactionNetwork(get_data($(QuoteNode(blockex)))...))
end

### Functions that process the input and rephrase it as a reaction system ###

"""
Recursively replace each interpolation `\$(x)` in `ex` with an `esc`aped reference, so a value from the
macro call site is spliced into the generated network unchanged. Mutates and returns `ex`.
"""
function esc_dollars!(ex)
    if ex isa Expr
        if ex.head == :$
            return esc(:($(ex.args[1])))
        else
            for i in 1:length(ex.args)
                ex.args[i] = esc_dollars!(ex.args[i])
            end
        end
    end
    return ex
end

# Turn a binary-op Expr `a op b` into the pair `a => b` (used for `key => value` attribute terms);
# a bare number passes through unchanged.
symbolize(pairex) = pairex isa Number ? pairex : (pairex.args[2] => pairex.args[3])

"""
Parse a `@reaction_network`/`@push` body `ex` into `(transitions, arcs, pcs, events)` — the tuple
the network constructor consumes. Splits a `begin … end` block line-by-line (or handles a single line),
dispatching each to [`get_data!`](@ref). `pcs` collects the lifted `@register` computed-value declarations.
"""
function get_data(ex)
    trans = []
    evs = []
    arcs = []
    pcs = []
    if isexpr(ex, :block)
        ex = striplines(ex)
        esc_dollars!(ex)
        foreach(
            l -> get_data!(trans, arcs, pcs, evs, isexpr(l, :tuple) ? l.args : [l]),
            ex.args,
        )
    elseif ex != :()
        get_data!(trans, arcs, pcs, evs, ex)
    end

    return trans, arcs, pcs, evs
end

"""
Route one parsed line into the right bucket: an `if`/conditional line becomes an event
([`get_events!`](@ref)), anything else a transition ([`get_transitions!`](@ref)). Accumulates into the
caller's `trans`/`arcs`/`pcs`/`evs` collections. Called per line by [`get_data`](@ref).
"""
function get_data!(trans, arcs, pcs, evs, exs)
    length(exs) == 0 && return

    return if exs[1] isa Expr && (exs[1].head ∈ ifs)
        get_events!(evs, normalize_pcs!(pcs, exs[1]))
    else
        get_transitions!(trans, arcs, pcs, exs)
    end
end

"""
Append the event(s) from a conditional line to `evs`: a `trigger && action` line becomes one
[`Event`](@ref) directly, while an `if/elseif/else` chain is expanded into one guarded `Event` per branch
(with the negated preceding conditions AND-ed in) via [`recursively_expand_actions!`](@ref).
"""
get_events!(evs, ex) =
if ex.head == :&&
    push!(evs, Event(ex.args...))
else
    recursively_expand_actions!(evs, Expr(:call, :&), ex)
end

"""
Expand an `if/elseif/else` chain `event` into one guarded [`Event`](@ref) per branch, threading the
accumulated condition `condex`: each branch fires under (all previous conditions negated) AND its own
condition, and the recursion carries the negation forward into the `else`. The event-side counterpart of
lowering a conditional. Appends to `evs`.
"""
function recursively_expand_actions!(evs, condex, event)
    return if isexpr(event, :if)
        condex_ = deepcopy(condex)
        push!(condex_.args, event.args[1])
        push!(evs, Event(condex_, event.args[2]))
        push!(condex.args, Expr(:call, :!, event.args[1]))
        length(event.args) >= 3 && recursively_expand_actions!(evs, condex, event.args[3])
    else
        push!(evs, Event(condex, event))
    end
end

"""
Lower a transition's authored rate into its genesis-intensity Expr: by default wrap it as a Poisson draw
from the state-owned RNG (`rand(state.rng, Poisson(state.dt * rate))`, §4), UNLESS it is `@deterministic(r)`
(then use `r` verbatim). Also rewrites an inline `@ct(x)` cycle-time macro to `1/x`. Called by
[`get_transitions!`](@ref).
"""
function expand_rate(rate)
    rate = if !(isexpr(rate, :macrocall) && (macroname(rate) == :deterministic))
        # Genesis intensity draws from the state-owned RNG (§4 D2/D5); `state` is in scope
        # because this Expr is compiled into a (state, transition) closure by wrap_expr.
        :(rand(state.rng, Poisson(max(state.dt * $rate, 0))))
    else
        rate.args[3]
    end

    return postwalk(rate) do ex
        if (isexpr(ex, :macrocall) && (macroname(ex) ∈ prettynames[:transCycleTime]))
            :(1 / $(ex.args[3]))
        else
            ex
        end
    end
end

"""
Parse one transition line `exs` (`rate, reaction_line, key => value…`) into `trans`: prune the reaction
line into arc terms ([`prune_reaction_line!`](@ref)), lower the rate ([`expand_rate`](@ref)), and
fold the trailing attributes (resolving their pretty-name aliases against `prettynames`, defaults from
`defargs[:T]`) into the transition's attribute dict. Appends `(rate, rxs) => args` entries.
"""
function get_transitions!(trans, arcs, pcs, exs)
    args = empty(defargs[:T])

    (rate, r_line) = exs[1:2]
    rxs = prune_reaction_line!(pcs, arcs, r_line)
    rate = expand_rate(rate)
    rxs = rxs isa Tuple ? tuple.(fill(rate, length(rxs)), rxs) : ((rate, rxs),)

    exs = exs[3:end]
    empty!(args)
    ix = 1
    while ix <= length(exs)
        (!isa(exs[ix], Expr) || (exs[ix].head != :call)) && (ix += 1; continue)
        karg = (
            xi = findfirst(k -> exs[ix].args[2] ∈ k, prettynames);
            isnothing(xi) && (ix += 1; continue);
            xi
        )
        push!(args, karg => normalize_pcs!(pcs, exs[ix].args[3]))
        deleteat!(exs, ix)
    end
    args = merge(defargs[:T], args)

    append!(trans, tuple.(rxs, Ref(args)))
    return trans
end

"""
Substitute sub-expressions throughout `expr` by the `old => new` `pairs`, replacing every occurrence
(via `prewalk`). Used to resolve a bidirectional-arrow reaction line into its forward/backward variants.
"""
function replace_in_expr(expr, pairs...)
    dict = Dict(pairs...)

    return prewalk(ex -> haskey(dict, ex) ? dict[ex] : ex, expr)
end

"""
Hoist inline `@register(expr)` computed-value declarations out of `expr` into the `pcs` accumulator,
rewriting each site to a `@take` of a generated name so the transition reads the registered value. The
per-expression worker behind [`register_observables`](@ref). Mutates `expr`/`pcs` and returns the walked
`expr`.
"""
function normalize_pcs!(pcs, expr)
    return postwalk(expr) do ex
        isexpr(ex, :macrocall) &&
            macroname(ex) == :register &&
            (
            push!(pcs, deepcopy(ex));
            ex.args[1] = Symbol("@", :take);
            ex.args = ex.args[1:3]
        )
        if isexpr(ex, :macrocall) && macroname(ex) == :register
            r_sym = gensym()
            (
                push!(pcs, (ex_ = deepcopy(ex); insert!(ex_.args, 3, r_sym); ex_));
                ex.args[1] = Symbol("@", :take);
                ex.args = [
                    r_sym
                    ex.args[1:2]
                ]
            )
        end
        return ex
    end
end

"""
Normalize a reaction `line` and collect its places into `arcs`: rewrite `-->` to the canonical
`→`, split a bidirectional `⟷` line into its forward/backward pair, and descend into the LHS/RHS to
register arc names ([`recursively_find_arcs!`](@ref)). Also threads the `pcs` computed-value
accumulator through. Returns the normalized line(s). Called by [`get_transitions!`](@ref).
"""
function prune_reaction_line!(pcs, arcs, line)
    line isa Expr &&
        (line.head == :-->) &&
        (line = Expr(:call, :→, line.args[1], line.args[2]))
    if isexpr(line, :macrocall)
        lines = copy(line.args[3:end])
        line.args = line.args[1:2]
        for l in lines
            biarrow = nothing
            prewalk(ex -> (ex ∈ double_arrows && (biarrow = ex); ex), l)
            append!(
                line.args,
                if isnothing(biarrow)
                    [l]
                else
                    [replace_in_expr(l, biarrow => :⟶), replace_in_expr(l, biarrow => :⟵)]
                end,
            )
        end
        for i in 3:length(line.args)
            line.args[i] = if isexpr(line.args[i], :tuple)
                Expr(
                    :tuple,
                    line.args[i].args[1],
                    prune_reaction_line!(pcs, arcs, line.args[i].args[2]),
                )
            else
                prune_reaction_line!(pcs, arcs, line.args[i])
            end
        end
    elseif line isa Expr && line.args[1] ∈ union(fwd_arrows, bwd_arrows)
        line.args[2:3] =
            recursively_find_arcs!.(Ref(arcs), Ref(pcs), line.args[2:3])
    elseif line isa Expr && line.args[1] ∈ double_arrows
        biarrow = nothing
        prewalk(ex -> (ex ∈ double_arrows && (biarrow = ex); ex), line)
        line = prune_reaction_line!.(
            Ref(pcs),
            Ref(arcs),
            (replace_in_expr(line, biarrow => :⟶), replace_in_expr(line, biarrow => :⟵)),
        )
    end

    return line
end

"""
Walk one side of a reaction line at AUTHORING time and register each place NAME into `arcs`
(distributing `*`/`+`, and descending into `@choose` alternatives). RHS macrocalls are handled specially:
`@structured` is validated to the named `(:Kind, field=value…)` form and left intact (the raw
constructor form is rejected so the IR stays eval-free), `@move`/`@advance` pass through, and `@select`
registers only its KIND as a place (its clause fields are token attributes, not place). The
authoring-time twin of the runtime [`recursive_find_arcs!`](@ref) in reaction_parser.jl.
"""
function recursively_find_arcs!(arcs, pcs, ex)
    if typeof(ex) != Expr || isexpr(ex, :.) || (ex.head == :escape)
        if (ex == 0 || in(ex, empty_set))
            return :∅
        else
            push!(arcs, recursively_expand_dots(ex))
        end
    elseif ex.args[1] == :*
        recursively_find_arcs!(arcs, pcs, ex.args[end])
        foreach(i -> ex.args[i] = normalize_pcs!(pcs, ex.args[i]), 2:(length(ex.args) - 1))
    elseif ex.args[1] == :+
        for i in 2:length(ex.args)
            recursively_find_arcs!(arcs, pcs, ex.args[i])
        end
    elseif isexpr(ex, :macrocall) && macroname(ex) == :choose
        for i in 3:length(ex.args)
            recursively_find_arcs!(
                arcs,
                pcs,
                isexpr(ex.args[i], :tuple) ? ex.args[i].args[2] : ex.args[i],
            )
        end
    elseif isexpr(ex, :macrocall) && macroname(ex) == :structured
        # @structured(:Kind, field = value, …) — the named, registry-resolved genesis product is
        # the ONLY supported form (ADR 0005 §39 / 0006 §C): args[3] is the quoted kind symbol, the
        # rest are `field = value` pairs. The raw `@structured(Ctor(…))` form (an inline host
        # constructor Expr) was REMOVED — it was the sole arc construct that could not
        # round-trip through the eval-free JSON IR, so forbidding it makes eval-free serialization a
        # TOTAL invariant (every genesis product is data). Reject the raw form at construction time.
        (length(ex.args) >= 3 && ex.args[3] isa QuoteNode) || error(
            "@structured: expected the named form `@structured(:Kind, field = value, …)` — a " *
                "quoted kind symbol whose constructor is resolved through the network registry by " *
                "name. The raw `@structured(Ctor(…))` constructor form is not supported (it cannot " *
                "serialize eval-free); register the kind's `(state, fields) -> token` constructor and " *
                "reference it by name instead (ADR 0006 §C). Got: $(ex.args[3])"
        )
        return ex
    elseif isexpr(ex, :macrocall) && macroname(ex) ∈ [:move, :advance]
        return ex
    elseif isexpr(ex, :macrocall) && macroname(ex) == :select
        # @select(Kind, clauses): register only the KIND as a place; the clause fields
        # (phase, npv, …) are token attributes, NOT place, so they must not be registered.
        push!(arcs, ex.args[3])
    elseif isexpr(ex, :macrocall)
        pass_value = ex.args[3] isa QuoteNode ? ex.args[3].value : ex.args[3]
        recursively_find_arcs!(arcs, pcs, pass_value)
    elseif isexpr(ex, :call)
        push!(arcs, ex.args[1])
    else
        push!(arcs, underscorize(ex))
    end

    return ex
end

"""
    @append_transitions net transitions

Append a runtime-built collection of reaction lines to an existing `net`. `transitions` evaluates to a collection of strings, each one reaction line in the [`@reaction_network`](@ref) surface syntax; they are joined into a single `begin…end` block, parsed, and handed to `@push`. Use this when the set of transitions is assembled programmatically (a vector built in a loop, read from a table) rather than written literally — the literal-authoring path is `@push`.

# Examples

```julia
lines = ["ν * I, I --> R, name => I2R", "γ, R --> S, name => R2S"]
@append_transitions net lines
```
"""
macro append_transitions(network, transitions)
    return quote
        transitions_expr = """
        begin
            $(join($(esc(transitions)), '\n'))
        end""" |> Meta.parseall |> striplines

        push_expr = quote
            @push $($(esc(network))) $(transitions_expr.args[1])
        end

        Base.eval(@__MODULE__, push_expr)
    end
end
