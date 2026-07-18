# Reaction-line parser: lower a `:trans` reaction line into its reactant terms. Parts of this file were
# originally adapted from Catalyst.jl.

using MacroTools: postwalk

"""
One parsed reactant term of a reaction line, before it is unfolded against the live state: a `species`
(a name, or an RHS macrocall Expr like `@structured`/`@move`), its `stoich` multiplicity, its `modality`
set (`:nonblock`/`:conserved`/`:rate` plus any custom tags), and an optional `predicate` — a
[`TokenPredicate`](@ref) when the term was written as `@select(Kind, clauses)`, else `nothing` (a plain
name/kind bind). Produced by [`recursive_find_reactants!`](@ref).
"""
struct FoldedReactant
    species::Union{Expr, Symbol}
    stoich::SampleableValues
    modality::Set{Symbol}
    predicate::Any   # nothing, or a TokenPredicate built from a @select(kind, clauses) LHS
end
FoldedReactant(species, stoich, modality) = FoldedReactant(species, stoich, modality, nothing)

"""
Resolve every `@choose(alts…)` in a reaction line to a concrete draw, walking `r_line` bottom-up and
replacing each `@choose` with a weighted `sample_range` pick from the state-owned RNG (§4). A `@choose`
alternative may itself be `(weight, value)` or a bare value, and may nest further `@choose`es (resolved
recursively). Called by [`extract_reactants`](@ref) before the line is unfolded into reactant terms.
"""
function recursively_choose(r_line, state)
    return postwalk(r_line) do ex
        if isexpr(ex, :macrocall) && (macroname(ex) == :choose)
            sample_range(
                [
                    (
                            if isexpr(r, :tuple)
                                (r.args[1], recursively_choose(r.args[2], state))
                        else
                                recursively_choose(r, state)
                        end
                        ) for r in ex.args[3:end]
                ],
                state,
            )
        else
            ex
        end
    end
end

"""
Parse a reaction line `r_line` into its `Vector{FoldedReactant}` against the live `state`: first resolve
any `@choose` ([`recursively_choose`](@ref)), escape indexed species names ([`escape_ref`](@ref)), then
recurse the additive/multiplicative structure ([`recursive_find_reactants!`](@ref)). The per-tick entry
point that turns one side of a `:trans` line into the reactant terms the step loop consumes.
"""
function extract_reactants(r_line, state::ReactionNetworkProblem)
    r_line = recursively_choose(r_line, state)

    return recursive_find_reactants!(
        escape_ref(r_line, state[:, :specName]),
        1.0,
        Set{Symbol}(),
        Vector{FoldedReactant}(undef, 0),
    )
end

"""
Parse a `@select(Kind, clause && clause …)` reactant into a `(Kind::Symbol, TokenPredicate)` pair (ADR
0008 §E). Each clause is `field op value` with `op ∈ PRED_OP_WHITELIST`; the value is kept as a raw
Expr/literal, evaluated later through the seeded closure in `matches`. An empty clause list means "any
token of `Kind`". A non-symbol kind is a hard `error`. Clause collection is done by
[`_collect_pred_clauses!`](@ref).
"""
function parse_token_predicate(ex)
    kind = ex.args[3]
    kind isa Symbol || error("@select: first argument must be a structured kind symbol, got $kind")
    clauses = Clause[]
    if length(ex.args) >= 4
        _collect_pred_clauses!(clauses, ex.args[4])
    end
    return kind, TokenPredicate(kind, clauses)
end

"""
Flatten an `@select` predicate expression into `clauses`, splitting a chain of `&&` into its individual
`field op value` binary clauses (`op ∈ PRED_OP_WHITELIST`). A chained `:comparison` (`a < b < c`) or any
other malformed shape is a hard `error`. The recursion behind [`parse_token_predicate`](@ref).
"""
function _collect_pred_clauses!(clauses, ex)
    if isexpr(ex, :(&&))
        _collect_pred_clauses!(clauses, ex.args[1])
        _collect_pred_clauses!(clauses, ex.args[2])
    elseif isexpr(ex, :call) && length(ex.args) == 3 && ex.args[1] ∈ PRED_OP_WHITELIST
        field = ex.args[2]
        field isa Symbol || error("@select clause LHS must be a token field name, got $field")
        push!(clauses, Clause(field, ex.args[1], ex.args[3]))
    elseif isexpr(ex, :comparison)
        error("@select: chained comparisons unsupported; use && between binary clauses")
    else
        error("@select: malformed clause `$ex` (expected `field op value` with op ∈ $(PRED_OP_WHITELIST))")
    end
    return clauses
end

"""
Recursively decompose one side of a reaction line into [`FoldedReactant`](@ref) terms, accumulating the
running multiplicity `mult` and modality set `mods` as it descends. Distributes `*` into the
multiplicity ([`multiplex`](@ref)) and `+` into separate terms; treats a `@select` as a predicate bind,
an RHS macrocall (`@structured`/`@move`/`@advance`) or plain call as an opaque term, and a modality
macrocall (`@nonblock`/…) as a modality tag wrapping its inner term. A zero or empty-set term is dropped;
a malformed line logs `@error`. Appends to and returns `reactants`.
"""
function recursive_find_reactants!(
        ex::SampleableValues,
        mult::SampleableValues,
        mods::Set{Symbol},
        reactants::Vector{FoldedReactant},
    )
    if typeof(ex) != Expr || isexpr(ex, :.) || (ex.head == :escape)
        if (ex == 0 || in(ex, empty_set))
            return reactants
        else
            push!(reactants, FoldedReactant(recursively_expand_dots(ex), mult, mods))
        end
    elseif ex.args[1] == :*
        recursive_find_reactants!(
            ex.args[end],
            multiplex(mult, ex.args[2:(end - 1)]...),
            mods,
            reactants,
        )
    elseif ex.args[1] == :+
        for i in 2:length(ex.args)
            recursive_find_reactants!(ex.args[i], mult, mods, reactants)
        end
    elseif ex.head == :macrocall && macroname(ex) == :select
        # @select(Kind, clause && clause …) — predicate-based token selection (ADR 0008 §E).
        kind, pred = parse_token_predicate(ex)
        push!(reactants, FoldedReactant(kind, mult, mods, pred))
    elseif isexpr(ex, :call) ||
            (ex.head == :macrocall && macroname(ex) ∈ [:structured, :move, :advance])
        push!(reactants, FoldedReactant(ex, mult, mods))
    elseif ex.head == :macrocall
        mods = copy(mods)
        macroname(ex) in species_modalities && push!(mods, macroname(ex))
        foreach(
            i -> push!(mods, ex.args[i] isa Symbol ? ex.args[i] : ex.args[i].value),
            4:length(ex.args),
        )
        recursive_find_reactants!(ex.args[3], mult, mods, reactants)
    else
        @error("malformed reaction")
    end

    return reactants
end
