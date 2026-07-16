# parts of the code were taken from Catalyst.jl and adapted

using MacroTools: postwalk

struct FoldedReactant
    species::Union{Expr, Symbol}
    stoich::SampleableValues
    modality::Set{Symbol}
    predicate::Any   # nothing, or a TokenPredicate built from a @select(kind, clauses) LHS
end
FoldedReactant(species, stoich, modality) = FoldedReactant(species, stoich, modality, nothing)

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

function extract_reactants(r_line, state::ReactionNetworkProblem)
    r_line = recursively_choose(r_line, state)

    return recursive_find_reactants!(
        escape_ref(r_line, state[:, :specName]),
        1.0,
        Set{Symbol}(),
        Vector{FoldedReactant}(undef, 0),
    )
end

# Parse `@select(Kind, clause && clause …)` into (Kind::Symbol, TokenPredicate). Each clause is
# `field op value` where op ∈ PRED_OP_WHITELIST; the value stays a raw Expr/literal (evaluated
# later through the seeded closure in `matches`). An empty clause list ⇒ "any token of Kind".
function parse_token_predicate(ex)
    kind = ex.args[3]
    kind isa Symbol || error("@select: first argument must be a structured kind symbol, got $kind")
    clauses = Clause[]
    if length(ex.args) >= 4
        _collect_pred_clauses!(clauses, ex.args[4])
    end
    return kind, TokenPredicate(kind, clauses)
end

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
