# model joins
export merge_networks!, @join

using MacroTools
using MacroTools: prewalk

"""
    merge_networks!(net1, net2, name = gensym("net"), eqs = []) -> ReactionNetwork

Merge `net2` into `net1` IN PLACE and return `net1`. `net2` is deep-copied and its place namespaced
under `name` (via [`prepend!`](@ref)) before merging, so the two fragments' private place stay
distinct; a place already present in `net1` is identified by name and its attribute cells overwritten
from `net2` (later fragment wins), with modality sets unioned. Transitions, params, metadata, events
(`:E`) and observables (`:obs`) are all carried across — `:E`/`:obs` STRUCTURALLY (appended, never
deduplicated). The `eqs` equation blocks drive place identification across fragments (see
[`normalize_name`](@ref)). The engine behind [`@join`](@ref).
"""
function merge_networks!(net1, net2, name = gensym("net"), eqs = [])
    net2 = deepcopy(net2)
    prepend!(net2, name, eqs)

    for i in row_ids(net2, :S)
        inc = find_rows(net1, net2[i, :placeName], :placeName)

        if isempty(inc)
            inc = add_row!(net1, :S; placeName = net2[i, :placeName])
            assign_defaults!(net1)
        end

        union!(net1[first(inc), :placeModality], net2[i, :placeModality])

        for attr in propertynames(net1.columns)
            !occursin("place", string(attr)) && continue
            !ismissing(net2[i, attr]) && (net1[first(inc), attr] = net2[i, attr])
        end
    end

    new_trans_ix = add_rows!(net1, :T, nrows(net2, :T))
    for attr in propertynames(net2.columns)
        !occursin("trans", string(attr)) && continue
        for (ix1, ix2) in enumerate(new_trans_ix)
            net1[ix2, attr] = net2[ix1, attr]
        end
    end

    foreach(
        i -> (
            net1[i, :transName] =
                normalize_name(Symbol(coalesce(net1[i, :transName], i)), name)
        ),
        new_trans_ix,
    )

    for i in row_ids(net2, :P)
        inc = find_rows(net1, net2[i, :prmName], :prmName)
        isempty(inc) && (inc = add_row!(net1, :P; prmName = net2[i, :prmName]))
        !ismissing(net2[i, :prmVal]) && (net1[first(inc), :prmVal] = net2[i, :prmVal])
    end

    for i in row_ids(net2, :M)
        inc = find_rows(net1, net2[i, :metaKeyword], :metaKeyword)
        isempty(inc) && (inc = add_row!(net1, :M; metaKeyword = net2[i, :metaKeyword]))
        !ismissing(net2[i, :metaVal]) && (net1[first(inc), :metaVal] = net2[i, :metaVal])
    end

    # Events (:E) and observables (:obs) are STRUCTURAL — appended, never deduplicated (like :T,
    # §7/J2). The historic gap (merge_networks! walked only :S/:T/:P/:M) silently DROPPED both on join
    # (determinism_composition_bugs.jl §"merge_networks! does NOT merge observables/events"). `prepend!`
    # (above) already namespaced the place referenced inside each event's trigger/action Expr and
    # inside each observable's option-Exprs (via prepend_obs), so both merges are pure structural
    # copies of already-namespaced rows.
    for i in row_ids(net2, :E)
        add_row!(
            net1,
            :E;
            eventTrigger = net2[i, :eventTrigger],
            eventAction = net2[i, :eventAction],
        )
    end

    for i in row_ids(net2, :obs)
        add_row!(net1, :obs; obsName = net2[i, :obsName], obsOpts = net2[i, :obsOpts])
    end

    return net1
end

# Deprecated ACSets-vocabulary alias (ADR 0015 Tier 2): `union_acs!` → `merge_networks!`.
@deprecate union_acs!(net1, net2, name = gensym("net"), eqs = []) merge_networks!(net1, net2, name, eqs)

"""
Namespace `net`'s place in place: rename each `X → name__X` and rewrite every reference to it across
all attribute columns (and, structurally, inside observable option Exprs via [`prepend_obs!`](@ref)), so
merging two fragments cannot conflate their private place. A `:shared`-role place (the first-class
`@catchall`, ADR 0009 §A / CONTRACT §11.1) is identified by BARE name and left un-namespaced; `:private`
(default) and the open `:input`/`:output` ports namespace here, with [`@compose`](@ref) re-identifying
the open ports afterwards by FK-repoint. `eqs` drives cross-fragment identification via
[`normalize_name`](@ref). Called by [`merge_networks!`](@ref) before it copies rows across.
"""
function prepend!(net::ReactionNetwork, name = gensym("net"), eqs = [])
    specmap = Dict()
    for i in row_ids(net, :S)
        # ADR 0009 §A / CONTRACT §11.1: a `shared`-role place is identified by BARE name across all
        # fragments (the first-class @catchall) — it is NOT namespaced. `private` (default) and the
        # open `input`/`output` ports namespace as usual here; @compose (§E) re-identifies the open
        # ports afterwards by FK-repoint. (A place carrying no role reads :private via port_role.)
        if port_role(net, i) === :shared
            continue
        end
        new_name = normalize_name(name, i, net[i, :placeName], eqs)
        push!(specmap, net[i, :placeName] => (net[i, :placeName] = new_name))
    end

    for attr in propertynames(net.columns)
        attr == :placeName && continue
        # Observable options live inside a FoldedObservable struct (the :obsOpts column), not as a
        # bare Expr the loop below rewrites — handle them structurally via prepend_obs! so place
        # referenced inside `on`/`range` exprs are namespaced consistently with every other attr.
        attr == :obsOpts && continue
        attr_ = net[:, attr]
        for i in eachindex(attr_)
            attr_[i] = escape_ref(attr_[i], collect(keys(specmap)))
            attr_[i] = recursively_substitute_vars!(specmap, attr_[i])
            net[i, attr] = attr_[i]
        end
    end

    for i in row_ids(net, :obs)
        prepend_obs!(net[i, :obsOpts], specmap)
    end

    return net
end

"""
Namespace the place referenced inside an observable's option expressions.

`prepend!` renames every place `X → parent__X` and records the map in `specmap`. An observable's
sampling triggers (`on`) and range endpoints (`range`) are stored as Exprs inside a `FoldedObservable`
(the `:obsOpts` column), which `prepend!`'s attribute loop skips — so without this the observable would
still reference the pre-namespaced place and silently read the wrong (or a missing) pool after a join.
This mirrors the per-attribute `escape_ref` + `recursively_substitute_vars!` rewrite `prepend!` applies
to every other spec-referencing attribute. The observable's own NAME (`obsName`) is intentionally left
un-namespaced — rate/guard exprs reference observables by bare name via `@obs(x)`.
"""
function prepend_obs!(opts::FoldedObservable, specmap)
    keys_ = collect(keys(specmap))
    subst(ex) = recursively_substitute_vars!(specmap, escape_ref(ex, keys_))
    opts.on = map(subst, opts.on)
    opts.range = map(
        r -> r isa Tuple ? (r[1], subst(r[2])) : subst(r),
        opts.range,
    )
    return opts
end

## place name normalization
normalize_name(name::Symbol, parent_name) = Symbol("$(parent_name)__$name")
normalize_name(name::String, parent_name) = "$(parent_name)__$name"
normalize_name(name, parent_name) = Symbol(parent_name, "__", name)

"""
The namespaced name for the `i`-th place (named `name`) of the fragment `parent`, honoring the
identification blocks in `eqs`: if the place is named by a block — by exact `:S` index, by a
`:catchall` name match, or by a `parent`-qualified name match — it collapses to that block's alias (its
`:alias` entry, else a generated `shared_place_N`); otherwise it namespaces to `parent__name`. This is
what lets `@equalize`/`@join` fuse place across fragments. Used by [`prepend!`](@ref).
"""
function normalize_name(parent, i::Int, name::Symbol, eqs = [])
    for (block_ix, block) in enumerate(eqs)
        block_alias = findfirst(e -> e[1] == :alias, block)
        block_alias = if !isnothing(block_alias)
            block[block_alias][2]
        else
            Symbol(:shared_place_, block_ix)
        end
        for e in block
            (
                (i == e[2]) ||
                    (
                    e[1] == :catchall &&
                        (normalize_name(e[2], parent) == normalize_name(name, parent))
                ) ||
                    (
                    e[1] == parent &&
                        (normalize_name(e[2], parent) == normalize_name(name, parent))
                )
            ) && return block_alias
        end
    end

    return normalize_name(name, parent)
end

# The two names a place may match after a join: its bare `name` and its namespaced `parent__name`.
matching_name(name::Symbol, parent_name) = [name, Symbol("$(parent_name)__$name")]

# Parse one side of a `@join`/`@equalize` equation into `(qualifier, name)` tuples: a dotted `net.X`
# reconstructs to `(net, X)`, an `@alias(X)` yields `(:alias, X)`, any other macrocall yields both a
# `:catchall` and an `:alias` entry, and a bare symbol is a `:catchall`.
expand_name(ex) =
if isexpr(ex, :.)
    reconstruct(ex)
elseif isexpr(ex, :macrocall)
    (
        if macroname(ex) == :alias
            [(:alias, ex.args[3])]
        else
            [(:catchall, ex.args[3]), (:alias, ex.args[3])]
        end
    )
else
    (:catchall, ex)
end

# Flatten a dotted access `a.b.c` into the symbol list `[a, b, c]` (the leaves of the `:.` Expr tree).
function recursively_get_syms(ex)
    return isexpr(ex, :.) ? [recursively_get_syms(ex.args[1]); ex.args[2].value] : ex
end
# Reconstruct a dotted name into a `(parent, name)` pair: `net.X` → `(net, :X)`, `net.A.B` →
# `(net, :A__B)`; a bare symbol stays a 1-tuple. The dotted-syntax counterpart of `expand_name`.
function reconstruct(ex)
    return if ex isa Symbol
        (ex,)
    else
        (syms = recursively_get_syms(ex); (syms[1], Symbol(join(syms[2:end], "__"))))
    end
end

"""
Flatten one `@join` equation (a possibly-chained `A = B = C`) into its full list of `(qualifier, name)`
pairs via [`expand_name`](@ref), recursing through the right-nested `:(=)` Exprs. Consumed by
[`merge_eqs!`](@ref) to build the identification blocks the [`@join`](@ref) macro passes on.
"""
function get_eqs(eq)
    return if isexpr(eq, :macrocall)
        expand_name(eq)
    elseif eq isa Expr
        [
            expand_name(eq.args[1])
            isexpr(eq.args[2], :(=)) ? get_eqs(eq.args[2]) : expand_name(eq.args[2])
        ]
    else
        [eq]
    end
end

# Merge `eqblock` into the accumulated identification blocks `eqs`, coalescing any existing blocks that
# share a member with it into one (so `A=B` then `B=C` fuse into a single `{A,B,C}` block). Mutates `eqs`.
function merge_eqs!(eqs, eqblock)
    eqs_ = []
    for s in eqblock
        ix = findfirst(e -> s in e, eqs)
        !isnothing(ix) && push!(eqs_, ix)
    end

    foreach(i -> append!(eqblock, eqs[i]), eqs_)
    foreach(i -> deleteat!(eqs, i), Iterators.reverse(sort!(eqs_)))

    push!(eqs, eqblock)
    return eqs_
end

"""
    @join models... [equalize...]

Performs join of models and identifies model variables, as specified.

Model variables / parameter values and metadata are propagated; the last model takes precedence.

# Examples

```julia
@join net1 net2 @catchall(A) = net2.Z @catchall(XY) @catchall(B)
```
"""
macro join(exs...)
    callex = :(
        begin
            merged = ReactionNetwork()
        end
    )
    exs = collect(exs)
    foreach(i -> (exs[i] = MacroTools.striplines(exs[i])), 1:length(exs))
    eqs = []
    ix = 1
    while ix <= length(exs)
        if exs[ix] isa Expr && MacroTools.isexpr(exs[ix], :macrocall, :(=))
            merge_eqs!(eqs, get_eqs(exs[ix]))
            deleteat!(exs, ix)
            continue
        end
        ix += 1
    end

    for netex in exs
        (netex, symex) = if isexpr(netex, :macrocall)
            str_inc = string(
                isexpr(netex.args[3], :(=)) ? netex.args[3].args[2] : netex.args[3],
            )
            if isexpr(netex.args[3], :(=))
                (:(include_model($str_inc)), netex.args[3].args[1])
            else
                (:(include_model($str_inc)), gensym(:net))
            end
        else
            (netex, netex)
        end
        push!(callex.args, :(merge_networks!(merged, $(esc(netex)), $(QuoteNode(symex)), $eqs)))
    end
    push!(callex.args, :(merged))

    return callex
end
