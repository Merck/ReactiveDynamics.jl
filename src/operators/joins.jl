# model joins
export merge_networks!, @join

using MacroTools
using MacroTools: prewalk

"""
Merge `acs2` into `acs1`, the attributes in `acs2` taking precedence. Identify respective species given `eqs`, renaming species in `acs2`.
"""
function merge_networks!(acs1, acs2, name = gensym("net"), eqs = [])
    acs2 = deepcopy(acs2)
    prepend!(acs2, name, eqs)

    for i in row_ids(acs2, :S)
        inc = find_rows(acs1, acs2[i, :specName], :specName)

        if isempty(inc)
            inc = add_row!(acs1, :S; specName = acs2[i, :specName])
            assign_defaults!(acs1)
        end

        union!(acs1[first(inc), :specModality], acs2[i, :specModality])

        for attr in propertynames(acs1.columns)
            !occursin("spec", string(attr)) && continue
            !ismissing(acs2[i, attr]) && (acs1[first(inc), attr] = acs2[i, attr])
        end
    end

    new_trans_ix = add_rows!(acs1, :T, nrows(acs2, :T))
    for attr in propertynames(acs2.columns)
        !occursin("trans", string(attr)) && continue
        for (ix1, ix2) in enumerate(new_trans_ix)
            acs1[ix2, attr] = acs2[ix1, attr]
        end
    end

    foreach(
        i -> (
            acs1[i, :transName] =
                normalize_name(Symbol(coalesce(acs1[i, :transName], i)), name)
        ),
        new_trans_ix,
    )

    for i in row_ids(acs2, :P)
        inc = find_rows(acs1, acs2[i, :prmName], :prmName)
        isempty(inc) && (inc = add_row!(acs1, :P; prmName = acs2[i, :prmName]))
        !ismissing(acs2[i, :prmVal]) && (acs1[first(inc), :prmVal] = acs2[i, :prmVal])
    end

    for i in row_ids(acs2, :M)
        inc = find_rows(acs1, acs2[i, :metaKeyword], :metaKeyword)
        isempty(inc) && (inc = add_row!(acs1, :M; metaKeyword = acs2[i, :metaKeyword]))
        !ismissing(acs2[i, :metaVal]) && (acs1[first(inc), :metaVal] = acs2[i, :metaVal])
    end

    # Events (:E) and observables (:obs) are STRUCTURAL — appended, never deduplicated (like :T,
    # §7/J2). The historic gap (merge_networks! walked only :S/:T/:P/:M) silently DROPPED both on join
    # (determinism_composition_bugs.jl §"merge_networks! does NOT merge observables/events"). `prepend!`
    # (above) already namespaced the species referenced inside each event's trigger/action Expr and
    # inside each observable's option-Exprs (via prepend_obs), so both merges are pure structural
    # copies of already-namespaced rows.
    for i in row_ids(acs2, :E)
        add_row!(
            acs1,
            :E;
            eventTrigger = acs2[i, :eventTrigger],
            eventAction = acs2[i, :eventAction],
        )
    end

    for i in row_ids(acs2, :obs)
        add_row!(acs1, :obs; obsName = acs2[i, :obsName], obsOpts = acs2[i, :obsOpts])
    end

    return acs1
end

# Deprecated ACSets-vocabulary alias (ADR 0015 Tier 2): `union_acs!` → `merge_networks!`.
@deprecate union_acs!(net1, net2, name = gensym("net"), eqs = []) merge_networks!(net1, net2, name, eqs)

"""
Prepend species names with a model identifier (unless a global species name).
"""
function prepend!(net::ReactionNetwork, name = gensym("net"), eqs = [])
    specmap = Dict()
    for i in row_ids(net, :S)
        # ADR 0009 §A / CONTRACT §11.1: a `shared`-role species is identified by BARE name across all
        # fragments (the first-class @catchall) — it is NOT namespaced. `private` (default) and the
        # open `input`/`output` ports namespace as usual here; @compose (§E) re-identifies the open
        # ports afterwards by FK-repoint. (A species carrying no role reads :private via port_role.)
        if port_role(net, i) === :shared
            continue
        end
        new_name = normalize_name(name, i, net[i, :specName], eqs)
        push!(specmap, net[i, :specName] => (net[i, :specName] = new_name))
    end

    for attr in propertynames(net.columns)
        attr == :specName && continue
        # Observable options live inside a FoldedObservable struct (the :obsOpts column), not as a
        # bare Expr the loop below rewrites — handle them structurally via prepend_obs! so species
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
Namespace the species referenced inside an observable's option expressions.

`prepend!` renames every species `X → parent__X` and records the map in `specmap`. An observable's
sampling triggers (`on`) and range endpoints (`range`) are stored as Exprs inside a `FoldedObservable`
(the `:obsOpts` column), which `prepend!`'s attribute loop skips — so without this the observable would
still reference the pre-namespaced species and silently read the wrong (or a missing) pool after a join.
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

## species name normalization
normalize_name(name::Symbol, parent_name) = Symbol("$(parent_name)__$name")
normalize_name(name::String, parent_name) = "$(parent_name)__$name"
normalize_name(name, parent_name) = Symbol(parent_name, "__", name)

function normalize_name(acs_name, i::Int, name::Symbol, eqs = [])
    for (block_ix, block) in enumerate(eqs)
        block_alias = findfirst(e -> e[1] == :alias, block)
        block_alias = if !isnothing(block_alias)
            block[block_alias][2]
        else
            Symbol(:shared_species_, block_ix)
        end
        for e in block
            (
                (i == e[2]) ||
                    (
                    e[1] == :catchall &&
                        (normalize_name(e[2], acs_name) == normalize_name(name, acs_name))
                ) ||
                    (
                    e[1] == acs_name &&
                        (normalize_name(e[2], acs_name) == normalize_name(name, acs_name))
                )
            ) && return block_alias
        end
    end

    return normalize_name(name, acs_name)
end

matching_name(name::Symbol, parent_name) = [name, Symbol("$(parent_name)__$name")]

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

function recursively_get_syms(ex)
    return isexpr(ex, :.) ? [recursively_get_syms(ex.args[1]); ex.args[2].value] : ex
end
function reconstruct(ex)
    return if ex isa Symbol
        (ex,)
    else
        (syms = recursively_get_syms(ex); (syms[1], Symbol(join(syms[2:end], "__"))))
    end
end

"""
Parse species equation blocks.
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
@join acs1 acs2 @catchall(A) = acs2.Z @catchall(XY) @catchall(B)
```
"""
macro join(exs...)
    callex = :(
        begin
            acs_new = ReactionNetwork()
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

    for acsex in exs
        (acsex, symex) = if isexpr(acsex, :macrocall)
            str_inc = string(
                isexpr(acsex.args[3], :(=)) ? acsex.args[3].args[2] : acsex.args[3],
            )
            if isexpr(acsex.args[3], :(=))
                (:(include_model($str_inc)), acsex.args[3].args[1])
            else
                (:(include_model($str_inc)), gensym(:net))
            end
        else
            (acsex, acsex)
        end
        push!(callex.args, :(merge_networks!(acs_new, $(esc(acsex)), $(QuoteNode(symex)), $eqs)))
    end
    push!(callex.args, :(acs_new))

    return callex
end
