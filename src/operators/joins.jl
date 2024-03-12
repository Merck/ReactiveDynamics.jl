# model joins
export union_acs!, @join

using MacroTools
using MacroTools: prewalk

"""
Merge `acs2` into `acs1`, the attributes in `acs2` taking precedence. Identify respective species given `eqs`, renaming species in `acs2`.
"""
function union_acs!(acs1, acs2, name = gensym("acs"), eqs = [])
    acs2 = deepcopy(acs2)
    prepend!(acs2, name, eqs)

    for i in parts(acs2, :S)
        inc = incident(acs1, acs2[i, :specName], :specName)

        if isempty(inc)
            inc = add_part!(acs1, :S; specName = acs2[i, :specName])
            assign_defaults!(acs1)
        end

        union!(acs1[first(inc), :specModality], acs2[i, :specModality])

        for attr in propertynames(acs1.subparts)
            !occursin("spec", string(attr)) && continue
            !ismissing(acs2[i, attr]) && (acs1[first(inc), attr] = acs2[i, attr])
        end
    end

    new_trans_ix = add_parts!(acs1, :T, nparts(acs2, :T))
    for attr in propertynames(acs2.subparts)
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

    for i in parts(acs2, :P)
        inc = incident(acs1, acs2[i, :prmName], :prmName)
        isempty(inc) && (inc = add_part!(acs1, :P; prmName = acs2[i, :prmName]))
        !ismissing(acs2[i, :prmVal]) && (acs1[first(inc), :prmVal] = acs2[i, :prmVal])
    end

    for i in parts(acs2, :M)
        inc = incident(acs1, acs2[i, :metaKeyword], :metaKeyword)
        isempty(inc) && (inc = add_part!(acs1, :M; metaKeyword = acs2[i, :metaKeyword]))
        !ismissing(acs2[i, :metaVal]) && (acs1[first(inc), :metaVal] = acs2[i, :metaVal])
    end

    return acs1
end

"""
Prepend species names with a model identifier (unless a global species name).
"""
function prepend!(acs::ReactionNetworkSchema, name = gensym("acs"), eqs = [])
    specmap = Dict()
    for i in parts(acs, :S)
        new_name = normalize_name(name, i, acs[i, :specName], eqs)
        push!(specmap, acs[i, :specName] => (acs[i, :specName] = new_name))
    end

    for attr in propertynames(acs.subparts)
        attr == :specName && continue
        attr_ = acs[:, attr]
        for i in eachindex(attr_)
            attr_[i] = escape_ref(attr_[i], collect(keys(specmap)))
            attr_[i] = recursively_substitute_vars!(specmap, attr_[i])
            acs[i, attr] = attr_[i]
        end
    end

    return acs
end

"""
Prepend identifier of an observable with a model identifier.
"""
function prepend_obs end

function prepend_obs(ex::Expr, name)
    return prewalk(ex -> if (isexpr(ex, :macrocall) && (macroname(ex) == :obs))
        (ex.args[end] = Symbol(name, "__", ex.args[end]); ex)
    else
        ex
    end, ex)
end

prepend_obs(ex, _) = ex

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
    if isexpr(eq, :macrocall)
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
            acs_new = ReactionNetworkSchema()
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
                (:(include_model($str_inc)), gensym(:acs))
            end
        else
            (acsex, acsex)
        end
        push!(callex.args, :(union_acs!(acs_new, $(esc(acsex)), $(QuoteNode(symex)), $eqs)))
    end
    push!(callex.args, :(acs_new))

    return callex
end
