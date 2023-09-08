# parts of the code were taken from Catalyst.jl and adapted

using MacroTools: postwalk

struct FoldedReactant
    species::Symbol
    stoich::SampleableValues
    modality::Set{Symbol}
end

function recursively_choose(r_line, state)
    postwalk(r_line) do ex
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

function extract_reactants(r_line, state::ReactiveNetwork)
    r_line = recursively_choose(r_line, state)

    return recursive_find_reactants!(
        escape_ref(r_line, state[:, :specName]),
        1.0,
        Set{Symbol}(),
        Vector{FoldedReactant}(undef, 0),
    )
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
            multiplex(mult, ex.args[2:(end-1)]...),
            mods,
            reactants,
        )
    elseif ex.args[1] == :+
        for i = 2:length(ex.args)
            recursive_find_reactants!(ex.args[i], mult, mods, reactants)
        end
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
