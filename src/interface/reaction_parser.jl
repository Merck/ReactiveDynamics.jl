# parts of the code were taken from Catalyst.jl

struct FoldedReactant
    species::Symbol
    stoich
    modality::Set{Symbol}
end

function recursively_choose(r_line, state)
    postwalk(r_line) do ex
        if isexpr(ex, :macrocall) && (macroname(ex) == :choose)
            sample_range([(isexpr(r, :tuple) ? (r.args[1], recursively_choose(r.args[2], state)) : recursively_choose(r, state)) for r in ex.args[3:end]], state)
        else ex end
    end
end

function extract_reactants(r_line, state::DyVEState)
    r_line = recursively_choose(r_line, state)

    recursive_find_reactants!(escape_ref(r_line, state[:, :specName]), 1., Set{Symbol}(), Vector{FoldedReactant}(undef, 0))
end

function recursive_find_reactants!(ex, mult, mods::Set{Symbol}, reactants::Vector{FoldedReactant})
    if typeof(ex) != Expr || isexpr(ex, :.) || (ex.head == :escape)   
        if (ex == 0 || in(ex, empty_set)) return reactants
        else push!(reactants, FoldedReactant(recursively_expand_dots(ex), mult, mods)) end
    elseif ex.args[1] == :*
        recursive_find_reactants!(ex.args[end], multiplex(mult, ex.args[2:end-1]...), mods, reactants)
    elseif ex.args[1] == :+
        for i = 2:length(ex.args) recursive_find_reactants!(ex.args[i], mult, mods, reactants) end
    elseif ex.head == :macrocall
        mods = copy(mods)
        macroname(ex) in species_modalities && push!(mods, macroname(ex))
        foreach(i -> push!(mods, ex.args[i] isa Symbol ? ex.args[i] : ex.args[i].value), 4:length(ex.args))
        recursive_find_reactants!(ex.args[3], mult, mods, reactants)
    else @error("malformed reaction") end

    reactants
end


function multiplex(mult, mults...)
    all(m -> isa(m, Number), [mult] ∪ mults) && return mult * prod(mults; init=1.)
    multarray = Any[]; recursively_find_mults!(multarray, mults...)
    mults_numeric = prod(filter(m -> isa(m, Number), multarray); init=1.) * (mult isa Number ? mult : 1.)
    mults_expr = filter(m -> !isa(m, Number), multarray)
    mult = mult isa Expr ? deepcopy(mult) : :(*())
    mults_numeric == 1 && length(mults_expr) == 1 && return mults_expr[1]
    mults_numeric != 1. && push!(mult.args, mults_numeric); append!(mult.args, mults_expr)
    
    mult
end

function recursively_find_mults!(multarray, mults...)
    for m in mults; isa(m, Expr) && m.args[1] == :* ? recursively_find_mults!(multarray, m.args[2:end]...) : push!(multarray, m) end
end