# transition line parsing

## first round: extract parameters and species
"Extract transitions and species from an expression."
function parse_transitions(expr)
    transitions = []; species = []
    if isexpr(expr, :block)
        expr = striplines(expr); esc_dollars!(expr)
        foreach(l -> get_transitions!(transitions, species, l), expr.args)
    elseif expr != :()
        get_transitions!(transitions, species, expr)
    end

    transitions, unique!(species)
end

"Expand a rate expression: turns `@ct(x)` to `1/x` (to a rate)."
function expand_rate(rate)
    rate = if !(isexpr(rate, :macrocall) && (macroname(rate) == :per_step))
        :(rand(Poisson(max($rate, 0))))
    else rate.args[3] end

    postwalk(rate) do ex
        if (isexpr(ex, :macrocall) && (macroname(ex) ∈ prettynames[:transition][:cycle_time]))
            :(1/$(ex.args[3]))
        else ex end
    end
end

"Parse a (rate, transition) tuple and append transitions and species."
function get_transitions!(transitions, species, ex)
    exs = isexpr(ex, :tuple) ? ex.args : [ex]

    rate, rline = exs[1:2]
    rline = parse_rline!(species, rline); rline = make_paths(rline, species)
    rate = expand_rate(rate)
    rxs = rline isa Tuple ? tuple.(rate.args, rline) : ((rate, rline),)

    args = Dict()
    for arg in exs[3:end]
        (arg isa Expr && isexpr(arg, :call)) || continue
        key = findfirst(k -> arg.args[2] ∈ k, prettynames[:transition]) # localize kwargs
        
        # extract params
        !isnothing(key) && push!(args, key => arg.args[3])
    end

    append!(transitions, tuple.(rxs, Ref(args)))
    
    transitions
end

"Parse a transition line and append transitions and species."
function parse_rline!(species, rline)
    rline isa Expr && (rline.head == :-->) && (rline = Expr(:call, :→, rline.args[1], rline.args[2]))

    if isexpr(rline, :macrocall)
        lines = copy(rline.args[3:end]); rline.args = rline.args[1:2]
        for l in lines
            biarrow = nothing; prewalk(ex -> (ex ∈ double_arrows && (biarrow = ex); ex), l)
            append!(rline.args, isnothing(biarrow) ? [l] : [replace_in_expr(l, biarrow => :⟶), replace_in_expr(l, biarrow => :⟵)])
        end
        for i in 3:length(rline.args)
            rline.args[i] = 
                if isexpr(rline.args[i], :tuple)
                    Expr(:tuple, rline.args[i].args[1], parse_rline!(species, rline.args[i].args[2]))
                else parse_rline!(species, rline.args[i]) end
        end
    elseif rline isa Expr && rline.args[1] ∈ union(fwd_arrows, bwd_arrows)
        rline.args[2:3] = find_species!.(Ref(species), rline.args[2:3])
    elseif rline isa Expr && rline.args[1] ∈ double_arrows
        biarrow = nothing; prewalk(ex -> (ex ∈ double_arrows && (biarrow = ex); ex), rline)
        rline = parse_rline!.(Ref(species), (replace_in_expr(rline, biarrow => :⟶), replace_in_expr(rline, biarrow => :⟵)))
    end

    rline
end

"Recursively extract species."
function find_species!(species, ex)
    if typeof(ex) != Expr || (ex.head == :escape)
        if (ex == 0 || in(ex, empty_set)) return :∅
        else push!(species, to_string(ex)) end
    elseif ex.args[1] == :*
        find_species!(species, ex.args[end])
        foreach(i -> ex.args[i] = ex.args[i], 2:(length(ex.args)-1))
    elseif ex.args[1] == :+
        for i = 2:length(ex.args) find_species!(species, ex.args[i]) end
    elseif isexpr(ex, :macrocall) && macroname(ex) == :choose
        for i in 3:length(ex.args)
            find_species!(species, isexpr(ex.args[i], :tuple) ? ex.args[i].args[2] : ex.args[i])
        end
    elseif ispath(ex)
        push!(species, to_string(ex))
    elseif isexpr(ex, :macrocall) && macroname(ex) != :action
        find_species!(species, ex.args[3])
    else return end

    ex
end

## second round: extract species w/ stoich coefficients
"Parametrized substrate or product in a transition."
struct Species
    agent # species identifier
    stoich # stoich coeff can be a general expression
    modalities::Set{Symbol} # local modalities
end

"Split transition line into (LHS, RHS) tuple."
function split_rline(rline, network, this)
    if rline isa Expr && rline.args[1] ∈ fwd_arrows; rline.args[2:3]
    elseif rline isa Expr && rline.args[1] ∈ bwd_arrows; rline.args[[3, 2]]
    elseif isexpr(rline, :macrocall) && (macroname(rline) == :choose)
        choices = map(rline.args[3:end]) do r
            if isexpr(r, :tuple)
                (r.args[1], split_rline(r.args[2], network, this))
            else split_rline(r, network, this) end
        end

        rand_polyrange(choices, network, this)
    end
end

"Walks a transition half-line (LHS, RHS) and takes choices, returning a sampled (structural) transition half-line."
function recursively_choose(rline, network, this)
    postwalk(rline) do ex
        if isexpr(ex, :macrocall) && (macroname(ex) == :choose)
            choices = map(ex.args[3:end]) do r
                if isexpr(r, :tuple)
                    (r.args[1], recursively_choose(r.args[2], network, this))
                else recursively_choose(r, network, this) end
            end

            rand_polyrange(choices, network, this)
        else ex end
    end
end

"Parse a transition half-line (LHS, RHS) and return species with stoich coefficients and modalities."
function extract_species(rline, network, this)
    rline = recursively_choose(rline, network, this)

    append_species!(rline, 1., Set{Symbol}(), Species[])
end

"Extract species with stoich coefficients and modalities."
function append_species!(species::Vector{Species}, mult, mods::Set{Symbol}, ex)
    if typeof(ex) != Expr || (ex.head == :escape)   
        if (ex == 0 || in(ex, empty_set)) return species
        else push!(species, Species(ex, mult, mods)) end
    elseif ex.args[1] == :*
        append_species!(species, multiplex(mult, ex.args[2:end-1]...), mods, ex.args[end])
    elseif ex.args[1] == :+
        for i = 2:length(ex.args) append_species!(species, mult, mods, ex.args[i]) end
    elseif isexpr(ex, :macrocall) && macroname(ex) == :action
        push!(species, Species(wrap_macro(ex) , 1., Set{Symbol}()))
    elseif ex.head == :macrocall
        mods = copy(mods)
        macroname(ex) in species_modalities && push!(mods, macroname(ex))
        foreach(i -> push!(mods, ex.args[i] isa Symbol ? ex.args[i] : ex.args[i].value), 4:length(ex.args))
        append_species!(species, mult, mods, ex.args[3])
    else @error("malformed transition") end

    species
end

"Wrap the macro's body into an expression that returns `nothing`."
wrap_macro(ex) = :(begin $(ex.args[end]); nothing end)

"Return compound stoichiometry coefficient (join expressions as stoich coefficients)."
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

"Recursively extract stoichiometry coefficients."
function recursively_find_mults!(multarray, mults...)
    for m in mults; isa(m, Expr) && m.args[1] == :* ? recursively_find_mults!(multarray, m.args[2:end]...) : push!(multarray, m) end
end