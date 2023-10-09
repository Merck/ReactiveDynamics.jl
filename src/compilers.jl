using MacroTools: prewalk

"""
Recursively find model variables in expressions.
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

function get_contained_params(ex, prms)
    r = []
    return recursively_find_vars!(r, ex) ∩ prms
end

"""
Recursively substitute model variables. Subsitution pairs are specified in `varmap`.
"""
function recursively_substitute_vars!(varmap, ex)
    ex isa Symbol && return (haskey(varmap, ex) ? varmap[ex] : ex)
    ex isa Expr && for i = 1:length(ex.args)
        if ex.args[i] isa Expr
            recursively_substitute_vars!(varmap, ex.args[i])
        else
            (
                ex.args[i] isa Symbol &&
                haskey(varmap, ex.args[i]) &&
                (ex.args[i] = varmap[ex.args[i]])
            )
        end
    end

    return ex
end

"""
Recursively normalize dotted notation in `ex` to species names whenever a name is contained in `vars`, else left unchanged.
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
    ex isa Expr && for i = 1:length(ex.args)
        ex.args[i] isa Union{Expr,Symbol} &&
            (ex.args[i] = recursively_expand_dots_in_ex!(ex.args[i], vars))
    end

    return ex
end

reserved_names =
    [:t, :state, :obs, :resample, :solverarg, :take, :log, :periodic, :set_params]
push!(reserved_names, :state)

function escape_ref(ex, species)
    return if ex isa Symbol
        ex
    else
        prewalk(
            ex ->
                isexpr(ex, :ref) && Symbol(string(ex)) ∈ species ? Symbol(string(ex)) : ex,
            ex,
        )
    end
end

function wrap_expr(fex, species_names, prm_names, varmap)
    !isa(fex, Union{Expr,Symbol}) && return fex
    # escape refs in species names: A[1] -> Symbol("A[1]")
    fex = escape_ref(fex, species_names)
    # escape dots in species' names: A.B -> Symbol("A.B")
    fex = deepcopy(fex)
    fex = recursively_expand_dots_in_ex!(fex, species_names)

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
        else
            x
        end
    end

    # substitute the species names with "pointers" into the state space: S -> state.u[1]
    fex = recursively_substitute_vars!(varmap, fex)
    # substitute the params names with "pointers" into the parameter space: β -> state.p[:β]
    # params can't be mutated!
    foreach(
        v -> push!(letex.args[1].args, :($v = state.p[$(QuoteNode(v))])),
        get_contained_params(fex, prm_names),
    )
    push!(letex.args[2].args, fex)

    # the function shall be a function of the dynamic ReactionNetworkSchema structure: letex -> :(state -> $letex)
    # eval the expression to a Julia function, save that function into the "compiled" acset
    return eval(:(state -> $letex))
end

function get_wrap_fun(acs::ReactionNetworkSchema)
    species_names = collect(acs[:, :specName])
    prm_names = collect(acs[:, :prmName])
    varmap = Dict([name => :(state.u[$i]) for (i, name) in enumerate(species_names)])
    for name in prm_names
        push!(varmap, name => :(state.p[$(QuoteNode(name))]))
    end

    return ex -> wrap_expr(ex, species_names, prm_names, varmap)
end

function skip_compile(attr)
    return any(contains.(Ref(string(attr)), ("Name", "obs", "meta"))) ||
           (string(attr) == "trans")
end

function compile_attrs(acs::ReactionNetworkSchema)
    species_names = collect(acs[:, :specName])
    prm_names = collect(acs[:, :prmName])
    varmap = Dict([name => :(state.u[$i]) for (i, name) in enumerate(species_names)])
    for name in prm_names
        push!(varmap, name => :(state.p[$(QuoteNode(name))]))
    end
    wrap_fun = ex -> wrap_expr(ex, species_names, prm_names, varmap)
    attrs = Dict{Symbol,Vector}()
    transitions = Dict{Symbol,Vector}()
    for attr in propertynames(acs.subparts)
        attrs_ = subpart(acs, attr)
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
    transitions[:transActivated] = fill(true, nparts(acs, :T))
    transitions[:transToSpawn] = zeros(nparts(acs, :T))
    transitions[:transHash] =
        [coalesce(acs[i, :transName], gensym()) for i in parts(acs, :T)]

    return attrs, transitions, wrap_fun
end

function remove_choose(acs::ReactionNetworkSchema)
    acs = deepcopy(acs)
    pcs = []
    for attr in propertynames(acs.subparts)
        attrs_ = subpart(acs, attr)
        foreach(
            i ->
                !isnothing(attrs_[i]) &&
                    attrs_[i] isa Expr &&
                    (attrs_[i] = normalize_pcs!(pcs, attrs_[i])),
            1:length(attrs_),
        )
    end

    add_obs!(acs, pcs)
    return acs
end
