@reexport using AlgebraicAgents
using DataFrames

struct UnfoldedReactant
    index::Int
    species::Symbol
    stoich::ActionableValues
    modality::Set{Symbol}
end

"""
Ongoing transition auxiliary structure.
"""
@aagent struct Transition
    i::Int

    trans::Dict{Symbol,Any}

    bound_structured_agents::Vector{AbstractAlgebraicAgent}
    nonblock_structured_agents::Vector{AbstractAlgebraicAgent}
    structured_to_agents::Vector

    t::Float64
    q::Float64
    state::Float64
end

Base.getindex(state::Transition, key) = state.trans[key]
Base.setindex!(state::Transition, val, key) = state.trans[key] = val

@aagent struct Observable
    last::Float64 # last sampling time
    range::Vector{Union{Tuple{Float64,SampleableValues},SampleableValues}}
    every::Float64
    on::Vector{ActionableValues}

    sampled::Any
end

@aagent struct ReactionNetworkProblem
    acs::ReactionNetworkSchema

    attrs::Dict{Symbol,Vector}
    transition_recipes::Dict{Symbol,Vector}

    u::Vector{Float64}
    p::Any
    t::Float64

    structured_token::Vector{Symbol}

    tspan::Tuple{Float64,Float64}
    dt::Float64

    transitions::Dict{Symbol,Vector}
    ongoing_transitions::Vector{Transition}
    log::Vector{Tuple}

    observables::Dict{Symbol,Observable}

    wrap_fun::Any
    sol::DataFrame
end

# get value of a numeric expression
# evaluate compiled numeric expression in context of (u, p, t)
function context_eval(state::ReactionNetworkProblem, transition, o)
    o = o isa Function ? Base.invokelatest(o, state, transition) : o

    return o isa Sampleable ? rand(o) : o
end

function Base.getindex(state::ReactionNetworkProblem, keys...)
    if any(occursin.(["transPreAction", "transPostAction"], Ref(string(keys[2]))))
        return state.acs[keys[1], keys[2]]
    else
        return context_eval(
            state,
            nothing,
            (contains(string(keys[2]), "trans") ? state.transitions : state.attrs)[keys[2]][keys[1]],
        )
    end
end

function init_u!(state::ReactionNetworkProblem)
    return (u = fill(0.0, nparts(state, :S));
    foreach(i -> u[i] = state[i, :specInitVal], parts(state, :S));
    state.u = u)
end
save!(state::ReactionNetworkProblem) = push!(state.sol, (state.t, state.u[:]...))

function compile_observables(acs::ReactionNetworkSchema)
    observables = Dict{Symbol,Observable}()
    species_names = collect(acs[:, :specName])
    prm_names = collect(acs[:, :prmName])
    varmap = Dict([name => :(state.u[$i]) for (i, name) in enumerate(species_names)])

    for (name, opts) in Iterators.zip(acs[:, :obsName], acs[:, :obsOpts])
        on = map(on -> wrap_expr(on, species_names, prm_names, varmap), opts.on)
        range = map(
            r -> begin
                r = r isa Tuple ? r : (1.0, r)
                (r[1], wrap_expr(r[2], species_names, prm_names, varmap))
            end,
            opts.range,
        )

        push!(
            observables,
            name => Observable(string(name), -Inf, range, opts.every, on, missing),
        )
    end

    return observables
end

compileval(ex, state) = !isa(ex, Expr) ? ex : eval(state.wrap_fun(ex))

function sample_range(rng, state)
    isempty(rng) && return missing
    r = rand() * sum(r -> r isa Tuple ? eval(r[1]) : 1, rng)
    ix = 0
    s = 0
    while s <= r && (ix < length(rng))
        ix += 1
        s += rng[ix] isa Tuple ? compileval(rng[ix][1], state) : 1
    end

    r = rng[ix] isa Tuple ? rng[ix][2] : rng[ix]
    return r isa Sampleable ? rand(r) : r
end

function resample!(state::ReactionNetworkProblem, o::Observable)
    o.last = state.t
    isempty(o.range) && (return o.val = missing)

    return o.sampled = context_eval(state, nothing, sample_range(o.range, state))
end

resample(state::ReactionNetworkProblem, o::Symbol) = resample!(state, state.observables[o])

function update_observables(state::ReactionNetworkProblem)
    return foreach(
        o -> (state.t - o.last) >= o.every && resample!(state, o),
        values(state.observables),
    )
end

function prune_r_line(r_line)
    return if r_line isa Expr && r_line.args[1] ∈ fwd_arrows
        r_line.args[[2, 3]]
    elseif r_line isa Expr && r_line.args[1] ∈ bwd_arrows
        r_line.args[[3, 2]]
    elseif isexpr(r_line, :macrocall) && (macroname(r_line) == :choose)
        sample_range(
            [(
                if isexpr(r, :tuple)
                    (r.args[1], prune_r_line(r.args[2]))
                else
                    prune_r_line(r)
                end
            ) for r in r_line.args[3:end]],
            state,
        )
    end
end

function find_index(species::Symbol, state::ReactionNetworkProblem)
    return findfirst(i -> state[i, :specName] == species, parts(state, :S))
end

function sample_transitions!(state::ReactionNetworkProblem)
    for (_, v) in state.transitions
        empty!(v)
    end
    for i = 1:length(state.transition_recipes[:trans])
        !(state.transition_recipes[:transActivated][i]) && continue
        l_line, r_line = prune_r_line(state.transition_recipes[:trans][i])

        for attr in keys(state.transition_recipes)
            (
                attr ∈
                [:trans, :transPreAction, :transPostAction, :transActivated, :transHash]
            ) && continue
            push!(
                state.transitions[attr],
                context_eval(state, nothing, state.transition_recipes[attr][i]),
            )
        end

        reactants = []
        for r in extract_reactants(l_line, state)
            j = find_index(r.species, state)
            push!(
                reactants,
                UnfoldedReactant(
                    j,
                    r.species,
                    context_eval(state, nothing, state.wrap_fun(r.stoich)),
                    r.modality ∪ state[j, :specModality],
                ),
            )
        end

        push!(state.transitions[:transLHS], reactants)
        push!(state.transitions[:transRHS], r_line)

        foreach(
            k -> push!(state.transitions[k], state.transition_recipes[k][i]),
            [:transPreAction, :transPostAction, :transToSpawn, :transHash],
        )

        state.transition_recipes[:transToSpawn] .= 0
    end
end

function as_state(u, t, state::ReactionNetworkProblem)
    return (state = deepcopy(state); state.u .= u; state.t = t; state)
end

function ACSets.ACSetInterface.nparts(state::ReactionNetworkProblem, obj::Symbol)
    return nparts(state.acs, obj)
end

function ACSets.ACSetInterface.parts(state::ReactionNetworkProblem, obj::Symbol)
    return parts(state.acs, obj)
end

## query the state

t(state::ReactionNetworkProblem) = state.t
solverarg(state::ReactionNetworkProblem, arg) = state.p[arg]
take(state::ReactionNetworkProblem, pcs::Symbol) = state.observables[pcs].sampled
log(state::ReactionNetworkProblem, msg) = (println(msg); push!(state.log, (:log, msg)))
state(state::ReactionNetworkProblem) = state

function periodic(state::ReactionNetworkProblem, period)
    return period == 0.0 || (
        length(state.sol.t) > 1 &&
        (fld(state.t, period) - fld(state.sol.t[end-1], period) > 0)
    )
end

set_params(state::ReactionNetworkProblem, vals...) =
    for (p, v) in vals
        state.p[p] = v
    end

function add_to_spawn!(state::ReactionNetworkProblem, hash, n)
    ix = findfirst(ix -> state.transition_recipes[:transHash][ix] == hash)
    return !isnothing(ix) && (state.transition_recipes[:transHash][ix] += n)
end
