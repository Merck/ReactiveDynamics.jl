using DiffEqBase: NullParameters

struct UnfoldedReactant
    index::Int
    species::Symbol
    stoich::ActionableValues
    modality::Set{Symbol}
end

"Ongoing transition auxiliary structure."
mutable struct Transition
    trans::Dict{Symbol, Any}

    t::Float64
    q::Float64
    state::Float64
end

Base.getindex(state::Transition, key) = state.trans[key]

mutable struct Observable
    last::Float64 # last sampling time
    range::Vector{Union{Tuple{Float64, SampleableValues}, SampleableValues}}
    every::Float64
    on::Vector{ActionableValues}

    sampled
end

mutable struct ReactiveDynamicsState
    acs::ReactionNetwork

    attrs::Dict{Symbol, Vector}
    transition_recipes::Dict{Symbol, Vector}

    u::Vector{Float64};; p; t::Float64

    transitions::Dict{Symbol, Vector}
    ongoing_transitions::Vector{Transition}
    log::Vector{Tuple}

    observables::Dict{Symbol, Observable}
    solverargs

    wrap_fun
    history_u::Vector{Vector{Float64}}; history_t::Vector{Float64}

    function ReactiveDynamicsState(acs::ReactionNetwork, attrs, transition_recipes, wrap_fun, t0=0; kwargs...)
        ongoing_transitions = Transition[]
        log = NamedTuple[]
        observables = compile_observables(acs)
        transitions_attrs = setdiff(filter(a -> contains(string(a), "trans"), keys(acs.attrs)), (:trans,)) ∪ 
            [:transLHS, :transRHS, :transToSpawn, :transHash]
        transitions = Dict{Symbol, Vector}(a => [] for a in transitions_attrs)

        new(acs, attrs, transition_recipes, zeros(nparts(acs, :S)), fetch_params(acs), t0,
            transitions, ongoing_transitions, log, observables, kwargs, wrap_fun, Vector{Float64}[], Float64[])
    end
end

# get value of a numeric expression
# evaluate compiled numeric expression in context of (u, p, t)
function context_eval(state::ReactiveDynamicsState, o)
    o = o isa Function ? Base.invokelatest(o, state) : o
    
    o isa Sampleable ? rand(o) : o
end

Base.getindex(state::ReactiveDynamicsState, keys...) = context_eval(state, (contains(string(keys[2]), "trans") ? state.transitions : state.attrs)[keys[2]][keys[1]])

init_u!(state::ReactiveDynamicsState) = (u = fill(.0, nparts(state, :S)); foreach(i -> u[i] = state[i, :specInitVal], 1:nparts(state, :S)); state.u = u)
save!(state::ReactiveDynamicsState) = (push!(state.history_u, state.u); push!(state.history_t, state.t))

function compile_observables(acs::ReactionNetwork)
    observables = Dict{Symbol, Observable}()
    species_names = collect(acs[:, :specName]); prm_names = collect(acs[:, :prmName])
    varmap = Dict([name => :(state.u[$i]) for (i, name) in enumerate(species_names)])

    for (name, opts) in Iterators.zip(acs[:, :obsName], acs[:, :obsOpts])
        on = map(on -> wrap_expr(on, species_names, prm_names, varmap), opts.on)
        range = map(r -> begin
            r = r isa Tuple ? r : (1., r)
            (r[1], wrap_expr(r[2], species_names, prm_names, varmap))
        end, opts.range)

        push!(observables, name => Observable(-Inf, range, opts.every, on, missing))
    end

    observables
end

compileval(ex, state) = !isa(ex, Expr) ? ex : eval(state.wrap_fun(ex))

function sample_range(rng, state)
    isempty(rng) && return missing
    r = rand() * sum(r -> r isa Tuple ? eval(r[1]) : 1, rng)
    ix = 0; s = 0; while s <= r && (ix < length(rng)); ix += 1; s += rng[ix] isa Tuple ? compileval(rng[ix][1], state) : 1; end

    r = rng[ix] isa Tuple ? rng[ix][2] : rng[ix]
    r isa Sampleable ? rand(r) : r
end

function resample!(state::ReactiveDynamicsState, o::Observable)
    o.last = state.t
    isempty(o.range) && (return o.val = missing)

    o.sampled = context_eval(state, sample_range(o.range, state))
end

resample(state::ReactiveDynamicsState, o::Symbol) = resample!(state, state.observables[o])

update_observables(state::ReactiveDynamicsState) = 
    foreach(o -> (state.t - o.last) >= o.every && resample!(state, o), values(state.observables))

function prune_r_line(r_line)
    return if r_line isa Expr && r_line.args[1] ∈ fwd_arrows; r_line.args[2:3]
    elseif r_line isa Expr && r_line.args[1] ∈ bwd_arrows; r_line.args[[3, 2]]
    elseif isexpr(r_line, :macrocall) && (macroname(r_line) == :choose)
        sample_range([(isexpr(r, :tuple) ? (r.args[1], prune_r_line(r.args[2])) : prune_r_line(r)) for r in r_line.args[3:end]], state)
    end
end

find_index(species::Symbol, state::ReactiveDynamicsState) = findfirst(i -> state[i, :specName] == species, 1:nparts(state, :S))

function sample_transitions!(state::ReactiveDynamicsState)
    for (_, v) in state.transitions; empty!(v) end
    for i in 1:length(state.transition_recipes[:trans])
        !(state.transition_recipes[:transActivated][i]) && continue
        l_line, r_line = prune_r_line(state.transition_recipes[:trans][i])

        for attr in keys(state.transition_recipes)
            attr ∈ [:trans, :transPostAction, :transActivated, :transHash] && continue
            push!(state.transitions[attr], context_eval(state, state.transition_recipes[attr][i]))
        end

        reactants = []; for r in extract_reactants(l_line, state)
            j = find_index(r.species, state)
            push!(reactants, UnfoldedReactant(j, r.species, context_eval(state, state.wrap_fun(r.stoich)), r.modality ∪ state[j, :specModality]))
        end
        push!(state.transitions[:transLHS], reactants); push!(state.transitions[:transRHS], r_line)
        foreach(k -> push!(state.transitions[k], state.transition_recipes[k][i]), [:transPostAction, :transToSpawn, :transHash])
        state.transition_recipes[:transToSpawn] .= 0
    end
end

## sync
update_u!(state::ReactiveDynamicsState, u) = (state.u .= u)
update_t!(state::ReactiveDynamicsState, t) = (state.t = t)
sync_p!(p, state::ReactiveDynamicsState) = merge!(p, state.p)

function sync!(state::ReactiveDynamicsState, u, p)
    state.u .= u
    for k in keys(state.p) haskey(p, k) && (state.p[k] = p[k]) end
end

as_state(u, t, state::ReactiveDynamicsState) = (state = deepcopy(state); state.u .= u; state.t = t; state)

Catlab.CategoricalAlgebra.nparts(state::ReactiveDynamicsState, obj::Symbol) = obj == :T ? length(state.transitions[:transLHS]) : nparts(state.acs, obj)

## query the state

t(state::ReactiveDynamicsState) = state.t
solverarg(state::ReactiveDynamicsState, arg) = state.solverargs[arg]
take(state::ReactiveDynamicsState, pcs::Symbol) = state.observables[pcs].sampled
log(state::ReactiveDynamicsState, msg) = (println(msg); push!(state.log, (:log, msg)))
state(state::ReactiveDynamicsState) = state

function periodic(state::ReactiveDynamicsState, period)
    period == .0 || 
        (length(state.history_t) > 1 && (fld(state.t, period)-fld(state.history_t[end-1], period) > 0))
end

function set_params(state::ReactiveDynamicsState, vals...)
    for (p, v) in vals state.p[p] = v end 
end

function add_to_spawn!(state::ReactiveDynamicsState, hash, n)
    ix = findfirst(ix -> state.transition_recipes[:transHash][ix] == hash)
    !isnothing(ix) && (state.transition_recipes[:transHash][ix] += n)
end