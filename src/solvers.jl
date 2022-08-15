# assortment of SciML-compatible problem solvers

export DiscreteProblem

using DiffEqBase, DifferentialEquations
using Distributions

function get_sampled_transition(state, i)
    transition = Dict{Symbol, Any}()
    foreach(k -> push!(transition, k => state[i, k]), keys(state.transitions))

    transition
end

"Compute resource requirements given transition quantities."
function get_reqs_init!(reqs, qs, state) 
    reqs .= .0
    for i in 1:size(reqs, 2)
        for tok in state[i, :transLHS]
            !any(m -> m in tok.modality, [:rate, :nonblock]) && (reqs[tok.index, i] += qs[i] * tok.stoich)
        end
    end

    reqs
end

"Compute resource requirements given transition quantities."
function get_reqs_ongoing!(reqs, qs, state) 
    reqs .= .0
    for i in 1:length(state.ongoing_transitions)
        for tok in state.ongoing_transitions[i][:transLHS]
            in(:rate, tok.modality) && (state.ongoing_transitions[i][:transCycleTime] > 0) &&
                (reqs[tok.index, i] += qs[i] * tok.stoich * state.solverargs[:tstep])
            in(:nonblock, tok.modality) && (reqs[tok.index, i] += qs[i] * tok.stoich)
        end
    end

    reqs
end

"Given requirements, return available allocation."
function get_allocs!(reqs, u, state, priorities, strategy=:weighted)
    strategy == :weighted ? alloc_weighted!(reqs, u, priorities, state) : alloc_greedy!(reqs, u, priorities, state)
end

function alloc_weighted!(reqs, u, priorities, state)
    allocs = zero(reqs)
    for i in 1:size(reqs, 1)
        s = sum(reqs[i, :]); u[i] >= s && (allocs[i, :] .= reqs[i, :]; continue)
        foreach(j -> allocs[i, j] = reqs[i, j] * priorities[j], 1:size(reqs, 2))
        s = sum(allocs[i, :]); allocs[i, :] .*= (s == 0) ? 0. : (max(0, u[i])/s)
    end

    allocs
end

function alloc_greedy!(reqs, u, priorities, state) 
    allocs = zero(reqs); sorted_trans = sort(1:size(reqs, 2), by=i -> -priorities[i])
    for i in 1:size(reqs, 1)
        s = sum(reqs[i, :]); u[i] >= s && (allocs[i, :] .= reqs[i, :]; continue)
        a = u[i]; j = 1
        while a > 0 && j <= size(reqs, 2)
            allocs[i, sorted_trans[j]] = min(reqs[i, sorted_trans[j]], a); a -= allocs[i, sorted_trans[j]]; j += 1
        end
    end

    allocs
end

"Given resource requirements and available allocations, output resulting shift size for each transition."
function get_frac_satisfied(allocs, reqs, state)
    for i in eachindex(allocs); allocs[i] = min(1, (reqs[i] == .0 ? 1 : (allocs[i] / reqs[i]))) end    
    qs = vec(minimum(allocs; dims=1)); foreach(i -> allocs[:, i] .= reqs[:, i] * qs[i], 1:size(reqs, 2))

    qs
end

"Given available allocations and qties of transitions requested to spawn, return number of spawned transitions. Update `alloc` to match actual allocation."
function get_init_satisfied(allocs, qs, state)
    reqs = zero(allocs)
    for i in 1:size(allocs, 2)
        all(allocs[:, i] .>= 0) || (allocs[:, i] .= 0.; qs[i] = 0)
        for tok in state[i, :transLHS]
            !any(m -> m in tok.modality, [:rate, :nonblock]) && (reqs[tok.index, i] += tok.stoich)
        end
    end
    for i in eachindex(allocs); allocs[i] = reqs[i] == .0 ? Inf : floor(allocs[i] / reqs[i]) end
    foreach(i -> qs[i] = min(qs[i], minimum(allocs[:, i])), 1:size(reqs, 2))
    foreach(i -> allocs[:, i] .= reqs[:, i] * qs[i], 1:size(reqs, 2))

    qs
end

"Evolve transitions, spawn new transitions."
function evolve!(u, state)
    update_u!(state, u); actual_allocs = zero(u)
    
    ## schedule new transitions
    reqs = zeros(nparts(state, :S), nparts(state, :T)); qs = zeros(nparts(state, :T))

    foreach(i -> qs[i] = state[i, :transRate] * state[i, :transMultiplier], 1:nparts(state, :T))
    for i in 1:nparts(state, :T)
        new_instances = rand(Poisson(max(state.solverargs[:tstep] * qs[i], 0))) + state[i, :transToSpawn]
        capacity = state[i, :transCapacity] - count(t -> t[:transHash] == state[i, :transHash], state.ongoing_transitions)
        (capacity < new_instances) && add_to_spawn!(state, state[i, :transHash], new_instances - capacity)
        qs[i] = min(capacity, new_instances)
    end

    reqs = get_reqs_init!(reqs, qs, state)
    allocs = get_allocs!(reqs, u, state, state[:, :transPriority], state.solverargs[:strategy])
    qs .= get_init_satisfied(allocs, qs, state)

    push!(state.log, (:new_transitions, state.t, [(hash, q) for (hash, q) in zip(state[:, :transHash], qs)]...))
    u .-= sum(allocs, dims=2); actual_allocs .+= sum(allocs, dims=2)

    # add spawned transitions to the heap
    for i in 1:nparts(state, :T)
        qs[i] != 0 && push!(state.ongoing_transitions, Transition(get_sampled_transition(state, i), state.t, qs[i], 0.))
    end
    
    update_u!(state, u)
    ## evolve ongoing transitions 
    reqs = zeros(nparts(state, :S), length(state.ongoing_transitions))
    qs = map(t -> t.q, state.ongoing_transitions)

    get_reqs_ongoing!(reqs, qs, state)
    allocs = get_allocs!(reqs, u, state, map(t -> t[:transPriority], state.ongoing_transitions), state.solverargs[:strategy])
    qs .= get_frac_satisfied(allocs, reqs, state)
    push!(state.log, (:saturation, state.t, [(state.ongoing_transitions[i][:transHash], qs[i]) for i in eachindex(state.ongoing_transitions)]...))
    u .-= sum(allocs, dims=2); actual_allocs .+= sum(allocs, dims=2)

    foreach(i -> state.ongoing_transitions[i].state += qs[i] * state.solverargs[:tstep], eachindex(state.ongoing_transitions))

    push!(state.log, (:allocation, state.t, actual_allocs))
    push!(state.log, (:valuation_cost, state.t, actual_allocs' * [state[i, :specCost] for i in 1:nparts(state, :S)]))
end

# execute callbacks
function event_action!(state)
    for i in 1:nparts(state, :E)
        isassigned(state.attrs[:eventTrigger], i) && isassigned(state.attrs[:eventAction], i) || continue
        v = state[i, :eventTrigger]
        q = v isa Bool ? (v ? 1 : 0) : (v isa Number ? rand(Poisson(v)) : 0)
        for _ in 1:q state[i, :eventAction] end
    end
end

# collect terminated transitions
function finish!(u, state)
    update_u!(state, u); val_reward = 0
    terminated_all = Dict{Symbol, Float64}(); terminated_success = Dict{Symbol, Float64}()
    
    ix = 1; while ix <= length(state.ongoing_transitions)
        trans_ = state.ongoing_transitions[ix]
        ((state.t-trans_.t) < trans_.trans[:transMaxLifeTime]) && trans_.state < trans_[:transCycleTime] && (ix += 1; continue)
        toks_rhs = []; for r in extract_reactants(trans_[:transRHS], state)
            i = find_index(r.species, state)
            push!(toks_rhs, UnfoldedReactant(i, r.species, context_eval(state, state.wrap_fun(r.stoich)), r.modality ∪ state[i, :specModality]))
        end
        for tok in trans_[:transLHS]
            in(:conserved, tok.modality) && (u[tok.index] += trans_.q * tok.stoich * (in(:rate, tok.modality) ? trans_[:transCycleTime] : 1))
        end
        q = trans_.state >= trans_[:transCycleTime] ? rand(Distributions.Binomial(Int(trans_.q), trans_[:transProbOfSuccess])) : 0
        foreach(tok -> (u[tok.index] += q * tok.stoich; val_reward += state[tok.index, :specReward] * q * tok.stoich), toks_rhs)

        update_u!(state, u); context_eval(state, trans_.trans[:transPostAction])
        terminated_all[trans_[:transHash]] = get(terminated_all, trans_[:transHash], 0) + trans_.q
        terminated_success[trans_[:transHash]] = get(terminated_success, trans_[:transHash], 0) + q

        ix += 1
    end

    filter!(s -> s.state < s[:transCycleTime], state.ongoing_transitions)

    push!(state.log, (:terminated_all, state.t, terminated_all...))
    push!(state.log, (:terminated_success, state.t, terminated_success...))
    push!(state.log, (:valuation_reward, state.t, val_reward))

    u
end

function free_blocked_species!(state)
    for trans in state.ongoing_transitions, tok in trans[:transLHS]
        in(:nonblock, tok.modality) && (state.u[tok.index] += q * tok.stoich)
    end
end

"""
Transform an `acs` to a `DiscreteProblem` instance, compatible with standard solvers.

# Examples
```julia
transform(DiscreteProblem, acs, schedule=schedule_weighted!)
```
"""
function transform(::Type{DiffEqBase.DiscreteProblem}, state::DyVEState; kwargs...)
    f = function (du, u, p, t)
        state = p[:__state__]
        free_blocked_species!(state); du .= state.u
        update_observables(state); sample_transitions!(state)
        evolve!(du, state); finish!(du, state); update_u!(state, du)
        event_action!(state)

        du .= state.u
        push!(state.log, (:valuation, t, du' * [state[i, :specValuation] for i in 1:nparts(state, :S)]))

        t = (state.t += state.solverargs[:tstep]); update_u!(state, du); save!(state)
        sync_p!(p, state)

        du
    end

    DiffEqBase.DiscreteProblem(f, state.u, (.0, Inf), Dict(state.p..., :__state__ => state, :__state0__ => deepcopy(state)); kwargs...)
end

## resolve tspan, tstep

function get_tcontrol(tspan, args)
    tspan isa Tuple && (tspan = tspan[2]-tspan[1])
    tunit = get(args, :tunit, oneunit(tspan))
    tspan = tspan / tunit

    tstep = get(args, :tstep, haskey(args, :tstops) ? tspan / args[:tstops] : tunit) / tunit

    ((.0, tspan), tstep)
end

"""
Transform an `acs` to a `DiscreteProblem` instance, compatible with standard solvers.

Optionally accepts initial values and parameters, which take precedence over specifications in `acs`.

# Examples
```julia
DiscreteProblem(acs, u0, p; tspan=(.0, 100.), schedule=schedule_weighted!)
```
"""
function DiffEqBase.DiscreteProblem(acs::ReactionNetwork, u0=Dict(), p=DiffEqBase.NullParameters(); kwargs...)
    assign_defaults!(acs)
    keywords = Dict{Symbol, Any}([acs[i, :metaKeyword] => acs[i, :metaVal] for i in 1:nparts(acs, :M) if isassigned(acs.attrs[:metaKeyword], i) && isassigned(acs.attrs[:metaVal], i)])
    merge!(keywords, Dict(collect(kwargs))); merge!(keywords,  Dict(:strategy => get(keywords, :alloc_strategy, :weighted)))       
    keywords[:tspan], keywords[:tstep] = get_tcontrol(keywords[:tspan], keywords)
    
    acs = remove_choose(acs); attrs, transitions, wrap_fun = compile_attrs(acs)
    state = DyVEState(acs, attrs, transitions, wrap_fun, keywords[:tspan][1]; keywords...); init_u!(state); save!(state)

    prob = transform(DiffEqBase.DiscreteProblem, state; kwargs...)

    u0 isa Dict &&
        foreach(i -> prob.u0[i] = isassigned(acs.attrs[:specName], i) && haskey(u0, acs[i, :specName]) ? u0[acs[i, :specName]] : prob.u0[i], 1:nparts(state, :S))
    p_ = p == DiffEqBase.NullParameters() ? Dict() : Dict(k => v for (k, v) in p)
    prob = remake(prob; u0=prob.u0, tspan=keywords[:tspan], dt=get(keywords, :tstep, 1),
        p=merge(prob.p, p_, Dict(:tstep => get(keywords, :tstep, 1), :strategy => get(keywords, :alloc_strategy, :weighted)))
    )

    prob
end

fetch_params(acs::ReactionNetwork) = Dict{Symbol, Float64}((acs[i, :prmName] => acs[i, :prmVal] for i in Iterators.filter(i -> isassigned(acs.attrs[:prmVal], i), 1:nparts(acs, :P))))

# EnsembleProblem's prob_func: sample initial values
function get_prob_func(prob)
    vars = prob.p[:__state__][:, :specInitUncertainty]

    prob_func = function (prob, _, _)
        prob.p[:__state__] = deepcopy(prob.p[:__state0__])
        for i in eachindex(prob.u0)
            rv = randn()*vars[i]
            prob.u0[i] = (sign(rv + prob.u0[i]) == sign(prob.u0[i])) ? rv + prob.u0[i] : prob.u0[i]
        end

        prob
    end

    prob_func
end