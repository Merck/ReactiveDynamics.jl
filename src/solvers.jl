using Distributions
using Random

function get_sampled_transition(state, i)
    transition = Dict{Symbol,Any}()
    foreach(k -> push!(transition, k => state[i, k]), keys(state.transitions))

    return transition
end

"""
Compute resource requirements given transition quantities.
"""
function get_reqs_init!(reqs, qs, state)
    reqs .= 0.0
    for i = 1:size(reqs, 2)
        for tok in state[i, :transLHS]
            !any(m -> m in tok.modality, [:rate, :nonblock]) &&
                (reqs[tok.index, i] += qs[i] * tok.stoich)
        end
    end

    return reqs
end

"""
Compute resource requirements given transition quantities.
"""
function get_reqs_ongoing!(reqs, qs, state)
    reqs .= 0.0
    for i = 1:length(state.ongoing_transitions)
        for tok in state.ongoing_transitions[i][:transLHS]
            in(:rate, tok.modality) &&
                (state.ongoing_transitions[i][:transCycleTime] > 0) &&
                (reqs[tok.index, i] += qs[i] * tok.stoich * state.dt)
            in(:nonblock, tok.modality) && (reqs[tok.index, i] += qs[i] * tok.stoich)
        end
    end

    return reqs
end

"""
Given requirements, return available allocation.
"""
function get_allocs!(reqs, u, state, priorities, strategy = :weighted)
    return if strategy == :weighted
        alloc_weighted!(reqs, u, priorities, state)
    else
        alloc_greedy!(reqs, u, priorities, state)
    end
end

function alloc_weighted!(reqs, u, priorities, state)
    allocs = zero(reqs)
    for i = 1:size(reqs, 1)
        s = sum(reqs[i, :])
        u[i] >= s && (allocs[i, :] .= reqs[i, :]; continue)
        foreach(j -> allocs[i, j] = reqs[i, j] * priorities[j], 1:size(reqs, 2))
        s = sum(allocs[i, :])
        allocs[i, :] .*= (s == 0) ? 0.0 : (max(0, u[i]) / s)
    end

    return allocs
end

function alloc_greedy!(reqs, u, priorities, state)
    allocs = zero(reqs)
    sorted_trans = sort(1:size(reqs, 2); by = i -> -priorities[i])
    for i = 1:size(reqs, 1)
        s = sum(reqs[i, :])
        u[i] >= s && (allocs[i, :] .= reqs[i, :]; continue)
        a = u[i]
        j = 1
        while a > 0 && j <= size(reqs, 2)
            allocs[i, sorted_trans[j]] = min(reqs[i, sorted_trans[j]], a)
            a -= allocs[i, sorted_trans[j]]
            j += 1
        end
    end

    return allocs
end

"""
Given resource requirements and available allocations, output resulting shift size for each transition.
"""
function get_frac_satisfied(allocs, reqs, state)
    for i in eachindex(allocs)
        allocs[i] = min(1, (reqs[i] == 0.0 ? 1 : (allocs[i] / reqs[i])))
    end
    qs = vec(minimum(allocs; dims = 1))
    foreach(i -> allocs[:, i] .= reqs[:, i] * qs[i], 1:size(reqs, 2))

    return qs
end

"""
Given available allocations and qties of transitions requested to spawn, return number of spawned transitions. Update `alloc` to match actual allocation.
"""
function get_init_satisfied(allocs, qs, state)
    reqs = zero(allocs)
    for i = 1:size(allocs, 2)
        all(allocs[:, i] .>= 0) || (allocs[:, i] .= 0.0; qs[i] = 0)
        for tok in state[i, :transLHS]
            !any(m -> m in tok.modality, [:rate, :nonblock]) &&
                (reqs[tok.index, i] += tok.stoich)
        end
    end
    for i in eachindex(allocs)
        allocs[i] = reqs[i] == 0.0 ? Inf : floor(allocs[i] / reqs[i])
    end
    foreach(i -> qs[i] = min(qs[i], minimum(allocs[:, i])), 1:size(reqs, 2))
    foreach(i -> allocs[:, i] .= reqs[:, i] * qs[i], 1:size(reqs, 2))

    return qs
end

"""
Evolve transitions, spawn new transitions.
"""
function evolve!(state)
    actual_allocs = zero(state.u)

    ## schedule new transitions
    reqs = zeros(nparts(state, :S), nparts(state, :T))
    qs = zeros(nparts(state, :T))

    foreach(
        i -> qs[i] = state[i, :transRate] * state[i, :transMultiplier],
        1:nparts(state, :T),
    )
    qs .= ceil.(Ref(Int), qs)
    @show qs
    for i = 1:nparts(state, :T)
        new_instances = state.dt * qs[i] + state[i, :transToSpawn]
        capacity =
            state[i, :transCapacity] -
            count(t -> t[:transHash] == state[i, :transHash], state.ongoing_transitions)
        (capacity < new_instances) &&
            add_to_spawn!(state, state[i, :transHash], new_instances - capacity)
        qs[i] = min(capacity, new_instances)
    end

    reqs = get_reqs_init!(reqs, qs, state)
    @show reqs
    allocs =
        get_allocs!(reqs, state.u, state, state[:, :transPriority], state.p[:strategy])
    @show allocs
    qs .= get_init_satisfied(allocs, qs, state)
    @show qs
    println("====")
    push!(
        state.log,
        (
            :new_transitions,
            state.t,
            [(hash, q) for (hash, q) in zip(state[:, :transHash], qs)]...,
        ),
    )
    state.u .-= sum(allocs; dims = 2)
    actual_allocs .+= sum(allocs; dims = 2)

    # add spawned transitions to the heap
    for i = 1:nparts(state, :T)
        qs[i] != 0 && push!(
            state.ongoing_transitions,
            Transition(state[i, :transName] * "_@$(state.t)", get_sampled_transition(state, i), state.t, qs[i], 0.0),
        )
    end

    ## evolve ongoing transitions 
    reqs = zeros(nparts(state, :S), length(state.ongoing_transitions))
    qs = map(t -> t.q, state.ongoing_transitions)

    get_reqs_ongoing!(reqs, qs, state)
    allocs = get_allocs!(
        reqs,
        state.u,
        state,
        map(t -> t[:transPriority], state.ongoing_transitions),
        state.p[:strategy],
    )
    qs .= get_frac_satisfied(allocs, reqs, state)
    push!(
        state.log,
        (
            :saturation,
            state.t,
            [
                (state.ongoing_transitions[i][:transHash], qs[i]) for
                i in eachindex(state.ongoing_transitions)
            ]...,
        ),
    )
    state.u .-= sum(allocs; dims = 2)
    actual_allocs .+= sum(allocs; dims = 2)

    foreach(
        i -> state.ongoing_transitions[i].state += qs[i] * state.dt,
        eachindex(state.ongoing_transitions),
    )

    push!(state.log, (:allocation, state.t, actual_allocs))
    return push!(
        state.log,
        (
            :valuation_cost,
            state.t,
            actual_allocs' * [state[i, :specCost] for i = 1:nparts(state, :S)],
        ),
    )
end

# execute callbacks
function event_action!(state)
    for i = 1:nparts(state, :E)
        !isnothing(state[i, :eventTrigger]) && !isnothing(state[i, :eventAction]) ||
            continue
        v = state[i, :eventTrigger]
        q = v isa Bool ? (v ? 1 : 0) : (v isa Number ? rand(Poisson(v)) : 0)
        for _ = 1:q
            state[i, :eventAction]
        end
    end
end

# collect terminated transitions
function finish!(state)
    val_reward = 0
    terminated_all = Dict{Symbol,Float64}()
    terminated_success = Dict{Symbol,Float64}()

    ix = 1
    while ix <= length(state.ongoing_transitions)
        trans_ = state.ongoing_transitions[ix]
        ((state.t - trans_.t) < trans_.trans[:transMaxLifeTime]) &&
            trans_.state < trans_[:transCycleTime] &&
            (ix += 1; continue)
        toks_rhs = []
        for r in extract_reactants(trans_[:transRHS], state)
            i = find_index(r.species, state)
            push!(
                toks_rhs,
                UnfoldedReactant(
                    i,
                    r.species,
                    context_eval(state, state.wrap_fun(r.stoich)),
                    r.modality ∪ state[i, :specModality],
                ),
            )
        end
        for tok in trans_[:transLHS]
            in(:conserved, tok.modality) && (
                state.u[tok.index] +=
                    trans_.q *
                    tok.stoich *
                    (in(:rate, tok.modality) ? trans_[:transCycleTime] : 1)
            )
        end
        q = if trans_.state >= trans_[:transCycleTime]
            rand(Distributions.Binomial(Int(trans_.q), trans_[:transProbOfSuccess]))
        else
            0
        end
        foreach(
            tok -> (state.u[tok.index] += q * tok.stoich;
            val_reward += state[tok.index, :specReward] * q * tok.stoich),
            toks_rhs,
        )

        context_eval(state, trans_.trans[:transPostAction])
        terminated_all[trans_[:transHash]] =
            get(terminated_all, trans_[:transHash], 0) + trans_.q
        terminated_success[trans_[:transHash]] =
            get(terminated_success, trans_[:transHash], 0) + q

        ix += 1
    end

    filter!(s -> s.state < s[:transCycleTime], state.ongoing_transitions)

    push!(state.log, (:terminated_all, state.t, terminated_all...))
    push!(state.log, (:terminated_success, state.t, terminated_success...))
    push!(state.log, (:valuation_reward, state.t, val_reward))

    return state.u
end

function free_blocked_species!(state)
    for trans in state.ongoing_transitions, tok in trans[:transLHS]
        in(:nonblock, tok.modality) && (state.u[tok.index] += q * tok.stoich)
    end
end


## resolve tspan, tstep

function get_tcontrol(tspan, args)
    tspan isa Tuple && (tspan = tspan[2] - tspan[1])
    tunit = get(args, :tunit, oneunit(tspan))
    tspan = tspan / tunit

    dt = get(args, :dt, haskey(args, :tstops) ? tspan / args[:tstops] : tunit) / tunit

    return ((0.0, tspan), dt)
end

function ReactiveNetwork(
    acs::ReactionNetwork,
    u0 = Dict(),
    p = Dict();
    name = "reactive_network",
    kwargs...,
)
    assign_defaults!(acs)
    keywords = Dict{Symbol,Any}([
        acs[i, :metaKeyword] => acs[i, :metaVal] for i = 1:nparts(acs, :M) if
        !isnothing(acs[i, :metaKeyword]) && !isnothing(acs[i, :metaVal])
    ])
    merge!(keywords, Dict(collect(kwargs)))
    merge!(keywords, Dict(:strategy => get(keywords, :alloc_strategy, :weighted)))

    keywords[:tspan], keywords[:tstep] = get_tcontrol(keywords[:tspan], keywords)

    acs = remove_choose(acs)
    attrs, transitions, wrap_fun = compile_attrs(acs)
    transition_recipes = transitions
    u0_init = zeros(nparts(acs, :S))

    for i in 1:nparts(acs, :S)
        if !isnothing(acs[i, :specName]) && haskey(u0, acs[i, :specName])
            u0_init[i] = u0[acs[i, :specName]]
        else
            u0_init[i] = acs[i, :specInitVal]
        end
    end

    ongoing_transitions = Transition[]
    log = NamedTuple[]
    observables = compile_observables(acs)
    transitions_attrs =
        setdiff(
            filter(a -> contains(string(a), "trans"), propertynames(acs.subparts)),
            (:trans,),
        ) ∪ [:transLHS, :transRHS, :transToSpawn, :transHash]
    transitions = Dict{Symbol,Vector}(a => [] for a in transitions_attrs)

    network =  ReactiveNetwork(
        name,
        acs,
        attrs,
        transition_recipes,
        u0_init, 
        merge(
            p,
            Dict(
                :tstep => get(keywords, :tstep, 1),
                :strategy => get(keywords, :alloc_strategy, :weighted),
            ),
        ),
        keywords[:tspan][1],
        keywords[:tspan],
        get(keywords, :tstep, 1),
        transitions,
        ongoing_transitions,
        log,
        observables,
        kwargs,
        wrap_fun,
        Vector{Float64}[],
        Float64[],
    )

    save!(network)

    return network
end

function AlgebraicAgents.step!(state::ReactiveNetwork)
    #du = copy(state.u)
    free_blocked_species!(state)
    #du .= state.u
    update_observables(state)
    sample_transitions!(state)
    evolve!(state)
    finish!(state)
    #update_u!(state, u)
    event_action!(state)

    push!(
        state.log,
        (:valuation, state.t, state.u' * [state[i, :specValuation] for i = 1:nparts(state, :S)]),
    )

    #update_u!(state, du)
    save!(state)

    #state.u .= du
    state.t += state.dt
end

function AlgebraicAgents._projected_to(state::ReactiveNetwork)
    if state.t >= state.tspan[2]
        true
    else
        state.t
    end
end

function fetch_params(acs::ReactionNetwork)
    return Dict{Symbol,Any}((
        acs[i, :prmName] => acs[i, :prmVal] for
        i in Iterators.filter(i -> !isnothing(acs[i, :prmVal]), 1:nparts(acs, :P))
    ))
end