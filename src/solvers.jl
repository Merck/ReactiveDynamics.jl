using Distributions
using Random

export ReactionNetworkProblem

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
    for i in axes(reqs, 2)
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
    for i in eachindex(state.ongoing_transitions)
        for tok in state.ongoing_transitions[i][:transLHS]
            in(:rate, tok.modality) &&
                (state.ongoing_transitions[i][:transCycleTime] > 0) &&
                (reqs[tok.index, i] += qs[i] * tok.stoich * state.dt)
            if in(:rate, tok.modality) && in(tok.species, state.structured_species)
                error("Modality `:rate` is not supported for structured species in transition $(trans[:transName]).")
            end
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
    for i in axes(reqs, 1)
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
    for i in axes(reqs, 1)
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

isinteger(x::Number) = x == trunc(x)

"""
Given available allocations and qties of transitions requested to spawn, return number of spawned transitions. Update `alloc` to match actual allocation.
"""
function get_init_satisfied(allocs, qs, state)
    reqs = zero(allocs)
    for i in axes(allocs, 2)
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
        parts(state, :T),
    )
    qs .= ceil.(Ref(Int), qs)

    for i in parts(state, :T)
        new_instances = qs[i] + state[i, :transToSpawn]
        capacity =
            state[i, :transCapacity] -
            count(t -> t[:transHash] == state[i, :transHash], state.ongoing_transitions)
        (capacity < new_instances) &&
            add_to_spawn!(state, state[i, :transHash], new_instances - capacity)
        qs[i] = min(capacity, new_instances)
    end

    reqs = get_reqs_init!(reqs, qs, state)

    allocs = get_allocs!(reqs, state.u, state, state[:, :transPriority], state.p[:strategy])

    qs .= get_init_satisfied(allocs, qs, state)

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
    for i in parts(state, :T)
        if qs[i] != 0
            transition = Transition(
                string(state[i, :transName]) * "_@$(state.t)",
                i,
                get_sampled_transition(state, i),
                AbstractAlgebraicAgent[],
                AbstractAlgebraicAgent[],
                state.t,
                qs[i],
                0.0,
            )
            push!(state.ongoing_transitions, transition)

            bound = transition.bound_structured_agents

            for (j, type) in enumerate(state.acs[:, :specName])
                if type ∈ state.structured_species
                    if !isinteger(allocs[j, i])
                        error("For structured species, stoichiometry coefficient must be integer in transition $i.")
                    end

                    all_of_type = collect(values(inners(getagent(state, "structured/$(string(type))"))))
                    filter!(!isblocked, all_of_type)
                    sort!(all_of_type, by=a->priority(a, state.acs[i, :transName]), rev=true)

                    ix = 1
                    while allocs[j, i] > 0 && ix <= length(all_of_type)
                        all_of_type[ix].bound_transition = transition
                        push!(bound, all_of_type[ix])
                        allocs[j, i] -= 1
                        ix += 1
                    end
                end
            end

            context_eval(state, transition, state.wrap_fun(state.acs[i, :transPreAction]))
        end
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

    for i in eachindex(state.ongoing_transitions)
        transition = state.ongoing_transitions[i]
        if qs[i] != 0
            transition.state += qs[i] * state.dt
            bound = transition.nonblock_structured_agents

            for (j, type) in enumerate(state.acs[:, :specName])
                if type ∈ state.structured_species
                    if !isinteger(allocs[j, i])
                        error("For structured species, stoichiometry coefficient must be integer in transition $i.")
                    end

                    all_of_type = collect(values(inners(getagent(state, "structured/$(string(type))"))))
                    filter!(!isblocked, all_of_type)
                    sort!(all_of_type, by=a -> priority(a, state.acs[i, :transName]), rev=true)

                    ix = 1
                    while allocs[j, i] > 0 && ix <= length(all_of_type)
                        all_of_type[ix].bound_transition = transition
                        push!(bound, all_of_type[ix])
                        allocs[j, i] -= 1
                        ix += 1
                    end
                end
            end
        end
    end

    push!(state.log, (:allocation, state.t, actual_allocs))
    return push!(
        state.log,
        (
            :valuation_cost,
            state.t,
            actual_allocs' * [state[i, :specCost] for i in parts(state, :S)],
        ),
    )
end

# execute callbacks
function event_action!(state)
    for i in parts(state, :E)
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
            (trans_.state < trans_[:transCycleTime]) &&
            (ix += 1; continue)
        
        q = if trans_.state >= trans_[:transCycleTime]
            rand(Distributions.Binomial(Int(trans_.q), trans_[:transProbOfSuccess]))
        else
            0
        end

        for r in extract_reactants(trans_[:transRHS], state)
            if r.species isa Expr
                if !(r.species.args[1] ∈ state.structured_species)
                    error("Unresolved structured species species $(r.species.args[1]).")
                end

                i = find_index(r.species.args[1], state)
                stoich = context_eval(state, trans_, state.wrap_fun(r.stoich))

                state.u[i] += q * stoich
                val_reward += state[i, :specReward] * q * stoich

                for _ in 1:q
                    a = context_eval(state, trans_, state.wrap_fun(r.species))
                    entangle!(getagent(state, "structured/$(r.species.args[1])"), a)
                end
            else
                i = find_index(r.species, state)
                stoich = context_eval(state, trans_, state.wrap_fun(r.stoich))

                state.u[i] += q * stoich
                val_reward += state[i, :specReward] * q * stoich
            end
        end

        for tok in trans_[:transLHS]
            if in(:conserved, tok.modality)
                state.u[tok.index] +=
                    trans_.q *
                    tok.stoich *
                    (in(:rate, tok.modality) ? trans_[:transCycleTime] : 1)
                if tok.species ∈ state.structured_species
                    for _ in 1:(trans_.q * tok.stoich)
                        trans_.bound_structured_agents[begin].bound_transition = nothing
                        deleteat!(trans_.bound_structured_agents, 1)
                    end
                end
            end

            if in(:nonblock, tok.modality)
                if in(:conserved, tok.modality)
                    error("Modalities `:conserved` and `:nonblock` cannot be specified at the same time.")
                end

                state.u[tok.index] += trans_.q * tok.stoich 
                if tok.species ∈ state.structured_species
                    for _ in 1:(trans_.q * tok.stoich)
                        trans_.nonblock_structured_agents[begin].bound_transition = nothing
                        deleteat!(trans_.nonblock_structured_agents, 1)
                    end
                end
            end
        end

        context_eval(state, trans_, state.wrap_fun(state.acs[trans_.i, :transPostAction]))

        terminated_all[Symbol(trans_[:transHash])] =
            get(terminated_all, Symbol(trans_[:transHash]), 0) + trans_.q
        
        terminated_success[Symbol(trans_[:transHash])] =
            get(terminated_success, Symbol(trans_[:transHash]), 0) + q

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

    for trans in state.ongoing_transitions
        for a in trans.nonblock_structured_agents
            a.bound_transition = nothing
        end

        empty!(trans.nonblock_structured_agents)
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

function ReactionNetworkProblem(
    acs::ReactionNetworkSchema,
    u0 = Dict(),
    p = Dict();
    name = "reaction_network",
    kwargs...,
)
    assign_defaults!(acs)
    keywords = Dict{Symbol,Any}([
        acs[i, :metaKeyword] => acs[i, :metaVal] for i in parts(acs, :M) if
        !isnothing(acs[i, :metaKeyword]) && !isnothing(acs[i, :metaVal])
    ])
    merge!(keywords, Dict(collect(kwargs)))
    merge!(keywords, Dict(:strategy => get(keywords, :alloc_strategy, :weighted)))

    keywords[:tspan], keywords[:tstep] = get_tcontrol(keywords[:tspan], keywords)

    acs = remove_choose(acs)

    structured_species_names = acs[filter(i -> acs[i, :specStructured], 1:nparts(acs, :S)), :specName]

    attrs, transitions, wrap_fun = compile_attrs(acs, structured_species_names)
    transition_recipes = transitions
    u0_init = zeros(nparts(acs, :S))

    for i in parts(acs, :S)
        if !isnothing(acs[i, :specName]) && haskey(u0, acs[i, :specName])
            u0_init[i] = u0[acs[i, :specName]]
        else
            u0_init[i] = acs[i, :specInitVal]
        end
    end

    prms = Dict{Symbol,Any}((
        acs[i, :prmName] => acs[i, :prmVal] for
        i in Iterators.filter(i -> !isnothing(acs[i, :prmVal]), 1:nparts(acs, :P))
    ))

    merge!(p, prms)

    ongoing_transitions = Transition[]
    log = NamedTuple[]
    observables = compile_observables(acs)
    transitions_attrs =
        setdiff(
            filter(a -> contains(string(a), "trans"), propertynames(acs.subparts)),
            (:trans,),
        ) ∪ [:transLHS, :transRHS, :transToSpawn, :transHash]
    transitions = Dict{Symbol,Vector}(a => [] for a in transitions_attrs)

    sol = DataFrame(
        "t" => Float64[],
        (string(name) => Float64[] for name in acs[:, :specName])...,
    )
    network = ReactionNetworkProblem(
        name,
        acs,
        attrs,
        transition_recipes,
        u0_init,
        merge(p, Dict(:strategy => get(keywords, :alloc_strategy, :weighted))),
        keywords[:tspan][1],
        Symbol[],
        keywords[:tspan],
        get(keywords, :tstep, 1),
        transitions,
        ongoing_transitions,
        log,
        observables,
        wrap_fun,
        sol,
    )

    entangle!(network, FreeAgent("structured"))

    structured_species = filter(i -> acs[i, :specStructured], 1:nparts(acs, :S))
    for i in structured_species
        push!(network.structured_species, acs[i, :specName])
        entangle!(getagent(network, "structured"), FreeAgent(string(acs[i, :specName])))
    end

    save!(network)

    return network
end

function AlgebraicAgents._reinit!(state::ReactionNetworkProblem)
    state.u .= isempty(state.sol) ? state.u : Vector(state.sol[1, 2:end])
    state.t = state.tspan[1]
    empty!(state.ongoing_transitions)
    empty!(state.log)
    state.observables = compile_observables(state.acs)
    empty!(state.sol)

    return state
end

function update_u_structured!(state)
    for (i, type) in enumerate(state.acs[:, :specName])
        structured_agents_type = values(inners(getagent(state, "structured/$type")))
        state.u[i] = count(!isblocked, structured_agents_type)
    end
end

function AlgebraicAgents._step!(state::ReactionNetworkProblem)
    free_blocked_species!(state)
    update_observables(state)
    sample_transitions!(state)
    evolve!(state)
    finish!(state)

    event_action!(state)

    push!(
        state.log,
        (
            :valuation,
            state.t,
            state.u' * [state[i, :specValuation] for i in parts(state, :S)],
        ),
    )

    save!(state)
    return state.t += state.dt
end

function AlgebraicAgents._projected_to(state::ReactionNetworkProblem)
    return state.t > state.tspan[2] ? true : state.t
end

function fetch_params(acs::ReactionNetworkSchema)
    return Dict{Symbol,Any}((
        acs[i, :prmName] => acs[i, :prmVal] for
        i in Iterators.filter(i -> !isnothing(acs[i, :prmVal]), parts(acs, :P))
    ))
end
