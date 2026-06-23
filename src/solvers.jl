using Distributions
using Random

export ReactionNetworkProblem

function get_sampled_transition(state, i)
    transition = Dict{Symbol,Any}()
    foreach(k -> push!(transition, k => state[i, k]), keys(state.transitions))

    return transition
end

# The token-selection predicate (ADR 0008) for structured species `type` in an LHS reactant
# list, or `nothing` (kind-only bind) when that reactant carries none.
function lhs_predicate(lhs, type::Symbol)
    for r in lhs
        r isa UnfoldedReactant && r.species == type && return r.predicate
    end
    return nothing
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
            if in(:rate, tok.modality) && in(tok.species, state.structured_token)
                error(
                    "Modality `:rate` is not supported for structured species in transition $(trans[:transName]).",
                )
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
    # A transition gated off this tick (deactivated or guard false, ADR 0010 §B) proposes no
    # new instances — zero its genesis quantity so it never competes for resources.
    foreach(i -> state.transitions[:transFiring][i] || (qs[i] = 0), parts(state, :T))

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

    structured_token = collect(values(inners(getagent(state, "structured"))))

    # add spawned transitions to the heap
    for i in parts(state, :T)
        if qs[i] != 0
            transition = Transition(
                string(state[i, :transName]) * "_@$(state.t)",
                i,
                get_sampled_transition(state, i),
                AbstractAlgebraicAgent[],
                AbstractAlgebraicAgent[],
                [],
                state.t,
                qs[i],
                0.0,
            )
            push!(state.ongoing_transitions, transition)

            bound = transition.bound_structured_agents
            structured_to_agents = transition.structured_to_agents

            for (j, type) in enumerate(state.acs[:, :specName])
                if type ∈ state.structured_token
                    if !isinteger(allocs[j, i])
                        error(
                            "For structured species, stoichiometry coefficient must be integer in transition $i.",
                        )
                    end

                    # ADR 0008 §B: narrow the candidate set by the LHS reactant's predicate
                    # (kind-only when none), then the unchanged priority sort + integer take.
                    pred = lhs_predicate(state.transitions[:transLHS][i], type)
                    available_species = filter(
                        a ->
                            get_species(a) == type &&
                                !isblocked(a) &&
                                matches(pred, a, state, transition),
                        structured_token,
                    )

                    # Total order (ADR 0008 inv 3): highest priority first, ties broken by the
                    # deterministic (species, creation_index) key — NOT the AA Dict / random-name
                    # order, which would make WHICH equal-priority token binds non-reproducible.
                    sort!(
                        available_species;
                        by = a -> (-priority(a, state.acs[i, :transName]), token_sortkey(state, a)),
                    )

                    ix = 1
                    while allocs[j, i] > 0 && ix <= length(available_species)
                        set_bound_transition!(available_species[ix], transition)

                        push!(bound, available_species[ix])
                        push!(structured_to_agents, type => available_species[ix])
                        add_to_log!(available_species[ix], type, state.t, transition)

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
            structured_to_agents = transition.structured_to_agents

            for (j, type) in enumerate(state.acs[:, :specName])
                if type ∈ state.structured_token
                    if !isinteger(allocs[j, i])
                        error(
                            "For structured species, stoichiometry coefficient must be integer in transition $i.",
                        )
                    end

                    # ADR 0008 §B: narrow by the in-flight transition's LHS predicate.
                    pred = lhs_predicate(transition[:transLHS], type)
                    available_species = filter(
                        a ->
                            get_species(a) == type &&
                                !isblocked(a) &&
                                matches(pred, a, state, transition),
                        structured_token,
                    )

                    # Total order (ADR 0008 inv 3): highest priority first, ties broken by the
                    # deterministic (species, creation_index) key — NOT the AA Dict / random-name
                    # order, which would make WHICH equal-priority token binds non-reproducible.
                    sort!(
                        available_species;
                        by = a -> (-priority(a, state.acs[i, :transName]), token_sortkey(state, a)),
                    )

                    ix = 1
                    while allocs[j, i] > 0 && ix <= length(available_species)
                        set_bound_transition!(available_species[ix], transition)

                        push!(bound, available_species[ix])
                        push!(structured_to_agents, type => available_species[ix])
                        add_to_log!(available_species[ix], type, state.t, transition)

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

# The legacy `event_action!` (a no-op fetch of :eventAction, the repaired CONTRACT §3.4 Inv 7
# defect) is superseded by the endogenous decision channel `fire_rules!` (ADR 0010, src/actions.jl),
# which evaluates each Rule's guard and runs its action at _step! step 10. :E rows are lifted to
# Rules at construction.

function allocate_for_move(t::Transition, s::Symbol)
    return t.bound_structured_agents ∩
           map(x -> x[2], filter(x -> x[1] == s, t.structured_to_agents))
end

function structured_rhs(expr::Expr, state, transition)
    if isexpr(expr, :macrocall) && macroname(expr) == :structured
        if length(expr.args) == 3
            expr = quote
                return $(expr.args[end])
            end

            token = context_eval(state, transition, state.wrap_fun(expr))

            entangle!(getagent(state, "structured"), token)

            return token, get_species(token)
        else
            expr = quote
                token = $(expr.args[end-1])
                species = $(expr.args[end])

                return token, species
            end

            token, species = context_eval(state, transition, state.wrap_fun(expr))
            set_species!(token, Symbol(species))

            entangle!(getagent(state, "structured"), token)

            return token, get_species(token)
        end
    elseif isexpr(expr, :macrocall) && macroname(expr) == :move
        expr = quote
            species_from = $(expr.args[end-1])
            species_to = $(expr.args[end])

            return species_from, species_to
        end

        species_from, species_to =
            Symbol.(context_eval(state, transition, state.wrap_fun(expr)))

        tokens =
            filter(x -> get_species(x) == species_from, transition.bound_structured_agents)

        if !isempty(tokens)
            token = first(tokens)
            entangle!(getagent(state, "structured"), token)

            set_species!(token, species_to)
            ix = findfirst(
                i -> transition.bound_structured_agents[i] == token,
                eachindex(transition.bound_structured_agents),
            )
            deleteat!(transition.bound_structured_agents, ix)
            set_bound_transition!(token, nothing)

            return token, species_to
        else
            @error "Not enough tokens to allocate for a move."
        end

    elseif isexpr(expr, :macrocall) && macroname(expr) == :advance
        # @advance(field, value): advance a bound token's lifecycle by writing one field, keeping
        # its identity/kind/uuid/creation_index/past_bonds (ADR 0008 §D). The phase-as-attribute
        # generalization of @move (which writes the `species` field). `value` may read the token's
        # own current fields via @field(name), and MAY draw (§F). The advanced token is then
        # released. It is consumed from the first bound token of this transition.
        field = expr.args[3]
        field isa Symbol || error("@advance: first argument must be a field name, got $field")
        valex = expr.args[4]
        # No bound token to advance (the @select predicate matched nothing this firing) — a
        # silent no-op: the instance produced no advance. finish! skips the nothing return.
        isempty(transition.bound_structured_agents) && return nothing, nothing
        token = first(transition.bound_structured_agents)
        # Evaluate the value with the bound token in scope so @field(name) reads its attributes.
        val = eval_with_token(state, transition, token, valex)
        if field === :species
            set_species!(token, Symbol(val))
        else
            setproperty!(token, field, val)
        end
        deleteat!(transition.bound_structured_agents, 1)
        set_bound_transition!(token, nothing)
        return token, get_species(token)

    else
        token = context_eval(state, transition, state.wrap_fun(expr))
        entangle!(getagent(state, "structured"), token)

        return token, get_species(token)
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
            rand(state.rng, Distributions.Binomial(Int(trans_.q), trans_[:transProbOfSuccess]))
        else
            0
        end

        for r in extract_reactants(trans_[:transRHS], state)
            if r.species isa Expr
                stoich = context_eval(state, trans_, state.wrap_fun(r.stoich))

                for _ = 1:(q*stoich)
                    token, species = structured_rhs(r.species, state, trans_)
                    # A structured-RHS op may legitimately produce nothing (e.g. @advance with no
                    # bound token to advance) — skip the count/reward in that case.
                    species === nothing && continue
                    i = find_index(species, state)
                    state.u[i] += 1
                    val_reward += state[i, :specReward]
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
                if tok.species ∈ state.structured_token
                    for _ = 1:(trans_.q*tok.stoich)
                        isempty(trans_.bound_structured_agents) && break
                        agent_ix = findfirst(
                            a -> get_species(a) == tok.species,
                            trans_.bound_structured_agents,
                        )

                        set_bound_transition!(
                            trans_.bound_structured_agents[agent_ix],
                            nothing,
                        )
                        deleteat!(trans_.bound_structured_agents, agent_ix)
                    end
                end
            end

            if in(:nonblock, tok.modality)
                if in(:conserved, tok.modality)
                    error(
                        "Modalities `:conserved` and `:nonblock` cannot be specified at the same time.",
                    )
                end

                state.u[tok.index] += trans_.q * tok.stoich
                if tok.species ∈ state.structured_token
                    for _ = 1:(trans_.q*tok.stoich)
                        agent_ix = findfirst(
                            a -> get_species(a) == tok.species,
                            trans_.nonblock_structured_agents,
                        )

                        set_bound_transition!(
                            trans_.nonblock_structured_agents[agent_ix],
                            nothing,
                        )
                        deleteat!(trans_.nonblock_structured_agents, agent_ix)
                    end
                end
            end
        end

        context_eval(state, trans_, state.wrap_fun(state.acs[trans_.i, :transPostAction]))

        for agent in trans_.bound_structured_agents
            set_species!(agent, :removed)
            set_bound_transition!(agent, nothing)
        end

        terminated_all[Symbol(trans_[:transHash])] =
            get(terminated_all, Symbol(trans_[:transHash]), 0) + trans_.q

        terminated_success[Symbol(trans_[:transHash])] =
            get(terminated_success, Symbol(trans_[:transHash]), 0) + q

        ix += 1
    end

    # Prune every instance that passed the terminal test above — i.e. keep only those that
    # have neither completed their cycle NOR aged out. This must mirror the skip condition at
    # the top of the loop; the old predicate ignored max-lifetime, so timed-out instances
    # (state < cycleTime but age ≥ maxLifeTime) were retained and re-emitted/​re-credited every
    # subsequent tick (violating conservation + termination-completeness, §3.4 INV2/INV6).
    filter!(
        s ->
            ((state.t - s.t) < s[:transMaxLifeTime]) && (s.state < s[:transCycleTime]),
        state.ongoing_transitions,
    )

    push!(state.log, (:terminated_all, state.t, terminated_all...))
    push!(state.log, (:terminated_success, state.t, terminated_success...))
    push!(state.log, (:valuation_reward, state.t, val_reward))

    return state.u
end

function free_blocked_species!(state)
    for trans in state.ongoing_transitions, tok in trans[:transLHS]
        in(:nonblock, tok.modality) && (state.u[tok.index] += trans.q * tok.stoich)
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

    # Determinism (§4 D5/D6): build the state-owned RNG. A `seed` kwarg fixes the stream;
    # absent it, draw a fresh seed from system entropy so a default run is still self-contained.
    # The REALIZED seed is stored on the struct (so an entropy-seeded run is replayable, D6)
    # and `initial_rng` snapshots the stream at t=0 so `_reinit!` restores it exactly (D7).
    # Any `Integer` seed is accepted verbatim (e.g. a `hash((root, k))` ensemble member key, §4 D8).
    seed = get(keywords, :seed, nothing)
    seed = isnothing(seed) ? rand(Random.RandomDevice(), UInt64) : seed
    rng = Random.Xoshiro(seed)
    initial_rng = copy(rng)

    acs = remove_choose(acs)

    structured_token_names =
        acs[filter(i -> acs[i, :specStructured], 1:nparts(acs, :S)), :specName]

    attrs, transitions, wrap_fun = compile_attrs(acs, structured_token_names)
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
        ) ∪ [:transLHS, :transRHS, :transToSpawn, :transHash, :transFiring]
    transitions = Dict{Symbol,Vector}(a => [] for a in transitions_attrs)

    sol = DataFrame(
        "t" => Float64[],
        (string(name) => Float64[] for name in acs[:, :specName])...,
    )

    # Endogenous decision channel (ADR 0010 §12). Per-network host registry for AddToken/Invoke
    # (ADR 0006 §C — by-name, never eval'd). Rules are built from the :E rows: a legacy event
    # `trigger && action` becomes a Rule{guard=trigger, action=RawExpr(action), every_tick}.
    # Typed Rules can also be supplied directly via the `rules=` kwarg / @rule authoring.
    registry = Dict{Symbol,Any}(get(keywords, :registry, Dict{Symbol,Any}()))
    rules = Any[
        Rule(Symbol("rule_", i), acs[i, :eventTrigger], RawExpr(acs[i, :eventAction]))
        for i in parts(acs, :E) if
        !isnothing(acs[i, :eventTrigger]) && !isnothing(acs[i, :eventAction])
    ]
    append!(rules, get(keywords, :rules, Any[]))

    network = ReactionNetworkProblem(
        name,
        acs,
        attrs,
        transition_recipes,
        u0_init,
        merge(p, Dict(:strategy => get(keywords, :alloc_strategy, :weighted))),
        keywords[:tspan][1],
        structured_token_names,
        keywords[:tspan],
        get(keywords, :tstep, 1),
        transitions,
        ongoing_transitions,
        log,
        observables,
        wrap_fun,
        sol,
        rng,
        seed,
        initial_rng,
        rules,
        registry,
        Dict{Symbol,Int}(),
        Dict{String,Int}(),
    )

    entangle!(network, FreeAgent("structured"))

    # save!(network)

    return network
end

function AlgebraicAgents._reinit!(state::ReactionNetworkProblem)
    state.u .= isempty(state.sol) ? state.u : Vector(state.sol[1, 2:end])
    state.t = state.tspan[1]
    empty!(state.ongoing_transitions)
    empty!(state.log)
    state.observables = compile_observables(state.acs)
    empty!(state.sol)
    # Restore the RNG to its construction state so the second run reproduces the first (§4 D7).
    state.rng = copy(state.initial_rng)
    # Reset every `once` rule's latch so a re-run from the same seed reproduces the lever (§4 D7).
    for r in state.rules
        r.fire_mode === :once && (r.enabled = true)
    end
    # Reset structured-token creation counters (the token population itself is rebuilt in Stage D).
    empty!(state.creation_counters)
    empty!(state.creation_index)

    return state
end

function update_u_structured!(state)
    structured_tokens = collect(values(inners(getagent(state, "structured"))))
    for (i, species) in enumerate(state.acs[:, :specName])
        if state.acs[i, :specStructured]
            state.u[i] =
                count(a -> get_species(a) == species && !isblocked(a), structured_tokens)
        end
    end

    return state.u
end

function AlgebraicAgents._step!(state::ReactionNetworkProblem)
    update_u_structured!(state)
    if isempty(state.sol)
        save!(state)
    end

    free_blocked_species!(state)
    update_u_structured!(state)
    update_observables(state)
    sample_transitions!(state)
    evolve!(state)
    update_u_structured!(state)
    finish!(state)
    update_u_structured!(state)

    # Step 10 (§3.3): fire the endogenous decision channel (ADR 0010). Rules see this tick's
    # post-finish state; their writes land on this tick's ledger row and the next tick's
    # genesis/guards. Replaces the old no-op event_action! slot.
    fire_rules!(state)
    update_u_structured!(state)

    push!(
        state.log,
        (
            :valuation,
            state.t,
            state.u' * [state[i, :specValuation] for i in parts(state, :S)],
        ),
    )

    state.t += state.dt

    save!(state)

    return state.t
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
