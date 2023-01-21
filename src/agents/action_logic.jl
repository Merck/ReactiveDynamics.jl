mutable struct DyVEState
    acs::ReactionNetwork

    attrs::Dict{Symbol, Vector}
    transition_recipes::Dict{Symbol, Vector}

    u::Vector{Float64}
    p::Any
    t::Float64

    transitions::Dict{Symbol, Vector}
    ongoing_transitions::Vector{Transition}
    log::Vector{Tuple}

    sampleables::Dict{Symbol, Sampleable}
    solverargs::Any

    wrap_fun::Any
    history_u::Vector{Vector{Float64}}
    history_t::Vector{Float64}

    function DyVEState(acs::ReactionNetwork, attrs, transition_recipes, wrap_fun, t0 = 0;
                       kwargs...)
        ongoing_transitions = Transition[]
        log = NamedTuple[]
        sampleables = compile_observables(acs)
        transitions_attrs = setdiff(filter(a -> contains(string(a), "trans"),
                                           keys(acs.attrs)), (:trans,)) ∪
                            [:transLHS, :transRHS, :transToSpawn, :transHash]
        transitions = Dict{Symbol, Vector}(a => [] for a in transitions_attrs)

        new(acs, attrs, transition_recipes, zeros(nparts(acs, :S)), fetch_params(acs), t0,
            transitions, ongoing_transitions, log, sampleables, kwargs, wrap_fun,
            Vector{Float64}[], Float64[])
    end
end

function get_modalities(path::AbstractString, network)
    get(getpath(network, path).p, :resModality, Set{Symbol}())
end

macro register_demand(demand, agent, uuid, q, priority)
    quote
        $(esc(demand))[($(esc(agent)), $(esc(uuid)))] = get($(esc(demand)),
                                                            ($(esc(agent)), $(esc(uuid))),
                                                            0.0) + esc($q)
    end
end

function get_requirements(t::ReactionNetworkTransition, network, dt, time)
    # requirements
    demand = Dict{String, Any}()
    empty!(t.allocs)
    ## new transitions
    transRate = t.params_interpreted(network, t)
    i = (t[:rate] * dt) |> Poisson |> rand
    if i > 0 ## should spawn new transitions? resample parameters and calculate requirements
        resample!(t, network)
        t[:rate] = transRate
        uuid = "new-initial"
        for r in t[:LHS]
            !isa(r.agent, AbstractString) && continue
            any(m -> m in r.modalities, [:rate, :nonblock]) && continue

            @register_demand demand r.agent uuid i*r_stoich t[:priority]
        end
        push!(t.ongoing_transitions,
              Transition("new-nominal", copy(t.params_sampled), time, 1.0, 0.0))
    else
        t[:LHS] = []
        push!(t.ongoing_transitions,
              Transition("new-nominal", copy(t.params_sampled), time, 0.0, 0.0))
    end

    ## ongoing transitions
    for o in t.ongoing_transitions
        for r in o.p[:LHS]
            if :rate ∈ r.modalities && (o.p[:cycle_time] > 0)
                @register_demand demand r.agent o.uuid o.q*r.stoich*dt
            end
            if :nonblock ∈ r.modality
                @register_demand demand r.agent o.uuid o.q*r.stoich o.p[:priority]
            end
        end
    end

    for (k, v) in demand
        register_demand!(getpath(network, k[1]), (k[2], t) => v)
    end
end

function evolve!(t::ReactionNetworkTransition, dt, tstep)
end

"Resample transition's parameters (incl. LHS). RHS remains not evaluated."
function resample!(t::ReactionNetworkTransition, network)
    for (k, v) in t.params_interpreted
        k ∈ [:trans] && continue
        t[k] = t.params_interpreted[v](network, t)
    end

    # sample the LHS, RHS of a reaction
    lhs, rhs = split_rline(t[:trans])

    reactants = extract_reactants(lhs, network)
    foreach(reactants) do r
        r.stoich = interpret_eval(r.stoich)(t, network)
        r.modalities = r.modalities ∪ get_modalities(r.agent, network)
    end
    t[:LHS] = reactants
    t[:RHS] = rhs
end

function resample!(state::DyVEState, o::Sampleable)
    o.last = state.t
    isempty(o.range) && (return o.val = missing)

    o.sampled = context_eval(state, sample_range(o.range, state))
end

resample(state::DyVEState, o::Symbol) = resample!(state, state.sampleables[o])

function update_observables(state::DyVEState)
    foreach(o -> (state.t - o.last) >= o.every && resample!(state, o),
            values(state.sampleables))
end
