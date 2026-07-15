export AbstractStructuredToken, BaseStructuredToken
export @structured_token
export register_structured_species!, add_structured_token!

# Abstract supertype of all structured species.
abstract type AbstractStructuredToken <: AbstractAlgebraicAgent end

# It comes handy to keep track of the transition the entity is assigned to (if).
# In general, we will probably assume that each "structured agent" type implements this field.
# Otherwise, it would be possible to implement getter and setter interface and use it from within ReaDyn.
@aagent FreeAgent struct BaseStructuredToken
    species::Union{Nothing,Symbol}
    bound_transition::Union{Nothing,ReactiveDynamics.Transition}
    past_bonds::Vector{Tuple{Symbol,Float64,Transition}}
end

# We use this to let the network know that the type is structured.
function register_structured_species!(reaction_network, type)
    if !(type ∈ reaction_network[:, :specName])
        add_row!(reaction_network, :S; specName = type)
    end

    i = first(find_rows(reaction_network, type, :specName))
    reaction_network[i, :specStructured] = true

    return nothing
end

# Convenience macro to define structured species.
macro structured_token(network, type)
    return quote
        $(AlgebraicAgents.aagent(
            BaseStructuredToken,
            AbstractStructuredToken,
            type,
            ReactiveDynamics,
        ))
    end
end

# Add a structured agent instance to an instance of a reaction network. Assigns the token a
# per-species monotonic creation index (ADR 0006 §E) so the deterministic (species,
# creation_index) selection order is well-defined before recording it in the network.
function add_structured_token!(problem::ReactionNetworkProblem, agent)
    entangle!(getagent(problem, "structured"), agent)
    sp = get_species(agent)
    if sp !== nothing
        k = get(problem.creation_counters, sp, 0) + 1
        problem.creation_counters[sp] = k
        problem.creation_index[AlgebraicAgents.getname(agent)] = k
    end
    return agent
end

export PopulationEntry

# A declarative initial-marking entry (ADR 0007 §B): `count` instances of a structured `kind`
# (a registry key → host constructor), each built from `attributes` (a Dict field => Expr/literal,
# sampled once at t=0 through state.rng — §4 D5). The structured analogue of specInitVal. The
# alternative authoring form — already-constructed host token agents — is just put in the
# `population` vector directly (no PopulationEntry needed; the registry isn't consulted).
struct PopulationEntry
    species::Symbol
    kind::Symbol
    count::Int
    attributes::Dict{Symbol,Any}
end
PopulationEntry(species, kind; count = 1, attributes = Dict{Symbol,Any}()) =
    PopulationEntry(species, kind, count, Dict{Symbol,Any}(attributes))

# Evaluate an initial-marking attribute through the seeded closure path; a bare QuoteNode (a
# literal `:Phase2`) is the symbol it wraps (wrap_fun/context_eval pass QuoteNodes through).
function _eval_attr(problem::ReactionNetworkProblem, v)
    v isa QuoteNode && return v.value
    r = context_eval(problem, nothing, problem.wrap_fun(v))
    return r isa QuoteNode ? r.value : r
end

# Instantiate the declarative initial marking into the network's structured container, in
# DECLARED ORDER (ADR 0007 §B instantiation contract). For a PopulationEntry the `kind` is
# resolved against the registry to a host constructor `(state, fields::Dict) -> token`, and each
# attribute expr is evaluated through the seeded closure; an already-constructed token agent is
# entangled as-is. Creation indices are assigned by add_structured_token! in this same order.
function instantiate_population!(problem::ReactionNetworkProblem)
    for entry in problem.population
        if entry isa PopulationEntry
            haskey(problem.registry, entry.kind) || error(
                "population: kind $(entry.kind) not in the network registry (ADR 0006 §C)",
            )
            ctor = problem.registry[entry.kind]
            for _ = 1:entry.count
                fields = Dict{Symbol,Any}(
                    f => _eval_attr(problem, v) for (f, v) in entry.attributes
                )
                add_structured_token!(problem, ctor(problem, fields))
            end
        else
            # an already-constructed host token agent
            add_structured_token!(problem, entry)
        end
    end
    return problem
end

# Record each live token's current field values keyed by token name, so the explicit-host-token
# population form can be restored to its t=0 attributes on reinit! (the PopulationEntry form
# rebuilds fresh tokens instead, so it needs no snapshot). Captures the modeling attributes plus
# the protocol `species` field (a soft-retired token has species==:removed and must be restored).
const _SNAPSHOT_SKIP = (:uuid, :name, :parent, :inners, :relpathrefs, :opera,
    :bound_transition, :past_bonds)
function snapshot_population!(problem::ReactionNetworkProblem)
    empty!(problem.init_snapshot)
    for tok in values(inners(getagent(problem, "structured")))
        fields = filter(f -> !(f in _SNAPSHOT_SKIP), fieldnames(typeof(tok)))
        problem.init_snapshot[AlgebraicAgents.getname(tok)] =
            Dict{Symbol,Any}(f => getproperty(tok, f) for f in fields)
    end
    return problem
end

# Restore an explicit-host-token's fields to its captured t=0 snapshot (used by reinit!).
function restore_token_snapshot!(problem::ReactionNetworkProblem, tok)
    snap = get(problem.init_snapshot, AlgebraicAgents.getname(tok), nothing)
    snap === nothing && return tok
    for (f, v) in snap
        setproperty!(tok, f, v)
    end
    return tok
end

import AlgebraicAgents

# By default, structured agents have no evolutionary rule.
AlgebraicAgents._projected_to(::AbstractStructuredToken) = nothing
AlgebraicAgents._step!(::AbstractStructuredToken) = nothing

# Tell if an agent is assigned to a transition, as a resource.
isblocked(a::AbstractStructuredToken) = !isnothing(get_bound_transition(a))

# Add a record that an agent was used as "species" in a "transition".
function add_to_log!(a::AbstractStructuredToken, species::Symbol, t, transition::Transition)
    return push!(a.past_bonds, (species, Float64(t), transition))
end

# Set the transition a token is bound to.
get_bound_transition(a::AbstractStructuredToken) = a.bound_transition
function set_bound_transition!(a::AbstractStructuredToken, t::Union{Nothing,Transition})
    return a.bound_transition = t
end

# Priority with which an unbound agent will be assigned to a transition.
priority(a::AbstractStructuredToken, transition) = 0.0

export log_token_fields

# Per-token trajectory-log hook (ADR 0013 §A1 / CONTRACT §14.1). A host KIND overrides this to
# declare WHICH fields the orchestrator records into `state.token_trajectory` each tick, e.g.
#
#     ReactiveDynamics.log_token_fields(t::ProjectToken) = (; t.phase, t.npv_peak, t.pos_remaining)
#
# The default logs NOTHING — the trajectory log is bounded by per-kind opt-in (Invariant 2), so it
# does not grow for kinds that don't opt in (load-bearing because retired tokens are KEPT under the
# Milestone-1 soft-`:removed` decision). The hook MUST be 𝓕ₜ-measurable: a pure field read, no RNG,
# no future (Invariant 1) — it returns the snapshot rather than holding it, keeping storage central.
log_token_fields(::AbstractStructuredToken) = NamedTuple()

# What species (place) is an agent currently assigned to.
get_species(a::AbstractStructuredToken) = a.species
set_species!(a::AbstractStructuredToken, species::Symbol) = a.species = species
