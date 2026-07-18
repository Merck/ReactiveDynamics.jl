export AbstractStructuredToken, BaseStructuredToken
export @structured_token
export register_structured_species!, add_structured_token!

"""
    AbstractStructuredToken <: AbstractAlgebraicAgent

Abstract supertype of every structured (agentic) token kind (ADR 0006/0008). A structured token is a first-class entity — it carries its own attributes, a stable identity, and a history, and moves through a lifecycle (created → bound to a transition as a consumed resource → advanced/retired) rather than being an anonymous unit of a plain-species count. Each concrete kind is an `AlgebraicAgents.@aagent`, so a token is a genuine node in the AA hierarchy under the problem's `"structured"` container. Define a kind with [`@structured_token`](@ref), register it with [`register_structured_species!`](@ref), and add instances with [`add_structured_token!`](@ref) or the declarative [`PopulationEntry`](@ref) initial marking. Kinds get their behavior for free (no evolution rule by default — the orchestrator advances them); override [`log_token_fields`](@ref) to record fields into the per-token trajectory log.
"""
abstract type AbstractStructuredToken <: AbstractAlgebraicAgent end

"""
    BaseStructuredToken <: AbstractStructuredToken

The base structured-token layout every `@structured_token` kind inherits (via `@aagent FreeAgent`). It supplies the protocol fields the engine relies on: `species` (the `:S` row the token currently occupies, or `:removed` once soft-retired), `bound_transition` (the [`Transition`](@ref) instance the token is currently committed to as a consumed resource, or `nothing` when free), and `past_bonds` (the token's audit history — the `(species, t, transition)` triples of every transition it has been used in). A concrete kind adds its own modeling attributes (e.g. `phase`, `npv`) on top of these; the base layout is what makes a token bindable, selectable, and auditable without per-kind boilerplate.
"""
@aagent FreeAgent struct BaseStructuredToken
    species::Union{Nothing, Symbol}
    bound_transition::Union{Nothing, ReactiveDynamics.Transition}
    past_bonds::Vector{Tuple{Symbol, Float64, Transition}}
end

"""
    register_structured_species!(net, type)

Register the structured-token kind `type` (a `Symbol`) as a species of the static network `net`, adding a `:S` row named `type` if one does not already exist and flagging it `specStructured = true`. This is what tells the engine that occupants of that place are first-class token agents (counted from the `"structured"` container), not a plain scalar count. Returns `nothing`. The `@register` sugar and [`@structured_token`](@ref) (which defines the host struct) are the usual companions; a kind must be registered before instances can be added to a `ReactionNetworkProblem`.
"""
function register_structured_species!(reaction_network, type)
    if !(type ∈ reaction_network[:, :specName])
        add_row!(reaction_network, :S; specName = type)
    end

    i = first(find_rows(reaction_network, type, :specName))
    reaction_network[i, :specStructured] = true

    return nothing
end

"""
    @structured_token net Kind

Define a new structured-token kind named `Kind` — a concrete subtype of [`AbstractStructuredToken`](@ref) built off [`BaseStructuredToken`](@ref)'s protocol fields (via `AlgebraicAgents.aagent`). This declares the host struct for the kind; add your own modeling attributes to the generated type, register it as a species with [`register_structured_species!`](@ref), and instantiate it with [`add_structured_token!`](@ref) or a [`PopulationEntry`](@ref).
"""
macro structured_token(network, type)
    return quote
        $(
            AlgebraicAgents.aagent(
                BaseStructuredToken,
                AbstractStructuredToken,
                type,
                ReactiveDynamics,
            )
        )
    end
end

"""
    add_structured_token!(problem, agent) -> agent

Add a structured-token instance `agent` to the live `problem`: entangle it under the `"structured"` AA container and assign it a per-species monotonic creation index (ADR 0006 §E), so the deterministic `(species, creation_index)` selection order every token operation relies on is well-defined. Returns the added `agent`. This is the imperative counterpart to declaring the token in the [`PopulationEntry`](@ref) initial marking; both funnel through here so creation indices are contiguous and reproducible under `(model, seed)`.
"""
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

"""
    PopulationEntry(species, kind; count = 1, attributes = Dict{Symbol, Any}())

A declarative initial-marking entry (ADR 0007 §B): `count` instances of the structured `kind` (a registry key resolved to a host constructor), each token placed on `species` and built from `attributes` — a `Dict` of `field => Expr`-or-literal evaluated once at t=0 through the state's seeded `rng` (§4 D5). It is the structured-token analogue of `specInitVal` for plain species: put `PopulationEntry`s in the `population=` vector of `ReactionNetworkProblem` to declare the t=0 marking declaratively, so it re-instantiates cleanly on `reinit!` / across ensemble members. The alternative form — an already-constructed host token agent — is placed in the `population` vector directly (no `PopulationEntry`, and the registry is not consulted).
"""
struct PopulationEntry
    species::Symbol
    kind::Symbol
    count::Int
    attributes::Dict{Symbol, Any}
end
PopulationEntry(species, kind; count = 1, attributes = Dict{Symbol, Any}()) =
    PopulationEntry(species, kind, count, Dict{Symbol, Any}(attributes))

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
            for _ in 1:entry.count
                fields = Dict{Symbol, Any}(
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
const _SNAPSHOT_SKIP = (
    :uuid, :name, :parent, :inners, :relpathrefs, :opera,
    :bound_transition, :past_bonds,
)
function snapshot_population!(problem::ReactionNetworkProblem)
    empty!(problem.init_snapshot)
    for tok in values(inners(getagent(problem, "structured")))
        fields = filter(f -> !(f in _SNAPSHOT_SKIP), fieldnames(typeof(tok)))
        problem.init_snapshot[AlgebraicAgents.getname(tok)] =
            Dict{Symbol, Any}(f => getproperty(tok, f) for f in fields)
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
function set_bound_transition!(a::AbstractStructuredToken, t::Union{Nothing, Transition})
    return a.bound_transition = t
end

# Priority with which an unbound agent will be assigned to a transition.
priority(a::AbstractStructuredToken, transition) = 0.0

export log_token_fields

"""
    log_token_fields(token) -> NamedTuple

Per-token trajectory-log hook (ADR 0013 §A1 / CONTRACT §14.1). A host token KIND overrides this to declare WHICH fields the orchestrator records into `state.token_trajectory` each tick, e.g.

    ReactiveDynamics.log_token_fields(t::ProjectToken) = (; t.phase, t.npv_peak, t.pos_remaining)

The default logs NOTHING — the trajectory log is bounded by per-kind opt-in (Invariant 2), so it does not grow for kinds that don't opt in (load-bearing because retired tokens are KEPT under the Milestone-1 soft-`:removed` decision). The hook MUST be 𝓕ₜ-measurable: a pure field read, no RNG, no future (Invariant 1) — it returns the snapshot rather than holding it, keeping storage central.
"""
log_token_fields(::AbstractStructuredToken) = NamedTuple()

# What species (place) is an agent currently assigned to.
get_species(a::AbstractStructuredToken) = a.species
set_species!(a::AbstractStructuredToken, species::Symbol) = a.species = species
