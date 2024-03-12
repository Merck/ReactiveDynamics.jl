export AbstractStructuredToken, BaseStructuredToken
export @structured_token
export add_structured_token!

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
        add_part!(reaction_network, :S; specName = type)
    end

    i = first(incident(reaction_network, type, :specName))
    reaction_network[i, :specStructured] = true

    return nothing
end

# Convenience macro to define structured species.
macro structured_token(network, type)
    quote
        $(AlgebraicAgents.aagent(
            BaseStructuredToken,
            AbstractStructuredToken,
            type,
            ReactiveDynamics,
        ))
    end
end

# Add a structured agent instance to an instance of a reaction network.
function add_structured_token!(problem::ReactionNetworkProblem, agent)
    return entangle!(getagent(problem, "structured"), agent)
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

# What species (place) is an agent currently assigned to.
get_species(a::AbstractStructuredToken) = a.species
set_species!(a::AbstractStructuredToken, species::Symbol) = a.species = species
