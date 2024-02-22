export AbstractStructuredSpecies, BaseStructuredSpecies
export @structured
export add_structured_species!

# Abstract supertype of all structured species.
abstract type AbstractStructuredSpecies <: AbstractAlgebraicAgent end

# It comes handy to keep track of the transition the entity is assigned to (if).
# In general, we will probably assume that each "structured agent" type implements this field.
# Otherwise, it would be possible to implement getter and setter interface and use it from within ReaDyn.
@aagent FreeAgent struct BaseStructuredSpecies
    bound_transition::Union{Nothing,ReactiveDynamics.Transition}
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
macro structured(network, type)
    name = Docs.namify(type.args[2])

    quote
        $(AlgebraicAgents.aagent(
            BaseStructuredSpecies,
            AbstractStructuredSpecies,
            type,
            ReactiveDynamics,
        ))
        register_structured_species!($(esc(network)), $(QuoteNode(name)))
    end
end

# Add a structured agent instance to an instance of a reaction network.
function add_structured_species!(problem::ReactionNetworkProblem, agent)
    return entangle!(getagent(problem, "structured/$(nameof(typeof(agent)))"), agent)
end

import AlgebraicAgents

# By default, structured agents have no evolutionary rule.
AlgebraicAgents._projected_to(::AbstractStructuredSpecies) = nothing
AlgebraicAgents._step!(::AbstractStructuredSpecies) = nothing

# Tell if an agent is assigned to a transition, as a resource.
isblocked(a) = !isnothing(a.bound_transition)

# Priority with which an unbound agent will be assigned to a transition.
priority(a, transition) = 0.0
